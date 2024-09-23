import logging
import random
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn import preprocessing
from sklearn.decomposition import PCA

from gsMap.GNN_VAE.adjacency_matrix import Construct_Adjacency_Matrix
from gsMap.GNN_VAE.train import Model_Train
from gsMap.config import FindLatentRepresentationsConfig

logger = logging.getLogger(__name__)


def set_seed(seed_value):
    """
    Set seed for reproducibility in PyTorch and other libraries.
    """
    torch.manual_seed(seed_value)
    np.random.seed(seed_value)
    random.seed(seed_value)
    if torch.cuda.is_available():
        logger.info('Using GPU for computations.')
        torch.cuda.manual_seed(seed_value)
        torch.cuda.manual_seed_all(seed_value)
    else:
        logger.info('Using CPU for computations.')


class LatentRepresentationFinder:
    def __init__(self, adata, args: FindLatentRepresentationsConfig):
        self.adata = adata
        self.Params = args

        # Preprocess data without copying the entire AnnData object
        if self.Params.data_layer in ['count', 'counts']:
            self.adata.X = self.adata.layers[self.Params.data_layer]
            sc.pp.highly_variable_genes(
                self.adata, flavor="seurat_v3", n_top_genes=self.Params.feat_cell
            )
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
            sc.pp.scale(self.adata)
        else:
            if self.Params.data_layer != 'X':
                self.adata.X = self.adata.layers[self.Params.data_layer]
            sc.pp.highly_variable_genes(
                self.adata, n_top_genes=self.Params.feat_cell
            )

        # Process the feature matrix
        self.node_X = self.adata[:, self.adata.var.highly_variable].X
        self.latent_pca = None  # Initialize latent_pca to None

    def compute_pca(self):
        """
        Compute PCA if it hasn't been computed yet.
        """
        if self.latent_pca is None:
            logger.info('Computing PCA...')
            pca = PCA(n_components=self.Params.n_comps)
            self.latent_pca = pca.fit_transform(self.node_X)
            logger.info('PCA computation completed.')
        return self.latent_pca

    def run_gnn_vae(self, label, verbose='whole ST data'):
        # Construct the neighboring graph
        graph_dict = Construct_Adjacency_Matrix(self.adata, self.Params)

        # Use PCA if specified
        if self.Params.input_pca:
            node_X = self.compute_pca()
        else:
            node_X = self.node_X

        # Update the input shape
        self.Params.n_nodes = node_X.shape[0]
        self.Params.feat_cell = node_X.shape[1]

        # Run GNN-VAE
        logger.info(f'------Finding latent representations for {verbose}...')
        gvae = Model_Train(node_X, graph_dict, self.Params, label)
        gvae.run_train()

        # Clean up to free memory
        del graph_dict

        return gvae.get_latent()

    def run_pca(self):
        return self.compute_pca()


def run_find_latent_representation(args: FindLatentRepresentationsConfig):
    set_seed(2024)
    num_features = args.feat_cell
    Path(args.hdf5_with_latent_path).parent.mkdir(
        parents=True, exist_ok=True, mode=0o755
    )

    # Load the ST data
    logger.info(f'------Loading ST data of {args.sample_name}...')
    adata = sc.read_h5ad(f'{args.input_hdf5_path}')
    adata.var_names_make_unique()
    adata.X = adata.layers[args.data_layer] if args.data_layer in adata.layers else adata.X
    logger.info(f'The ST data contains {adata.shape[0]} cells, {adata.shape[1]} genes.')

    # Load the cell type annotation
    if args.annotation is not None:
        # Remove cells without enough annotations
        adata = adata[~adata.obs[args.annotation].isnull(), :]
        num = adata.obs[args.annotation].value_counts()
        valid_annotations = num[num >= 30].index.to_list()
        adata = adata[adata.obs[args.annotation].isin(valid_annotations)]

        le = preprocessing.LabelEncoder()
        adata.obs['categorical_label'] = le.fit_transform(adata.obs[args.annotation])
        label = adata.obs['categorical_label'].to_numpy()
    else:
        label = None

    # Find latent representations
    latent_rep = LatentRepresentationFinder(adata, args)
    latent_gvae = latent_rep.run_gnn_vae(label)
    latent_pca = latent_rep.run_pca()

    # Add latent representations to the AnnData object
    logger.info('------Adding latent representations...')
    adata.obsm["latent_GVAE"] = latent_gvae
    adata.obsm["latent_PCA"] = latent_pca

    # Run UMAP based on latent representations
    for name in ['latent_GVAE', 'latent_PCA']:
        sc.pp.neighbors(adata, n_neighbors=10, use_rep=name)
        sc.tl.umap(adata)
        adata.obsm['X_umap_' + name] = adata.obsm['X_umap']

    # Find the latent representations hierarchically (optionally)
    if args.annotation is not None and args.hierarchically:
        logger.info('------Finding latent representations hierarchically...')
        num_cells = adata.shape[0]
        n_comps = args.n_comps

        latent_gvae_hierarchy = np.zeros((num_cells, n_comps))
        latent_pca_hierarchy = np.zeros((num_cells, n_comps))

        # Mapping from cell names to indices
        cell_indices = {cell: idx for idx, cell in enumerate(adata.obs_names)}

        for ct in adata.obs[args.annotation].unique():
            adata_part = adata[adata.obs[args.annotation] == ct, :]
            logger.info(f'Processing {ct} with shape {adata_part.shape}')

            # Find latent representations for the selected cell type
            latent_rep_part = LatentRepresentationFinder(adata_part, args)

            # Use the already computed PCA if possible
            latent_pca_part = latent_rep_part.compute_pca()

            if adata_part.shape[0] <= args.n_comps:
                latent_gvae_part = latent_pca_part
            else:
                latent_gvae_part = latent_rep_part.run_gnn_vae(label=None, verbose=ct)

            # Map indices and store the representations
            indices = [cell_indices[cell] for cell in adata_part.obs_names]
            latent_gvae_hierarchy[indices, :] = latent_gvae_part
            latent_pca_hierarchy[indices, :] = latent_pca_part

            # Clean up to free memory
            del adata_part, latent_rep_part, latent_pca_part, latent_gvae_part

            # Reset feat_cell to original value
            args.feat_cell = num_features

        adata.obsm["latent_GVAE_hierarchy"] = latent_gvae_hierarchy
        adata.obsm["latent_PCA_hierarchy"] = latent_pca_hierarchy

    logger.info('------Saving ST data...')
    adata.write(args.hdf5_with_latent_path)

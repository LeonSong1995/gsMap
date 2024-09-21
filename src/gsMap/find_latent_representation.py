import logging
import random

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn import preprocessing

from gsMap.GNN_VAE.adjacency_matrix import Construct_Adjacency_Matrix
from gsMap.GNN_VAE.train import Model_Train
from gsMap.config import FindLatentRepresentationsConfig

logger = logging.getLogger(__name__)

def set_seed(seed_value):
    """
    Set seed for reproducibility in PyTorch.
    """
    torch.manual_seed(seed_value)  # Set the seed for PyTorch
    np.random.seed(seed_value)  # Set the seed for NumPy
    random.seed(seed_value)  # Set the seed for Python random module
    if torch.cuda.is_available():
        logger.info('Running use GPU')
        torch.cuda.manual_seed(seed_value)  # Set seed for all CUDA devices
        torch.cuda.manual_seed_all(seed_value)  # Set seed for all CUDA devices
    else:
        logger.info('Running use CPU')


# The class for finding latent representations
class Latent_Representation_Finder:

    def __init__(self, adata, args:FindLatentRepresentationsConfig):
        self.adata = adata.copy()
        self.Params = args

        # Standard process
        if self.Params.data_layer == 'count' or self.Params.data_layer == 'counts':
            self.adata.X = self.adata.layers[self.Params.data_layer]
            sc.pp.highly_variable_genes(self.adata, flavor="seurat_v3", n_top_genes=self.Params.feat_cell)
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
            sc.pp.scale(self.adata)
        else:
            if self.Params.data_layer != 'X':
                self.adata.X = self.adata.layers[self.Params.data_layer]
            sc.pp.highly_variable_genes(self.adata, n_top_genes=self.Params.feat_cell)

    def Run_GNN_VAE(self, label, verbose='whole ST data'):

        # Construct the neighbouring graph
        graph_dict = Construct_Adjacency_Matrix(self.adata, self.Params)

        # Process the feature matrix
        node_X = self.adata[:, self.adata.var.highly_variable].X
        logger.info(f'The shape of feature matrix is {node_X.shape}.')
        if self.Params.input_pca:
            node_X = sc.pp.pca(node_X, n_comps=self.Params.n_comps)

        # Update the input shape
        self.Params.n_nodes = node_X.shape[0]
        self.Params.feat_cell = node_X.shape[1]

        # Run GNN-VAE
        logger.info(f'------Finding latent representations for {verbose}...')
        gvae = Model_Train(node_X, graph_dict, self.Params, label)
        gvae.run_train()

        return gvae.get_latent()

    def Run_PCA(self):
        sc.tl.pca(self.adata)
        return self.adata.obsm['X_pca'][:, 0:self.Params.n_comps]


def run_find_latent_representation(args:FindLatentRepresentationsConfig):
    set_seed(2024)
    num_features = args.feat_cell
    args.hdf5_with_latent_path.parent.mkdir(parents=True, exist_ok=True,mode=0o755)
    # Load the ST data
    logger.info(f'------Loading ST data of {args.sample_name}...')
    adata = sc.read_h5ad(f'{args.input_hdf5_path}')
    adata.var_names_make_unique()
    adata.X = adata.layers[args.data_layer] if args.data_layer in adata.layers.keys() else adata.X
    logger.info('The ST data contains %d cells, %d genes.' % (adata.shape[0], adata.shape[1]))
    # Load the cell type annotation
    if not args.annotation is None:
        # remove cells without enough annotations
        adata = adata[~pd.isnull(adata.obs[args.annotation]), :]
        num = adata.obs[args.annotation].value_counts()
        adata = adata[adata.obs[args.annotation].isin(num[num >= 30].index.to_list())]

        le = preprocessing.LabelEncoder()
        le.fit(adata.obs[args.annotation])
        adata.obs['categorical_label'] = le.transform(adata.obs[args.annotation])
        label = adata.obs['categorical_label'].to_list()
    else:
        label = None
    # Find latent representations
    latent_rep = Latent_Representation_Finder(adata, args)
    latent_GVAE = latent_rep.Run_GNN_VAE(label)
    latent_PCA = latent_rep.Run_PCA()
    # Add latent representations to the spe data
    logger.info(f'------Adding latent representations...')
    adata.obsm["latent_GVAE"] = latent_GVAE
    adata.obsm["latent_PCA"] = latent_PCA
    # Run umap based on latent representations
    for name in ['latent_GVAE', 'latent_PCA']:
        sc.pp.neighbors(adata, n_neighbors=10, use_rep=name)
        sc.tl.umap(adata)
        adata.obsm['X_umap_' + name] = adata.obsm['X_umap']

        # Find the latent representations hierarchically (optionally)
    if not args.annotation is None and args.hierarchically:
        logger.info(f'------Finding latent representations hierarchically...')
        PCA_all = pd.DataFrame()
        GVAE_all = pd.DataFrame()

        for ct in adata.obs[args.annotation].unique():
            adata_part = adata[adata.obs[args.annotation] == ct, :]
            logger.info(adata_part.shape)

            # Find latent representations for the selected ct
            latent_rep = Latent_Representation_Finder(adata_part, args)

            latent_PCA_part = pd.DataFrame(latent_rep.Run_PCA())
            if adata_part.shape[0] <= args.n_comps:
                latent_GVAE_part = latent_PCA_part
            else:
                latent_GVAE_part = pd.DataFrame(latent_rep.Run_GNN_VAE(label=None, verbose=ct))

            latent_GVAE_part.index = adata_part.obs_names
            latent_PCA_part.index = adata_part.obs_names

            GVAE_all = pd.concat((GVAE_all, latent_GVAE_part), axis=0)
            PCA_all = pd.concat((PCA_all, latent_PCA_part), axis=0)

            args.feat_cell = num_features

            adata.obsm["latent_GVAE_hierarchy"] = np.array(GVAE_all.loc[adata.obs_names,])
            adata.obsm["latent_PCA_hierarchy"] = np.array(PCA_all.loc[adata.obs_names,])
    logger.info(f'------Saving ST data...')
    adata.write(args.hdf5_with_latent_path)


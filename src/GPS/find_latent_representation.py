import argparse
import logging
import pprint
import random
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn import preprocessing

from GPS.GNN_VAE.adjacency_matrix import Construct_Adjacency_Matrix
from GPS.GNN_VAE.train import Model_Train
from GPS.config import add_find_latent_representations_args, FindLatentRepresentationsConfig

# seed all

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def set_seed(seed_value):
    """
    Set seed for reproducibility in PyTorch.
    """
    torch.manual_seed(seed_value)  # Set the seed for PyTorch
    np.random.seed(seed_value)  # Set the seed for NumPy
    random.seed(seed_value)  # Set the seed for Python random module
    if torch.cuda.is_available():
        print('Running use GPU')
        torch.cuda.manual_seed(seed_value)  # Set seed for all CUDA devices
        torch.cuda.manual_seed_all(seed_value)  # Set seed for all CUDA devices
    else:
        print('Running use CPU')

set_seed(2024)

# The class for finding latent representations
class Latent_Representation_Finder:

    def __init__(self, adata, Params):
        self.adata = adata.copy()
        self.Params = Params

        # Standard process
        if self.Params.type == 'count' or self.Params.type == 'counts':
            self.adata.X = self.adata.layers[self.Params.type]
            sc.pp.highly_variable_genes(self.adata, flavor="seurat_v3", n_top_genes=self.Params.feat_cell)
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
            sc.pp.scale(self.adata)
        else:
            self.adata.X = self.adata.layers[self.Params.type]
            sc.pp.highly_variable_genes(self.adata, n_top_genes=self.Params.feat_cell)

    def Run_GNN_VAE(self, label, verbose='whole ST data'):

        # Construct the neighbouring graph
        graph_dict = Construct_Adjacency_Matrix(self.adata, self.Params)

        # Process the feature matrix
        node_X = self.adata[:, self.adata.var.highly_variable].X
        print(f'The shape of feature matrix is {node_X.shape}.')
        if self.Params.input_pca:
            node_X = sc.pp.pca(node_X, n_comps=self.Params.n_comps)

        # Update the input shape
        self.Params.n_nodes = node_X.shape[0]
        self.Params.feat_cell = node_X.shape[1]

        # Run GNN-VAE
        print(f'------Finding latent representations for {verbose}...')
        gvae = Model_Train(node_X, graph_dict, self.Params, label)
        gvae.run_train()

        return gvae.get_latent()

    def Run_PCA(self):
        sc.tl.pca(self.adata)
        return self.adata.obsm['X_pca'][:, 0:self.Params.n_comps]


def run_find_latent_representation(args:FindLatentRepresentationsConfig):
    num_features = args.feat_cell
    args.output_dir = Path(args.output_hdf5_path).parent
    args.output_dir.mkdir(parents=True, exist_ok=True,mode=0o755)
    # Load the ST data
    print(f'------Loading ST data of {args.sample_name}...')
    adata = sc.read_h5ad(f'{args.input_hdf5_path}')
    adata.var_names_make_unique()
    print('The ST data contains %d cells, %d genes.' % (adata.shape[0], adata.shape[1]))
    # Load the cell type annotation
    if not args.annotation is None:
        # remove cells without enough annotations
        adata = adata[~pd.isnull(adata.obs[args.annotation]), :]
        num = adata.obs[args.annotation].value_counts()
        adata = adata[adata.obs[args.annotation].isin(num[num >= 30].index.to_list()),]

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
    print(f'------Adding latent representations...')
    adata.obsm["latent_GVAE"] = latent_GVAE
    adata.obsm["latent_PCA"] = latent_PCA
    # Run umap based on latent representations
    for name in ['latent_GVAE', 'latent_PCA']:
        sc.pp.neighbors(adata, n_neighbors=10, use_rep=name)
        sc.tl.umap(adata)
        adata.obsm['X_umap_' + name] = adata.obsm['X_umap']

        # Find the latent representations hierarchically (optionally)
    if not args.annotation is None and args.hierarchically:
        print(f'------Finding latent representations hierarchically...')
        PCA_all = pd.DataFrame()
        GVAE_all = pd.DataFrame()

        for ct in adata.obs[args.annotation].unique():
            adata_part = adata[adata.obs[args.annotation] == ct, :]
            print(adata_part.shape)

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
    print(f'------Saving ST data...')
    adata.write(args.output_hdf5_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script is designed to find latent representations in spatial transcriptomics data using a Graph Neural Network Variational Autoencoder (GNN-VAE). It processes input data, constructs a neighboring graph, and runs GNN-VAE to output latent representations.")
    add_find_latent_representations_args(parser)
    TEST=True
    if TEST:
        test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021'
        name = 'Cortex_151507'


        args = parser.parse_args(
            [
                '--input_hdf5_path','/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/processed/h5ad/Cortex_151507.h5ad',
                '--output_hdf5_path',f'{test_dir}/{name}/hdf5/{name}_add_latent.h5ad',
                '--sample_name', name,
                '--annotation','layer_guess',
                '--type','count',
            ]

        )

    else:
        args = parser.parse_args()
    config=FindLatentRepresentationsConfig(**{'annotation': 'SubClass',
 'convergence_threshold': 0.0001,
 'epochs': 300,
 'feat_cell': 3000,
 'feat_hidden1': 256,
 'feat_hidden2': 128,
 'gcn_decay': 0.01,
 'gcn_hidden1': 64,
 'gcn_hidden2': 30,
 'gcn_lr': 0.001,
 'hierarchically': False,
 'input_hdf5_path': '/storage/yangjianLab/songliyang/SpatialData/Data/Brain/macaque/Cell/processed/h5ad/T862_macaque3.h5ad',
 'label_w': 1.0,
 'n_comps': 300,
 'n_neighbors': 11,
 'nheads': 3,
 'output_hdf5_path': 'T862_macaque3/find_latent_representations/T862_macaque3_add_latent.h5ad',
 'p_drop': 0.1,
 'rec_w': 1.0,
 'sample_name': 'T862_macaque3',
 'type': 'SCT',
 'var': False,
 'weighted_adj': False})
    # config=FindLatentRepresentationsConfig(**vars(args))
    start_time = time.time()
    logger.info(f'Find latent representations for {config.sample_name}...')
    pprint.pprint(config)
    run_find_latent_representation(config)
    end_time = time.time()
    logger.info(f'Find latent representations for {config.sample_name} finished. Time spent: {(end_time - start_time) / 60:.2f} min.')

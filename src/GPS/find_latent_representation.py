
import os
import random

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import preprocessing

from GPS.GNN_VAE.Adjacency_Matrix import Construct_Adjacency_Matrix
from GPS.GNN_VAE.Train import Model_Train

random.seed(20230609)

# Set args
import argparse
def add_find_latent_representations_args(parser):

    parser.add_argument('--epochs', default=300, type=int, help="Number of training epochs for the GNN-VAE model. Default is 300.")
    parser.add_argument('--feat_hidden1', default=256, type=int, help="Number of neurons in the first hidden layer of the feature extraction network. Default is 256.")
    parser.add_argument('--feat_hidden2', default=128, type=int, help="Number of neurons in the second hidden layer of the feature extraction network. Default is 128.")
    parser.add_argument('--feat_cell', default=3000, type=int, help="Number of top variable genes to select. Default is 3000.")
    parser.add_argument('--gcn_hidden1', default=64, type=int, help="Number of units in the first hidden layer of the GCN. Default is 64.")
    parser.add_argument('--gcn_hidden2', default=30, type=int, help="Number of units in the second hidden layer of the GCN. Default is 30.")
    parser.add_argument('--p_drop', default=0.1, type=float, help="Dropout rate used in the GNN-VAE model. Default is 0.1.")
    parser.add_argument('--gcn_lr', default=0.001, type=float, help="Learning rate for the GCN network. Default is 0.001.")
    parser.add_argument('--gcn_decay', default=0.01, type=float, help="Weight decay (L2 penalty) for the GCN network. Default is 0.01.")
    parser.add_argument('--n_neighbors', default=11, type=int, help="Number of neighbors to consider for graph construction in GCN. Default is 11.")
    parser.add_argument('--label_w', default=1, type=float, help="Weight of the label loss in the loss function. Default is 1.")
    parser.add_argument('--rec_w', default=1, type=float, help="Weight of the reconstruction loss in the loss function. Default is 1.")
    parser.add_argument('--input_pca', default=True, type=bool, help="Whether to perform PCA on input features. Default is True.")
    parser.add_argument('--n_comps', default=300, type=int, help="Number of principal components to keep if PCA is performed. Default is 300.")
    parser.add_argument('--weighted_adj', default=False, type=bool, help="Whether to use a weighted adjacency matrix in GCN. Default is False.")
    parser.add_argument('--var', default=False, type=bool, )
    parser.add_argument('--nheads', default=3, type=int, help="Number of heads in the attention mechanism of the GNN. Default is 3.")
    parser.add_argument('--convergence_threshold', default=1e-4, type=float, help="Threshold for convergence during training. Training stops if the loss change is below this threshold. Default is 1e-4.")
    parser.add_argument('--hierarchically', default=False, type=bool, help="Whether to find latent representations hierarchically. Default is False.")
    parser.add_argument('--output_dir', default=None, type=str, help="Path to the output directory. Default is None.")
def add_sample_info_args(parser):
    parser.add_argument('--sample_hdf5', default=None, type=str, help='Path to the sample hdf5 file', required=True)
    parser.add_argument('--sample_name', type=str, help='Name of the sample', required=True)
    parser.add_argument('--annotation_layer_name', default=None, type=str, help='Name of the annotation layer',dest='annotation')
    parser.add_argument('--is_count', action='store_true', help='Whether the data is count data')


# The class for finding latent representations
class Latent_Representation_Finder:

    def __init__(self, adata, Params):
        self.adata = adata.copy()
        self.Params = Params

        # Standard process
        if self.Params.is_count:
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


def run_find_latent_representation(args):
    num_features = args.feat_cell
    args.output_dir = os.path.join(args.output_dir,'find_latent_representations')
    # Load the ST data
    print(f'------Loading ST data of {args.sample_name}...')
    adata = sc.read_h5ad(f'{args.sample_hdf5}')
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

        # TODO : Don't know the meaning of the following code
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
    # Save the data
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, mode=0o777, exist_ok=True)
    data_name = args.sample_name
    adata.write(f'{args.output_dir}/{data_name}_add_latent.h5ad')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script is designed to find latent representations in spatial transcriptomics data using a Graph Neural Network Variational Autoencoder (GNN-VAE). It processes input data, constructs a neighboring graph, and runs GNN-VAE to output latent representations.")
    add_sample_info_args(parser)
    add_find_latent_representations_args(parser)
    TEST=True
    if TEST:
        test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021'
        name = 'Cortex_151507'


        args = parser.parse_args(
            [
                '--sample_hdf5','/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/processed/h5ad/Cortex_151507.h5ad',
                '--sample_name', name,
                '--annotation_layer_name','layer_guess',
                '--is_count',
                '--output_dir',f'{test_dir}/{name}',
            ]

        )
    else:
        args = parser.parse_args()
    run_find_latent_representation(args)

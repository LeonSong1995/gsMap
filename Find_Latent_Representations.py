#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:34:35 2023

@author: songliyang
"""

import scanpy as sc
import umap
import torch
import os
import anndata
import pandas as pd
import numpy as np
import plotnine as pn
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from skmisc.loess import loess
from sklearn import metrics
import sklearn.neighbors
from scipy.spatial.distance import pdist, squareform
from progress.bar import IncrementalBar
from scipy.sparse import csr_matrix
import scipy.sparse as sp
from scipy.spatial import distance_matrix
from sklearn import preprocessing
from torch.nn.modules.module import Module
from torch.nn.parameter import Parameter
import torch.nn as nn
import random
import argparse

# Import the In-house GNN_VAE 
import sys
sys.path.insert(0,'/storage/yangjianLab/songliyang/SpatialData/spatial_ldsc_v1/GNN_VAE/')
from Adjacency_Matrix import Construct_Adjacency_Matrix
from Train import Model_Train
random.seed(20230609)

# Set args
parser = argparse.ArgumentParser()
parser.add_argument('--epochs', default=300, type=int)
parser.add_argument('--feat_hidden1', default=256, type=int)
parser.add_argument('--feat_hidden2', default=128, type=int)
parser.add_argument('--feat_cell', default=3000, type=int)
parser.add_argument('--gcn_hidden1', default=64, type=int)
parser.add_argument('--gcn_hidden2', default=30, type=int)
parser.add_argument('--p_drop', default=0.1, type=float)
parser.add_argument('--gcn_lr', default=0.001, type=float)
parser.add_argument('--gcn_decay', default=0.01, type=float)
parser.add_argument('--n_neighbors', default=11, type=int)
parser.add_argument('--label_w', default=1, type=float)
parser.add_argument('--rec_w', default=1, type=float)
parser.add_argument('--input_pca', default=True,type=bool)
parser.add_argument('--n_comps', default=300, type=int)
parser.add_argument('--weighted_adj', default=False,type=bool)
parser.add_argument('--var', default=False,type=bool)
parser.add_argument('--type', default=None, type=str)
parser.add_argument('--nheads', default=3, type=int)
parser.add_argument('--convergence_threshold', default=1e-4, type=float)
parser.add_argument('--annotation', default=None, type=str)
parser.add_argument('--hierarchically', default=False, type=bool)
parser.add_argument('--spe_path', default=None, type=str)
parser.add_argument('--spe_name', default=None, type=str)
parser.add_argument('--spe_out', default=None, type=str)


# The class for finding latent representations
class Latent_Representation_Finder:
    
    def __init__(self, adata, Params):
        self.adata = adata.copy()
        self.Params = Params

        # Standard process
        if self.Params.type == 'count' or self.Params.type == 'counts':
            self.adata.X = self.adata.layers[self.Params.type]
            sc.pp.highly_variable_genes(self.adata, flavor="seurat_v3",n_top_genes = self.Params.feat_cell)
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
            sc.pp.scale(self.adata)
        else:
            self.adata.X = self.adata.layers[self.Params.type]
            sc.pp.highly_variable_genes(self.adata,n_top_genes= self.Params.feat_cell)
    
    
    def Run_GNN_VAE(self,label,verbose='whole ST data'):

        # Construct the neighbouring graph 
        graph_dict = Construct_Adjacency_Matrix(self.adata,self.Params)
        
        # Process the feature matrix
        node_X =  self.adata[:, self.adata.var.highly_variable].X
        print(f'The shape of feature matrix is {node_X.shape}.')
        if self.Params.input_pca:
            node_X = sc.pp.pca(node_X, n_comps = self.Params.n_comps)

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
        return self.adata.obsm['X_pca'][:,0:self.Params.n_comps] 


    
if __name__ == '__main__':

    args = parser.parse_args()
    num_features = args.feat_cell
    
    # Load the ST data
    print(f'------Loading ST data of {args.spe_name}...')
    adata = sc.read_h5ad(f'{args.spe_path}/{args.spe_name}')
    adata.var_names_make_unique()
    print('The ST data contains %d cells, %d genes.' %(adata.shape[0], adata.shape[1]))

    # Load the cell type annotation
    if not args.annotation is None:
        # remove cells without enough annotations
        adata = adata[~pd.isnull(adata.obs[args.annotation]),:]
        num = adata.obs[args.annotation].value_counts()
        adata = adata[adata.obs[args.annotation].isin(num[num >= 30].index.to_list()),]
        
        le = preprocessing.LabelEncoder()
        le.fit(adata.obs[args.annotation])
        adata.obs['categorical_label'] = le.transform(adata.obs[args.annotation])
        label = adata.obs['categorical_label'].to_list()
    else:
        label = None
        
    # Find latent representations
    latent_rep = Latent_Representation_Finder(adata,args)
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
            adata_part = adata[adata.obs[args.annotation]==ct,:]
            print(adata_part.shape)
            
            # Find latent representations for the selected ct
            latent_rep = Latent_Representation_Finder(adata_part,args)
            
            latent_PCA_part = pd.DataFrame(latent_rep.Run_PCA())
            if adata_part.shape[0] <= args.n_comps:
                latent_GVAE_part = latent_PCA_part
            else:
                latent_GVAE_part = pd.DataFrame(latent_rep.Run_GNN_VAE(label=None,verbose=ct))
            
            latent_GVAE_part.index = adata_part.obs_names
            latent_PCA_part.index = adata_part.obs_names
            
            GVAE_all = pd.concat((GVAE_all,latent_GVAE_part),axis=0)    
            PCA_all = pd.concat((PCA_all,latent_PCA_part),axis=0)
            
            args.feat_cell = num_features
            
            adata.obsm["latent_GVAE_hierarchy"] = np.array(GVAE_all.loc[adata.obs_names,])
            adata.obsm["latent_PCA_hierarchy"] = np.array(PCA_all.loc[adata.obs_names,])
    
    
    print(f'------Saving ST data...')
    # Save the data
    if not os.path.exists(args.spe_out):
        os.makedirs(args.spe_out, mode=0o777, exist_ok=True)
        
    data_name = args.spe_name.split('.h5ad')[0]    
    adata.write(f'{args.spe_out}/{data_name}_add_latent.h5ad')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 21:31:27 2023

@author: songliyang
"""
import pandas as pd
import numpy as np
import sklearn.neighbors
import torch
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import distance_matrix


def Cal_Spatial_Net(adata, n_neighbors=5, verbose=True):
    """\
    Construct the spatial neighbor networks.
    """
    #- 
    if verbose:
        print('------Calculating spatial graph...')
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    #- 
    nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors).fit(coor)
    #- 
    distances, indices = nbrs.kneighbors(coor, return_distance=True)
    KNN_list = []
    for it in range(indices.shape[0]):
        KNN_list.append(pd.DataFrame(zip([it]*indices[it].shape[0], indices[it], distances[it])))
    #- 
    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']
    #- 
    Spatial_Net = KNN_df.copy()
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance']>0,]
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index), ))
    Spatial_Net['Cell1'] = Spatial_Net['Cell1'].map(id_cell_trans)
    Spatial_Net['Cell2'] = Spatial_Net['Cell2'].map(id_cell_trans)
    #- 
    return Spatial_Net


def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse.FloatTensor(indices, values, shape)


def preprocess_graph(adj):
    adj = sp.coo_matrix(adj)
    adj_ = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj_.sum(1))
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5).flatten())
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
    return sparse_mx_to_torch_sparse_tensor(adj_normalized)



def Construct_Adjacency_Matrix(adata,Params, verbose=True):
    # Construct the neighbor graph 
    Spatial_Net = Cal_Spatial_Net(adata, n_neighbors=Params.n_neighbors)
    #- 
    if verbose:
        print('The graph contains %d edges, %d cells.' %(Spatial_Net.shape[0], adata.n_obs))
        print('%.4f neighbors per cell on average.' %(Spatial_Net.shape[0]/adata.n_obs))
    #-  
    cells = np.array(adata.obs.index)
    cells_id_tran = dict(zip(cells, range(cells.shape[0])))
    #- 
    G_df = Spatial_Net.copy()
    G_df['Cell1'] = G_df['Cell1'].map(cells_id_tran)
    G_df['Cell2'] = G_df['Cell2'].map(cells_id_tran)
    #- 
    if Params.weighted_adj:
        distance_normalized = G_df.Distance/(max(G_df.Distance)+1)
        adj_org = sp.coo_matrix((np.exp(-distance_normalized**2/(2)), (G_df['Cell1'], G_df['Cell2'])), shape=(adata.n_obs, adata.n_obs))
    else:
        adj_org = sp.coo_matrix((np.ones(G_df.shape[0]), (G_df['Cell1'], G_df['Cell2'])), shape=(adata.n_obs, adata.n_obs))
    #- 
    adj_m1 = adj_org
    adj_norm_m1 = preprocess_graph(adj_m1)
    adj_label_m1 = adj_m1 + sp.eye(adj_m1.shape[0])
    norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - adj_m1.sum()) * 2)
    #- 
    graph_dict = {
        "adj_org": adj_org,
        "adj_norm": adj_norm_m1,
        "norm_value": norm_m1
    }
    #- 
    return graph_dict

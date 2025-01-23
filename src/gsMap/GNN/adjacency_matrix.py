import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from sklearn.neighbors import NearestNeighbors


def cal_spatial_net(adata, n_neighbors=5, verbose=True):
    """Construct the spatial neighbor network."""
    if verbose:
        print("------Calculating spatial graph...")
    coor = pd.DataFrame(adata.obsm["spatial"], index=adata.obs.index)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coor)
    distances, indices = nbrs.kneighbors(coor)
    n_cells, n_neighbors = indices.shape
    cell_indices = np.arange(n_cells)
    cell1 = np.repeat(cell_indices, n_neighbors)
    cell2 = indices.flatten()
    distance = distances.flatten()
    knn_df = pd.DataFrame({"Cell1": cell1, "Cell2": cell2, "Distance": distance})
    knn_df = knn_df[knn_df["Distance"] > 0].copy()
    cell_id_map = dict(zip(cell_indices, coor.index, strict=False))
    knn_df["Cell1"] = knn_df["Cell1"].map(cell_id_map)
    knn_df["Cell2"] = knn_df["Cell2"].map(cell_id_map)
    return knn_df


def sparse_mx_to_torch_sparse_tensor(sparse_mx):
    """Convert a scipy sparse matrix to a torch sparse tensor."""
    sparse_mx = sparse_mx.tocoo().astype(np.float32)
    indices = torch.from_numpy(np.vstack((sparse_mx.row, sparse_mx.col)).astype(np.int64))
    values = torch.from_numpy(sparse_mx.data)
    shape = torch.Size(sparse_mx.shape)
    return torch.sparse_coo_tensor(indices, values, shape)


def preprocess_graph(adj):
    """Symmetrically normalize the adjacency matrix."""
    adj = sp.coo_matrix(adj)
    adj_ = adj + sp.eye(adj.shape[0])
    rowsum = np.array(adj_.sum(1)).flatten()
    degree_mat_inv_sqrt = sp.diags(np.power(rowsum, -0.5))
    adj_normalized = adj_.dot(degree_mat_inv_sqrt).transpose().dot(degree_mat_inv_sqrt).tocoo()
    return sparse_mx_to_torch_sparse_tensor(adj_normalized)


def construct_adjacency_matrix(adata, params, verbose=True):
    """Construct the adjacency matrix from spatial data."""
    spatial_net = cal_spatial_net(adata, n_neighbors=params.n_neighbors, verbose=verbose)
    if verbose:
        num_edges = spatial_net.shape[0]
        num_cells = adata.n_obs
        print(f"The graph contains {num_edges} edges, {num_cells} cells.")
        print(f"{num_edges / num_cells:.2f} neighbors per cell on average.")
    cell_ids = {cell: idx for idx, cell in enumerate(adata.obs.index)}
    spatial_net["Cell1"] = spatial_net["Cell1"].map(cell_ids)
    spatial_net["Cell2"] = spatial_net["Cell2"].map(cell_ids)
    if params.weighted_adj:
        distance_normalized = spatial_net["Distance"] / (spatial_net["Distance"].max() + 1)
        weights = np.exp(-0.5 * distance_normalized**2)
        adj_org = sp.coo_matrix(
            (weights, (spatial_net["Cell1"], spatial_net["Cell2"])),
            shape=(adata.n_obs, adata.n_obs),
        )
    else:
        adj_org = sp.coo_matrix(
            (np.ones(spatial_net.shape[0]), (spatial_net["Cell1"], spatial_net["Cell2"])),
            shape=(adata.n_obs, adata.n_obs),
        )
    adj_norm = preprocess_graph(adj_org)
    norm_value = adj_org.shape[0] ** 2 / ((adj_org.shape[0] ** 2 - adj_org.sum()) * 2)
    graph_dict = {"adj_org": adj_org, "adj_norm": adj_norm, "norm_value": norm_value}
    return graph_dict

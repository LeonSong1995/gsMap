#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:34:35 2023

@author: songliyang
"""

import scanpy as sc
import os
import anndata
import pandas as pd
import numpy as np
import multiprocessing
import plotnine as pn
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from sklearn import metrics
import sklearn
import sklearn.neighbors
from scipy.spatial.distance import pdist, squareform
from progress.bar import IncrementalBar
from math import floor
from scipy.sparse import issparse, vstack
from scipy import stats
from multiprocessing import Pool
from tqdm import tqdm
import tqdm.contrib.concurrent
from scipy.stats import rankdata
from sklearn.neighbors import NearestNeighbors
import argparse
from scipy.stats import gmean

from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import cosine_similarity


def find_Neighbors(coor, num_neighbour):
    """
    find Neighbors of each cell (based on spatial coordinates)
    """
    nbrs = NearestNeighbors(n_neighbors=num_neighbour).fit(coor)
    distances, indices = nbrs.kneighbors(coor, return_distance=True)

    KNN_list = [pd.DataFrame(zip([it] * indices[it].shape[0], indices[it], distances[it])) for it in
                range(indices.shape[0])]
    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']

    spatial_net = KNN_df.copy()
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index)))

    spatial_net['Cell1'] = spatial_net['Cell1'].map(id_cell_trans)
    spatial_net['Cell2'] = spatial_net['Cell2'].map(id_cell_trans)

    return spatial_net


def _build_spatial_net(adata, annotation, num_neighbour):
    """
    1 Build spatial neighbourhood matrix for each spot (cell) based on the spatial coord
    """
    print(f'------Building spatial graph based on spatial coordinates...')

    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index

    if not annotation is None:
        print(f'Cell annotations are provided...')
        spatial_net = pd.DataFrame()
        # Cells with annotations
        for ct in adata.obs[annotation].dropna().unique():
            coor_temp = coor.loc[adata.obs[annotation] == ct, :]
            spatial_net_temp = find_Neighbors(coor_temp, min(num_neighbour, coor_temp.shape[0]))
            spatial_net = pd.concat((spatial_net, spatial_net_temp), axis=0)
            print(f'{ct}: {coor_temp.shape[0]} cells')

        # Cells labeled as nan
        if pd.isnull(adata.obs[annotation]).any():
            cell_nan = adata.obs.index[np.where(pd.isnull(adata.obs[annotation]))[0]]
            print(f'Nan: {len(cell_nan)} cells')

            spatial_net_temp = find_Neighbors(coor, num_neighbour)
            spatial_net_temp = spatial_net_temp.loc[spatial_net_temp.Cell1.isin(cell_nan), :]
            spatial_net = pd.concat((spatial_net, spatial_net_temp), axis=0)
    else:
        print(f'Cell annotations are not provided...')
        spatial_net = find_Neighbors(coor, num_neighbour)

    return spatial_net


def find_Neighbors_Regional(cell):
    cell_use = spatial_net.loc[spatial_net.Cell1 == cell, 'Cell2'].to_list()
    # TODO: To modify the wrongly usage of global variable `coor_latent`
    similarity = cosine_similarity(coor_latent.loc[cell].values.reshape(1, -1),
                                   coor_latent.loc[cell_use].values).tolist()[0]
    if not args.annotation is None:
        annotation = adata.obs[args.annotation]
        df = pd.DataFrame({'Cell2': cell_use, 'Similarity': similarity, 'Annotation': annotation[cell_use]})
        df = df.loc[df.loc[cell_use, 'Annotation'] == df.loc[cell, 'Annotation']]
    else:
        df = pd.DataFrame({'Cell2': cell_use, 'Similarity': similarity})

    df = df.sort_values(by='Similarity', ascending=False)
    cell_select = df.Cell2[0:args.num_neighbour].to_list()

    return cell_select


def _ranks(X, mask=None, mask_rest=None):
    """
    calculate Wilcoxon rank
    """
    CONST_MAX_SIZE = 10000000
    n_genes = X.shape[1]

    if issparse(X):
        merge = lambda tpl: vstack(tpl).toarray()
        adapt = lambda X: X.toarray()
    else:
        merge = np.vstack
        adapt = lambda X: X

    masked = mask is not None and mask_rest is not None
    if masked:
        n_cells = np.count_nonzero(mask) + np.count_nonzero(mask_rest)
        get_chunk = lambda X, left, right: merge(
            (X[mask, left:right], X[mask_rest, left:right])
        )
    else:
        n_cells = X.shape[0]
        get_chunk = lambda X, left, right: adapt(X[:, left:right])

    max_chunk = floor(CONST_MAX_SIZE / n_cells)

    for left in range(0, n_genes, max_chunk):
        right = min(left + max_chunk, n_genes)
        df = pd.DataFrame(data=get_chunk(X, left, right))
        ranks = df.rank()
        yield ranks, left, right


def compute_z_score(adata, cell_select, fold, pst):
    """
    calculate z score of the Wilcoxon test
    """

    mask = adata.obs_names.isin(cell_select)
    mask_rest = ~mask

    rest_avg_exp = adata.X[mask_rest].mean(axis=0).A1
    target_avg_exp = adata.X[mask].mean(axis=0).A1

    mask_gene_expression_fraction = (adata.X[mask].toarray() > 0).sum(axis=0) / mask.sum()
    select = (target_avg_exp > fold * (rest_avg_exp + 1e-8)) & (mask_gene_expression_fraction >= pst)

    n_active = np.count_nonzero(mask)
    m_active = np.count_nonzero(mask_rest)

    X = adata.X[:, select]
    n_genes = X.shape[1]
    scores = np.zeros(n_genes)
    gene_indices = adata[:, select].var.index.tolist()

    T = 1

    if issparse(X):
        X.eliminate_zeros()

    std_dev = np.sqrt(T * n_active * m_active * (n_active + m_active + 1) / 12.0)
    mean_rank = n_active * ((n_active + m_active + 1) / 2.0)

    for ranks, left, right in _ranks(X, mask, mask_rest):
        ranks_sum = np.sum(ranks.iloc[0:n_active, :])
        scores[left:right] = (ranks_sum - mean_rank) / std_dev

    scores[np.isnan(scores)] = 0

    df = pd.DataFrame({'scores': scores}, index=gene_indices)

    return df


def _compute_dge(args):
    """
    calculate z score for one spatial spot or cell
    """
    cell_tg, fold, pst = args
    cell_select = find_Neighbors_Regional(cell_tg)

    scores = np.zeros(adata.shape[1])
    gene_indices = adata.var.index.tolist()
    scores_df = pd.DataFrame({'scores': scores}, index=gene_indices)

    scores_temp = compute_z_score(adata, cell_select, fold, pst)
    scores_temp[scores_temp < 1] = 0

    scores_df.loc[scores_temp.index, 'scores'] = scores_temp.scores

    scores_df.columns = [cell_tg]

    return scores_df


def _compute_regional_ranks(cell_tg):
    """
    compute gmean ranks of a region
    """
    cell_select = find_Neighbors_Regional(cell_tg)

    # Ratio of expression fractions
    frac_focal = np.array((adata[cell_select,].X > 0).sum(axis=0))[0] / (adata[cell_select,].shape[0])
    frac_region = frac_focal / frac_whole
    frac_region[frac_region <= 1] = 0
    frac_region[frac_region > 1] = 1

    # Ratio of expression ranks
    ranks_tg = ranks.loc[cell_select]
    gene_ranks_region = gmean(ranks_tg, axis=0)
    gene_ranks_region[gene_ranks_region <= 1] = 0

    # Simultaneously consider the ratio of expression fractions and ranks
    gene_ranks_region = pd.DataFrame(gene_ranks_region * frac_region)
    gene_ranks_region.columns = [cell_tg]

    return gene_ranks_region


parser = argparse.ArgumentParser()
def add_laten_to_gene_args(parser):
    parser.add_argument('--spe_path', default=None, type=str)
    parser.add_argument('--spe_name', default=None, type=str)
    parser.add_argument('--spe_out', default=None, type=str)
    parser.add_argument('--method', default='rank', type=str)
    parser.add_argument('--latent_representation', default='latent_GVAE', choices=['latent_GVAE', 'latent_PCA'], type=str,help='Name of the latent representation to be used.')
    parser.add_argument('--num_neighbour', default=21, type=int)
    parser.add_argument('--num_neighbour_spatial', default=101, type=int)
    parser.add_argument('--num_processes', default=4, type=int)
    parser.add_argument('--fold', default=1.0, type=float)
    parser.add_argument('--pst', default=0.2, type=float)
    parser.add_argument('--species', default=None, type=str)
    parser.add_argument('--gs_species', default=None, type=str)
    parser.add_argument('--gM_slices', default=None, type=str)

def add_sample_info_args(parser):
    parser.add_argument('--sample_hdf5', default=None, type=str, help='Path to the sample hdf5 file', required=True)
    parser.add_argument('--sample_name', type=str, help='Name of the sample', required=True)
    parser.add_argument('--annotation_layer_name', default=None, type=str, help='Name of the annotation layer',dest='annotation')
    parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts'). This specifies the data layer to be used.",)


if __name__ == '__main__':
    TEST = True
    if TEST:
        name = 'Cortex_151507'
        args = parser.parse_args([
            '--latent_representation', 'latent_GVAE',
            '--spe_path',
            f'/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/annotation/{name}/h5ad',
            '--spe_name', f'{name}_add_latent.h5ad',
            '--num_processes', '4',
            '--type', 'count',
            '--annotation', 'layer_guess',
            '--num_neighbour', '51',
            '--spe_out',
            f'/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/annotation/{name}/gene_markers'
        ])
    else:
        args = parser.parse_args()

    # Load and process the spatial data
    print('------Loading the spatial data...')
    adata = sc.read_h5ad(f'{args.spe_path}/{args.spe_name}')
    num_cpus = min(multiprocessing.cpu_count(), args.num_processes)

    if not args.annotation is None:
        print(f'------Cell annotations are provided as {args.annotation}...')
        adata = adata[~pd.isnull(adata.obs[args.annotation]), :]

    # Homologs transformation
    if not args.species is None:
        print(f'------Transforming the {args.species} to HUMAN_GENE_SYM...')
        homologs = pd.read_csv(args.gs_species, sep='\t')
        homologs.index = homologs[args.species]
        adata = adata[:, adata.var_names.isin(homologs[args.species])]
        print(f'{adata.shape[1]} genes left after homologs transformation.')
        adata.var_names = homologs.loc[adata.var_names, 'HUMAN_GENE_SYM']

    # Process the data
    if args.type == 'count':
        adata.X = adata.layers[args.type]
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        adata.X = adata.layers[args.type]

        # Remove cells that do not express any genes after transformation, and genes that are not expressed in any cells.
    print(f'Number of cells, genes of the input data: {adata.shape[0]},{adata.shape[1]}')
    adata = adata[adata.X.sum(axis=1) > 0, adata.X.sum(axis=0) > 0]
    print(f'Number of cells, genes after transformation: {adata.shape[0]},{adata.shape[1]}')

    # Buid the spatial graph
    spatial_net = _build_spatial_net(adata, args.annotation, args.num_neighbour_spatial)

    # Extract the latent representation
    coor_latent = pd.DataFrame(adata.obsm[args.latent_representation])
    coor_latent.index = adata.obs.index

    # Find marker genes
    mk_score = []
    cell_list = adata.obs.index.tolist()

    if args.method == 'rank':
        prefix = 'rank.feather'

        # Load the geometrical mean across slices
        if not args.gM_slices is None:
            print('Geometrical mean across multiple slices are provided.')
            gM = pd.read_parquet(args.gM_slices)
            # Select the common gene
            common_gene = np.intersect1d(adata.var_names, gM.index)
            gM = gM.loc[common_gene]
            gM = gM['G_Mean'].to_list()
            print('------Ranking the spatial data...')
            adata = adata[:, common_gene]
            ranks = np.apply_along_axis(rankdata, 1, adata.X.toarray())
        else:
            print('------Ranking the spatial data...')
            ranks = np.apply_along_axis(rankdata, 1, adata.X.toarray())
            gM = gmean(ranks, axis=0)

        # Compute the fraction of each gene across cells
        frac_whole = np.array((adata.X > 0).sum(axis=0))[0] / (adata.shape[0])

        # Normalize the geometrical mean
        ranks = ranks / gM
        ranks = pd.DataFrame(ranks, index=adata.obs_names)
        ranks.columns = adata.var.index

        with Pool(num_cpus) as p:
            with tqdm(total=len(cell_list), desc="Finding markers (Rank-based approach) | cells") as progress_bar:
                for mk_cell in p.imap(_compute_regional_ranks, [cell_tg for cell_tg in cell_list]):
                    mk_score.append(np.exp(mk_cell ** 2) - 1)
                    progress_bar.update()
        # use tqdm.progress_map instead of Pool.imap
        # mk_score = list(tqdm.contrib.concurrent.process_map(_compute_regional_ranks, [cell_tg for cell_tg in cell_list],
        #                                                     max_workers=num_cpus,
        #                                                     chunksize=len(cell_list)//num_cpus//10,
        #                                                     ))
    elif args.method == 'wilcox':
        prefix = 'wilcox.feather'

        with Pool(num_cpus) as p:
            with tqdm(total=len(cell_list), desc="Finding markers (Wilcoxon-based approach) | cells") as progress_bar:
                for mk_cell in p.imap(_compute_dge, [(cell_tg, args.fold, args.pst) for cell_tg in cell_list]):
                    mk_score.append(mk_cell ** 2)
                    progress_bar.update()

    # Normalize the marker scores
    mk_score = pd.concat(mk_score, axis=1)
    mk_score.index = adata.var_names
    # mk_score_normalized = mk_score.div(mk_score.sum())*1e+2

    # Remove the mitochondrial genes
    mt_genes = [gene for gene in mk_score.index if gene.startswith('MT-') or gene.startswith('mt-')]
    mk_score_normalized = mk_score.loc[~mk_score.index.isin(mt_genes), :]
    print(mk_score_normalized.shape)

    # Save the marker scores
    print(f'------Saving marker scores ...')
    if not os.path.exists(args.spe_out):
        os.makedirs(args.spe_out, mode=0o777, exist_ok=True)

    data_name = args.spe_name.split('_add_latent.h5ad')[0]
    mk_score_normalized = mk_score_normalized.reset_index()
    mk_score_normalized.rename(columns={mk_score_normalized.columns[0]: 'HUMAN_GENE_SYM'}, inplace=True)
    mk_score_normalized.to_feather(f'{args.spe_out}/{data_name}_{prefix}')









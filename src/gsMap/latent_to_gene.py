import argparse
import logging
import multiprocessing
import pprint
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import gmean
from scipy.stats import rankdata
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

from gsMap.config import add_latent_to_gene_args, LatentToGeneConfig

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


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
    cell_use = spatial_net_dict[cell]
    similarity = cosine_similarity(coor_latent.loc[cell].values.reshape(1, -1),
                                   coor_latent.loc[cell_use].values).reshape(-1)
    if not args.annotation is None:
        annotation = adata.obs[args.annotation]
        df = pd.DataFrame({'Cell2': cell_use, 'Similarity': similarity, 'Annotation': annotation[cell_use]})
        df = df.loc[df.loc[cell_use, 'Annotation'] == df.loc[cell, 'Annotation']]
    else:
        df = pd.DataFrame({'Cell2': cell_use, 'Similarity': similarity})

    df = df.sort_values(by='Similarity', ascending=False)
    cell_select = df.Cell2[0:args.num_neighbour].to_list()

    return cell_select


def _compute_regional_mkscore(cell_tg, ):
    """
    compute gmean ranks of a region
    """
    cell_select = find_Neighbors_Regional(cell_tg)

    # Ratio of expression ranks
    ranks_tg = ranks.loc[cell_select]
    gene_ranks_region = gmean(ranks_tg, axis=0)
    gene_ranks_region[gene_ranks_region <= 1] = 0

    if not args.no_expression_fraction:
        # Ratio of expression fractions
        frac_focal = expressed_mask.loc[cell_select].sum(0) / len(cell_select)
        frac_region = frac_focal / frac_whole
        frac_region[frac_region <= 1] = 0
        frac_region[frac_region > 1] = 1

        # Simultaneously consider the ratio of expression fractions and ranks
        gene_ranks_region = (gene_ranks_region * frac_region).values

    mkscore = np.exp(gene_ranks_region ** 1.5) - 1
    return mkscore.astype(np.float16, copy=False)


def run_latent_to_gene(config: LatentToGeneConfig):
    global adata, coor_latent, spatial_net, ranks, frac_whole, args, spatial_net_dict, expressed_mask
    args = config
    # Load and process the spatial data
    print('------Loading the spatial data...')
    adata = sc.read_h5ad(config.hdf5_with_latent_path)

    # Process the data
    adata.X = adata.layers[config.type]

    print('------Ranking the spatial data...')
    adata.layers['rank'] = rankdata(adata.X.toarray().astype(np.float32), axis=1).astype(np.float32)

    if not config.annotation is None:
        print(f'------Cell annotations are provided as {config.annotation}...')
        adata = adata[~pd.isnull(adata.obs[config.annotation]), :]

    # Homologs transformation
    if not config.homolog_file is None:
        print(f'------Transforming the {config.species} to HUMAN_GENE_SYM...')
        homologs = pd.read_csv(config.homolog_file, sep='\t')
        if homologs.shape[1] == 2:
            raise ValueError(
                "Homologs file must have two columns: one for the species and one for the human gene symbol.")

        homologs.columns = [config.species, 'HUMAN_GENE_SYM']
        homologs.set_index(config.species, inplace=True)
        adata = adata[:, adata.var_names.isin(homologs.index)]
        # Log the number of genes left after homolog transformation
        print(f"{adata.shape[1]} genes retained after homolog transformation.")
        adata.var_names = homologs.loc[adata.var_names, 'HUMAN_GENE_SYM'].values

    # Remove cells that do not express any genes after transformation, and genes that are not expressed in any cells.
    print(f'Number of cells, genes of the input data: {adata.shape[0]},{adata.shape[1]}')
    adata = adata[adata.X.sum(axis=1) > 0, adata.X.sum(axis=0) > 0]
    print(f'Number of cells, genes after transformation: {adata.shape[0]},{adata.shape[1]}')
    # Buid the spatial graph
    spatial_net = _build_spatial_net(adata, config.annotation, config.num_neighbour_spatial)
    spatial_net.set_index('Cell1', inplace=True)
    # convert the spatial graph to a dictionary cell1 to cells in the neighbourhood
    spatial_net_dict = spatial_net.groupby(spatial_net.index).Cell2.apply(list).to_dict()

    # Extract the latent representation
    coor_latent = pd.DataFrame(adata.obsm[config.latent_representation])
    coor_latent.index = adata.obs.index
    # Find marker genes
    cell_list = adata.obs.index.tolist()

    # Load the geometrical mean across slices
    if config.gM_slices is not None:
        print('Geometrical mean across multiple slices is provided.')
        gM = pd.read_parquet(config.gM_slices)
        if config.species is not None:
            homologs = pd.read_csv(config.homolog_file, sep='\t', header=None)
            if homologs.shape[1] < 2:
                raise ValueError(
                    "Homologs file must have at least two columns: one for the species and one for the human gene symbol.")
            homologs.columns = [config.species, 'HUMAN_GENE_SYM']
            homologs.set_index(config.species, inplace=True)
            gM = gM.loc[gM.index.isin(homologs.index)]
            gM.index = homologs.loc[gM.index, 'HUMAN_GENE_SYM'].values
        common_gene = np.intersect1d(adata.var_names, gM.index)
        gM = gM.loc[common_gene]
        gM = gM['G_Mean'].to_numpy()
        adata = adata[:, common_gene]
    else:
        gM = gmean(adata.layers['rank'], axis=0)

    # Compute the fraction of each gene across cells
    expressed_mask = pd.DataFrame((adata.X > 0).toarray(), index=adata.obs.index, columns=adata.var.index)
    # frac_whole = np.array((adata_layer > 0).sum(axis=0))[0] / (adata.shape[0])
    frac_whole = np.array(expressed_mask.sum(axis=0)) / (adata.shape[0])
    # Normalize the geometrical mean
    ranks = adata.layers['rank'] / gM
    ranks = pd.DataFrame(ranks, index=adata.obs_names)
    ranks.columns = adata.var.index
    mk_score = [
        _compute_regional_mkscore(cell_tg)
        for cell_tg in tqdm(cell_list,
                            desc="Finding markers (Rank-based approach) | cells")
    ]
    # Normalize the marker scores
    mk_score = pd.DataFrame(np.vstack(mk_score).T, index=adata.var.index, columns=cell_list)
    # mk_score_normalized = mk_score.div(mk_score.sum())*1e+2
    # Remove the mitochondrial genes
    mt_genes = [gene for gene in mk_score.index if gene.startswith('MT-') or gene.startswith('mt-')]
    mask = ~mk_score.index.isin(set(mt_genes))
    mk_score = mk_score[mask]  # Apply the mask to mk_score
    print(mk_score.shape)

    # save the modified adata
    adata.write(config.hdf5_with_latent_path)

    # Save the marker scores
    print(f'------Saving marker scores ...')
    output_file_path = Path(config.mkscore_feather_path)
    output_file_path.parent.mkdir(parents=True, exist_ok=True, mode=0o755)
    mk_score.reset_index(inplace=True)
    mk_score.rename(columns={mk_score.columns[0]: 'HUMAN_GENE_SYM'}, inplace=True)
    mk_score.to_feather(output_file_path)


# %%
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process latent to gene data.")
    add_latent_to_gene_args(parser)
    TEST = True
    if TEST:
        name = 'Cortex_151507'
        test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_gsMap/data/gsMap_test/Nature_Neuroscience_2021'

        args = parser.parse_args([
            '--input_hdf5_with_latent_path', f'{test_dir}/{name}/hdf5/{name}_add_latent.h5ad',
            '--sample_name', f'{name}',
            '--output_feather_path', f'{test_dir}/{name}/gene_markers/{name}_rank.feather',
            '--method', 'rank',
            '--latent_representation', 'latent_GVAE',
            '--type', 'count',
            '--annotation', 'layer_guess',
            '--num_neighbour', '51',
            # '--no_expression_fraction',

        ])

        # config = LatentToGeneConfig(
        #     **{'annotation': 'SubClass',
        #        'fold': 1.0,
        #        'gM_slices': None,
        #        'homolog_file': '/storage/yangjianLab/songliyang/SpatialData/homologs/macaque_human_homologs.txt',
        #        'input_hdf5_with_latent_path': '/storage/yangjianLab/chenwenhao/projects/202312_gsMap/data/gsMap_test/macaque/T121_macaque1/find_latent_representations/T121_macaque1_add_latent.h5ad',
        #        'latent_representation': 'latent_GVAE',
        #        'method': 'rank',
        #        'num_neighbour': 51,
        #        'num_neighbour_spatial': 201,
        #        'output_feather_path': '/storage/yangjianLab/chenwenhao/projects/202312_gsMap/data/gsMap_test/macaque/T121_macaque1/latent_to_gene/T121_macaque1_gene_marker_score.feather',
        #        'pst': 0.2,
        #        'sample_name': 'T121_macaque1',
        #        'species': 'MACAQUE_GENE_SYM',
        #        'type': 'SCT'}
        # )
    else:
        args = parser.parse_args()
        config = LatentToGeneConfig(**vars(args))
    logger.info(f'Latent to gene for {args.sample_name}...')
    pprint.pprint(config)
    start_time = time.time()
    run_latent_to_gene(config)
    end_time = time.time()
    logger.info(
        f'Latent to gene for {config.sample_name} finished. Time spent: {(end_time - start_time) / 60:.2f} min.')

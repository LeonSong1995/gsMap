import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from plotly.subplots import make_subplots
import plotly.io as pio
from scipy.stats import norm

from gsMap.config import DiagnosisConfig
from gsMap.utils.manhattan_plot import ManhattanPlot
from gsMap.visualize import draw_scatter, load_st_coord, estimate_point_size_for_plot

# %%
warnings.filterwarnings("ignore", category=FutureWarning)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {name} {levelname:6s} {message}', style='{'))
logger.addHandler(handler)


def convert_z_to_p(gwas_data):
    gwas_data['P'] = norm.sf(abs(gwas_data['Z'])) * 2
    return gwas_data


def load_ldsc(ldsc_input_file):
    ldsc = pd.read_csv(ldsc_input_file, compression='gzip')
    ldsc['spot'] = ldsc['spot'].astype(str).replace('\.0', '', regex=True)
    ldsc.set_index('spot', inplace=True)
    ldsc['logp'] = -np.log10(ldsc['p'])
    return ldsc


def generate_manhattan_plot(config: DiagnosisConfig):
    report_save_dir = config.get_report_dir(config.trait_name)

    # Load and process GWAS data
    logger.info('Loading and processing GWAS data...')
    gwas_data = pd.read_csv(config.sumstats_file, compression='gzip', sep='\t')
    gwas_data = convert_z_to_p(gwas_data)

    gene_diagnostic_info = load_gene_diagnostic_info(config)

    # Select top correlated genes if not provided in the configuration
    if not config.selected_genes:
        plot_genes = gene_diagnostic_info.Gene.iloc[:config.top_corr_genes].tolist()
    else:
        plot_genes = config.selected_genes

    # get SNP-gene pairs from generate_ldscore dir
    ldscore_save_dir = Path(config.ldscore_save_dir)
    SNP_gene_pair = pd.concat([
        pd.read_feather(ldscore_save_dir / f'SNP_gene_pair/{config.sample_name}_chr{chrom}.feather')
        for chrom in range(1, 23)
    ])

    # Merge GWAS data with SNP-gene pairs
    gwas_data_with_gene = SNP_gene_pair.merge(gwas_data, on='SNP', how='inner')
    gwas_data_with_gene.rename(columns={'gene_name': 'GENE'}, inplace=True)

    # Merge GWAS data with gene annotations
    gwas_data_with_gene_annotation = gwas_data_with_gene.merge(
        gene_diagnostic_info,
        left_on='GENE',
        right_on='Gene',
        how='left'
    )

    # subsample the GWAS data
    SUBSAMPLE_SNP_NUMBER = 100_000
    wo_gene_annotation_mask = gwas_data_with_gene_annotation['Annotation'].isna()
    logger.info('Only keeping SNPs with gene annotation...')
    gwas_data_with_gene_annotation = gwas_data_with_gene_annotation[~wo_gene_annotation_mask]

    gwas_data_with_gene_annotation_sort = gwas_data_with_gene_annotation.sort_values('P', )

    pass_suggestive_line_mask = gwas_data_with_gene_annotation_sort['P'] < 1e-5
    pass_suggestive_line_number = pass_suggestive_line_mask.sum()

    if pass_suggestive_line_number > SUBSAMPLE_SNP_NUMBER:
        snps2plot = gwas_data_with_gene_annotation_sort[pass_suggestive_line_mask].SNP
        logger.info(f'To reduce the number of SNPs to plot, only {snps2plot.shape[0]} SNPs with P < 1e-5 are plotted.')
    else:
        snps2plot = gwas_data_with_gene_annotation_sort.head(SUBSAMPLE_SNP_NUMBER).SNP
        logger.info(
            f'To reduce the number of SNPs to plot, only {SUBSAMPLE_SNP_NUMBER} SNPs with the smallest P-values are plotted.')

    gwas_data_to_plot = gwas_data_with_gene_annotation[
        gwas_data_with_gene_annotation['SNP'].isin(snps2plot)].reset_index(drop=True)
    gwas_data_to_plot['Annotation_text'] = 'PCC: ' + gwas_data_to_plot['PCC'].round(2).astype(
        str) + '<br>' + 'Annotation: ' + gwas_data_to_plot['Annotation']

    gwas_data_to_plot.value_counts('Annotation', dropna=False)
    # Create the Manhattan plot
    logger.info('Generating Diagnostic Manhattan Plot...')
    fig = ManhattanPlot(
        dataframe=gwas_data_to_plot,
        title='gsMap Diagnosis Manhattan Plot',
        point_size=3,
        highlight_gene_list=plot_genes,
        suggestiveline_value=-np.log10(1e-5),
        annotation='Annotation_text',
    )

    save_manhattan_plot_path = config.get_manhattan_plot_path(config.trait_name)
    fig.write_html(save_manhattan_plot_path)
    logger.info(f'Diagnostic Manhattan Plot saved to {save_manhattan_plot_path}.')


def load_gene_diagnostic_info(config):
    gene_diagnostic_info_save_path = config.get_gene_diagnostic_info_save_path(config.trait_name)
    if gene_diagnostic_info_save_path.exists():
        logger.info(f'Loading gene diagnostic information from {gene_diagnostic_info_save_path}...')
        gene_diagnostic_info = pd.read_csv(gene_diagnostic_info_save_path)
    else:
        logger.info('Gene diagnostic information not found. Calculating gene diagnostic information...')
        gene_diagnostic_info = get_gene_diagnostic_info(config)
    return gene_diagnostic_info


def get_gene_diagnostic_info(config):
    # Load the H5AD file and gene marker scores
    logger.info('Loading ST data...')
    adata = sc.read_h5ad(config.input_hdf5_path, backed='r')
    mk_score = pd.read_feather(config.mkscore_feather_file)
    mk_score.set_index('HUMAN_GENE_SYM', inplace=True)
    mk_score = mk_score.T
    # Load LDSC results
    logger.info('Loading Spatial LDSC results...')
    trait_ldsc_result = load_ldsc(f"{config.input_ldsc_dir}/{config.sample_name}_{config.trait_name}.csv.gz")

    mk_score_aligned = mk_score.loc[trait_ldsc_result.index]

    # drop genes mk score always 0, which will raise error in correlation calculation
    zero_mk_score_mask = mk_score_aligned.sum(axis=0) == 0
    mk_score_aligned = mk_score_aligned.loc[:, ~zero_mk_score_mask]

    logger.info('Calculating correlation between gene marker scores and trait logp-values...')
    corr = mk_score_aligned.corrwith(trait_ldsc_result['logp'].astype(np.float32))
    corr.name = 'PCC'

    # Group the marker score by annotation and determine the highest GSS annotation
    grouped_mk_score = mk_score_aligned.groupby(adata.obs[config.annotation]).median()
    max_annotations = grouped_mk_score.idxmax()
    high_GSS_Gene_annotation_pair = pd.DataFrame({
        'Gene': max_annotations.index,
        'Annotation': max_annotations.values,
        'Median_GSS': grouped_mk_score.max().values
    })
    # Filter out genes with Max Score < 1.0
    high_GSS_Gene_annotation_pair = high_GSS_Gene_annotation_pair[high_GSS_Gene_annotation_pair['Median_GSS'] >= 1.0]
    # Merge with PCC
    high_GSS_Gene_annotation_pair = high_GSS_Gene_annotation_pair.merge(corr, left_on='Gene', right_index=True, )
    # Save the gene diagnostic information
    gene_diagnostic_info_cols = ['Gene', 'Annotation', 'Median_GSS', 'PCC']
    gene_diagnostic_info = high_GSS_Gene_annotation_pair[gene_diagnostic_info_cols].drop_duplicates().dropna(
        subset=['Gene'])
    gene_diagnostic_info.sort_values('PCC', ascending=False, inplace=True)
    gene_diagnostic_info_save_path = config.get_gene_diagnostic_info_save_path(config.trait_name)
    gene_diagnostic_info.to_csv(gene_diagnostic_info_save_path, index=False)
    logger.info(f'Gene diagnostic information saved to {gene_diagnostic_info_save_path}.')

    # add the gene_diagnostic_info info to adata
    adata
    return gene_diagnostic_info


def generate_GSS_distribution(config: DiagnosisConfig):

    # Load the H5AD file and gene marker scores
    logger.info('Loading ST data...')
    adata = sc.read_h5ad(config.hdf5_with_latent_path)

    (pixel_width, pixel_height), point_size = estimate_point_size_for_plot(adata.obsm['spatial'])

    logger.info('Loading marker score...')
    mk_score = pd.read_feather(config.mkscore_feather_path)
    mk_score.set_index('HUMAN_GENE_SYM', inplace=True)
    mk_score = mk_score.T

    if config.selected_genes is not None:
        logger.info('Selected genes provided. Checking if they are present in the dataset...')

        common_genes = set(adata.var_names) & set(config.selected_genes)
        if len(common_genes) == 0:
            raise Warning(
                "Selcted genes don't contain any gene in the dataset: "
            )
        elif len(common_genes) < len(config.selected_genes):
            warnings.warn(
                f"Some genes don't contain in the dataset: {set(config.selected_genes) - common_genes}"
            )

        plot_genes = common_genes
    else:
        logger.info(f'Selected genes not provided. Using the top {config.top_corr_genes} correlated genes...')
        gene_diagnostic_info = load_gene_diagnostic_info(config)
        plot_genes = gene_diagnostic_info.Gene.iloc[:config.top_corr_genes].tolist()

    sub_fig_save_dir = config.get_GSS_plot_dir(config.trait_name)
    logger.info(f'Saving GSS distribution plots to {sub_fig_save_dir}...')

    for selected_gene in plot_genes:
        logger.info(f'Generating GSS & Expression distribution plot for {selected_gene}...')

        # Gene Expression Distribution
        expression = adata[:, selected_gene].X.toarray().flatten()
        expression_series = pd.Series(expression, index=adata.obs.index, name='Expression')
        select_gene_expression_with_space_coord = load_st_coord(adata, expression_series, 'annotation')
        sub_fig_1 = draw_scatter(select_gene_expression_with_space_coord,
                                 title=f'{selected_gene} (Expression)',
                                 annotation='annotation',
                                 color_by='Expression',
                                 point_size=point_size,
                                 width=pixel_width,
                                 height=pixel_height
                                 )
        save_sub_fig_1_path = sub_fig_save_dir / f'{config.sample_name}_{selected_gene}_Expression_Distribution.html'
        sub_fig_1.write_html(str(save_sub_fig_1_path))

        sub_fig_1.update_layout(showlegend=False)
        sub_fig_1.write_image(str(save_sub_fig_1_path).replace('.html', '.png'))

        # GSS distribution
        select_gene_GSS_with_space_coord = load_st_coord(adata, mk_score[selected_gene].rename('GSS'), 'annotation')
        sub_fig_2 = draw_scatter(select_gene_GSS_with_space_coord,
                                 title=f'{selected_gene} (GSS)',
                                 annotation='annotation',
                                 color_by='GSS',
                                 point_size=point_size,
                                 width=pixel_width,
                                 height=pixel_height
                                 )
        save_sub_fig_2_path = sub_fig_save_dir / f'{config.sample_name}_{selected_gene}_GSS_Distribution.html'
        sub_fig_2.write_html(str(save_sub_fig_2_path))

        # disable legend
        sub_fig_2.update_layout(showlegend=False)
        # output png
        sub_fig_2.write_image(str(save_sub_fig_2_path).replace('.html', '.png'))

        # Combine figures into a subplot for PNG output
        combined_fig = make_subplots(rows=1, cols=2,
                                     subplot_titles=(f'{selected_gene} (Expression)', f'{selected_gene} (GSS)'))

        # Add both figures, preserving their individual color scales
        for trace in sub_fig_1.data:
            combined_fig.add_trace(trace, row=1, col=1)

        for trace in sub_fig_2.data:
            combined_fig.add_trace(trace, row=1, col=2)

def run_Diagnosis(config: DiagnosisConfig):
    if config.plot_type == 'manhattan':
        generate_manhattan_plot(config)
    elif config.plot_type == 'GSS':
        generate_GSS_distribution(config)
    elif config.plot_type == 'both':
        generate_manhattan_plot(config)
        generate_GSS_distribution(config)
    else:
        raise ValueError(f'Invalid plot type: {config.plot_type}')


# %%
if __name__ == '__main__':
    config = DiagnosisConfig(
        sample_name='E16.5_E1S1.MOSTA',
        annotation='annotation',
        # trait_name='GIANT_EUR_Height_2022_Nature',
        trait_name='Depression_2023_NatureMed',
        top_corr_genes=5,
        # selected_genes=['COL11A1', 'MECOM'],
        sumstats_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/example_data/GWAS/Depression_2023_NatureMed.sumstats.gz',
        plot_type='manhattan',
    )

    # generate_manhattan_plot(config)

    generate_GSS_distribution(config)

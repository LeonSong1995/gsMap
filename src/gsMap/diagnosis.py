from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
import warnings
import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pandas.api.types import is_numeric_dtype
from plotly.subplots import make_subplots
from scipy.stats import norm
import scanpy as sc
import pyranges as pr
from gsMap import generate_ldscore
import gsMap.config
from gsMap.utils.manhattan_plot import ManhattanPlot
from gsMap.visualize import draw_scatter, load_st_coord

# %%
warnings.filterwarnings("ignore", category=FutureWarning)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {name} {levelname:6s} {message}', style='{'))
logger.addHandler(handler)


# %%
@dataclass
class DiagnosisConfig:
    sample_name: str
    input_hdf5_path: str
    annotation: str
    mkscore_feather_file: str
    bfile_root: str
    input_ldsc_dir: str
    trait_name: str
    gtf_annotation_file: str
    sumstats_file: str
    diagnosis_save_dir: str

    gene_window_size: int = 50000
    top_corr_genes: int = 50
    selected_genes: Optional[List[str]] = None


def _get_hover_text(df, snpname=None, genename=None, annotationname=None):
    hover_text = ''
    if snpname is not None and snpname in df.columns:
        hover_text = 'SNP: ' + df[snpname].astype(str)
    if genename is not None and genename in df.columns:
        hover_text += '<br>GENE: ' + df[genename].astype(str)
    if annotationname is not None and annotationname in df.columns:
        hover_text += '<br>' + df[annotationname].astype(str)
    return hover_text


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
    # Load and process GWAS data
    logger.info('Loading and processing GWAS data...')
    gwas_data = pd.read_csv(config.sumstats_file, compression='gzip', sep='\t')
    gwas_data = convert_z_to_p(gwas_data)

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

    logger.info('Calculating correlation between gene marker scores and trait logp-values...')
    corr = mk_score_aligned.corrwith(trait_ldsc_result['logp'])
    corr.name = 'PCC'

    # Select top correlated genes if not provided in the configuration
    if not config.selected_genes:
        # Calculate the correlation between gene marker scores and trait logp-values
        top_corr_genes = corr.sort_values(ascending=False).head(config.top_corr_genes)
        plot_genes = top_corr_genes.index.tolist()
    else:
        plot_genes = config.selected_genes

    # Load and prepare SNP-gene pairs
    gsMap.config.logger.setLevel(logging.ERROR)  # Suppress logging from generate_ldscore
    s_ldsc_boost = generate_ldscore.S_LDSC_Boost(
        config=gsMap.config.GenerateLDScoreConfig(
            sample_name=config.sample_name,
            chrom='all',
            ldscore_save_dir="",
            mkscore_feather_file=config.mkscore_feather_file,
            bfile_root=config.bfile_root,
            keep_snp_root=None,
            gtf_annotation_file=config.gtf_annotation_file,
            spots_per_chunk=config.gene_window_size,
        )
    )

    # Concatenate SNP position from BIM files
    bim = pd.concat([pd.read_csv(f'{config.bfile_root}.{chrom}.bim', sep='\t', header=None) for chrom in range(1, 23)],
                    axis=0)
    bim.columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]

    # Convert BIM file to PyRanges
    bim_pr = bim.copy()
    bim_pr.columns = ["Chromosome", "SNP", "CM", "Start", "A1", "A2"]
    bim_pr['End'] = bim_pr['Start'].copy()
    bim_pr['Start'] = bim_pr['Start'] - 1  # Due to bim file being 1-based
    bim_pr.Chromosome = 'chr' + bim_pr.Chromosome.astype(str)
    bim_pr = pr.PyRanges(bim_pr)

    # Get SNP-gene pairs using the loaded LDSC Boost instance
    logger.info('Creating SNP-gene pairs...')
    SNP_gene_pair = s_ldsc_boost.get_SNP_gene_pair_from_gtf(bim, bim_pr)

    # Merge GWAS data with SNP-gene pairs
    gwas_data_with_gene = gwas_data.merge(SNP_gene_pair, on='SNP', how='inner')
    gwas_data_with_gene.rename(columns={'gene_name': 'GENE'}, inplace=True)

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

    # Merge GWAS data with gene annotations
    gwas_data_with_gene_annotation = gwas_data_with_gene.merge(
        high_GSS_Gene_annotation_pair,
        left_on='GENE',
        right_on='Gene',
        how='left'
    )

    # Save the gene diagnostic information
    save_dir = Path(config.diagnosis_save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    gene_diagnostic_info_cols = ['Gene', 'Annotation', 'Median_GSS', 'PCC']
    gene_diagnostic_info = gwas_data_with_gene_annotation[gene_diagnostic_info_cols].drop_duplicates().dropna(
        subset=['Gene'])
    gene_diagnostic_info.sort_values('PCC', ascending=False, inplace=True)
    gene_diagnostic_info_save_path = save_dir / f'{config.sample_name}_{config.trait_name}_Gene_Diagnostic_Info.csv'
    gene_diagnostic_info.to_csv(gene_diagnostic_info_save_path, index=False)
    logger.info(f'Gene diagnostic information saved to {gene_diagnostic_info_save_path}.')

    gwas_data_with_gene_annotation['Annotation_text'] = 'PCC: ' + gwas_data_with_gene_annotation['PCC'].round(2).astype(
        str) + '<br>' + 'Annotation: ' + gwas_data_with_gene_annotation['Annotation']

    # Create the Manhattan plot
    logger.info('Generating Diagnostic Manhattan Plot...')
    fig = ManhattanPlot(
        dataframe=gwas_data_with_gene_annotation,
        title='gsMap Diagnosis Manhattan Plot',
        point_size=1,
        highlight_gene_list=plot_genes,
        suggestiveline_value=-np.log10(1e-5),
        annotation='Annotation_text',
    )

    fig.show()
    save_manhattan_plot_path = save_dir / f'{config.sample_name}_{config.trait_name}_Diagnostic_Manhattan_Plot.html'
    fig.write_html(save_manhattan_plot_path)
    logger.info(f'Diagnostic Manhattan Plot saved to {save_manhattan_plot_path}. Open in browser to view.')


def generate_GSS_distribution(config: DiagnosisConfig):
    save_dir = Path(config.diagnosis_save_dir)

    # Load the H5AD file and gene marker scores
    logger.info('Loading ST data...')
    adata = sc.read_h5ad(config.input_hdf5_path, )

    logger.info('Loading marker score...')
    mk_score = pd.read_feather(config.mkscore_feather_file)
    mk_score.set_index('HUMAN_GENE_SYM', inplace=True)
    mk_score = mk_score.T

    assert config.selected_genes is not None, 'Please provide selected genes for GSS distribution plot'

    common_genes = set(adata.var_names) & set(config.selected_genes)
    if len(common_genes) == 0:
        raise Warning(
            "Selcted genes don't contain any gene in the dataset: "
        )
    elif len(common_genes) < len(config.selected_genes):
        warnings.warn(
            f"Some genes don't contain in the dataset: {set(config.selected_genes) - common_genes}"
        )

    plot_genes = config.selected_genes

    sub_fig_save_dir = save_dir / 'sub_figures'
    sub_fig_save_dir.mkdir(parents=True, exist_ok=True)

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
                                 point_size=2
                                 )
        save_sub_fig_1_path = sub_fig_save_dir / f'{config.sample_name}_{selected_gene}_Expression_Distribution.html'
        sub_fig_1.write_html(str(save_sub_fig_1_path))

        # GSS distribution
        select_gene_GSS_with_space_coord = load_st_coord(adata, mk_score[selected_gene].rename('GSS'), 'annotation')
        sub_fig_2 = draw_scatter(select_gene_GSS_with_space_coord,
                                 title=f'{selected_gene} (GSS)',
                                 annotation='annotation',
                                 color_by='GSS',
                                 point_size=2
                                 )
        save_sub_fig_2_path = sub_fig_save_dir / f'{config.sample_name}_{selected_gene}_GSS_Distribution.html'
        sub_fig_2.write_html(str(save_sub_fig_2_path))



# %%
if __name__ == '__main__':
    config = DiagnosisConfig(
        sample_name='E16.5_E1S1.MOSTA',
        input_hdf5_path='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/find_latent_representations/E16.5_E1S1.MOSTA_add_latent.h5ad',
        annotation='annotation',
        mkscore_feather_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/latent_to_gene/E16.5_E1S1.MOSTA_gene_marker_score.feather',
        bfile_root='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/GPS_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC',
        input_ldsc_dir='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/spatial_ldsc',
        trait_name='GIANT_EUR_Height_2022_Nature',
        gtf_annotation_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/GPS_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf',
        gene_window_size=50000,
        top_corr_genes=10,
        selected_genes=['COL11A1', 'MECOM'],
        sumstats_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/example_data/GWAS/GIANT_EUR_Height_2022_Nature.sumstats.gz',
        diagnosis_save_dir='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/diagnosis'
    )

    # generate_manhattan_plot(config)

    generate_GSS_distribution(config)

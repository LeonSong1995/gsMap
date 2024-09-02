from dataclasses import dataclass
from typing import Optional, List
import warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pandas.api.types import is_numeric_dtype
from scipy.stats import norm
import scanpy as sc
from gsMap.utils.manhattan_plot import ManhattanPlot

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
    gene_window_size: int = 50000
    top_corr_genes: int = 10
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
    print('Loading and processing GWAS data...')
    gwas_data = pd.read_csv(config.sumstats_file, compression='gzip', sep='\t')
    gwas_data = convert_z_to_p(gwas_data)

    # Load the H5AD file and gene marker scores
    adata = sc.read_h5ad(config.input_hdf5_path, backed='r')
    mk_score = pd.read_feather(config.mkscore_feather_file)
    mk_score.set_index('HUMAN_GENE_SYM', inplace=True)
    mk_score = mk_score.T

    # Load LDSC results
    trait_ldsc_result = load_ldsc(f"{config.input_ldsc_dir}/{config.sample_name}_{config.trait_name}.csv.gz")
    mk_score_aligned = mk_score.loc[trait_ldsc_result.index]

    # Select top correlated genes if not provided in the configuration
    if not config.selected_genes:
        # Calculate the correlation between gene marker scores and trait logp-values
        print('Calculating correlation between gene marker scores and trait logp-values...')
        corr = mk_score_aligned.corrwith(trait_ldsc_result['logp'])
        top_corr_genes = corr.sort_values(ascending=False).head(config.top_corr_genes)
        plot_genes = top_corr_genes.index.tolist()
    else:
        plot_genes = config.selected_genes

    # Group the marker score by annotation and determine the highest GSS annotation
    grouped_mk_score = mk_score_aligned.groupby(adata.obs[config.annotation]).median()
    max_annotations = grouped_mk_score.idxmax()

    high_GSS_Gene_annotation_pair = pd.DataFrame({
        'Gene': max_annotations.index,
        'Max_Score_Annotation': max_annotations.values,
        'Max_Score': grouped_mk_score.max().values
    })

    # Filter out genes with Max Score < 1.0
    high_GSS_Gene_annotation_pair = high_GSS_Gene_annotation_pair[high_GSS_Gene_annotation_pair['Max_Score'] >= 1.0]

    # Merge GWAS data with gene annotations
    gwas_data_with_gene_annotation = gwas_data.merge(high_GSS_Gene_annotation_pair, left_on='GENE', right_on='Gene',
                                                     how='left')
    gwas_data_with_gene_annotation['High GSS Annotation'] = gwas_data_with_gene_annotation['Max_Score_Annotation']

    # Create the Manhattan plot
    fig = ManhattanPlot(
        dataframe=gwas_data_with_gene_annotation,
        title='gsMap Diagnosis Manhattan Plot',
        point_size=1,
        highlight_gene_list=plot_genes,
        suggestiveline_value=-np.log10(1e-5),
        annotation='High GSS Annotation',
    )

    fig.show()
    fig.write_html('gsMap_Diagnosis_Manhattan_Plot.html')

#%%
if __name__ == '__main__':
    config = DiagnosisConfig(
        sample_name='E16.5_E1S1.MOSTA',
        input_hdf5_path='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/example_data/ST/E16.5_E1S1.MOSTA.h5ad',
        annotation='annotation',
        mkscore_feather_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/latent_to_gene/E16.5_E1S1.MOSTA_gene_marker_score.feather',
        bfile_root='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/GPS_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC',
        input_ldsc_dir='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/spatial_ldsc',
        trait_name='GIANT_EUR_Height_2022_Nature',
        gtf_annotation_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/GPS_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf',
        gene_window_size=50000,
        top_corr_genes=10,
        selected_genes=None,
        sumstats_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/example_data/GWAS/GIANT_EUR_Height_2022_Nature.sumstats.gz'
    )

    generate_manhattan_plot(config)
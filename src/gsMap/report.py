import logging
import os
import shutil
from pathlib import Path

from jinja2 import Environment, FileSystemLoader
import pandas as pd

import gsMap
from gsMap.cauchy_combination_test import run_Cauchy_combination
from gsMap.config import CauchyCombinationConfig, ReportConfig
from gsMap.diagnosis import run_Diagnosis

# logging.basicConfig(level=logging.INFO,
#                     format='%(asctime)s - %(levelname)s - %(message)s',
#                     handlers=[
#                         logging.FileHandler("pipeline.log"),
#                         logging.StreamHandler()
#                     ])
logger = logging.getLogger(__name__)

# Load the Jinja2 environment
try:
    from importlib.resources import files

    template_dir = files('gsMap').joinpath('templates')
except (ImportError, FileNotFoundError):
    # Fallback to a relative path if running in development mode
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')

# Set up Jinja2 environment
env = Environment(loader=FileSystemLoader(template_dir))

# Load the template
template = env.get_template('report_template.html')

def copy_files_to_report_dir(result_dir, report_dir, files_to_copy):
    """Copy specified files (HTML or PNG) to the report directory."""
    os.makedirs(report_dir, exist_ok=True)
    for file in files_to_copy:
        shutil.copy2(file, os.path.join(report_dir, os.path.basename(file)))


def load_cauchy_table(csv_file):
    """Load the Cauchy combination table from a compressed CSV file using Pandas."""
    df = pd.read_csv(csv_file, compression='gzip')
    table_data = df[['annotation', 'p_cauchy', 'p_median']].to_dict(orient='records')
    return table_data


def load_gene_diagnostic_info(csv_file):
    """Load the Gene Diagnostic Info CSV file and return the top 50 rows."""
    df = pd.read_csv(csv_file)
    top_50 = df.head(50).to_dict(orient='records')
    return top_50


def embed_html_content(file_path):
    """Read the content of an HTML file and return it as a string."""
    with open(file_path, 'r') as f:
        return f.read()

def check_and_run_cauchy_combination(config):
    cauchy_result_file = config.get_cauchy_result_file(config.trait_name)
    if cauchy_result_file.exists():
        logger.info(
            f"Cauchy combination already done for trait {config.trait_name}. Results saved at {cauchy_result_file}. Skipping...")
    else:
        logger.info(f"Running Cauchy combination for trait {config.trait_name}...")
        cauchy_config = CauchyCombinationConfig(
            workdir=config.workdir,
            sample_name=config.sample_name,
            annotation=config.annotation,
            trait_name=config.trait_name,
        )
        run_Cauchy_combination(cauchy_config)

    df = pd.read_csv(cauchy_result_file, compression='gzip')
    table_data = df[['annotation', 'p_cauchy', 'p_median']].to_dict(orient='records')

    return table_data

def run_Report(config: ReportConfig, run_parameters=None):

    logger.info('Running gsMap Diagnosis Module')
    run_Diagnosis(config)
    logger.info('gsMap Diagnosis running successfully')

    report_dir = config.get_report_dir(config.trait_name)
    gene_diagnostic_info_file = config.get_gene_diagnostic_info_save_path(config.trait_name)
    gene_diagnostic_info = load_gene_diagnostic_info(gene_diagnostic_info_file)

    # Load data (Cauchy table and gene diagnostic info)
    cauchy_table = check_and_run_cauchy_combination(config)

    # Paths to PNGs for gene expression and GSS distribution
    gss_distribution_dir = config.get_GSS_plot_dir(config.trait_name)

    gene_plots = []
    plot_select_gene_list = config.get_GSS_plot_select_gene_file(config.trait_name).read_text().splitlines()
    for gene_name in plot_select_gene_list:
        expression_png = gss_distribution_dir / f"{config.sample_name}_{gene_name}_Expression_Distribution.png"
        gss_png = gss_distribution_dir / f"{config.sample_name}_{gene_name}_GSS_Distribution.png"
        # check if expression and GSS plots exist
        if not os.path.exists(expression_png) or not os.path.exists(gss_png):
            print(f"Skipping gene {gene_name} as expression or GSS plot is missing.")
            continue
        gene_plots.append({
            'name': gene_name,
            'expression_plot': expression_png.relative_to(report_dir),  # Path for gene expression plot
            'gss_plot': gss_png.relative_to(report_dir)  # Path for GSS distribution plot
        })

    # # Copy PNG files to the report directory
    # copy_files_to_report_dir(result_dir, report_dir, [gene['expression_plot'] for gene in gene_plots] + [gene['gss_plot'] for gene in gene_plots])

    # Update paths to point to copied images inside the report folder
    # for gene in gene_plots:
    #     gene['expression_plot'] = os.path.join(os.path.basename(gene['expression_plot']))
    #     gene['gss_plot'] = os.path.join(os.path.basename(gene['gss_plot']))

    # Sample data for other report components
    title = f"{config.sample_name} Genetic Spatial Mapping Report"

    genetic_mapping_plot = embed_html_content(config.get_gsMap_html_plot_save_path(config.trait_name))
    manhattan_plot = embed_html_content(config.get_manhattan_html_plot_path(config.trait_name))

    gsmap_version = gsMap.__version__
    # Render the template with dynamic content, including the run parameters
    output_html = template.render(
        title=title,
        genetic_mapping_plot=genetic_mapping_plot,  # Inlined genetic mapping plot
        manhattan_plot=manhattan_plot,  # Inlined Manhattan plot
        cauchy_table=cauchy_table,
        gene_plots=gene_plots,  # List of PNG paths for gene plots
        gsmap_version=gsmap_version,
        parameters=run_parameters,  # Pass the run parameters to the template
        gene_diagnostic_info=gene_diagnostic_info  # Include top 50 gene diagnostic info rows
    )

    # Save the generated HTML report in the 'report' directory
    report_file = config.get_gsMap_report_file(config.trait_name)
    with open(report_file, "w") as f:
        f.write(output_html)

    logger.info(f"Report generated successfully! Saved at {report_file}.")
    logger.info(f"Copy the report directory to your local PC and open the HTML report file in a web browser to view the report.")


#
#
# if __name__ == '__main__':
#
#     # Example usage
#     result_dir = "/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/"
#     sample_name = "E16.5_E1S1.MOSTA"
#     trait_name = "Depression_2023_NatureMed"
#
#
#     run_parameter_dict = {
#         "sample_name": 'config.SAMPLE_NAME',
#         "trait_name": trait_name,
#         "ldscore_dir": 'ldscore_config.ldscore_save_dir',
#         "w_file": 'config.W_FILE',
#         "annotation": 'config.ANNOTATION',
#         "gtf_annotation_file": 'config.GTFFILE',
#         "bfile_root": 'config.BFILE_ROOT',
#         "keep_snp_root": 'config.KEEP_SNP_ROOT',
#         "mkscore_feather_file": 'latent_to_gene_config.output_feather_path',
#         "spatial_ldsc_save_dir": 'spatial_ldsc_config.ldsc_save_dir',
#         "sumstats_file": 'sumstats_config[trait_name]',
#     }
#     run_Report(result_dir, sample_name, trait_name, run_parameter_dict)
#
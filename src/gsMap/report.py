import csv
import json
import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader

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


def generate_report(result_dir, sample_name, trait_name):
    """
    Generate the report by dynamically loading data based on the result directory,
    sample name, and trait name.
    """

    # Paths to different directories and files based on the provided result directory and sample/trait name
    cauchy_file = os.path.join(result_dir, 'cauchy_combination', f"{sample_name}_{trait_name}.Cauchy.csv.gz")
    diagnosis_dir = os.path.join(result_dir, 'diagnosis')
    gene_diagnostic_info_file = os.path.join(diagnosis_dir, f"{sample_name}_{trait_name}_Gene_Diagnostic_Info.csv")

    # Load data
    cauchy_table = load_cauchy_table(cauchy_file)

    # Load gene expression and GSS plots based on files in the diagnosis directory
    gss_distribution_dir = os.path.join(diagnosis_dir, 'GSS_distribution')
    gene_plots = [
        {'name': gene_name,
         'expression_plot': os.path.join(gss_distribution_dir,
                                         f"{sample_name}_{gene_name}_Expression_Distribution.html"),
         'gss_plot': os.path.join(gss_distribution_dir, f"{sample_name}_{gene_name}_GSS_Distribution.html")}
        for gene_name in ['CELF4', 'COL11A1', 'INA', 'MAP2', 'MAPT', 'MECOM', 'RAB3C']  # Add more gene names as needed
    ]

    # Load the top 50 rows of gene diagnostic info
    gene_diagnostic_info = load_gene_diagnostic_info(gene_diagnostic_info_file)

    # Sample data for other report components
    title = f"{sample_name} Genetic Spatial Mapping Report"
    genetic_mapping_plot = "trait.html"  # Update this with the actual path if needed
    manhattan_plot = "manhattan_plot.html"  # Update this with the actual path if needed
    gsmap_version = "1.0.0"
    parameters = "param1=value1, param2=value2"

    # Render the template with dynamic content
    output_html = template.render(
        title=title,
        genetic_mapping_plot=genetic_mapping_plot,
        manhattan_plot=manhattan_plot,
        cauchy_table=cauchy_table,
        gene_plots=gene_plots,
        gsmap_version=gsmap_version,
        parameters=parameters,
        gene_diagnostic_info=gene_diagnostic_info,  # Include top 50 gene diagnostic info rows
        gene_plots_json=json.dumps(gene_plots)  # Pass gene plots as JSON for JavaScript
    )

    # Save the generated HTML report
    report_file = os.path.join(result_dir, f"{sample_name}_{trait_name}_gsMap_report.html")
    with open(report_file, "w") as f:
        f.write(output_html)

    print(f"Report generated successfully! Saved at {report_file}")


# Example usage
result_dir = "/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/"
sample_name = "E16.5_E1S1.MOSTA"
trait_name = "Depression_2023_NatureMed"

generate_report(result_dir, sample_name, trait_name)

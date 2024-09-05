import os
import shutil

from jinja2 import Environment, FileSystemLoader
import pandas as pd

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


def copy_png_files(result_dir, report_dir, gene_plots):
    """Copy PNG files to the report directory."""
    os.makedirs(report_dir, exist_ok=True)
    for gene in gene_plots:
        expression_png = gene['expression_plot']
        gss_png = gene['gss_plot']

        # Copy expression plot
        shutil.copy2(expression_png, os.path.join(report_dir, os.path.basename(expression_png)))

        # Copy GSS plot
        shutil.copy2(gss_png, os.path.join(report_dir, os.path.basename(gss_png)))


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
    sample name, and trait name. Use PNG images for gene plots and copy them to the report folder.
    """

    # Paths to different directories and files based on the provided result directory and sample/trait name
    cauchy_file = os.path.join(result_dir, 'cauchy_combination', f"{sample_name}_{trait_name}.Cauchy.csv.gz")
    diagnosis_dir = os.path.join(result_dir, 'diagnosis')
    gene_diagnostic_info_file = os.path.join(diagnosis_dir, f"{sample_name}_{trait_name}_Gene_Diagnostic_Info.csv")
    report_dir = os.path.join(result_dir, 'report')

    # Load data (Cauchy table and gene diagnostic info)
    cauchy_table = load_cauchy_table(cauchy_file)
    gene_diagnostic_info = load_gene_diagnostic_info(gene_diagnostic_info_file)

    # Paths to PNGs for gene expression and GSS distribution
    gss_distribution_dir = os.path.join(diagnosis_dir, 'GSS_distribution')
    gene_plots = []
    for gene_name in ['CELF4', 'INA', 'MAP2', 'MAPT', 'MECOM', 'RAB3C']:  # Add more gene names as needed
        expression_png = os.path.join(gss_distribution_dir, f"{sample_name}_{gene_name}_Expression_Distribution.png")
        gss_png = os.path.join(gss_distribution_dir, f"{sample_name}_{gene_name}_GSS_Distribution.png")
        # check if expression and GSS plots exist
        if not os.path.exists(expression_png) or not os.path.exists(gss_png):
            print(f"Skipping gene {gene_name} as expression or GSS plot is missing.")
            continue
        gene_plots.append({
            'name': gene_name,
            'expression_plot': expression_png,  # Path for expression plot
            'gss_plot': gss_png  # Path for GSS plot
        })

    # Copy PNG files to the report directory
    copy_png_files(result_dir, report_dir, gene_plots)

    # Update paths to point to copied images inside the report folder
    for gene in gene_plots:
        gene['expression_plot'] = os.path.join('report', os.path.basename(gene['expression_plot']))
        gene['gss_plot'] = os.path.join('report', os.path.basename(gene['gss_plot']))

    # Sample data for other report components
    title = f"{sample_name} Genetic Spatial Mapping Report"

    # Embed the genetic mapping plot and Manhattan plot as HTML
    genetic_mapping_plot = os.path.join(result_dir, 'visualize', f'{sample_name}_{trait_name}.html')  # Embed HTML
    manhattan_plot = os.path.join(diagnosis_dir, f'{sample_name}_{trait_name}_Diagnostic_Manhattan_Plot.html')  # Embed HTML

    gsmap_version = "1.0.0"
    parameters = "param1=value1, param2=value2"

    # Render the template with dynamic content
    output_html = template.render(
        title=title,
        genetic_mapping_plot=genetic_mapping_plot,  # Inlined genetic mapping plot (HTML)
        manhattan_plot=manhattan_plot,  # Inlined Manhattan plot (HTML)
        cauchy_table=cauchy_table,
        gene_plots=gene_plots,  # List of PNG paths for gene plots
        gsmap_version=gsmap_version,
        parameters=parameters,
        gene_diagnostic_info=gene_diagnostic_info  # Include top 50 gene diagnostic info rows
    )

    # Save the generated HTML report in the 'report' directory
    report_file = os.path.join(report_dir, f"{sample_name}_{trait_name}_gsMap_report.html")
    with open(report_file, "w") as f:
        f.write(output_html)

    print(f"Report generated successfully! Saved at {report_file}")


# Example usage
result_dir = "/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/"
sample_name = "E16.5_E1S1.MOSTA"
trait_name = "Depression_2023_NatureMed"

generate_report(result_dir, sample_name, trait_name)
import csv
import json
from jinja2 import Environment, FileSystemLoader

# Load the Jinja2 environment
env = Environment(loader=FileSystemLoader('path/to/your/templates'))  # Update the path

# Define template
template = env.get_template('report_template.html')

# Load the Cauchy combination result from CSV file
def load_cauchy_table(csv_file):
    table_data = []
    with open(csv_file, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            table_data.append({
                'annotation': row['annotation'],
                'p_cauchy': row['p_cauchy'],
                'p_median': row['p_median']
            })
    return table_data

# Sample data
title = "gsMap Genetic Spatial Mapping Report"
genetic_mapping_plot = "trait.html"
manhattan_plot = "manhattan_plot.html"
gsmap_version = "1.0.0"
parameters = "param1=value1, param2=value2"

# Load the Cauchy combination table from a CSV file
cauchy_table = load_cauchy_table("cauchy_combination_results.csv")

# Sample gene expression and GSS plots
gene_plots = [
    {'name': 'geneA', 'expression_plot': 'geneA.expression.html', 'gss_plot': 'geneA.GSS.html'},
    {'name': 'geneB', 'expression_plot': 'geneB.expression.html', 'gss_plot': 'geneB.GSS.html'},
    {'name': 'geneC', 'expression_plot': 'geneC.expression.html', 'gss_plot': 'geneC.GSS.html'}
]

# Render the template with dynamic content
output_html = template.render(
    title=title,
    genetic_mapping_plot=genetic_mapping_plot,
    manhattan_plot=manhattan_plot,
    cauchy_table=cauchy_table,
    gene_plots=gene_plots,
    gsmap_version=gsmap_version,
    parameters=parameters,
    gene_plots_json=json.dumps(gene_plots)  # Pass gene plots as JSON for JavaScript
)

# Write the rendered HTML to a file
with open("gsMap_report.html", "w") as f:
    f.write(output_html)

print("Report generated successfully!")

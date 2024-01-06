# GPS: Genetics-informed Pathogenic Spatial Mapping Program

## Overview

GPS (Genetics-informed Pathogenic Spatial Mapping) is Python command-line tool designed for ....

## Features

....

## Installation

install use pip:

```bash
pip install GPS-mapping
```

install from source:

```bash
git clone
cd GPS-mapping
pip install -e .
```

## Usage

To use GPS, navigate to the command line and enter `GPS` followed by the subcommand that corresponds to the desired operation. Each subcommand may require specific arguments to run.

### Basic Command Structure

```bash
GPS subcommand [arguments...]
```

- `subcommand`: The specific operation you wish to perform.
- `arguments`: The arguments and options required for the subcommand.

### Available Subcommands

(Provide a list and brief description of each available subcommand. For example:)

- `run_find_latent_representations`: Finds latent representations using a GNN-VAE model.
- `run_latent_to_gene`: Maps latent representations to gene markers.
- `run_generate_ldscore`: Generates LD scores for genomic spots.
- `run_spatial_ldsc`: Conducts spatial LDSC analysis.
- `run_cauchy_combination`: Performs Cauchy combination tests for annotations.
- `run_all_mode`: Executes a comprehensive pipeline covering all steps.

### Examples

To run a specific functionality, you need to provide the appropriate subcommand and arguments. For example:

```bash
GPS run_find_latent_representations --input_hdf5_path <path> --output_hdf5_path <path> --sample_name <name>
```

This command initiates the process of finding latent representations based on the given HDF5 input and output paths and sample name.

## Contributing

...

## License

....

---

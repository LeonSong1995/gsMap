# gsMap (genetically informed spatial mapping of cells for complex traits)
[![stars-badge](https://img.shields.io/github/stars/LeonSong1995/MeDuSA?logo=GitHub&color=yellow)](https://github.com/LeonSong1995/gsMap/stargazers)
[![pypi-badge](https://img.shields.io/pypi/v/scglue)](https://pypi.org/project/scglue)
[![docs-badge](https://readthedocs.org/projects/gps-mapping/badge/?version=latest](https://gps-mapping.readthedocs.io/en/latest/)
[![license-badge](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## Features

....

## Installation

install use pip:

```bash
pip install gsMap
```

install from source:

```bash
git clone
cd gsMap
pip install -e .
```

## Usage

To use gsMap, navigate to the command line and enter `gsMap` followed by the subcommand that corresponds to the desired operation. Each subcommand may require specific arguments to run.

### Basic Command Structure

```bash
gsmap subcommand [arguments...]
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
### Running Requirement


```bash
gsmap run_find_latent_representations --input_hdf5_path <path> --output_hdf5_path <path> --sample_name <name>
```

This command initiates the process of finding latent representations based on the given HDF5 input and output paths and sample name.

## Contributing

...

## License

....

---

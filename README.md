[![GitHub Stars](https://img.shields.io/github/stars/JianYang-Lab/gsMap?logo=GitHub&color=yellow)](https://github.com/LeonSong1995/gsMap/stargazers)
[![PyPI Version](https://img.shields.io/pypi/v/gsMap)](https://pypi.org/project/gsMap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# gsMap

## Introduction

`gsMap` (genetically informed spatial mapping of cells for complex traits)
integrates spatial transcriptomics (ST) data with genome-wide association study (GWAS)
summary statistics to map cells to human complex traits, including diseases,
in a spatially resolved manner.

## Key Features

- **Spatially-aware High-Resolution Trait Mapping**
- **Spatial Region Identification**
- **Putative Causal Genes Identification**

![Model Architecture](schematic.png)

## Installation

Install using pip:

```bash
conda create -n gsMap python>=3.10
conda activate gsMap
pip install gsMap
```

Install from source:

```bash
git clone https://github.com/LeonSong1995/gsMap.git
cd gsMap
pip install -e .
```

Verify the installation by running the following command:

```bash
gsmap --help
```

## Usage

Please check out the documentation and tutorials at [gsMap Documentation](https://yanglab.westlake.edu.cn/gsmap/document/software).

## Online Visualization

To visualize the traits-cell association spatial maps,
please refer to [gsMap Visualization](https://yanglab.westlake.edu.cn/gsmap/visualize).

## Citation

Song, L., Chen, W., Hou, J., Guo, M. & Yang, J.
[Spatially reso lved mapping of cells associated with human complex traits.](https://www.medrxiv.org/content/10.1101/2024.10.31.24316538v1)
medRxiv, 2024.2010.2031.24316538 (2024) (Nature in press).

Please cite the paper and give us a STAR if you find gsMap useful for your research.

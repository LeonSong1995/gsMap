[![GitHub Stars](https://img.shields.io/github/stars/LeonSong1995/gsMap?logo=GitHub&color=yellow)](https://github.com/LeonSong1995/gsMap/stargazers)
[![PyPI Version](https://img.shields.io/pypi/v/gsMap)](https://pypi.org/project/gsMap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# gsMap

## Introduction

`gsMap` (genetically informed spatial mapping of cells for complex traits) integrates spatial transcriptomics (ST) data with genome-wide association study (GWAS) summary statistics to map cells to human complex traits, including diseases, in a spatially resolved manner.


## Key Features
- **Spatially-aware High-Resolution Trait Mapping**: Maps trait-associated cells at single-cell resolution, offering insights into their spatial distributions.
- **Spatial Region Identification**: Aggregates trait-cell association p-values into trait-tissue region association p-values, prioritizing tissue regions relevant to traits of interest.
- **Putative Causal Genes Identification**: Prioritizes putative causal genes by associating gene expression levels with cell-trait relevance.

![Model Architecture](schematic.png)

## Installation

Install using pip:

```bash
conda create -n gsMap python>3.8
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

To visualize the traits-embryo or traits-brain association spatial maps, please refer to [gsMap Visualization](https://yanglab.westlake.edu.cn/gsmap/visualize).

## How to Cite

If you use `gsMap` in your studies, please cite:

    Liyang Song, Wenhao Chen, Junren Hou, Minmin Guo, Jian Yang (2024) Spatially resolved mapping of cells associated with human complex traits. (Under review)

[![GitHub Stars](https://img.shields.io/github/stars/LeonSong1995/gsMap?logo=GitHub&color=yellow)](https://github.com/LeonSong1995/gsMap/stargazers)
[![PyPI Version](https://img.shields.io/pypi/v/gsMap)](https://pypi.org/project/gsMap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Introduction

`gsMap` (genetically informed spatial mapping of cells for complex traits) integrates spatial transcriptomics (ST) data with genome-wide association study (GWAS) summary statistics to map cells to human complex traits, including diseases, in a spatially resolved manner.


## Key Features
- **Spatially-aware High-Resolution Trait Mapping**
- **Spatial Region Identification**
- **Putative Causal Genes Identification**

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

To visualize the traits-cell association spatial maps, please refer to [gsMap Visualization](https://yanglab.westlake.edu.cn/gsmap/visualize).

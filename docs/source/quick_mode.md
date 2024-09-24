# Mouse Embryo (Quick Mode)

The `Quick Mode` option allows you to execute the entire `gsMap` pipeline in a simplified and efficient way, utilizing pre-calculated resources to minimize running time and configuration complexity. This mode is ideal for users who prefer a streamlined approach. For a more customizable experience with adjustable parameters, please refer to the {doc}`Step by Step <mouse_example>` guide.

## Preparation

Make sure you have {doc}`installed <install>` the `gsMap` package before proceeding.

### 1. Download Dependencies

The `gsMap` package in quick mode requires the following resources:
- **Gene transfer format (GTF) file**, for gene coordinates on the genome.
- **LD reference panel**, in quick mode, we provide a pre-built LD score snp-by-gene matrix based on 1000G_EUR_Phase3.
- **SNP weight file**, to adjust correlations between SNP-trait association statistics.
- **Homologous gene transformations file** (optional), to map genes between species.

To download all the required files:
```bash
wget https://yanglab.westlake.edu.cn/data/gsMap/gsMap_resource.tar.gz
tar -xvzf gsMap_resource.tar.gz
```

Directory structure:
```bash
tree -L 2

gsMap_resource
    ├── genome_annotation
    │   ├── enhancer
    │   └── gtf
    ├── homologs
    │   ├── macaque_human_homologs.txt
    │   └── mouse_human_homologs.txt
    ├── LD_Reference_Panel
    │   └── 1000G_EUR_Phase3_plink
    ├── LDSC_resource
    │   ├── hapmap3_snps
    │   └── weights_hm3_no_hla
    └── quick_mode
        ├── baseline
        ├── SNP_gene_pair
        └── snp_gene_weight_matrix.h5ad
```

### 2. Download Example Data

To run the quick mode example, you can download the example data as follows:

```bash
wget https://yanglab.westlake.edu.cn/data/gsMap/gsMap_example_data.tar.gz
tar -xvzf gsMap_example_data.tar.gz
```

Directory structure:
```bash
tree -L 2

gsMap_example_data/
├── GWAS
│   ├── BCX2_MCHC_EA_GWAMA.sumstats.gz
│   ├── GIANT_EUR_Height_2022_Nature.sumstats.gz
│   ├── gwas_config.yaml
│   └── IQ_NG_2018.sumstats.gz
└── ST
    └── E16.5_E1S1.MOSTA.h5ad
```

## Running `gsMap` in Quick Mode

<span style="color:#31a354"> Required memory: 80G (120K cells) </span>

```bash
gsmap quick_mode \
    --workdir './example_quick_mode/Mouse_Embryo' \
    --homolog_file 'gsMap_resource/homologs/mouse_human_homologs.txt' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --gsMap_resource_dir 'gsMap_resource' \
    --hdf5_path 'gsMap_example_data/ST/E16.5_E1S1.MOSTA.h5ad' \
    --annotation 'annotation' \
    --data_layer 'count' \
    --sumstats_file 'gsMap_example_data/GWAS/IQ_NG_2018.sumstats.gz' \
    --trait_name 'IQ'
```

### Parameters:

- `--workdir`: The working directory where output files will be saved.
- `--homolog_file`: The homologous gene file for converting gene names from different species to human.
- `--sample_name`: The name of the sample (e.g., `E16.5_E1S1.MOSTA`).
- `--gsMap_resource_dir`: Path to the directory containing the `gsMap` resources.
- `--hdf5_path`: Path to the input HDF5 file with spatial transcriptomics (ST) data.
- `--annotation`: The name of the annotation column in the `adata.obs` of the input HDF5 file.
- `--data_layer`: The layer of the gene expression matrix (e.g., `count`).
- `--sumstats_file`: Path to the GWAS summary statistics file.
- `--trait_name`: Name of the trait (e.g., `IQ`).

### Additional Options:

- If you want to analyze multiple traits at once, provide a configuration file (`--sumstats_config_file`) instead of a single summary statistics file:

```bash
gsmap quick_mode \
    --workdir './example_quick_mode/Mouse_Embryo' \
    --homolog_file 'gsMap_resource/homologs/mouse_human_homologs.txt' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --gsMap_resource_dir 'gsMap_resource' \
    --hdf5_path 'gsMap_example_data/ST/E16.5_E1S1.MOSTA.h5ad' \
    --annotation 'annotation' \
    --data_layer 'count' \
    --sumstats_config_file 'gsMap_example_data/GWAS/gwas_config.yaml'
```

The `gwas_config.yaml` file includes the following:

```yaml
Height: gsMap_example_data/GWAS/GIANT_EUR_Height_2022_Nature.sumstats.gz
IQ: gsMap_example_data/GWAS/IQ_NG_2018.sumstats.gz
SCZ: gsMap_example_data/GWAS/PGC3_SCZ_wave3_public_INFO80.sumstats.gz
```

### Output

- The output will be saved in the `--workdir` and will include intermediate files: latent representations, gene marker scores, LD scores. These intermediate files will be reused when you're running another GWAS trait of the same sample in this `--workdir`.
- A web report will be presented in the `report` folder, including visualizations of cell-trait associations and model diagnostic plots. Copy this folder to your local machine to open the HTML report in a web browser. This is a [example report](https://yanglab.westlake.edu.cn/data/gsMap/IQ/E16.5_E1S1.MOSTA_IQ_gsMap_Report.html) for the `IQ` trait.

```bash

### Example Output Structure

After running in quick mode, the following directory will be generated:

```bash
tree -L 3

example_quick_mode/Mouse_Embryo
├── E16.5_E1S1.MOSTA
│   ├── find_latent_representations # h5ad which contains latent representations
│   ├── latent_to_gene # gene marker scores
│   ├── generate_ldscore # LD scores
│   ├── spatial_ldsc # spatial cell-trait association results
│   ├── cauchy_combination # region-level or cell type-level association results
│   └── report # web report

```
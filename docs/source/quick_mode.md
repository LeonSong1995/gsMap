Here is the updated documentation for running `gsMap` in **quick mode** using the new `config.py` structure:

```markdown
# Quick Mode Example

## Preparation

Make sure you have installed the `gsMap` package before proceeding.

### 1. Download Dependencies

The `gsMap` package in quick mode requires the following resources:
- **Gene transfer format (GTF) file** for gene coordinates on the genome.
- **LD reference panel (PLINK bfile)** for computing LD scores.
- **SNP weight file** to adjust correlations between SNP statistics.
- **Homologous gene transformations file** (optional) for cross-species mapping.
- **Enhancer-gene mapping file** (optional) for SNP-to-gene enhancer linkages.

To download all the required files:
```bash
wget http://cnsgenomics.com/data/gsMap/gsMap_running_dependencies.tar.gz
tar -xvzf gsMap_running_dependencies.tar.gz
```

Directory structure:
```bash
tree -L 2

gsMap_resource
├── genome_annotation
│   ├── enhancer
│   └── gtf
├── LD_Reference_Panel
│   └── 1000G_EUR_Phase3_plink
└── LDSC_resource
    ├── hapmap3_snps
    └── weights_hm3_no_hla
```

### 2. Download Example Data

To run the quick mode example, you can download the example data as follows:

```bash
wget http://cnsgenomics.com/data/gsMap/gsMap_example_data.tar.gz
tar -xvzf gsMap_example_data.tar.gz
```

Directory structure:
```bash
tree -L 2

example_data
├── GWAS
│   ├── GIANT_EUR_Height_2022_Nature.sumstats.gz
│   ├── gwas_config.yaml
│   ├── IQ_NG_2018.sumstats.gz
│   └── BCX2_MCHC_EA_GWAMA.sumstats.gz
└── ST
    └── E16.5_E1S1.MOSTA.h5ad
```

## Running `gsMap` in Quick Mode

Quick mode allows you to run the entire `gsMap` pipeline in a streamlined manner, using pre-defined resources and configuration.

### Execution

The following command will execute the entire pipeline in quick mode:

```bash
gsmap quick_mode \
    --workdir './example/Mouse_Embryo' \
    --homolog_file 'gsMap_resource/homologs/mouse_human_homologs.txt' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --gsMap_resource_dir 'gsMap_resource' \
    --hdf5_path 'example_data/ST/E16.5_E1S1.MOSTA.h5ad' \
    --annotation 'annotation' \
    --data_layer 'count' \
    --sumstats_file 'example_data/GWAS/IQ_NG_2018.sumstats.gz' \
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

- If you want to analyze multiple traits in batch mode, provide a configuration file (`sumstats_config_file`) instead of a single summary statistics file:

```bash
gsmap quick_mode \
    --workdir './example/Mouse_Embryo' \
    --homolog_file 'gsMap_resource/homologs/mouse_human_homologs.txt' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --gsMap_resource_dir 'gsMap_resource' \
    --hdf5_path 'example_data/ST/E16.5_E1S1.MOSTA.h5ad' \
    --annotation 'annotation' \
    --data_layer 'count' \
    --sumstats_config_file 'example_data/GWAS/gwas_config.yaml'
```

The `gwas_config.yaml` file includes the following:

```yaml
Height: example_data/GWAS/GIANT_EUR_Height_2022_Nature.sumstats.gz
IQ: example_data/GWAS/IQ_NG_2018.sumstats.gz
SCZ: example_data/GWAS/PGC3_SCZ_wave3_public_INFO80.sumstats.gz
```

### Output

- The output will be saved in the `--workdir` directory and will include all the intermediate files, such as latent representations, gene marker scores, LD scores, and LDSC results.
- A final report will also be generated with spatial trait associations.

### Example Output Structure

After running in quick mode, the following directory structure will be created:

```bash
tree -L 2

example/Mouse_Embryo
├── E16.5_E1S1.MOSTA
│   ├── find_latent_representations
│   ├── latent_to_gene
│   ├── generate_ldscore
│   ├── spatial_ldsc
│   ├── report
│   └── cauchy_combination
```

## Summary

The `quick_mode`  lets you run the entire pipeline with a single command, making the process easier and reducing the need for manual steps.

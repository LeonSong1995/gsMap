# Mouse Embryo Example

## Preparation

Please ensure you have installed the `gsMap` package before proceeding.

### 1. Download dependencies

`gsMap` requires specific reference files:
- **Gene transfer format (GTF) file**, for gene coordinates on the genome.
- **LD reference panel (PLINK bfile)**, for computing LD scores.
- **SNP weight file**, to adjust correlations between SNP statistics.
- **Homologous gene transformations file** (optional), to map genes between species.
- **Enhancer-gene mapping file** (optional), for linking SNPs to genes based on enhancer regions.

To download the resources:
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

### 2. Download example data

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

## Running `gsMap`

### 1. find_latent_representations

**Objective**: Use GNN model to handle technical noise and find latent representations of each spot.

**Execution**:
```bash
gsmap run_find_latent_representations \
    --workdir './example/Mouse_Embryo' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --input_hdf5_path 'example_data/ST/E16.5_E1S1.MOSTA.h5ad' \
    --annotation 'annotation' \
    --data_layer 'count'
```

### 2. latent_to_gene

**Objective**: Generate gene marker scores by aggregating information from latent representations of neighboring spots.

**Execution**:
```bash
gsmap run_latent_to_gene \
    --workdir './example/Mouse_Embryo' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --annotation 'annotation' \
    --latent_representation 'latent_GVAE' \
    --num_neighbour 51 \
    --num_neighbour_spatial 201 \
    --homolog_file 'gsMap_resource/homologs/mouse_human_homologs.txt'
```

### 3. generate_ldscore

**Objective**: Assign gene specificity scores to SNPs using LD scores based on the GTF data.

**Execution**:

#### Method 1: Use TSS Only

```bash
for CHROM in {1..22}; do
    gsmap run_generate_ldscore \
        --workdir './example/Mouse_Embryo' \
        --sample_name 'E16.5_E1S1.MOSTA' \
        --chrom $CHROM \
        --bfile_root 'gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC' \
        --keep_snp_root 'gsMap_resource/LDSC_resource/hapmap3_snps/hm' \
        --gtf_annotation_file 'gsMap_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf' \
        --gene_window_size 50000
done
```

#### Method 2: Use Enhancer-Gene Linking Only

```bash
for CHROM in {1..22}; do
    gsmap run_generate_ldscore \
        --workdir './example/Mouse_Embryo' \
        --sample_name 'E16.5_E1S1.MOSTA' \
        --chrom $CHROM \
        --bfile_root 'gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC' \
        --keep_snp_root 'gsMap_resource/LDSC_resource/hapmap3_snps/hm' \
        --gtf_annotation_file 'gsMap_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf' \
        --enhancer_annotation_file 'gsMap_resource/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed' \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'enhancer_only'
done
```

#### Method 3: Use Both TSS and Enhancer-Gene Linking

```bash
for CHROM in {1..22}; do
    gsmap run_generate_ldscore \
        --workdir './example/Mouse_Embryo' \
        --sample_name 'E16.5_E1S1.MOSTA' \
        --chrom $CHROM \
        --bfile_root 'gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC' \
        --keep_snp_root 'gsMap_resource/LDSC_resource/hapmap3_snps/hm' \
        --gtf_annotation_file 'gsMap_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf' \
        --gene_window_size 50000 \
        --enhancer_annotation_file 'gsMap_resource/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed' \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'gene_window_first'
done
```

### 4. spatial_ldsc

**Objective**: Run spatial LDSC to map traits-associated spots.


```bash
gsmap run_spatial_ldsc \
    --workdir './example/Mouse_Embryo' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --trait_name 'IQ' \
    --sumstats_file 'example_data/GWAS/IQ_NG_2018.sumstats.gz' \
    --w_file 'gsMap_resource/LDSC_resource/weights_hm3_no_hla/weights.' \
    --num_processes 4
```


### 5. cauchy_combination (Optional)

**Objective**: Aggregate P values for spatial regions using the Cauchy Combination Test.

```bash
gsmap run_cauchy_combination \
    --workdir './example/Mouse_Embryo' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --trait_name 'IQ' \
    --annotation 'annotation'
```

### 6. Report Generation

Generate diagnostic reports and plots:

```bash
gsmap run_report \
    --workdir './example/Mouse_Embryo' \
    --sample_name 'E16.5_E1S1.MOSTA' \
    --trait_name 'IQ' \
    --sumstats_file 'example_data/GWAS/IQ_NG_2018.sumstats.gz' \
    --top_corr_genes 50
```

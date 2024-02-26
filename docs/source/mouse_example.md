
# Mouse Embryo Example

## 1. Preparation

First you should [install](install) the GPS package.

### 1.1 Download running dependencies

```bash

wget http://cnsgenomics.com/data/GPS/GPS_running_dependencies.tar.gz
tar -xvzf GPS_running_dependencies.tar.gz
```


### 1.2 Download example data

```bash
wget http://cnsgenomics.com/data/GPS/GPS_example_data.tar.gz
tar -xvzf GPS_example_data.tar.gz
```

Directory structure
```bash
tree -L 2
```

```
example_data
├── GWAS
│   ├── GIANT_EUR_Height_2022_Nature.sumstats.gz
│   ├── gwas_config.yaml
│   ├── IQ_NG_2018.sumstats.gz
│   └── PGC3_SCZ_wave3_public_INFO80.sumstats.gz
└── ST
    ├── Cortex_151507.h5ad
    └── E16.5_E1S1.MOSTA.h5ad
GPS_resource
├── genome_annotation
│   ├── enhancer
│   └── gtf
├── LD_Reference_Panel
│   └── 1000G_EUR_Phase3_plink
└── LDSC_resource
    ├── hapmap3_snps
    └── weights_hm3_no_hla

10 directories, 6 files
```



## 2. Run GPS

```shell
# Constants and configuration
WORKDIR='./example/Mouse_Embryo' # This should be the directory where the GPS output will be saved
SAMPLE_NAME="E16.5_E1S1" # This should be the name of the sample

# Input data
HDF5_PATH="example_data/ST/E16.5_E1S1.MOSTA.h5ad"
ANNOTATION="annotation" # This should be the cell type annotation of each spot
DATA_TYPE='count' # This should be the type of ST data, either 'count' or 'log1p'

# Running Dependencies and Resources
GTFFILE="GPS_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
ALL_ENHANCER_FILE="GPS_resource/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed"
BFILE_ROOT="GPS_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
KEEP_SNP_ROOT="GPS_resource/LDSC_resource/hapmap3_snps/hm"
W_FILE="GPS_resource/LDSC_resource/weights_hm3_no_hla/weights."
```


After setting up the constants and configuration, execute the following steps in order to run the GPS analysis on the human cortex data:

### 2.1 find_latent_representations

**Objective**: This initial step is for addressing technical noise and sparsity in spatial transcriptomics data. By employing the GNN model, the GPS finds latent representations for each spot.


**Input**: 
- HDF5 file containing ST data (`HDF5_PATH`) with cell type annotation (`ANNOTATION`).

**Output**: A .h5ad file containing the ST data with the refined latent representations for each spot.

**Execution**:
```shell
HDF5_WITH_LATENT_PATH="$WORKDIR/$SAMPLE_NAME/find_latent_representations/${SAMPLE_NAME}_add_latent.h5ad"
GPS run_find_latent_representations \
    --input_hdf5_path $HDF5_PATH \
    --sample_name $SAMPLE_NAME \
    --output_hdf5_path $HDF5_OUTPUT \
    --annotation $ANNOTATION \
    --type $DATA_TYPE
```
(latent_to_gene_mouse)=
### 2.2 latent_to_gene (Mouse)

**Objective**: Building on the latent representations, this step leverages the latent representations of the spot and its neighbors to generate gene marker scores.


**Input**:
- HDF5 file with latent representations from the previous step (`HDF5_WITH_LATENT_PATH`).

**Output**: 
- A .feather file with gene marker scores (`MKSCORE_FEATHER_PATH`).


**Execution**:
```shell
MKSCORE_FEATHER_PATH="$WORKDIR/$SAMPLE_NAME/latent_to_gene/${SAMPLE_NAME}_gene_marker_score.feather"
GPS run_latent_to_gene \
    --input_hdf5_with_latent_path $HDF5_WITH_LATENT_PATH \
    --sample_name $SAMPLE_NAME \
    --output_feather_path $MKSCORE_FEATHER_PATH \
    --latent_representation "latent_GVAE" \
    --num_neighbour 51 \
    --num_neighbour_spatial 201 \
    --annotation $ANNOTATION \
    --type $DATA_TYPE \
    --species MOUSE_GENE_SYM \
    --gs_species "homologs/mouse_human_homologs.txt"
```


### 2.3 generate_ldscore

**Objective**: By using the gene marker scores from the previous step, GPS first assigns gene specificity scores to each SNP by referencing the gene annotation data. Specifically, GPS assigns gene specificity scores by evaluating each SNP's closeness to the gene transcription start sites (TSS) and, optionally, enhancer-gene linking maps. Then use the gene specificity scores to generates LD scores for each SNP.



**Input**:
- Feather file with gene marker scores from the previous step (`MKSCORE_FEATHER_FILE`).
- Reference panel data (`BFILE_ROOT` and `KEEP_SNP_ROOT`).
- Gene annotation file.

**Output**: 
- A set of LD score chunks stored in the (`LDScoreDir`).


**Three SNP to gene linking methods are available:**
````{tab} 1. Use TSS Only
This will use TSS only to link SNPs to gene specificity.

`--gene_window_size = 50000` is the window size around the gene body to consider for gene specificity scores.
If a SNP is within the gene window, it will be assigned the gene specificity score of that gene. If a SNP is within the window of multiple genes, the nearest gene will be used.

**Execution**:


```shell
LDScoreDir="$WORKDIR/$SAMPLE_NAME/generate_ldscore"
for CHROM in {1..22}; do
    GPS run_generate_ldscore \
        --sample_name $SAMPLE_NAME \
        --chrom $CHROM \
        --ldscore_save_dir $LDScoreDir \
        --mkscore_feather_file $MKSCORE_FEATHER_PATH \
        --bfile_root $BFILE_ROOT \
        --keep_snp_root $KEEP_SNP_ROOT \
        --gtf_annotation_file $GTFFILE \
        --gene_window_size 50000 \
        --spots_per_chunk 1000 \
        --ld_wind 1 \
        --ld_unit "CM"
done
```
````

`````{tab} 2. Use Enhancer-Gene Linking Only


When a SNP mapped to multiple enhancers, the gene specificity score of the SNP will be set by the `--snp_multiple_enhancer_strategy` parameter. The default value is `max_mkscore`, which means the gene specificity score of the SNP will be set to the maximum gene specificity score of the enhancers that the SNP mapped to. Possible choices: max_mkscore, nearest_TSS

In this example we choose the all tissue enhancer annotation file (`ENHANCER_ANNOTATION_FILE`).

**Execution**:

```shell
LDScoreDir="$WORKDIR/$SAMPLE_NAME/generate_ldscore"
for CHROM in {1..22}; do
    GPS run_generate_ldscore \
        --sample_name $SAMPLE_NAME \
        --chrom $CHROM \
        --ldscore_save_dir $LDScoreDir \
        --mkscore_feather_file $MKSCORE_FEATHER_PATH \
        --bfile_root $BFILE_ROOT \
        --keep_snp_root $KEEP_SNP_ROOT \
        --gtf_annotation_file $GTFFILE \
        --enhancer_annotation_file $ENHANCER_ANNOTATION_FILE \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'enhancer_only' \
        --spots_per_chunk 1000 \
        --ld_wind 1 \
        --ld_unit "CM"
done

```


`````


`````{tab} 3. Use Both TSS and Enhancer-Gene Linking
This will use both TSS and enhancer-gene linking to link SNPs to gene specificity.

In the scenario where a single nucleotide polymorphism (SNP) falls within the proximity window of both a gene body window and an enhancer linked to a different gene, the determination of the SNP's gene specificity score is governed by the `--gene_window_enhancer_priority` parameter. Possible choices is `gene_window_first` or `enhancer_first`.


**Execution**:

```shell
LDScoreDir="$WORKDIR/$SAMPLE_NAME/generate_ldscore"
for CHROM in {1..22}; do
    GPS run_generate_ldscore \
        --sample_name $SAMPLE_NAME \
        --chrom $CHROM \
        --ldscore_save_dir $LDScoreDir \
        --mkscore_feather_file $MKSCORE_FEATHER_PATH \
        --bfile_root $BFILE_ROOT \
        --keep_snp_root $KEEP_SNP_ROOT \
        --gtf_annotation_file $GTFFILE \
        --gene_window_size 50000 \
        --enhancer_annotation_file ENHANCER_ANNOTATION_FILE \
        --snp_multiple_enhancer_strategy 'max_mkscore' \
        --gene_window_enhancer_priority 'gene_window_first' \
        --spots_per_chunk 1000 \
        --ld_wind 1 \
        --ld_unit "CM"
done
```


`````

```{caution}
If you runout of memory in this step or the next step,
you can reduce the `--spots_per_chunk` parameter to a smaller value.

In general, 40GB memory is required when `--spots_per_chunk` is set to 1000.
```

### 2.4 spatial_ldsc

**Objective**: Identify spots associated with complex traits using spatial S-LDSC method.


**Input**:
- Directory containing LD scores generated in the previous step (`LDScoreDir`).

**Output**: 
- Spatial LDSC results saved in the specified directory (`LDSC_SAVE_DIR`) for each trait.


`````{tab} Run single trait

**Execution**:

```shell
LDSC_DIR="$WORKDIR/$SAMPLE_NAME/spatial_ldsc"
SUMSTATS_FILE="example_data/GWAS/IQ_NG_2018.sumstats.gz"
TRAIT_NAME="IQ"

GPS run_spatial_ldsc \
    --sumstats_file $SUMSTATS_FILE \
    --trait_name $TRAIT_NAME \
    --w_file $W_FILE \
    --sample_name $SAMPLE_NAME \
    --num_processes 4 \
    --ldscore_input_dir $LDScoreDir \
    --ldsc_save_dir $LDSC_DIR \
 ```

`````

`````{tab} Run multiple traits in batch

**Execution**:

You could run multiple traits in batch by providing a YAML file with the `--sumstats_config_file` parameter. An example of `sumstats_config_file` is like this:

```yaml
Height: example_data/GWAS/GIANT_EUR_Height_2022_Nature.sumstats.gz
IQ: example_data/GWAS/IQ_NG_2018.sumstats.gz
SCZ: example_data/GWAS/PGC3_SCZ_wave3_public_INFO80.sumstats.gz
```

```shell
LDSC_DIR="$WORKDIR/$SAMPLE_NAME/spatial_ldsc"
SUMSTATS_CONFIG_FILE="example_data/GWAS/gwas_config.yaml"
GPS run_spatial_ldsc \
    --w_file $W_FILE \
    --sample_name $SAMPLE_NAME \
    --num_processes 4 \
    --ldscore_input_dir $LDScoreDir \
    --ldsc_save_dir $LDSC_DIR \
    --sumstats_config_file $SUMSTATS_CONFIG_FILE
```
`````


## 3. Post-processing


### 2.5 cauchy_combination (optional)

**Objective**: Use the Cauchy Combination Test to aggregate P values of individual spots within specific spatial regions or cell types to evaluate the association of these regions with complex traits.

**Execution**:

```shell
CAUCHY_SAVE_DIR="$WORKDIR/$SAMPLE_NAME/cauchy_combination"
TRAIT_NAME="IQ"
GPS run_cauchy_combination \
    --input_hdf5_path $HDF5_PATH \
    --input_ldsc_dir $LDSC_DIR \
    --sample_name $SAMPLE_NAME \
    --output_cauchy_dir $CAUCHY_SAVE_DIR \
    --trait_name $TRAIT_NAME \
    --annotation $ANNOTATION
```

You will get a csv file with the aggregated P values for each region or cell type to find out the most significant regions or cell types associated with the complex trait.

`````{tab} Height

| annotation           | p_cauchy               | p_median            |
|----------------------|------------------------|---------------------|
| Adipose tissue       | 7.99e-15               | 1.11e-08            |
| Adrenal gland        | 2.69e-13               | 2.20e-08            |
| Brain                | 1.00                   | 0.999               |
| Cartilage            | 5.49e-20               | 6.49e-09            |
| Connective tissue    | 1.74e-16               | 6.94e-09            |
| ...                  | ...                    | ...                 |
`````



`````{tab} IQ

| annotation           | p_cauchy               | p_median            |
|----------------------|------------------------|---------------------|
| Adipose tissue       | 0.0277                 | 0.4878              |
| Adrenal gland        | 0.0171                 | 0.1874              |
| Brain                | 1.86e-19               | 1.60e-12            |
| Cartilage            | 0.00016                | 0.1277              |
| Connective tissue    | 8.97e-05               | 0.201               |
| ...                  | ...                    | ...                 |

`````

`````{tab} MCHC

| annotation           | p_cauchy               | p_median            |
|----------------------|------------------------|---------------------|
| Adipose tissue       | 2.73e-08               | 0.0361              |
| Adrenal gland        | 1.54e-06               | 4.85e-05            |
| Brain                | 1.00                   | 1.00                |
| Cartilage            | 0.00155                | 0.158               |
| Connective tissue    | 4.11e-07               | 0.0143              |
| ...                  | ...                    | ...                 |


`````

### 3.2 visualization

```shell
GPS run_visualize \
    --input_hdf5_path $HDF5_PATH \
    --input_ldsc_dir $LDSC_DIR \
    --output_figure_dir $WORKDIR/$SAMPLE_NAME/figures \
    --sample_name $SAMPLE_NAME \
    --trait_name $TRAIT_NAME \
    --fig_title $TRAIT_NAME \
    --annotation $ANNOTATION \
    --point_size 7
    
```
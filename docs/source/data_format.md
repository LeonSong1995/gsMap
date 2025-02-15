(data-format)=
# Data Format


## ST Data

The input ST data must be an h5ad file containing at least the gene expression matrix and spatial coordinates. The gene expression matrix should be stored in the `layers` attribute, and the spatial coordinates should be in the obsm `attribute` with the key `spatial`. Optionally, the h5ad file may include spot (cell) annotations in the obs attribute.
```python
import scanpy as sc

adata = sc.read_h5ad("gsMap_example_data/ST/E16.5_E1S1.MOSTA.h5ad")

print(adata.layers["count"].shape)
print(adata.obsm["spatial"].shape)
print(adata.obs["annotation"].value_counts().head())
```



## GWAS Data
The input GWAS data is a text file containing at least the columns for SNP (rs number), Z (Z-statistics), and N (sample size).  Column headers are keywords used by gsMap.

```shell
zcat gsMap_example_data/GWAS/IQ_NG_2018.sumstats.gz | head -n 5

SNP		A1	A2	Z	N
rs12184267	T	C	0.916	225955
rs12184277	G	A	0.656	226215
rs12184279	A	C	1.050	226224
rs116801199	T	G	0.300	226626
```

###  How to format the GWAS data
You can convert GWAS summary data into the required format using custom code. For convenience, gsMap provides a command to do this. Below is an example of how to use the command.

Download the human height GWAS data and decompress it.
```bash
wget https://portals.broadinstitute.org/collaboration/giant/images/4/4e/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz

gzip -d GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz
```

Convert the summary statistics to the required format.
```bash
gsmap format_sumstats \
--sumstats 'GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL' \
--out 'HEIGHT'
```
You will obtain a file named HEIGHT.sumstats.gz

```bash
zcat HEIGHT.sumstats.gz | head -n 5

SNP		A1	A2	Z	N
rs3131969	G	A	0.328	1494218.000
rs3131967	C	T	0.386	1488150.000
rs12562034	A	G	1.714	1554976.000
rs4040617	G	A	-0.463	1602016.000
```

For more usage options, please refer to:
```bash
gsMap format_sumstats -h
```

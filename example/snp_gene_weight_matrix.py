
# %%
# set debug

# %%

# %%
import pandas as pd
import numpy as np
import pyranges as pr
# %%
from GPS.generate_r2_matrix import PlinkBEDFileWithR2Cache, getBlockLefts, ID_List_Factory


# %%
# load gtf
def load_gtf(gtf_file, mk_score, window_size):
    """
    Load the gene annotation file (gtf).
    """
    print("Loading gtf data")
    #
    # Load GTF file
    gtf = pr.read_gtf(gtf_file)
    gtf = gtf.df
    #
    # Select the common genes
    gtf = gtf[gtf['Feature'] == 'gene']
    common_gene = np.intersect1d(mk_score.index, gtf.gene_name)
    #
    gtf = gtf[gtf.gene_name.isin(common_gene)]
    mk_score = mk_score[mk_score.index.isin(common_gene)]
    #
    # Remove duplicated lines
    gtf = gtf.drop_duplicates(subset='gene_name', keep="first")
    #
    # Process the GTF (open 100-KB window: Tss - Ted)
    gtf_bed = gtf[['Chromosome', 'Start', 'End', 'gene_name', 'Strand']].copy()
    gtf_bed.loc[:, 'TSS'] = gtf_bed['Start']
    gtf_bed.loc[:, 'TED'] = gtf_bed['End']

    gtf_bed.loc[:, 'Start'] = gtf_bed['TSS'] - window_size
    gtf_bed.loc[:, 'End'] = gtf_bed['TED'] + window_size
    gtf_bed.loc[gtf_bed['Start'] < 0, 'Start'] = 0
    #
    # Correct the negative strand
    tss_neg = gtf_bed.loc[gtf_bed['Strand'] == '-', 'TSS']
    ted_neg = gtf_bed.loc[gtf_bed['Strand'] == '-', 'TED']
    gtf_bed.loc[gtf_bed['Strand'] == '-', 'TSS'] = ted_neg
    gtf_bed.loc[gtf_bed['Strand'] == '-', 'TED'] = tss_neg
    gtf_bed = gtf_bed.drop('Strand', axis=1)
    #
    # Transform the GTF to PyRanges
    gtf_pr = pr.PyRanges(gtf_bed)
    return gtf_pr, mk_score


# %%
# load gtf
def load_gtf_all_gene(gtf_file, window_size):
    """
    Load the gene annotation file (gtf).
    """
    print("Loading gtf data")
    #
    # Load GTF file
    gtf = pr.read_gtf(gtf_file)
    gtf = gtf.df
    #
    # Select the common genes
    gtf = gtf[gtf['Feature'] == 'gene']

    # Remove duplicated lines
    gtf = gtf.drop_duplicates(subset='gene_name', keep="first")
    #
    # Process the GTF (open 100-KB window: Tss - Ted)
    gtf_bed = gtf[['Chromosome', 'Start', 'End', 'gene_name', 'Strand']].copy()
    gtf_bed.loc[:, 'TSS'] = gtf_bed['Start']
    gtf_bed.loc[:, 'TED'] = gtf_bed['End']

    gtf_bed.loc[:, 'Start'] = gtf_bed['TSS'] - window_size
    gtf_bed.loc[:, 'End'] = gtf_bed['TED'] + window_size
    gtf_bed.loc[gtf_bed['Start'] < 0, 'Start'] = 0
    #
    # Correct the negative strand
    tss_neg = gtf_bed.loc[gtf_bed['Strand'] == '-', 'TSS']
    ted_neg = gtf_bed.loc[gtf_bed['Strand'] == '-', 'TED']
    gtf_bed.loc[gtf_bed['Strand'] == '-', 'TSS'] = ted_neg
    gtf_bed.loc[gtf_bed['Strand'] == '-', 'TED'] = tss_neg
    gtf_bed = gtf_bed.drop('Strand', axis=1)
    #
    # Transform the GTF to PyRanges
    gtf_pr = pr.PyRanges(gtf_bed)
    return gtf_pr


# %%
def load_marker_score(mk_score_file):
    """
    Load marker scores of each cell.
    """
    mk_score = pd.read_feather(mk_score_file).set_index('HUMAN_GENE_SYM').rename_axis('gene_name')
    mk_score.insert(0, 'all_gene', 1)
    return mk_score


# %%
# config = MakeAnnotationConfig(
#     input_feather_file=f'{test_dir}/{name}/gene_markers/{name}_rank.feather',
#     sample_name=name,
#     output_dir=f'{test_dir}/{name}/snp_annotation',
#     gtf_file='/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf',
#     bfile_root='/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC',
#     baseline_annotation=None,
#     keep_snp_root='/storage/yangjianLab/sharedata/LDSC_resource/hapmap3_snps/hm',
#     chr=TASK_ID,
#     window_size=50000,
#     cells_per_chunk=500,
#     ld_wind=1,
#     ld_wind_unit='CM',
#     r2_cache_dir='/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/r2_matrix',
#     use_gpu=False,
#     snps_per_chunk=100_000
# )

# %%
name = 'Cortex_151507'
TASK_ID = 2
test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021'
# %%
gtf_file = '/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf'
feather_file = f'{test_dir}/{name}/gene_markers/{name}_rank.feather'
bfile_root = '/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC'
# %%
mk_score = load_marker_score(feather_file)
# gtf_pr,mk_score_common = load_gtf(gtf_file, mk_score, window_size=50000)
gtf_pr = load_gtf_all_gene(gtf_file, window_size=50000)


# %%
# load mkscore get common gene
# %%
# load bim
def load_bim(bfile_root, chrom):
    """
    Load the bim file.
    """
    print("Loading bim data")
    bim = pd.read_csv(f'{bfile_root}.{chrom}.bim', sep='\t', header=None)
    bim.columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
    #
    # Transform bim to PyRanges
    bim_pr = bim.copy()
    bim_pr.columns = ["Chromosome", "SNP", "CM", "Start", "A1", "A2"]
    bim_pr['End'] = bim_pr['Start']
    bim_pr = pr.PyRanges(bim_pr)
    bim_pr.Chromosome = f'chr{chrom}'
    return bim, bim_pr


# %%
# get SNP gene pair
bim, bim_pr = load_bim(bfile_root, TASK_ID)


# %%
def Overlaps_gtf_bim(gtf_pr, bim_pr):
    """
    Find overlaps between gtf and bim file.
    """
    # Select the overlapped regions (SNPs in gene windows)
    overlaps = gtf_pr.join(bim_pr)
    overlaps = overlaps.df
    overlaps['Distance'] = np.abs(overlaps['Start_b'] - overlaps['TSS'])
    overlaps_small = overlaps.copy()
    overlaps_small = overlaps_small.loc[overlaps_small.groupby('SNP').Distance.idxmin()]
    return overlaps_small


# %%
overlaps_small = Overlaps_gtf_bim(gtf_pr, bim_pr)
# %%
overlaps_small
# %%
annot = bim[["CHR", "BP", "SNP", "CM"]]
SNP_gene_pair = overlaps_small[['SNP', 'gene_name']].set_index('SNP').join(annot.set_index('SNP'), how='right')
# %%
SNP_gene_pair
# %%
bim
# %%
SNP_gene_pair_dummy = pd.get_dummies(SNP_gene_pair['gene_name'])
# %%
SNP_gene_pair_dummy.shape
# %%
SNP_gene_pair_dummy.memory_usage()


# %%
def get_snp_gene_dummy(chrom, bfile_root, gtf_pr, ):
    """
    Get the dummy matrix of SNP-gene pairs.
    """
    # Load the bim file
    bim, bim_pr = load_bim(bfile_root, chrom)
    # Find overlaps between gtf and bim file
    overlaps_small = Overlaps_gtf_bim(gtf_pr, bim_pr)
    # Get the SNP-gene pair
    annot = bim[["CHR", "BP", "SNP", "CM"]]
    SNP_gene_pair = overlaps_small[['SNP', 'gene_name']].set_index('SNP').join(annot.set_index('SNP'), how='right')
    # Get the dummy matrix
    SNP_gene_pair_dummy = pd.get_dummies(SNP_gene_pair['gene_name'])
    return SNP_gene_pair_dummy


# %%
# snp_file, snp_obj = bfile + '.bim', PlinkBIMFile
# array_snps = snp_obj(snp_file)
# m = len(array_snps.IDList)
# print(f'Read list of {m} SNPs from {snp_file}')
# #
# # Load fam
# ind_file, ind_obj = bfile + '.fam', PlinkFAMFile
# array_indivs = ind_obj(ind_file)
# n = len(array_indivs.IDList)
# print(f'Read list of {n} individuals from {ind_file}')
#         array_file, array_obj = bfile + '.bed', PlinkBEDFileWithR2Cache
#         geno_array = array_obj(array_file, n, array_snps, keep_snps=None, keep_indivs=None, mafMin=None)
# # # Load the annotations of the baseline
# # if self.ld_wind_unit == 'SNP':
# #     max_dist = self.ld_wind
# #     coords = np.array(range(geno_array.m))
# # elif self.ld_wind_unit == 'BP':
# #     max_dist = self.ld_wind * 1000
# #     coords = np.array(array_snps.df['BP'])[geno_array.kept_snps]
# # elif self.ld_wind_unit == 'CM':
# #     max_dist = self.ld_wind
# #     coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]
# # block_left = getBlockLefts(coords, max_dist)
# %%
def get_ldscore(bfile_root, chrom, annot_matrix, ld_wind, ld_unit):
    PlinkBIMFile = ID_List_Factory(['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], 1, '.bim', usecols=[0, 1, 2, 3, 4, 5])
    PlinkFAMFile = ID_List_Factory(['IID'], 0, '.fam', usecols=[1])

    bfile = f'{bfile_root}.{chrom}'
    snp_file, snp_obj = bfile + '.bim', PlinkBIMFile
    array_snps = snp_obj(snp_file)
    m = len(array_snps.IDList)
    print(f'Read list of {m} SNPs from {snp_file}')

    # Load fam
    ind_file, ind_obj = bfile + '.fam', PlinkFAMFile
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    print(f'Read list of {n} individuals from {ind_file}')
    array_file, array_obj = bfile + '.bed', PlinkBEDFileWithR2Cache
    geno_array = array_obj(array_file, n, array_snps, keep_snps=None, keep_indivs=None, mafMin=None)
    # Load the annotations of the baseline
    if ld_unit == 'SNP':
        max_dist = ld_wind
        coords = np.array(range(geno_array.m))
    elif ld_unit == 'BP':
        max_dist = ld_wind * 1000
        coords = np.array(array_snps.df['BP'])[geno_array.kept_snps]
    elif ld_unit == 'CM':
        max_dist = ld_wind
        coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]
    else:
        raise ValueError(f'Invalid ld_wind_unit: {ld_unit}')
    block_left = getBlockLefts(coords, max_dist)
    # Calculate the LD score
    lN_df = pd.DataFrame(geno_array.ldScoreVarBlocks(block_left, 100, annot=annot_matrix))
    return lN_df


# %%
def calculate_SNP_Gene_weight_matrix(chrom, bfile_root, gtf_pr, ):
    """
    Calculate the SNP-gene weight matrix.
    """
    # Get the dummy matrix
    SNP_gene_pair_dummy = get_snp_gene_dummy(chrom, bfile_root, gtf_pr)
    # Get the SNP-gene weight matrix
    lN_df = get_ldscore(bfile_root, chrom, SNP_gene_pair_dummy, 1, 'SNP')
    return lN_df


# %%
# Get the SNP-gene weight matrix
lN_df = get_ldscore(bfile_root, chrom=TASK_ID, annot_matrix=SNP_gene_pair_dummy.values, ld_wind=1, ld_unit='CM')
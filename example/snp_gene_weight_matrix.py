from pathlib import Path

# %%
import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import trange, tqdm

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
    mk_score = mk_score.astype(np.float32, copy=False)
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
CHROM = 2
test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021'
# %%
gtf_file = '/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf'
feather_file = f'{test_dir}/{name}/gene_markers/{name}_rank.feather'
bfile_root = '/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC'
# %%
mk_score = load_marker_score(feather_file)
gtf_pr, mk_score_common = load_gtf(gtf_file, mk_score, window_size=50000)


# gtf_pr = load_gtf_all_gene(gtf_file, window_size=50000)


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
def get_snp_gene_dummy(chrom, bfile_root, gtf_pr, dummy_na=True):
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
    SNP_gene_pair_dummy = pd.get_dummies(SNP_gene_pair['gene_name'], dummy_na=dummy_na)
    return SNP_gene_pair_dummy


# %%
def get_snp_pass_maf(bfile_root, chrom, maf_min=0.05):
    """
    Get the dummy matrix of SNP-gene pairs.
    """
    # Load the bim file
    PlinkBIMFile = ID_List_Factory(['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], 1, '.bim', usecols=[0, 1, 2, 3, 4, 5])
    PlinkFAMFile = ID_List_Factory(['IID'], 0, '.fam', usecols=[1])

    bfile = f'{bfile_root}.{chrom}'
    snp_file, snp_obj = bfile + '.bim', PlinkBIMFile
    array_snps = snp_obj(snp_file)
    m = len(array_snps.IDList)

    # Load fam
    ind_file, ind_obj = bfile + '.fam', PlinkFAMFile
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    array_file, array_obj = bfile + '.bed', PlinkBEDFileWithR2Cache
    geno_array = array_obj(array_file, n, array_snps, keep_snps=None, keep_indivs=None, mafMin=None)
    ii = geno_array.maf > maf_min
    snp_pass_maf = array_snps.IDList[ii]
    print(f'After filtering SNPs with MAF < {maf_min}, {len(snp_pass_maf)} SNPs remain.')
    return snp_pass_maf.SNP.to_list()


def get_ldscore(bfile_root, chrom, annot_matrix, ld_wind, ld_unit='CM'):
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
    elif ld_unit == 'KB':
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
def calculate_SNP_Gene_weight_matrix(SNP_gene_pair_dummy, chrom, bfile_root, ld_wind=1, ld_unit='CM'):
    """
    Calculate the SNP-gene weight matrix.
    """
    # Get the dummy matrix
    # Get the SNP-gene weight matrix
    snp_gene_weight_matrix = get_ldscore(bfile_root, chrom, SNP_gene_pair_dummy.values, ld_wind=ld_wind,
                                         ld_unit=ld_unit)
    snp_gene_weight_matrix = snp_gene_weight_matrix.astype(np.float32, copy=False)
    snp_gene_weight_matrix.index = SNP_gene_pair_dummy.index
    snp_gene_weight_matrix.columns = SNP_gene_pair_dummy.columns
    return snp_gene_weight_matrix


# %%
def calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(snp_gene_weight_matrix,
                                                          mk_score_chunk,
                                                          chrom,
                                                          save_file_name,
                                                          keep_snp_root=None,
                                                          ):
    save_dir = Path(save_file_name).parent
    save_dir.mkdir(parents=True, exist_ok=True)
    if keep_snp_root is not None:
        keep_snp = pd.read_csv(f'{keep_snp_root}.{chrom}.snp', header=None)[0].to_list()
        keep_snp = np.intersect1d(keep_snp, snp_gene_weight_matrix.index)
    ldscore_chr_chunk = snp_gene_weight_matrix @ mk_score_chunk
    ldscore_chr_chunk = ldscore_chr_chunk.astype(np.float16, copy=False)
    ldscore_chr_chunk = ldscore_chr_chunk if keep_snp_root is None else ldscore_chr_chunk.loc[keep_snp]
    # save for each chunk
    ldscore_chr_chunk.reset_index().to_feather(save_file_name)


def calculate_M_use_SNP_gene_pair_dummy_by_chunk(SNP_gene_pair_dummy: pd.DataFrame,
                                                 mk_score_chunk,
                                                 snp_pass_maf,
                                                 M_file_path, M_5_file_path,
                                                 ):
    '''
    calculate M use SNP_gene_pair_dummy_sumed_along_snp_axis and mk_score_chunk
    '''
    SNP_gene_pair_dummy_sumed_along_snp_axis = SNP_gene_pair_dummy.values.sum(axis=0, keepdims=True)
    SNP_gene_pair_dummy_sumed_along_snp_axis_pass_maf = SNP_gene_pair_dummy.loc[snp_pass_maf].values.sum(axis=0,
                                                                                                         keepdims=True)
    save_dir = Path(M_file_path).parent
    save_dir.mkdir(parents=True, exist_ok=True)
    M_chr_chunk = SNP_gene_pair_dummy_sumed_along_snp_axis @ mk_score_chunk
    M_5_chr_chunk = SNP_gene_pair_dummy_sumed_along_snp_axis_pass_maf @ mk_score_chunk
    np.savetxt(M_file_path, M_chr_chunk, delimiter='\t', )
    np.savetxt(M_5_file_path, M_5_chr_chunk, delimiter='\t', )


def calculate_ldscore_use_SNP_Gene_weight_matrix_by_chr(snp_gene_weight_matrix, mk_score_common, spots_per_chunk,
                                                        save_dir, chrom,
                                                        snp_pass_maf,
                                                        keep_snp_root=None):
    """
    Calculate the LD score using the SNP-gene weight matrix.
    """
    # Calculate the LD score

    for i in trange(0, mk_score_common.shape[1], spots_per_chunk, desc=f'Calculating LD score by chunk for chr{chrom}'):
        mk_score_chunk = mk_score_common.iloc[:, i:i + spots_per_chunk]
        save_file_name = f'{save_dir}/chr{chrom}_ldscore_{i}.feather'
        calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(snp_gene_weight_matrix,
                                                              mk_score_chunk,
                                                              chrom,
                                                              save_file_name,
                                                              keep_snp_root=keep_snp_root)
        calculate_M_use_SNP_gene_pair_dummy_by_chunk(snp_gene_weight_matrix,
                                                     mk_score_chunk,
                                                     snp_pass_maf,
                                                     f'{save_dir}/chr{chrom}_M_{i}',
                                                     f'{save_dir}/chr{chrom}_M_5_{i}',
                                                     )


# %%
keep_snp_root = '/storage/yangjianLab/sharedata/LDSC_resource/hapmap3_snps/hm'
keep_snp = pd.read_csv(f'{keep_snp_root}.{CHROM}.snp', header=None)[0].to_list()
num_snp = len(keep_snp)
# %%
SNP_gene_pair_dummy = get_snp_gene_dummy(CHROM, bfile_root, gtf_pr, dummy_na=True)
snp_gene_weight_matrix = calculate_SNP_Gene_weight_matrix(SNP_gene_pair_dummy, CHROM, bfile_root, ld_wind=1,
                                                          ld_unit='CM')
# %%
# save baseline ld score
baseline_mk_score = np.ones((snp_gene_weight_matrix.shape[1], 2))
baseline_mk_score[-1, 0] = 0  # all_gene
baseline_mk_score_df = pd.DataFrame(baseline_mk_score, index=snp_gene_weight_matrix.columns,
                                    columns=['all_gene', 'base'])
baseline_ldscore_save_path = f'{test_dir}/{name}/snp_annotation/test/baseline/baseline.{CHROM}.feather'

calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(snp_gene_weight_matrix,
                                                      baseline_mk_score_df,
                                                      CHROM,
                                                      baseline_ldscore_save_path,
                                                      keep_snp_root=keep_snp_root)
# %%
# save baseline M
# %%
snp_pass_maf = get_snp_pass_maf(bfile_root, CHROM, maf_min=0.05)
calculate_M_use_SNP_gene_pair_dummy_by_chunk(SNP_gene_pair_dummy,
                                             baseline_mk_score_df,
                                             snp_pass_maf,
                                             f'{test_dir}/{name}/snp_annotation/test/baseline/baseline.{CHROM}.M',
                                             f'{test_dir}/{name}/snp_annotation/test/baseline/baseline.{CHROM}.M_5', )
# %%
common_gene_chr = SNP_gene_pair_dummy.columns[:-1]
save_ldscore = calculate_ldscore_use_SNP_Gene_weight_matrix_by_chr(snp_gene_weight_matrix.iloc[:, :-1],
                                                                   mk_score_common.loc[common_gene_chr],
                                                                   spots_per_chunk=10_000,
                                                                   save_dir=f'{test_dir}/{name}/snp_annotation/test/ldscore',
                                                                   chrom=CHROM,
                                                                   snp_pass_maf=snp_pass_maf,
                                                                   keep_snp_root=keep_snp_root)

# %%
# calcualte LDscore
mk_score_chr = mk_score_common.loc[SNP_gene_pair_dummy.columns].values
# # %%
# mk_score_chr_large= mk_score_chr.repeat(2, axis=1).astype(np.float16, )
# ldscore_chr_large= lN_df @ mk_score_chr_large

# %%
true_set_path = '/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/annotation/Cortex_151507/snp_annotation/baseline/baseline.2.l2.ldscore.feather'
true_set_df = pd.read_feather(true_set_path).set_index('SNP')

# %%
now = pd.read_feather(baseline_ldscore_save_path)
now_df = now.set_index('SNP')
# %%
before_path = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/Cortex_151507/snp_annotation/baseline/baseline.2.l2.ldscore.feather'
before_df = pd.read_feather(before_path).set_index('SNP')

# %%
# calcualte correlaton coefficient before and true
before_df.loc[true_set_df.index].corrwith(true_set_df)
# %%
# calcualte correlaton coefficient now and true
now_df.loc[true_set_df.index].corrwith(true_set_df)

# %%
true_chunk1_path = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/Cortex_151507/snp_annotation/Cortex_151507_chunk1/Cortex_151507.2.l2.ldscore.feather'
true_chunk1_df = pd.read_feather(true_chunk1_path).set_index('SNP')
# %%
now_chunk1_path = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/Cortex_151507/snp_annotation/test/ldscore/chr2_ldscore_0.feather'
now_chunk1_df = pd.read_feather(now_chunk1_path).set_index('SNP')
# %%
now_chunk1_df.loc[true_chunk1_df.index, true_chunk1_df.columns[6:]].corrwith(
    true_chunk1_df.loc[:, true_chunk1_df.columns[6:]])

# %%
now_new = pd.concat([true_chunk1_df.iloc[:, :6], now_chunk1_df.loc[true_chunk1_df.index, true_chunk1_df.columns[5:6]]],
                    axis=1)

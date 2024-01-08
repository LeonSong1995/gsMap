import argparse
from pathlib import Path

import numpy as np
# %%
import pandas as pd
import pyranges as pr
from scipy.sparse import csr_matrix
from tqdm import trange

from GPS.config import add_generate_ldscore_args, GenerateLDScoreConfig
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
    mk_score = mk_score.astype(np.float32, copy=False)
    return mk_score


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
    # %timeit snp_gene_weight_matrix @ mk_score_chunk
    # snp_gene_weight_matrix_sparse_float16 = csr_matrix(snp_gene_weight_matrix)
    # %timeit snp_gene_weight_matrix_sparse_float16 @ mk_score_chunk
    # snp_gene_weight_matrix_sparse_float16.data = snp_gene_weight_matrix_sparse_float16.data.astype(np.float16)
    # %time snp_gene_weight_matrix_sparse_float16 @ mk_score_chunk.values.astype(np.float16, )
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
                                                        save_dir, chrom, snp_pass_maf, sample_name, keep_snp_root=None):
    """
    Calculate the LD score using the SNP-gene weight matrix.
    :param sample_name:
    """
    # Calculate the LD score
    chunk_index = 1
    for i in trange(0, mk_score_common.shape[1], spots_per_chunk, desc=f'Calculating LD score by chunk for chr{chrom}'):
        mk_score_chunk = mk_score_common.iloc[:, i:i + spots_per_chunk]

        ld_score_file = f'{save_dir}/{sample_name}_chunk{chunk_index}/{sample_name}.{chrom}.l2.ldscore.feather'
        M_file = f'{save_dir}/{sample_name}_chunk{chunk_index}/{sample_name}.{chrom}.l2.M'
        M_5_file = f'{save_dir}/{sample_name}_chunk{chunk_index}/{sample_name}.{chrom}.l2.M_5_50'

        calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(snp_gene_weight_matrix,
                                                              mk_score_chunk,
                                                              chrom,
                                                              save_file_name=ld_score_file,
                                                              keep_snp_root=keep_snp_root)
        calculate_M_use_SNP_gene_pair_dummy_by_chunk(snp_gene_weight_matrix,
                                                     mk_score_chunk,
                                                     snp_pass_maf,
                                                     M_file,
                                                     M_5_file,
                                                     )

        chunk_index += 1


def calculate_ldscore_for_base_line(snp_gene_weight_matrix, SNP_gene_pair_dummy, snp_pass_maf, save_dir, chrom,
                                    sample_name, keep_snp_root=None):
    # save baseline ld score
    baseline_mk_score = np.ones((snp_gene_weight_matrix.shape[1], 2))
    baseline_mk_score[-1, 0] = 0  # all_gene
    baseline_mk_score_df = pd.DataFrame(baseline_mk_score, index=snp_gene_weight_matrix.columns,
                                        columns=['all_gene', 'base'])
    ld_score_file = f'{save_dir}/baseline/baseline.{chrom}.l2.ldscore.feather'
    M_file = f'{save_dir}/baseline/baseline.{chrom}.l2.M'
    M_5_file = f'{save_dir}/baseline/baseline.{chrom}.l2.M_5_50'

    calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(snp_gene_weight_matrix,
                                                          baseline_mk_score_df,
                                                          chrom,
                                                          ld_score_file,
                                                          keep_snp_root=keep_snp_root)
    # save baseline M
    calculate_M_use_SNP_gene_pair_dummy_by_chunk(SNP_gene_pair_dummy,
                                                 baseline_mk_score_df,
                                                 snp_pass_maf,
                                                 M_file,
                                                 M_5_file,
                                                 )






def run_generate_ldscore(config: GenerateLDScoreConfig):
    # Load marker score
    mk_score = load_marker_score(config.mkscore_feather_file)

    # Load GTF and get common markers
    gtf_pr, mk_score_common = load_gtf(config.gtf_file, mk_score, window_size=config.window_size)

    def process_chromosome(chrom: int):
        # Process SNPs
        snp_pass_maf = get_snp_pass_maf(config.bfile_root, chrom, maf_min=0.05)

        # Get SNP-Gene dummy pairs
        SNP_gene_pair_dummy = get_snp_gene_dummy(chrom, config.bfile_root, gtf_pr, dummy_na=True)

        # Calculate SNP-Gene weight matrix
        snp_gene_weight_matrix = calculate_SNP_Gene_weight_matrix(SNP_gene_pair_dummy, chrom, config.bfile_root,
                                                                  ld_wind=config.ld_wind, ld_unit=config.ld_unit)

        # Calculate baseline LD score
        calculate_ldscore_for_base_line(snp_gene_weight_matrix, SNP_gene_pair_dummy, snp_pass_maf, config.ldscore_save_dir,
                                        chrom, config.sample_name, keep_snp_root=config.keep_snp_root)

        # Process common genes and calculate LD score
        common_gene_chr = SNP_gene_pair_dummy.columns[:-1]
        calculate_ldscore_use_SNP_Gene_weight_matrix_by_chr(snp_gene_weight_matrix.iloc[:, :-1],
                                                            mk_score_common.loc[common_gene_chr],
                                                            spots_per_chunk=config.spots_per_chunk,
                                                            save_dir=config.ldscore_save_dir,
                                                            chrom=chrom, snp_pass_maf=snp_pass_maf,
                                                            sample_name=config.sample_name,
                                                            keep_snp_root=config.keep_snp_root)

    # Handle 'all' case
    if config.chrom == 'all':
        for chrom in range(1, 23):
            process_chromosome(chrom)
    else:
        process_chromosome(config.chrom)




# %%
if __name__ == '__main__':
    TEST = True
    if TEST:
        # %%
        sample_name = 'Cortex_151507'
        chrom = 'all'
        save_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/Cortex_151507/snp_annotation/test/0101'
        # %%
        gtf_file = '/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf'
        mkscore_feather_file = f'/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/{sample_name}/gene_markers/{sample_name}_rank.feather'
        bfile_root = '/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC'
        window_size = 50000
        keep_snp_root = '/storage/yangjianLab/sharedata/LDSC_resource/hapmap3_snps/hm'
        spots_per_chunk = 10_000

        # %%
        config1 = GenerateLDScoreConfig(
            sample_name=sample_name,
            chrom=chrom,
            ldscore_save_dir=save_dir,
            gtf_file=gtf_file,
            mkscore_feather_file=mkscore_feather_file,
            bfile_root=bfile_root,
            keep_snp_root=keep_snp_root,
            window_size=window_size,
            spots_per_chunk=spots_per_chunk,
        )
        config1 = GenerateLDScoreConfig(
            **{'bfile_root': '/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC',
               'chrom': 21,
               'gtf_file': '/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf',
               'keep_snp_root': '/storage/yangjianLab/sharedata/LDSC_resource/hapmap3_snps/hm',
               'ld_unit': 'CM',
               'ld_wind': 1,
               'ldscore_save_dir': '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/macaque/representative_slices2/T101_macaque1/generate_ldscore',
               'mkscore_feather_file': '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/macaque/representative_slices2/T101_macaque1/latent_to_gene/T101_macaque1_gene_marker_score.feather',
               'sample_name': 'T101_macaque1',
               'spots_per_chunk': 5000,
               'window_size': 50000})

        run_generate_ldscore(config1)
    else:
        parser = argparse.ArgumentParser(description="Configuration for the application.")
        add_generate_ldscore_args(parser)
        args = parser.parse_args()
        config = GenerateLDScoreConfig(**vars(args))
        run_generate_ldscore(config)

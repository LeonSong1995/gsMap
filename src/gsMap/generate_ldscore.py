import argparse
import logging
from pathlib import Path
import zarr
import numpy as np
# %%
import pandas as pd
import pyranges as pr
from scipy.sparse import csr_matrix
from tqdm import trange
import warnings
from gsMap.config import GenerateLDScoreConfig, add_generate_ldscore_args
# %%
from gsMap.generate_r2_matrix import PlinkBEDFileWithR2Cache, getBlockLefts, ID_List_Factory

warnings.filterwarnings("ignore", category=FutureWarning)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:6s} {message}', style='{'))
logger.addHandler(handler)


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
    bim = pd.read_csv(f'{bfile_root}.{chrom}.bim', sep='\t', header=None)
    bim.columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
    #
    # Transform bim to PyRanges
    bim_pr = bim.copy()
    bim_pr.columns = ["Chromosome", "SNP", "CM", "Start", "A1", "A2"]

    bim_pr['End'] = bim_pr['Start'].copy()
    bim_pr['Start'] = bim_pr['Start'] - 1 # Due to bim file is 1-based

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
def filter_snps_by_keep_snp(bim_df, keep_snp_file):
    # Load the keep_snp file and filter the BIM DataFrame
    keep_snp = pd.read_csv(keep_snp_file, header=None)[0].to_list()
    filtered_bim_df = bim_df[bim_df['SNP'].isin(keep_snp)]
    return filtered_bim_df


def get_snp_counts(config):
    snp_counts = {}
    total_snp = 0

    for chrom in range(1, 23):
        bim_df, _ = load_bim(config.bfile_root, chrom)

        if config.keep_snp_root:
            keep_snp_file = f'{config.keep_snp_root}.{chrom}.snp'
            filtered_bim_df = filter_snps_by_keep_snp(bim_df, keep_snp_file)
        else:
            filtered_bim_df = bim_df

        snp_counts[chrom] = filtered_bim_df.shape[0]
        total_snp += snp_counts[chrom]

    snp_counts['total'] = total_snp

    chrom_snp_length_array = np.array([snp_counts[chrom] for chrom in range(1, 23)]).cumsum()

    snp_counts['chrom_snp_start_point'] = [0] + chrom_snp_length_array.tolist()

    return snp_counts




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
def calculate_ldscore_from_annotation(SNP_annotation_df, chrom, bfile_root, ld_wind=1, ld_unit='CM'):
    """
    Calculate the SNP-gene weight matrix.
    """
    # Get the dummy matrix
    # Get the SNP-gene weight matrix
    snp_gene_weight_matrix = get_ldscore(bfile_root, chrom, SNP_annotation_df.values, ld_wind=ld_wind,
                                         ld_unit=ld_unit)
    snp_gene_weight_matrix = snp_gene_weight_matrix.astype(np.float32, copy=False)
    snp_gene_weight_matrix.index = SNP_annotation_df.index
    snp_gene_weight_matrix.columns = SNP_annotation_df.columns
    return snp_gene_weight_matrix


def calculate_ldscore_from_multiple_annotation(SNP_annotation_df_list, chrom, bfile_root, ld_wind=1, ld_unit='CM'):
    SNP_annotation_df = pd.concat(SNP_annotation_df_list, axis=1).astype(np.float32, copy=False)

    snp_gene_weight_matrix = get_ldscore(bfile_root, chrom, SNP_annotation_df.values, ld_wind=ld_wind,
                                         ld_unit=ld_unit)
    snp_gene_weight_matrix = snp_gene_weight_matrix.astype(np.float32, copy=False)
    snp_gene_weight_matrix.index = SNP_annotation_df.index
    snp_gene_weight_matrix.columns = SNP_annotation_df.columns

    # split to each annotation
    snp_annotation_len_list = [len(df.columns) for df in SNP_annotation_df_list]
    snp_gene_weight_matrix_list = []
    start = 0
    for snp_annotation_len in snp_annotation_len_list:
        snp_gene_weight_matrix_list.append(snp_gene_weight_matrix.iloc[:, start:start + snp_annotation_len])
        start += snp_annotation_len
    return snp_gene_weight_matrix_list


# %%
class S_LDSC_Boost:
    def __init__(self, config: GenerateLDScoreConfig):
        self.config = config

        self.mk_score = load_marker_score(config.mkscore_feather_file)

        # Load GTF and get common markers
        self.gtf_pr, self.mk_score_common = load_gtf(config.gtf_annotation_file, self.mk_score,
                                                     window_size=config.gene_window_size)

        # Load enhancer
        if config.enhancer_annotation_file is not None:
            enhancer_df = pr.read_bed(config.enhancer_annotation_file, as_df=True)
            enhancer_df.set_index('Name', inplace=True)
            enhancer_df.index.name = 'gene_name'

            # keep the common genes and add the enhancer score
            avg_mkscore = pd.DataFrame(self.mk_score_common.mean(axis=1), columns=['avg_mkscore'])
            enhancer_df = enhancer_df.join(avg_mkscore, how='inner', on='gene_name', )

            # add distance to TSS
            enhancer_df['TSS'] = self.gtf_pr.df.set_index('gene_name').reindex(enhancer_df.index)['TSS']

            # convert to pyranges
            self.enhancer_pr = pr.PyRanges(enhancer_df.reset_index())

        else:
            self.enhancer_pr = None

        # create tha zarr file
        if config.ldscore_save_format == 'zarr':

            chrom_snp_length_dict = get_snp_counts(config)
            self.chrom_snp_start_point = chrom_snp_length_dict['chrom_snp_start_point']

            zarr_path = Path(config.ldscore_save_dir)/f'{config.sample_name}.ldscore.zarr'
            if not zarr_path.exists():
                self.zarr_file = zarr.open(zarr_path.as_posix(), mode='a', dtype=np.float16,
                                           chunks=config.zarr_chunk_size,
                                           shape=(chrom_snp_length_dict['total'], self.mk_score_common.shape[1]))
                zarr_path.mkdir(parents=True, exist_ok=True)
                # save spot names
                self.zarr_file.attrs['spot_names'] = self.mk_score_common.columns.to_list()
                # save chrom_snp_length_dict
                self.zarr_file.attrs['chrom_snp_start_point'] = self.chrom_snp_start_point
            else:
                self.zarr_file = zarr.open(zarr_path.as_posix(), mode='a')


    def process_chromosome(self, chrom: int):
        self.snp_pass_maf = get_snp_pass_maf(self.config.bfile_root, chrom, maf_min=0.05)

        # Get SNP-Gene dummy pairs
        self.snp_gene_pair_dummy = self.get_snp_gene_dummy(chrom, )

        if self.config.keep_snp_root is not None:
            keep_snp = pd.read_csv(f'{self.config.keep_snp_root}.{chrom}.snp', header=None)[0].to_list()
            self.keep_snp_mask = self.snp_gene_pair_dummy.index.isin(keep_snp)
            # the SNP name of keeped
            self.snp_name = self.snp_gene_pair_dummy.index[self.keep_snp_mask].to_list()
        else:
            self.keep_snp_mask = None
            self.snp_name = self.snp_gene_pair_dummy.index.to_list()

        if self.config.additional_baseline_annotation_dir_path is not None:
            additional_baseline_annotation_dir_path = Path(self.config.additional_baseline_annotation_dir_path)
            additional_baseline_annotation_file_path = additional_baseline_annotation_dir_path / f'baseline.{chrom}.annot.gz'
            assert additional_baseline_annotation_file_path.exists(), f'additional_baseline_annotation_file_path not exists: {additional_baseline_annotation_file_path}'
            additional_baseline_annotation_df = pd.read_csv(additional_baseline_annotation_file_path, sep='\t')
            additional_baseline_annotation_df.set_index('SNP', inplace=True)

            # drop these columns if exists CHR         BP       CM]
            additional_baseline_annotation_df.drop(['CHR', 'BP', 'CM'], axis=1, inplace=True, errors='ignore')

            # reindex, for those SNPs not in additional_baseline_annotation_df, set to 0
            num_of_not_exist_snp = (~self.snp_gene_pair_dummy.index.isin(additional_baseline_annotation_df.index)).sum()
            if num_of_not_exist_snp > 0:
                logger.warning(
                    f'{num_of_not_exist_snp} SNPs not in additional_baseline_annotation_df but in the reference panel, so the additional baseline annotation of these SNP will set to 0')
                additional_baseline_annotation_df = additional_baseline_annotation_df.reindex(
                    self.snp_gene_pair_dummy.index,
                    fill_value=0)
            else:
                additional_baseline_annotation_df = additional_baseline_annotation_df.reindex(
                    self.snp_gene_pair_dummy.index)

            # do this for saving the cpu time, by only calculate r2 once
            self.snp_gene_weight_matrix, additional_baseline_annotation_ldscore = (
                calculate_ldscore_from_multiple_annotation(
                    [self.snp_gene_pair_dummy, additional_baseline_annotation_df],
                    chrom,
                    self.config.bfile_root,
                    ld_wind=self.config.ld_wind,
                    ld_unit=self.config.ld_unit))

            ld_score_file = f'{self.config.ldscore_save_dir}/additional_baseline/baseline.{chrom}.l2.ldscore.feather'
            M_file_path = f'{self.config.ldscore_save_dir}/additional_baseline/baseline.{chrom}.l2.M'
            M_5_file_path = f'{self.config.ldscore_save_dir}/additional_baseline/baseline.{chrom}.l2.M_5_50'

            # save additional baseline annotation ldscore
            self.save_ldscore_to_feather(additional_baseline_annotation_ldscore.values,
                                         column_names=additional_baseline_annotation_ldscore.columns,
                                         save_file_name=ld_score_file,
                                         )

            # caculate the M and save
            save_dir = Path(M_file_path).parent
            save_dir.mkdir(parents=True, exist_ok=True)
            M_chr_chunk = additional_baseline_annotation_df.values.sum(axis=0, keepdims=True)
            M_5_chr_chunk = additional_baseline_annotation_df.loc[self.snp_pass_maf].values.sum(axis=0, keepdims=True)
            np.savetxt(M_file_path, M_chr_chunk, delimiter='\t', )
            np.savetxt(M_5_file_path, M_5_chr_chunk, delimiter='\t', )

        else:
            # Calculate SNP-Gene weight matrix
            self.snp_gene_weight_matrix = calculate_ldscore_from_annotation(self.snp_gene_pair_dummy, chrom,
                                                                            self.config.bfile_root,
                                                                            ld_wind=self.config.ld_wind,
                                                                            ld_unit=self.config.ld_unit)
        # only keep the snp in keep_snp_root
        if self.keep_snp_mask is not None:
            self.snp_gene_weight_matrix = self.snp_gene_weight_matrix[self.keep_snp_mask]

        # convert to sparse
        self.snp_gene_weight_matrix = csr_matrix(self.snp_gene_weight_matrix)

        # calculate baseline ld score
        self.calculate_ldscore_for_base_line(chrom, self.config.sample_name, self.config.ldscore_save_dir)

        # calculate ld score for annotation
        self.calculate_ldscore_use_SNP_Gene_weight_matrix_by_chr(
            self.mk_score_common.loc[self.snp_gene_pair_dummy.columns[:-1]],
            chrom,
            self.config.sample_name,
            self.config.ldscore_save_dir,
        )

    def calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(self,
                                                              mk_score_chunk,
                                                              drop_dummy_na=True,
                                                              ):

        if drop_dummy_na:
            ldscore_chr_chunk = self.snp_gene_weight_matrix[:, :-1] @ mk_score_chunk
        else:
            ldscore_chr_chunk = self.snp_gene_weight_matrix @ mk_score_chunk

        return ldscore_chr_chunk

    def save_ldscore_to_feather(self, ldscore_chr_chunk: np.ndarray, column_names, save_file_name):
        save_dir = Path(save_file_name).parent
        save_dir.mkdir(parents=True, exist_ok=True)

        ldscore_chr_chunk = ldscore_chr_chunk.astype(np.float16, copy=False)
        # avoid overflow of float16, if inf, set to max of float16
        ldscore_chr_chunk[np.isinf(ldscore_chr_chunk)] = np.finfo(np.float16).max
        # ldscore_chr_chunk = ldscore_chr_chunk if self.config.keep_snp_root is None else ldscore_chr_chunk[
        #     self.keep_snp_mask]

        # save for each chunk
        df = pd.DataFrame(ldscore_chr_chunk,
                          index=self.snp_name,
                          columns=column_names,
                          )
        df.index.name = 'SNP'
        df.reset_index().to_feather(save_file_name)

    def save_ldscore_chunk_to_zarr(self, ldscore_chr_chunk: np.ndarray,
                                   chrom:int, start_col_index,
                                   ):
        ldscore_chr_chunk = ldscore_chr_chunk.astype(np.float16, copy=False)
        # avoid overflow of float16, if inf, set to max of float16
        ldscore_chr_chunk[np.isinf(ldscore_chr_chunk)] = np.finfo(np.float16).max

        # save for each chunk
        chrom_snp_start_point = self.chrom_snp_start_point[chrom - 1]
        chrom_snp_end_point = self.chrom_snp_start_point[chrom]

        self.zarr_file[chrom_snp_start_point:chrom_snp_end_point, start_col_index:start_col_index + ldscore_chr_chunk.shape[1]] = ldscore_chr_chunk

    def calculate_M_use_SNP_gene_pair_dummy_by_chunk(self,
                                                     mk_score_chunk,
                                                     M_file_path, M_5_file_path,
                                                     drop_dummy_na=True,
                                                     ):
        '''
        calculate M use SNP_gene_pair_dummy_sumed_along_snp_axis and mk_score_chunk
        '''
        SNP_gene_pair_dummy_sumed_along_snp_axis = self.snp_gene_pair_dummy.values.sum(axis=0, keepdims=True)
        SNP_gene_pair_dummy_sumed_along_snp_axis_pass_maf = self.snp_gene_pair_dummy.loc[self.snp_pass_maf].values.sum(
            axis=0,
            keepdims=True)
        if drop_dummy_na:
            SNP_gene_pair_dummy_sumed_along_snp_axis = SNP_gene_pair_dummy_sumed_along_snp_axis[:, :-1]
            SNP_gene_pair_dummy_sumed_along_snp_axis_pass_maf = SNP_gene_pair_dummy_sumed_along_snp_axis_pass_maf[:,
                                                                :-1]
        save_dir = Path(M_file_path).parent
        save_dir.mkdir(parents=True, exist_ok=True)
        M_chr_chunk = SNP_gene_pair_dummy_sumed_along_snp_axis @ mk_score_chunk
        M_5_chr_chunk = SNP_gene_pair_dummy_sumed_along_snp_axis_pass_maf @ mk_score_chunk
        np.savetxt(M_file_path, M_chr_chunk, delimiter='\t', )
        np.savetxt(M_5_file_path, M_5_chr_chunk, delimiter='\t', )

    def calculate_ldscore_use_SNP_Gene_weight_matrix_by_chr(self, mk_score_common, chrom, sample_name, save_dir):
        """
        Calculate the LD score using the SNP-gene weight matrix.
        :param sample_name:
        """
        # Calculate the LD score
        chunk_index = 1
        for i in trange(0, mk_score_common.shape[1], self.config.spots_per_chunk,
                        desc=f'Calculating LD score by chunk for chr{chrom}'):
            mk_score_chunk = mk_score_common.iloc[:, i:i + self.config.spots_per_chunk]

            ld_score_file = f'{save_dir}/{sample_name}_chunk{chunk_index}/{sample_name}.{chrom}.l2.ldscore.feather'
            M_file = f'{save_dir}/{sample_name}_chunk{chunk_index}/{sample_name}.{chrom}.l2.M'
            M_5_file = f'{save_dir}/{sample_name}_chunk{chunk_index}/{sample_name}.{chrom}.l2.M_5_50'

            ldscore_chr_chunk = self.calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(
                mk_score_chunk,
                drop_dummy_na=True,
            )
            if self.config.ldscore_save_format == 'feather':
                self.save_ldscore_to_feather(ldscore_chr_chunk,
                                             column_names=mk_score_chunk.columns,
                                             save_file_name=ld_score_file,
                                             )
            elif self.config.ldscore_save_format == 'zarr':
                self.save_ldscore_chunk_to_zarr(ldscore_chr_chunk,
                                                chrom=chrom,
                                                start_col_index=i,
                                                )
            else:
                raise ValueError(f'Invalid ldscore_save_format: {self.config.ldscore_save_format}')

            self.calculate_M_use_SNP_gene_pair_dummy_by_chunk(
                mk_score_chunk,
                M_file,
                M_5_file,
                drop_dummy_na=True,
            )

            chunk_index += 1

    def calculate_ldscore_for_base_line(self, chrom, sample_name, save_dir):
        # save baseline ld score
        baseline_mk_score = np.ones((self.snp_gene_pair_dummy.shape[1], 2))
        baseline_mk_score[-1, 0] = 0  # all_gene
        baseline_mk_score_df = pd.DataFrame(baseline_mk_score, index=self.snp_gene_pair_dummy.columns,
                                            columns=['all_gene', 'base'])
        ld_score_file = f'{save_dir}/baseline/baseline.{chrom}.l2.ldscore.feather'
        M_file = f'{save_dir}/baseline/baseline.{chrom}.l2.M'
        M_5_file = f'{save_dir}/baseline/baseline.{chrom}.l2.M_5_50'

        ldscore_chr_chunk = self.calculate_ldscore_use_SNP_Gene_weight_matrix_by_chunk(
            baseline_mk_score_df,
            drop_dummy_na=False,
        )

        self.save_ldscore_to_feather(ldscore_chr_chunk,
                                     column_names=baseline_mk_score_df.columns,
                                     save_file_name=ld_score_file,
                                     )
        # save baseline M
        self.calculate_M_use_SNP_gene_pair_dummy_by_chunk(
            baseline_mk_score_df,
            M_file,
            M_5_file,
            drop_dummy_na=False,
        )

    def get_snp_gene_dummy(self, chrom, ):
        """
        Get the dummy matrix of SNP-gene pairs.
        """
        # Load the bim file
        print("Loading bim data")
        bim, bim_pr = load_bim(self.config.bfile_root, chrom)

        if self.config.gene_window_enhancer_priority in ['gene_window_first', 'enhancer_first']:

            SNP_gene_pair_gtf = self.get_SNP_gene_pair_from_gtf(bim, bim_pr, )
            SNP_gene_pair_enhancer = self.get_SNP_gene_pair_from_enhancer(bim, bim_pr, )
            # total_SNP_gene_pair = SNP_gene_pair_gtf.join(SNP_gene_pair_enhancer, how='outer', lsuffix='_gtf', )

            mask_of_nan_gtf = SNP_gene_pair_gtf.gene_name.isna()
            mask_of_nan_enhancer = SNP_gene_pair_enhancer.gene_name.isna()

            if self.config.gene_window_enhancer_priority == 'gene_window_first':
                SNP_gene_pair = SNP_gene_pair_gtf
                SNP_gene_pair.loc[mask_of_nan_gtf, 'gene_name'] = SNP_gene_pair_enhancer.loc[
                    mask_of_nan_gtf, 'gene_name']
            elif self.config.gene_window_enhancer_priority == 'enhancer_first':
                SNP_gene_pair = SNP_gene_pair_enhancer
                SNP_gene_pair.loc[mask_of_nan_enhancer, 'gene_name'] = SNP_gene_pair_gtf.loc[
                    mask_of_nan_enhancer, 'gene_name']
            else:
                raise ValueError(
                    f'Invalid self.config.gene_window_enhancer_priority: {self.config.gene_window_enhancer_priority}')

        elif self.config.gene_window_enhancer_priority is None:  # use gtf only
            SNP_gene_pair_gtf = self.get_SNP_gene_pair_from_gtf(bim, bim_pr, )
            SNP_gene_pair = SNP_gene_pair_gtf

        elif self.config.gene_window_enhancer_priority == 'enhancer_only':
            SNP_gene_pair_enhancer = self.get_SNP_gene_pair_from_enhancer(bim, bim_pr, )
            SNP_gene_pair = SNP_gene_pair_enhancer
        else:
            raise ValueError('gtf_pr and enhancer_pr cannot be None at the same time')

        # save the SNP_gene_pair to feather
        SNP_gene_pair_save_path=Path(self.config.ldscore_save_dir)/f'SNP_gene_pair/{self.config.sample_name}_chr{chrom}.feather'
        SNP_gene_pair_save_path.parent.mkdir(parents=True, exist_ok=True)
        SNP_gene_pair.to_feather(SNP_gene_pair_save_path)

        # Get the dummy matrix
        SNP_gene_pair_dummy = pd.get_dummies(SNP_gene_pair['gene_name'], dummy_na=True)
        return SNP_gene_pair_dummy

    def get_SNP_gene_pair_from_gtf(self, bim, bim_pr):
        logger.info(
            "Get SNP-gene pair from gtf, if a SNP is in multiple genes, it will be assigned to the most nearby gene (TSS)")
        overlaps_small = Overlaps_gtf_bim(self.gtf_pr, bim_pr)
        # Get the SNP-gene pair
        annot = bim[["CHR", "BP", "SNP", "CM"]]
        SNP_gene_pair = overlaps_small[['SNP', 'gene_name']].set_index('SNP').join(annot.set_index('SNP'), how='right')
        return SNP_gene_pair

    def get_SNP_gene_pair_from_enhancer(self, bim, bim_pr, ):
        logger.info(
            "Get SNP-gene pair from enhancer, if a SNP is in multiple genes, it will be assigned to the gene with highest marker score")
        # Get the SNP-gene pair
        overlaps_small = self.enhancer_pr.join(bim_pr).df
        annot = bim[["CHR", "BP", "SNP", "CM"]]
        if self.config.snp_multiple_enhancer_strategy == 'max_mkscore':
            logger.debug('select the gene with highest marker score')
            overlaps_small = overlaps_small.loc[overlaps_small.groupby('SNP').avg_mkscore.idxmax()]

        elif self.config.snp_multiple_enhancer_strategy == 'nearest_TSS':
            logger.debug('select the gene with nearest TSS')
            overlaps_small['Distance'] = np.abs(overlaps_small['Start_b'] - overlaps_small['TSS'])
            overlaps_small = overlaps_small.loc[overlaps_small.groupby('SNP').Distance.idxmin()]

        SNP_gene_pair = overlaps_small[['SNP', 'gene_name']].set_index('SNP').join(annot.set_index('SNP'), how='right')

        return SNP_gene_pair


def run_generate_ldscore(config: GenerateLDScoreConfig):
    s_ldsc_boost = S_LDSC_Boost(config)
    if config.chrom == 'all':
        for chrom in range(1, 23):
            s_ldsc_boost.process_chromosome(chrom)
    else:
        s_ldsc_boost.process_chromosome(config.chrom)


# %%
if __name__ == '__main__':
    TEST = True
    if TEST:
        # %%
        sample_name = 'Cortex_151507'
        chrom = 'all'
        save_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/Cortex_151507/snp_annotation/test/0101/sparse'
        # %%
        gtf_file = '/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf'
        mkscore_feather_file = f'/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/{sample_name}/gene_markers/{sample_name}_rank.feather'
        bfile_root = '/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC'
        window_size = 50000
        keep_snp_root = '/storage/yangjianLab/sharedata/LDSC_resource/hapmap3_snps/hm'
        spots_per_chunk = 3_000
        enhancer_annotation = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/resource/epigenome/cleaned_data/by_tissue/BRN/ABC_roadmap_merged.bed'
        # %%
        import submitit
        log_folder = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/20240605_without_denoise/mouse_embryo/submitit/log%j"
        executor = submitit.AutoExecutor(folder=log_folder)
        # executor = submitit.DebugExecutor(folder=log_folder)
        executor.update_parameters(slurm_partition='intel-sc3,amd-ep2,amd-ep2-short',
                                   slurm_cpus_per_task=5,
                                   slurm_mem='30G',
                                   slurm_qos='normal',
                                   slurm_time='1-00:00:00',
                                   )
        jobs = []
        for chrom in range(1, 23):
            config = GenerateLDScoreConfig(
                sample_name=sample_name,
                chrom=chrom,
                ldscore_save_dir=save_dir,
                gtf_annotation_file=gtf_file,
                mkscore_feather_file=mkscore_feather_file,
                bfile_root=bfile_root,
                keep_snp_root=keep_snp_root,
                gene_window_size=window_size,
                spots_per_chunk=spots_per_chunk,
                ldscore_save_format='feather',
                # enhancer_annotation_file=enhancer_annotation,
                # gene_window_enhancer_priority='enhancer_first',
                # additional_baseline_annotation_dir_path='/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/resource/ldsc/baseline_v1.2/remove_base'
            )
            executor.update_parameters(job_name=f'ldsc_{sample_name}_{chrom}')
            job = executor.submit(run_generate_ldscore, config)
            jobs.append(job)
        for job in submitit.helpers.as_completed(jobs):
            try:
                print(job.stdout())
            except Exception as e:
                print(e)
        # %%
        # run_generate_ldscore(config)
    else:
        parser = argparse.ArgumentParser(description="Configuration for the application.")
        add_generate_ldscore_args(parser)
        args = parser.parse_args()
        config = GenerateLDScoreConfig(**vars(args))
        run_generate_ldscore(config)

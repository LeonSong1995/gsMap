import argparse
import logging
import os
import time
from pathlib import Path

import cupy as cp
import numpy as np
import pandas as pd
import pyranges as pr
from progress.bar import IncrementalBar

import GPS.config
from GPS.generate_r2_matrix import PlinkBEDFileWithR2Cache, getBlockLefts, ID_List_Factory

pool = cp.cuda.MemoryPool(cp.cuda.malloc_async)
cp.cuda.set_allocator(pool.malloc)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


class Snp_Annotator:
    """
    1. Annotate SNPs based on score of genes.
    2. Add baseline annotations.
    """

    def __init__(self, mk_score_file, gtf_file, bfile_root, annot_root, annot_name, chr=None, base_root=None,
                 window_size=50000, const_max_size=100):
        #
        # marker score
        self.mk_score_file = mk_score_file
        self.mk_score = self.load_marker_score()
        #
        # chunk cells
        # self.const_max_size = const_max_size
        self.n_cells = len(self.mk_score.columns)
        self.max_chunk = const_max_size
        # self.max_chunk = floor(self.n_cells / self.const_max_size)
        #
        # gtf data
        self.gtf_file = gtf_file
        self.window_size = window_size
        self.gtf_pr = self.load_gtf(mk_score=self.mk_score)
        #
        self.bfile_root = bfile_root
        self.annot_root = annot_root
        self.base_root = base_root
        self.chr = chr

        self.data_name = annot_name

    #
    def load_marker_score(self):
        """
        Load marker scores of each cell.
        """
        mk_score = pd.read_feather(self.mk_score_file).set_index('HUMAN_GENE_SYM').rename_axis('gene_name')
        mk_score.insert(0, 'all_gene', 1)
        return mk_score

    #
    def load_gtf(self, mk_score):
        """
        Load the gene annotation file (gtf).
        """
        print("Loading gtf data")
        #
        # Load GTF file
        gtf = pr.read_gtf(self.gtf_file)
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

        gtf_bed.loc[:, 'Start'] = gtf_bed['TSS'] - self.window_size
        gtf_bed.loc[:, 'End'] = gtf_bed['TED'] + self.window_size
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

    #
    def load_baseline(self, chr):
        """
        Load baseline annotations.
        """
        baseline = pd.read_csv(f'{self.base_root}.{chr}.annot.gz', sep='\t')
        baseline.drop(['CHR', 'BP', 'CM'], axis=1, inplace=True)
        return baseline

    # -
    def Load_bim(self, chr):
        """
        Load bim files.
        """
        bim_file = f'{self.bfile_root}.{chr}.bim'
        bim = pd.read_csv(bim_file, sep='\t', header=None)
        bim.columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
        #
        # Transform bim to PyRanges
        bim_pr = bim.copy()
        bim_pr.columns = ["Chromosome", "SNP", "CM", "Start", "A1", "A2"]
        bim_pr['End'] = bim_pr['Start']
        bim_pr = pr.PyRanges(bim_pr)
        bim_pr.Chromosome = f'chr{chr}'
        return bim_pr, bim

    # -
    def Overlaps_gtf_bim(self, bim_pr):
        """
        Find overlaps between gtf and bim file.
        """
        # Select the overlapped regions (SNPs in gene windows)
        overlaps = self.gtf_pr.join(bim_pr)
        overlaps = overlaps.df
        overlaps['Distance'] = np.abs(overlaps['Start_b'] - overlaps['TSS'])
        overlaps_small = overlaps.copy()
        overlaps_small = overlaps_small.loc[overlaps_small.groupby('SNP').Distance.idxmin()]
        return overlaps_small

    # -
    def map_baseline(self, snp_score, baseline, chr):
        """
        Generate the baseline annotations for SNPs.
        """

        header = snp_score.columns[0:6].to_list()

        if baseline is None:
            print(f'Baseline annotations of chr{chr} are not provided, using uniform annotations for genes and SNPs')
            baseline_score = snp_score[header + ['all_gene']].copy()
            baseline_score.loc[:, 'base'] = 1

        else:
            print(f'Mapping baseline annotations of chr{chr}')
            snp_score_baseline = pd.merge(snp_score, baseline, how='left', on='SNP').fillna(0).copy()

            baseline_score = snp_score_baseline[header + ['all_gene'] + baseline.columns.to_list()]
            baseline_score = baseline_score.loc[:, ~baseline_score.columns.duplicated()].copy()

        # Create the folder (for baseline annotation)
        file_base_root = f'{self.annot_root}/baseline'
        if not os.path.exists(file_base_root):
            os.makedirs(file_base_root, mode=0o777, exist_ok=True)

            # Save baseline annotations (in parquet format)
        file_base = f'{file_base_root}/baseline.{chr}.feather'
        baseline_score.to_feather(file_base)

        return 0

    # -
    def annotate_chr(self, chr):
        """
        Annotate SNPs of each chr.
        """
        # Load the baseline file
        baseline = None
        if self.base_root is not None:
            baseline = self.load_baseline(chr)

        # Load the bim file
        bim_pr, bim = self.Load_bim(chr)

        # Find overlapping
        overlaps_small = self.Overlaps_gtf_bim(bim_pr)

        # Do annotations
        all_chunks = int(np.ceil(self.n_cells / self.max_chunk))
        bar = IncrementalBar(f'Mapping the gene marker scores to SNPs in chr{chr}', max=all_chunks)
        bar.check_tty = False

        # Preprocess bim outside the loop as it doesn't change
        anno_template = bim[["CHR", "BP", "SNP", "CM"]]

        for chunk_index, left in enumerate(range(0, self.n_cells, self.max_chunk), start=1):
            right = min(left + self.max_chunk, self.n_cells)
            mk_score_current = self.mk_score.iloc[:, left:right]

            # Process marker scores for SNPs
            anno = anno_template.copy()
            merged_data = overlaps_small[['SNP', 'gene_name', 'TSS']].merge(mk_score_current, on='gene_name',
                                                                            how='left')
            snp_score = pd.merge(anno, merged_data, how='left', on='SNP').fillna(0)
            snp_score = snp_score.rename(columns={'gene_name': 'Gene'})
            snp_score.loc[snp_score.Gene == 0, 'Gene'] = 'None'

            # Process baseline annotations for the first chunk
            if chunk_index == 1:
                self.map_baseline(snp_score, baseline, chr)
                snp_score = snp_score.drop('all_gene', axis=1)

            # Create the folder and save SNP annotations
            file_root = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}'
            os.makedirs(file_root, mode=0o777, exist_ok=True)
            file_anno = f'{file_root}/{self.data_name}.{chr}.feather'
            snp_score.to_feather(file_anno)

            bar.next()

        bar.finish()

        return all_chunks

    #
    def annotate(self):
        """
        Perform SNP annotations for each chromosome.
        """
        if self.chr == None:
            for chr in range(1, 23):
                const_max_size = self.annotate_chr(chr=chr)
        else:
            const_max_size = self.annotate_chr(chr=self.chr)

        return const_max_size


class LDscore_Generator:
    def __init__(self, make_annotation_config: GPS.config.MAKE_ANNOTATION_Conifg, const_max_size):
        self.bfile_root = make_annotation_config.bfile_root
        self.annot_root = Path(make_annotation_config.output_dir) / 'snp_annotation'
        self.const_max_size = const_max_size
        self.data_name = make_annotation_config.sample_info.sample_name
        self.chr = make_annotation_config.chr
        self.ld_wind = make_annotation_config.ld_wind
        self.ld_wind_unit = make_annotation_config.ld_wind_unit
        self.keep_snp = make_annotation_config.keep_snp_root
        self.r2_cache_dir = make_annotation_config.r2_cache_dir
        self.use_gpu = make_annotation_config.use_gpu
        self.config = make_annotation_config
        self.generate_r2_cache = False

        # Set the r2 cache
        if self.r2_cache_dir is None:
            logger.info('No r2 cache directory specified, will not use r2 cache')
            self.chr_r2_cache_dir = None
        else:
            assert self.chr is not None, 'Must specify chr when using r2 cache'
            chr_r2_cache_dir = os.path.join(self.r2_cache_dir, f'chr{self.chr}')
            self.chr_r2_cache_dir = chr_r2_cache_dir
            if not os.path.exists(os.path.join(chr_r2_cache_dir, 'combined_r2_matrix.npz')):
                logger.warning(
                    f'No r2 cache found for chr{self.chr}, will generate r2 cache for this chromosome, first time may take a while')
                os.makedirs(chr_r2_cache_dir, exist_ok=True, mode=0o777, )
                self.generate_r2_cache = True
            else:
                logger.info(f'Found r2 cache for chr{self.chr}, will use r2 cache for this chromosome')

    def compute_ldscore(self):
        """
        Compute LD scores.
        """
        start_time = time.time()
        if self.chr == None:
            for chr in range(1, 23):
                logger.info(f'Computing LD scores for chr{chr}')
                self.compute_ldscore_chr(chr=chr)
                logger.info(f'Finished computing LD scores for chr{chr}')
        else:
            logger.info(f'Computing LD scores for chr{self.chr}')
            self.compute_ldscore_chr(chr=self.chr)
            logger.info(f'Finished computing LD scores for chr{self.chr}')
        end_time = time.time()
        logger.info(f'Finished computing LD scores, time elapsed: {(end_time - start_time) / 60} minutes')

    def compute_ldscore_chunk(self, annot_file, ld_score_file, M_file, M_5_file, geno_array: PlinkBEDFileWithR2Cache,
                              block_left, snp):
        """
        Compute and save LD scores for each chunk
        :param annot_file: Path to the annotation file
        :param ld_score_file: Path to the LD score file
        :param M_file: Path to the M file
        :param M_5_file: Path to the M_5_50 file
        :param geno_array: Genotype array
        :param block_left: Block left
        :param snp: SNP to be kept
        :return: None
        """
        annot_df = pd.read_feather(annot_file)
        n_annot, ma = len(annot_df.columns) - 6, len(annot_df)

        # print("Read {A} annotations for {M} SNPs from {f}".format(f=annot_file, A=n_annot, M=ma))
        annot_matrix = np.array(annot_df.iloc[:, 6:])
        annot_colnames = annot_df.columns[6:]

        # Reset the SNP point
        geno_array.__restart__()

        # Compute annotated LD score
        if self.chr_r2_cache_dir is None:
            lN_df = pd.DataFrame(geno_array.ldScoreVarBlocks(block_left, 50, annot=annot_matrix))
        else:
            lN_df = pd.DataFrame(self.get_ldscore_use_cache(annot_matrix))

        ldscore = pd.concat([annot_df.iloc[:, 0:6], lN_df], axis=1)
        ldscore.columns = annot_df.columns

        # Keep the targeted SNPs
        if not snp is None:
            ldscore = ldscore.loc[ldscore.SNP.isin(snp)]

        # Save the LD score annotations
        ldscore = ldscore.reset_index()
        ldscore.drop(columns=['index'], inplace=True)
        # ldscore.to_feather(ld_score_file)

        # Compute the .M (.M_5_50) file
        M = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix, axis=0))))
        ii = geno_array.maf > 0.05
        M_5_50 = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix[ii, :], axis=0))))

        # Save the sum of score annotations (all and maf > 0.05)
        np.savetxt(M_file, M, delimiter='\t')
        np.savetxt(M_5_file, M_5_50, delimiter='\t')

    def get_ldscore_use_cache(self, annot_matrix, ):
        if self.use_gpu:
            logger.debug('Using GPU to compute LD score')
            annot_matrix = cp.asarray(annot_matrix, dtype=cp.float32)
            for i,r2_matrix_chunk in enumerate(self.r2_matrix_chunk_list):
                r2_matrix_chunk = cp.sparse.csr_matrix(r2_matrix_chunk, dtype=cp.float32)
                lN_chunk = cp.asnumpy(r2_matrix_chunk @ annot_matrix)
                # convert to float16
                lN_chunk = lN_chunk.astype(np.float16)
                if i == 0:
                    lN = lN_chunk
                else:
                    lN = np.concatenate([lN, lN_chunk], axis=0)
        else:
            logger.debug('Using CPU to compute LD score')
            for i,r2_matrix_chunk in enumerate(self.r2_matrix_chunk_list):
                lN_chunk = r2_matrix_chunk @ annot_matrix
                # convert to float16
                lN_chunk = lN_chunk.astype(np.float16)
                if i == 0:
                    lN = lN_chunk
                else:
                    lN = np.concatenate([lN, lN_chunk], axis=0)
        return lN
    def compute_ldscore_chr(self, chr):
        PlinkBIMFile = ID_List_Factory(['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], 1, '.bim', usecols=[0, 1, 2, 3, 4, 5])
        PlinkFAMFile = ID_List_Factory(['IID'], 0, '.fam', usecols=[1])

        bfile = f"{self.bfile_root}.{chr}"
        #
        # Load bim file
        snp_file, snp_obj = bfile + '.bim', PlinkBIMFile
        array_snps = snp_obj(snp_file)
        m = len(array_snps.IDList)
        print(f'Read list of {m} SNPs from {snp_file}')
        #
        # Load fam
        ind_file, ind_obj = bfile + '.fam', PlinkFAMFile
        array_indivs = ind_obj(ind_file)
        n = len(array_indivs.IDList)
        print(f'Read list of {n} individuals from {ind_file}')
        #
        # Load genotype array
        array_file, array_obj = bfile + '.bed', PlinkBEDFileWithR2Cache
        geno_array = array_obj(array_file, n, array_snps, keep_snps=None, keep_indivs=None, mafMin=None)

        # Load the snp to be print
        if not self.keep_snp is None:
            snp = pd.read_csv(f'{self.keep_snp}.{chr}.snp', header=None)[0].to_list()
            num_snp = len(snp)
            print(f'Loading {num_snp} SNPs')
        else:
            snp = None

        # Load the annotations of the baseline
        if self.ld_wind_unit == 'SNP':
            max_dist = self.ld_wind
            coords = np.array(range(geno_array.m))
        elif self.ld_wind_unit == 'BP':
            max_dist = self.ld_wind * 1000
            coords = np.array(array_snps.df['BP'])[geno_array.kept_snps]
        elif self.ld_wind_unit == 'CM':
            max_dist = self.ld_wind
            coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]
        block_left = getBlockLefts(coords, max_dist)
        if self.generate_r2_cache:
            logger.info(f'Generating r2 cache for chr{chr}, this may take a while')
            geno_array.compute_r2_cache(block_left,
                                        Path(self.chr_r2_cache_dir))
            logger.info(f'Finished generating r2 cache for chr{chr}')
        if self.chr_r2_cache_dir is not None:
            logger.info('Loading r2 cache')
            r2_matrix = geno_array.load_combined_r2_matrix(cached_r2_matrix_dir=self.chr_r2_cache_dir)
            self.r2_matrix_chunk_list = [r2_matrix[i:i + self.config.snps_per_chunk, :] for i in
                                    range(0, r2_matrix.shape[0], self.config.snps_per_chunk)]
            logger.info('Finished loading r2 cache')
        # Set the baseline root
        annot_file = f'{self.annot_root}/baseline/baseline.{chr}.feather'
        ld_score_file = f'{self.annot_root}/baseline/baseline.{chr}.l2.ldscore.feather'
        M_file = f'{self.annot_root}/baseline/baseline.{chr}.l2.M'
        M_5_file = f'{self.annot_root}/baseline/baseline.{chr}.l2.M_5_50'

        # Compute annotations of the baseline
        print(f"Computing LD score for baseline annotations of chr{chr}")
        self.compute_ldscore_chunk(annot_file, ld_score_file, M_file, M_5_file, geno_array, block_left, snp)

        # Load annotations of chunks
        bar = IncrementalBar(f"Computing LD scores for spatial data annotations of chr{chr}", max=self.const_max_size)
        bar.check_tty = False
        for chunk_index in range(1, self.const_max_size + 1):
            # Set the file root
            annot_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.feather'
            ld_score_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.l2.ldscore.feather'
            M_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.l2.M'
            M_5_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.l2.M_5_50'

            # Compute annotations of the current chunk
            self.compute_ldscore_chunk(annot_file, ld_score_file, M_file, M_5_file, geno_array, block_left, snp)

            bar.next()

        bar.finish()


def get_sample_info_parser(parser):
    parser.add_argument('--sample_hdf5', default=None, type=str, help='Path to the sample hdf5 file', )
    parser.add_argument('--sample_name', type=str, help='Name of the sample', required=True)
    parser.add_argument('--annotation_layer_name', default=None, type=str, help='Name of the annotation layer')
    parser.add_argument('--is_count', action='store_true', help='Whether the data is count data')


def get_make_annotation_parser(parser):
    parser.add_argument('--gtf_file', default=None, type=str, help='Path to the GTF file', required=True)
    parser.add_argument('--bfile_root', default=None, type=str, help='Bfile root for LD score', required=True)
    parser.add_argument('--baseline_annotation', default=None, type=str, help='Baseline annotation')
    parser.add_argument('--keep_snp_root', default=None, type=str,
                        help='Only keep these SNP file after calculating LD score')
    parser.add_argument('--chr', default=None, type=int, help='Chromosome ID', )
    parser.add_argument('--window_size', default=50000, type=int,
                        help='Window size for SNP annotation')
    parser.add_argument('--cells_per_chunk_size', default=500, type=int,
                        help='Chunk size for number of cells for batch processing')
    parser.add_argument('--ld_wind', default=1, type=float)
    parser.add_argument('--ld_wind_unit', default='CM', type=str, choices=['CM', 'BP', 'SNP'],
                        help='LD window size unit')
    parser.add_argument('--r2_cache_dir', default=None, type=str, help='Directory for r2 cache')
    parser.add_argument('--use_gpu', action='store_true', help='Whether to use GPU to compute LD score')
    parser.add_argument('--snps_per_chunk', default=50_000, type=int,
                        help='Chunk size for number of SNPs for batch processing')
    parser.add_argument('--output_dir', default=None, type=str, help='Output directory', required=True)


# Defin the Container for plink files

def run_make_annotation(make_annotation_config: GPS.config.MAKE_ANNOTATION_Conifg):
    mk_score_file = Path(
        make_annotation_config.output_dir) / f'gene_markers/{make_annotation_config.sample_info.sample_name}_rank.feather'
    # snp_annotate = Snp_Annotator(mk_score_file=mk_score_file,
    #                              gtf_file=make_annotation_config.gtf_file,
    #                              bfile_root=make_annotation_config.bfile_root,
    #                              annot_root=Path(make_annotation_config.output_dir)/'snp_annotation',
    #                              annot_name=make_annotation_config.sample_info.sample_name,
    #                              chr=make_annotation_config.chr,
    #                              base_root=make_annotation_config.baseline_annotation,
    #                              window_size=make_annotation_config.window_size,
    #                              const_max_size=make_annotation_config.chunk_size,
    #                              )
    # const_max_size = snp_annotate.annotate()
    const_max_size = 9
    ldscore_generate = LDscore_Generator(
        make_annotation_config, const_max_size
    )
    ldscore_generate.compute_ldscore()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make_annotations.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    get_sample_info_parser(parser)
    get_make_annotation_parser(parser)
    parser.print_help()

    # Store the Params
    TEST = True
    if TEST:
        name = 'Cortex_151507'
        TASK_ID = 1
        test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021'

        sample_config = GPS.config.ST_SAMPLE_INFO(
            sample_hdf5=f'/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/processed/h5ad/{name}.hdf5',
            sample_name=name,
            annotation_layer_name='spatial',
            is_count=True
        )
        make_annotation_config = GPS.config.MAKE_ANNOTATION_Conifg(
            sample_info=sample_config,
            gtf_file='/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf',
            bfile_root='/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC',
            baseline_annotation=None,
            keep_snp_root='/storage/yangjianLab/sharedata/LDSC_resource/hapmap3_snps/hm',
            chr=TASK_ID,
            window_size=50000,
            cells_per_chunk=500,
            ld_wind=1,
            ld_wind_unit='CM',
            r2_cache_dir='/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/r2_matrix',
            output_dir=f'{test_dir}/{name}/',
            use_gpu=True,
            snps_per_chunk=100_000
        )

    else:
        args = parser.parse_args()
        sample_config = GPS.config.ST_SAMPLE_INFO(
            sample_hdf5=args.sample_hdf5,
            sample_name=args.sample_name,
            annotation_layer_name=args.annotation_layer_name,
            is_count=args.is_count
        )
        make_annotation_config = GPS.config.MAKE_ANNOTATION_Conifg(
            sample_info=sample_config,
            gtf_file=args.gtf_file,
            bfile_root=args.bfile_root,
            baseline_annotation=args.baseline_annotation,
            keep_snp_root=args.keep_snp_root,
            chr=args.chr,
            window_size=args.window_size,
            chunk_size=args.chunk_size,
            ld_wind=args.ld_wind,
            ld_wind_unit=args.ld_wind_unit,
            r2_cache_dir=args.r2_cache_dir,
            output_dir=args.output_dir,
        )

    run_make_annotation(make_annotation_config)

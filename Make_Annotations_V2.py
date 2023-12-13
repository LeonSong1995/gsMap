#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:34:35 2023

@author: songliyang
"""

import pandas as pd
import numpy as np
import pyranges as pr
import sys
import argparse
import pandas as pd
import os
from progress.bar import IncrementalBar
sys.path.append('/storage/yangjianLab/songliyang/SpatialData/spatial_ldsc_v1')
from Build_LD_Score_old import *
from math import floor

class Snp_Annotator:
    """
    1. Annotate SNPs based on score of genes. 
    2. Add baseline annotations. 
    """
    def __init__(self, mk_score_file, gtf_file, bfile_root,annot_root,annot_name, chr = None,base_root = None,
                 window_size = 50000, const_max_size = 100):
        #
        # marker score
        self.mk_score_file = mk_score_file
        self.mk_score = self.load_marker_score()
        #
        # chunk cells
        #self.const_max_size = const_max_size
        self.n_cells = len(self.mk_score.columns)
        self.max_chunk = const_max_size
        #self.max_chunk = floor(self.n_cells / self.const_max_size)
        #
        # gtf data
        self.gtf_file = gtf_file
        self.window_size = window_size
        self.gtf_pr = self.load_gtf(mk_score = self.mk_score)
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
        # mk_score = pd.read_parquet(self.mk_score_file)
        mk_score = pd.read_feather(self.mk_score_file)
        mk_score.index = mk_score.HUMAN_GENE_SYM
        mk_score.drop(columns=['HUMAN_GENE_SYM'], inplace=True)
        mk_score.index.name = 'gene_name'
        mk_score['all_gene'] = 1
        mk_score = mk_score[[mk_score.columns[-1]] + mk_score.columns[0:len(mk_score.columns)-1].to_list()]        
        #
        return mk_score
    
    #
    def load_gtf(self,mk_score):
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
        gtf_bed = gtf[['Chromosome', 'Start', 'End', 'gene_name']].copy()
        gtf_bed.loc[:,'TSS'] = gtf_bed['Start']
        gtf_bed.loc[:,'TED'] = gtf_bed['End']
        
        gtf_bed.loc[:,'Start'] = gtf_bed['TSS'] - self.window_size
        gtf_bed.loc[:,'End'] = gtf_bed['TED'] + self.window_size
        gtf_bed.loc[gtf_bed['Start'] < 0, 'Start'] = 0
        #
        # Transform the GTF to PyRanges
        gtf_pr = pr.PyRanges(gtf_bed)
        return gtf_pr
    
    #
    def load_baseline(self,chr):
        """
        Load baseline annotations.
        """
        baseline = pd.read_csv(f'{self.base_root}.{chr}.annot.gz',sep='\t')
        baseline.drop(['CHR','BP','CM'], axis=1, inplace=True)
        return baseline
    
    #-
    def Load_bim(self,chr):
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
   
    #-
    def Overlaps_gtf_bim(self,bim_pr):
        """
        Find overlaps between gtf and bim file. 
        """
        # Select the overlapped regions (SNPs in gene windows)
        overlaps = self.gtf_pr.join(bim_pr)
        overlaps = overlaps.df
        overlaps['Distance'] = np.abs(overlaps['Start_b'] - overlaps['TSS'])
        # 
        # For SNPs in multiple gene windows, assign them to the nearest genes (snp pos - gene tss)
        overlaps_small = overlaps.copy()
        overlaps_small = overlaps_small.loc[overlaps_small.groupby('SNP').Distance.idxmin()]
        return overlaps_small
    
    #-
    def map_baseline(self,snp_score,baseline,chr):
        """
        Generate the baseline annotations for SNPs.
        """
        
        header = snp_score.columns[0:6].to_list()
        
        if baseline is None:
            print(f'Baseline annotations of chr{chr} are not provided, using uniform annotations for genes and SNPs')
            baseline_score = snp_score[header + ['all_gene']].copy()
            baseline_score.loc[:,'base'] = 1
            
        else:    
            print(f'Mapping baseline annotations of chr{chr}')
            snp_score_baseline = pd.merge(snp_score,baseline,how='left', on='SNP').fillna(0).copy()

            baseline_score = snp_score_baseline[header + ['all_gene'] + baseline.columns.to_list()]
            baseline_score = baseline_score.loc[:,~baseline_score.columns.duplicated()].copy()

        # Create the folder (for baseline annotation)     
        file_base_root = f'{self.annot_root}/baseline'
        if not os.path.exists(file_base_root):
            os.makedirs(file_base_root, mode=0o777, exist_ok=True)    

        # Save baseline annotations (in parquet format)
        file_base = f'{file_base_root}/baseline.{chr}.feather'
        baseline_score.to_feather(file_base) 
        
        return 0
    
    #-   
    def annotate_chr(self,chr):
        """
        Annotate SNPs of each chr. 
        """
        chunk_index = 1
        
        # Load the baseline file
        baseline = None
        if self.base_root is not None:
            baseline = self.load_baseline(chr)
        
        # Load the bim file    
        bim_pr, bim = self.Load_bim(chr)
        
        # Find overlapping
        overlaps_small = self.Overlaps_gtf_bim(bim_pr)
        
        # Do annotations
        all_chunks = len(range(0, self.n_cells, self.max_chunk))
        
        bar = IncrementalBar(f'Mapping the gene marker scores to SNPs in chr{chr}', max = all_chunks)
        bar.check_tty = False 
        for left in range(0, self.n_cells, self.max_chunk):
            
            right = min(left + self.max_chunk, self.n_cells)
            mk_score_current = self.mk_score.iloc[:,left:right]
             
            # Assign marker score for SNPs
            anno = bim.copy()
            anno = anno[["CHR", "BP", "SNP", "CM"]]
            temp = overlaps_small[['SNP', 'gene_name','TSS']].merge(mk_score_current, on='gene_name', how='left')
            snp_score = pd.merge(anno, temp, how='left', on='SNP').fillna(0)
            snp_score = snp_score.rename(columns={'gene_name': 'Gene'})
            snp_score.loc[snp_score.Gene == 0, 'Gene'] = 'None'
            
            # Process the baseline annotations (only for the fisrt chunk, as this will be shared by all chunks)
            if chunk_index == 1:       
                self.map_baseline(snp_score,baseline,chr)
                snp_score = snp_score.drop('all_gene', axis=1)
                
            # Create the folder (for each chunk)
            file_root = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}'
            if not os.path.exists(file_root):
                os.makedirs(file_root, mode=0o777, exist_ok=True)
                
            # Save SNP annotations (in parquet format)
            file_anno = f'{file_root}/{self.data_name}.{chr}.feather'
            snp_score.to_feather(file_anno)
            
            # Update the chunk index
            chunk_index = chunk_index + 1
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
    def __init__(self, bfile_root, annot_root, const_max_size,annot_name,
                 chr=None,ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1, keep_snp=None):
        self.bfile_root = bfile_root
        self.annot_root = annot_root
        self.ld_wind_snps = ld_wind_snps
        self.ld_wind_kb = ld_wind_kb
        self.ld_wind_cm = ld_wind_cm
        self.keep_snp = keep_snp
        self.chr = chr
        
        self.data_name = annot_name
        self.const_max_size = const_max_size
        
    #    
    def compute_ldscore(self):
        """
        Compute LD scores.
        """
        if self.chr == None:
            for chr in range(1, 23):
                self.compute_ldscore_chr(chr=chr)
        else:
            self.compute_ldscore_chr(chr=self.chr) 
            
            
    def compute_ldscore_chunk(self,annot_file,ld_score_file,M_file,M_5_file,geno_array, block_left, snp):
        """
        Compute and save LD scores for each chunk
        """
        annot_df = pd.read_feather(annot_file)
        n_annot, ma = len(annot_df.columns) - 6, len(annot_df)

        # print("Read {A} annotations for {M} SNPs from {f}".format(f=annot_file, A=n_annot, M=ma))
        annot_matrix = np.array(annot_df.iloc[:, 6:])
        annot_colnames = annot_df.columns[6:]

        # Reset the SNP point
        geno_array.__restart__()

        # Compute annotated LD score
        lN_df = pd.DataFrame(geno_array.ldScoreVarBlocks(block_left, 50, annot=annot_matrix))
        ldscore = pd.concat([annot_df.iloc[:, 0:6], lN_df], axis=1)
        ldscore.columns = annot_df.columns
        
        # Keep the targeted SNPs
        if not snp is None:
            ldscore = ldscore.loc[ldscore.SNP.isin(snp)]
        
        # Save the LD score annotations
        ldscore = ldscore.reset_index()
        ldscore.drop(columns=['index'], inplace=True)
        ldscore.to_feather(ld_score_file)
        
        # Compute the .M (.M_5_50) file
        M = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix, axis=0))))
        ii = geno_array.maf > 0.05
        M_5_50 = np.atleast_1d(np.squeeze(np.asarray(np.sum(annot_matrix[ii, :], axis=0))))
        
        # Save the sum of score annotations (all and maf > 0.05)
        np.savetxt(M_file, M, delimiter='\t')
        np.savetxt(M_5_file, M_5_50, delimiter='\t')

    
    def compute_ldscore_chr(self, chr):
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
        array_file, array_obj = bfile + '.bed', PlinkBEDFile
        geno_array = array_obj(array_file, n, array_snps, keep_snps=None, keep_indivs=None, mafMin=None)
        
        # Load the snp to be print
        if not self.keep_snp is None:
            snp = pd.read_csv(f'{self.keep_snp}.{chr}.snp',header=None)[0].to_list()
            num_snp = len(snp)
            print(f'Loading {num_snp} SNPs')
        else:
            snp = None
        
        #Determin LD blocks
        x = np.array((self.ld_wind_snps, self.ld_wind_kb, self.ld_wind_cm), dtype=bool)
        if np.sum(x) != 1:
            raise ValueError('Must specify exactly one ld-wind option')
        # 
        if self.ld_wind_snps:
            max_dist = self.ld_wind_snps
            coords = np.array(range(geno_array.m))
        elif self.ld_wind_kb:
            max_dist = self.ld_wind_kb * 1000
            coords = np.array(array_snps.df['BP'])[geno_array.kept_snps]
        elif self.ld_wind_cm:
            max_dist = self.ld_wind_cm
            coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]
        block_left = getBlockLefts(coords, max_dist)
        
        # Set the baseline root   
        annot_file = f'{self.annot_root}/baseline/baseline.{chr}.feather'
        ld_score_file = f'{self.annot_root}/baseline/baseline.{chr}.l2.ldscore.feather'
        M_file = f'{self.annot_root}/baseline/baseline.{chr}.l2.M'
        M_5_file = f'{self.annot_root}/baseline/baseline.{chr}.l2.M_5_50'

        # Compute annotations of the baseline
        print(f"Computing LD score for baseline annotations of chr{chr}")
        self.compute_ldscore_chunk(annot_file,ld_score_file,M_file,M_5_file,geno_array, block_left, snp)

        # Load annotations of chunks
        bar = IncrementalBar(f"Computing LD scores for spatial data annotations of chr{chr}", max = self.const_max_size)
        bar.check_tty = False 
        for chunk_index in range(1,self.const_max_size+1):
            
            # Set the file root
            annot_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.feather'
            ld_score_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.l2.ldscore.feather'
            M_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.l2.M'
            M_5_file = f'{self.annot_root}/{self.data_name}_chunk{chunk_index}/{self.data_name}.{chr}.l2.M_5_50'
            
            # Compute annotations of the current chunk
            self.compute_ldscore_chunk(annot_file,ld_score_file,M_file,M_5_file,geno_array, block_left, snp)
            
            bar.next()
            
        bar.finish()


parser = argparse.ArgumentParser()
parser.add_argument('--mk_score_file', default=None, type=str)
parser.add_argument('--gtf_file', default=None, type=str)
parser.add_argument('--bfile_root', default=None, type=str)
parser.add_argument('--annot_root', default=None, type=str)
parser.add_argument('--annot_name', default=None, type=str)
parser.add_argument('--base_root', default=None, type=str)
parser.add_argument('--keep_snp', default=None, type=str)

parser.add_argument('--chr', default=None, type=int)
parser.add_argument('--window_size', default=50000, type=int)
parser.add_argument('--const_max_size', default=100, type=int)
parser.add_argument('--ld_wind_snps', default=None, type=float)
parser.add_argument('--ld_wind_kb', default=None, type=float)
parser.add_argument('--ld_wind_cm', default=None, type=float)


# Defin the Container for plink files   
PlinkBIMFile = ID_List_Factory(['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], 1, '.bim', usecols=[0, 1, 2, 3, 4, 5])
PlinkFAMFile = ID_List_Factory(['IID'], 0, '.fam', usecols=[1])
FilterFile = ID_List_Factory(['ID'], 0, None, usecols=[0])


if __name__ == '__main__':

    # Store the Params
    args = parser.parse_args()

    # Mapping gene score to SNPs
    snp_annotate = Snp_Annotator(mk_score_file=args.mk_score_file, gtf_file=args.gtf_file,
                                 bfile_root=args.bfile_root, annot_root=args.annot_root, 
                                 base_root=args.base_root,annot_name=args.annot_name,
                                 window_size=args.window_size, chr=args.chr, const_max_size=args.const_max_size)
    const_max_size = snp_annotate.annotate()


    # Generate LD scores annotations
    ldscore_generate = LDscore_Generator(bfile_root=args.bfile_root, annot_root=args.annot_root, const_max_size=const_max_size,
                                         annot_name=args.annot_name, chr=args.chr, ld_wind_snps=args.ld_wind_snps, ld_wind_kb=args.ld_wind_kb, 
                                         ld_wind_cm=args.ld_wind_cm, keep_snp = args.keep_snp)
    ldscore_generate.compute_ldscore()


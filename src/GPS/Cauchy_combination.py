#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 11:34:35 2023

@author: songliyang
"""

import pandas as pd
import numpy as np
import scipy as sp
import scanpy as sc
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--ldsc_path', default=None, type=str)
parser.add_argument('--ldsc_name', default=None, type=str)

parser.add_argument('--spe_path', default=None, type=str)
parser.add_argument('--spe_name', default=None, type=str)

parser.add_argument('--meta', default=None, type=str)
parser.add_argument('--slide', default=None, type=str)

parser.add_argument('--annotation', default=None, type=str)


# The fun of cauchy combination 
def acat_test(pvalues,weights=None):
    '''acat_test()
    Aggregated Cauchy Assocaition Test
    A p-value combination method using the Cauchy distribution.
    
    Inspired by: https://github.com/yaowuliu/ACAT/blob/master/R/ACAT.R
    Inputs:
        pvalues: <list or numpy array>
            The p-values you want to combine.
        weights: <list or numpy array>, default=None
            The weights for each of the p-values. If None, equal weights are used.
    
    Returns:
        pval: <float>
            The ACAT combined p-value.
    '''
    if any(np.isnan(pvalues)):
        raise Exception("Cannot have NAs in the p-values.")
    if any([(i>1)|(i<0) for i in pvalues]):
        raise Exception("P-values must be between 0 and 1.")
    if any([i==1 for i in pvalues])&any([i==0 for i in pvalues]):
        raise Exception("Cannot have both 0 and 1 p-values.")
    if any([i==0 for i in pvalues]):
        print("Warn: p-values are exactly 0.")
        return 0
    if any([i==1 for i in pvalues]):
        print("Warn: p-values are exactly 1.")
        return 1
    if weights==None:
        weights = [1/len(pvalues) for i in pvalues]
    elif len(weights)!=len(pvalues):
        raise Exception("Length of weights and p-values differs.")
    elif any([i<0 for i in weights]):
        raise Exception("All weights must be positive.")
    else:
        weights = [i/len(weights) for i in weights]
    
    pvalues = np.array(pvalues)
    weights = np.array(weights)
    
    if any([i<1e-16 for i in pvalues])==False:
        cct_stat = sum(weights*np.tan((0.5-pvalues)*np.pi))
    else:
        is_small = [i<(1e-16) for i in pvalues]
        is_large = [i>=(1e-16) for i in pvalues]
        cct_stat = sum((weights[is_small]/pvalues[is_small])/np.pi)
        cct_stat += sum(weights[is_large]*np.tan((0.5-pvalues[is_large])*np.pi))
    
    if cct_stat>1e15:
        pval = (1/cct_stat)/np.pi
    else:
        pval = 1 - sp.stats.cauchy.cdf(cct_stat)
    
    return pval





if __name__ == '__main__':

    args = parser.parse_args()

    # Load the ldsc results
    print(f'------Loading LDSC results of {args.ldsc_name}...')
    ldsc = pd.read_csv(f'{args.ldsc_path}/{args.ldsc_name}')
    ldsc.spot = ldsc.spot.astype(str).replace('\.0', '', regex=True)
    ldsc.index = ldsc.spot
    
    
    if args.meta is None:
        # Load the spatial data
        print(f'------Loading ST data of {args.spe_name}...')
        spe = sc.read_h5ad(f'{args.spe_path}/{args.spe_name}')
        
        common_cell = np.intersect1d(ldsc.index,spe.obs_names)
        spe = spe[common_cell,]
        ldsc = ldsc.loc[common_cell]
        
        # Add the annotation
        ldsc['annotation'] = spe.obs.loc[ldsc.spot][args.annotation].to_list()
    
    elif args.meta is not None:
        # Or Load the additional annotation (just for the macaque data at this stage: 2023Nov25)
        print(f'------Loading additional annotation...')
        meta = pd.read_csv(args.meta,index_col=0)
        meta = meta.loc[meta.slide==args.slide]
        meta.index = meta.cell_id.astype(str).replace('\.0', '', regex=True)
        
        common_cell = np.intersect1d(ldsc.index,meta.index)
        meta = meta.loc[common_cell]
        ldsc = ldsc.loc[common_cell]
        
        # Add the annotation
        ldsc['annotation'] = meta.loc[ldsc.spot][args.annotation].to_list()

    # Perform the Cauchy combination based on the given annotations
    p_cauchy = []
    p_median = []
    for ct in np.unique(ldsc.annotation):
        p_temp  = ldsc.loc[ldsc['annotation'] == ct,'p']
        p_cauchy_temp = acat_test(p_temp)
        p_median_temp = np.median(p_temp)
        
        p_cauchy.append(p_cauchy_temp)
        p_median.append(p_median_temp)

#     p_tissue = pd.DataFrame(p_cauchy,p_median,np.unique(ldsc.annotation))
    data = {'p_cauchy': p_cauchy, 'p_median': p_median, 'annotation': np.unique(ldsc.annotation)}
    p_tissue = pd.DataFrame(data)
    p_tissue.columns = ['p_cauchy','p_median','annotation']

    # Save the results
    name = args.ldsc_name.split('.gz')[0]
    p_tissue.to_csv(f'{args.ldsc_path}/{name}.{args.annotation}.Cauchy.gz',compression='gzip',index=False) 


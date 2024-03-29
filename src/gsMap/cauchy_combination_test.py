import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp

from gsMap.config import CauchyCombinationConfig, add_Cauchy_combination_args

# The fun of cauchy combination
def acat_test(pvalues, weights=None):
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
    if any([(i > 1) | (i < 0) for i in pvalues]):
        raise Exception("P-values must be between 0 and 1.")
    if any([i == 1 for i in pvalues]) & any([i == 0 for i in pvalues]):
        raise Exception("Cannot have both 0 and 1 p-values.")
    if any([i == 0 for i in pvalues]):
        print("Warn: p-values are exactly 0.")
        return 0
    if any([i == 1 for i in pvalues]):
        print("Warn: p-values are exactly 1.")
        return 1
    if weights == None:
        weights = [1 / len(pvalues) for i in pvalues]
    elif len(weights) != len(pvalues):
        raise Exception("Length of weights and p-values differs.")
    elif any([i < 0 for i in weights]):
        raise Exception("All weights must be positive.")
    else:
        weights = [i / len(weights) for i in weights]

    pvalues = np.array(pvalues)
    weights = np.array(weights)

    if any([i < 1e-16 for i in pvalues]) == False:
        cct_stat = sum(weights * np.tan((0.5 - pvalues) * np.pi))
    else:
        is_small = [i < (1e-16) for i in pvalues]
        is_large = [i >= (1e-16) for i in pvalues]
        cct_stat = sum((weights[is_small] / pvalues[is_small]) / np.pi)
        cct_stat += sum(weights[is_large] * np.tan((0.5 - pvalues[is_large]) * np.pi))

    if cct_stat > 1e15:
        pval = (1 / cct_stat) / np.pi
    else:
        pval = 1 - sp.stats.cauchy.cdf(cct_stat)

    return pval


def run_Cauchy_combination(config:CauchyCombinationConfig):
    # Load the ldsc results
    print(f'------Loading LDSC results of {config.input_ldsc_dir}...')
    ldsc_input_file= Path(config.input_ldsc_dir)/f'{config.sample_name}_{config.trait_name}.csv.gz'
    ldsc = pd.read_csv(ldsc_input_file, compression='gzip')
    ldsc.spot = ldsc.spot.astype(str).replace('\.0', '', regex=True)
    ldsc.index = ldsc.spot
    if config.meta is None:
        # Load the spatial data
        print(f'------Loading ST data of {config.input_hdf5_path}...')
        spe = sc.read_h5ad(f'{config.input_hdf5_path}')

        common_cell = np.intersect1d(ldsc.index, spe.obs_names)
        spe = spe[common_cell,]
        ldsc = ldsc.loc[common_cell]

        # Add the annotation
        ldsc['annotation'] = spe.obs.loc[ldsc.spot][config.annotation].to_list()

    elif config.meta is not None:
        # Or Load the additional annotation (just for the macaque data at this stage: 2023Nov25)
        print(f'------Loading additional annotation...')
        meta = pd.read_csv(config.meta, index_col=0)
        meta = meta.loc[meta.slide == config.slide]
        meta.index = meta.cell_id.astype(str).replace('\.0', '', regex=True)

        common_cell = np.intersect1d(ldsc.index, meta.index)
        meta = meta.loc[common_cell]
        ldsc = ldsc.loc[common_cell]

        # Add the annotation
        ldsc['annotation'] = meta.loc[ldsc.spot][config.annotation].to_list()
    # Perform the Cauchy combination based on the given annotations
    p_cauchy = []
    p_median = []
    for ct in np.unique(ldsc.annotation):
        p_temp = ldsc.loc[ldsc['annotation'] == ct, 'p']
        
        # The Cauchy test is sensitive to very small p-values, so extreme outliers should be considered for removal...
        # to enhance robustness, particularly in cases where spot annotations may be incorrect. 
        # p_cauchy_temp = acat_test(p_temp[p_temp != np.min(p_temp)])
        p_temp_log = -np.log10(p_temp)
        median_log = np.median(p_temp_log)
        IQR_log = np.percentile(p_temp_log, 75) - np.percentile(p_temp_log, 25)
        
        p_use = p_temp[p_temp_log < median_log + 3*IQR_log]
        n_remove = len(p_temp) - len(p_use)
        
        # Outlier: -log10(p) < median + 3IQR && len(outlier set) < 20
        if (0 < n_remove < 20):
            print(f'Remove {n_remove}/{len(p_temp)} outliers (median + 3IQR) for {ct}.')
            p_cauchy_temp = acat_test(p_use)
        else:
             p_cauchy_temp = acat_test(p_temp)
                
        p_median_temp = np.median(p_temp)

        p_cauchy.append(p_cauchy_temp)
        p_median.append(p_median_temp)
    #     p_tissue = pd.DataFrame(p_cauchy,p_median,np.unique(ldsc.annotation))
    data = {'p_cauchy': p_cauchy, 'p_median': p_median, 'annotation': np.unique(ldsc.annotation)}
    p_tissue = pd.DataFrame(data)
    p_tissue.columns = ['p_cauchy', 'p_median', 'annotation']
    # Save the results
    output_dir = Path(config.output_cauchy_dir)
    output_dir.mkdir(parents=True, exist_ok=True, mode=0o755)
    output_file = output_dir / f'{config.sample_name}_{config.trait_name}.Cauchy.csv.gz'
    p_tissue.to_csv(
        output_file,
        compression='gzip',
        index=False,
    )


if __name__ == '__main__':
    TEST = True
    if TEST:
        test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_gsMap/data/gsMap_test/Nature_Neuroscience_2021'
        name = 'Cortex_151507'

        config = CauchyCombinationConfig(
            input_hdf5_path= f'{test_dir}/{name}/hdf5/{name}_add_latent.h5ad',
            input_ldsc_dir=
            f'/storage/yangjianLab/chenwenhao/projects/202312_gsMap/data/gsMap_test/Nature_Neuroscience_2021/snake_workdir/Cortex_151507/ldsc/',
            sample_name=name,
            annotation='layer_guess',
            output_cauchy_dir='/storage/yangjianLab/chenwenhao/projects/202312_gsMap/data/gsMap_test/Nature_Neuroscience_2021/snake_workdir/Cortex_151507/cauchy/',
            trait_name='adult1_adult2_onset_asthma',
        )
    else:

        parser = argparse.ArgumentParser(description="Run Cauchy Combination Analysis")
        add_Cauchy_combination_args(parser)
        args = parser.parse_args()
        config = CauchyCombinationConfig(**vars(args))
    run_Cauchy_combination(config)

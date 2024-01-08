import argparse
import logging
import multiprocessing
import os
from pathlib import Path

from scipy.stats import norm
from tqdm import trange
from tqdm.contrib.concurrent import process_map

import GPS.jackknife as jk
from GPS.config import add_spatial_ldsc_args, SpatialLDSCConfig
from GPS.regression_read import _read_sumstats, _read_w_ld, _read_ref_ld_v2, _read_M_v2, _check_variance_v2

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd


# %%
# Helper Functions
def warn_length(sumstats):
    if len(sumstats) < 200000:
        print('WARNING: number of SNPs less than 200k; this is almost always bad.')


def _coef_new(jknife):
    # return coef[0], coef_se[0], z[0]]
    est_ = jknife.est[0, 0]
    se_ = jknife.jknife_se[0, 0]
    return est_, se_


def append_intercept(x):
    n_row = x.shape[0]
    intercept = np.ones((n_row, 1))
    x_new = np.concatenate((x, intercept), axis=1)
    return x_new


def apply_weights(x, w):
    if np.any(w <= 0):
        raise ValueError('Weights must be > 0')
    (n, p) = x.shape
    if w.shape != (n, 1):
        raise ValueError('w must have shape (n, 1).')
    w = w / float(np.sum(w))
    x_new = np.multiply(x, w)
    return x_new


def filter_sumstats_by_chisq(sumstats, chisq_max):
    before_len = len(sumstats)
    if chisq_max is None:
        chisq_max = max(0.001 * sumstats.N.max(), 80)
        logger.info(f'No chi^2 threshold provided, using {chisq_max} as default')
    sumstats['chisq'] = sumstats.Z ** 2
    sumstats = sumstats[sumstats.chisq < chisq_max]
    after_len = len(sumstats)
    if after_len < before_len:
        logger.info(f'Removed {before_len - after_len} SNPs with chi^2 > {chisq_max} ({after_len} SNPs remain)')
    else:
        logger.info(f'No SNPs removed with chi^2 > {chisq_max} ({after_len} SNPs remain)')
    return sumstats


# Core Functionalities
def get_weight(sumstats, ref_ld_baseline, ref_ld_spatial, w_ld, M_annot, intercept=1):
    M_tot = np.sum(M_annot)
    x_tot = ref_ld_baseline.sum(axis=1).values + ref_ld_spatial.sum(axis=1).values
    tot_agg = aggregate(sumstats.chisq, x_tot, sumstats.N, M_tot, intercept)
    initial_w = weights(x_tot, w_ld.LD_weights.values, sumstats.N.values, M_tot, tot_agg, intercept)
    initial_w = np.sqrt(initial_w)
    return initial_w


def aggregate(y, x, N, M, intercept=1):
    num = M * (np.mean(y) - intercept)
    denom = np.mean(np.multiply(x, N))
    return num / denom


def weights(ld, w_ld, N, M, hsq, intercept=1):
    M = float(M)
    hsq = np.clip(hsq, 0.0, 1.0)
    ld = np.maximum(ld, 1.0)
    w_ld = np.maximum(w_ld, 1.0)
    c = hsq * N / M
    het_w = 1.0 / (2 * np.square(intercept + np.multiply(c, ld)))
    oc_w = 1.0 / w_ld
    w = np.multiply(het_w, oc_w)
    return w


def jackknife_for_processmap(spot_id):
    x_focal = np.concatenate((np.reshape(ref_ld_spatial[:, spot_id], (-1, 1)),
                              baseline_annotation), axis=1)
    # random shuffle x and y
    jknife = jk.LstsqJackknifeFast(x_focal, y, n_blocks)
    return _coef_new(jknife)



def run_spatial_ldsc(config: SpatialLDSCConfig):

    global ref_ld_spatial, baseline_annotation, y, n_blocks
    n_blocks = config.n_blocks
    trait_name = config.trait_name
    print(f'------Running Spatial LDSC for {trait_name}...')
    data_name = config.sample_name
    num_cpus = min(multiprocessing.cpu_count(), config.num_processes)
    # Load the gwas summary statistics
    sumstats = _read_sumstats(fh=config.h2, alleles=False, dropna=False)
    sumstats.set_index('SNP', inplace=True)
    sumstats = sumstats.astype(np.float32)
    sumstats = filter_sumstats_by_chisq(sumstats, config.chisq_max)

    # Load the regression weights
    w_ld = _read_w_ld(config.w_file)
    w_ld_cname = w_ld.columns[1]
    w_ld.set_index('SNP', inplace=True)
    # Load the baseline annotations
    ld_file_baseline = f'{config.ldscore_input_dir}/baseline/baseline.'
    ref_ld_baseline = _read_ref_ld_v2(ld_file_baseline)
    common_snp = ref_ld_baseline.index.intersection(w_ld.index).intersection(sumstats.index)
    logger.info(f'Find {len(common_snp)} common SNPs between GWAS and baseline annotations')

    filter_by_common_snp = lambda df: df.loc[common_snp]
    sumstats = filter_by_common_snp(sumstats)
    ref_ld_baseline = filter_by_common_snp(ref_ld_baseline)
    w_ld = filter_by_common_snp(w_ld)
    n_annot_baseline = len(ref_ld_baseline.columns)
    M_annot_baseline = _read_M_v2(ld_file_baseline, n_annot_baseline, config.not_M_5_50)
    # Detect chunk files
    all_file = os.listdir(config.ldscore_input_dir)
    if config.all_chunk is None:
        all_chunk = sum('chunk' in name for name in all_file)
        print(f'\t')
        print(f'Find {all_chunk} chunked files')
    else:
        all_chunk = config.all_chunk
        print(f'\t')
        print(f'Input {all_chunk} chunked files')
    # Process each chunk
    out_all = pd.DataFrame()
    for chunk_index in range(1, all_chunk + 1):
        print(f'------Processing chunk-{chunk_index}')
        # Load the spatial ldscore annotations
        ld_file_spatial = f'{config.ldscore_input_dir}/{data_name}_chunk{chunk_index}/{data_name}.'
        ref_ld_spatial = _read_ref_ld_v2(ld_file_spatial)
        ref_ld_spatial = filter_by_common_snp(ref_ld_spatial).astype(np.float32)
        ref_ld_spatial_cnames = ref_ld_spatial.columns

        n_annot_spatial = len(ref_ld_spatial.columns)
        M_annot_spatial = _read_M_v2(ld_file_spatial, n_annot_spatial, config.not_M_5_50)

        # Merge the spatial annotations and baseline annotations
        n_annot = n_annot_baseline + n_annot_spatial
        M_annot = np.concatenate((M_annot_baseline, M_annot_spatial), axis=1)

        # Check the variance of the design matrix
        # M_annot, ref_ld_spatial = _check_variance_v2(M_annot, ref_ld_spatial)

        initial_w = get_weight(sumstats, ref_ld_baseline, ref_ld_spatial, w_ld, M_annot).astype(np.float32).reshape(
            (-1, 1))

        assert np.any(initial_w > 0), 'Weights must be > 0'

        baseline_annotation = ref_ld_baseline.copy()
        baseline_annotation = baseline_annotation * sumstats.N.values.reshape((-1, 1)) / sumstats.N.mean()
        # append intercept
        baseline_annotation = append_intercept(baseline_annotation)

        # apply weight
        initial_w = initial_w / np.sum(initial_w)
        baseline_annotation *= initial_w
        ref_ld_spatial *= initial_w
        ref_ld_spatial = ref_ld_spatial.to_numpy(dtype=np.float32)

        y = sumstats.chisq.to_numpy(dtype=np.float32).reshape((-1, 1))
        y *= initial_w

        chunk_size = ref_ld_spatial.shape[1]
        out_chunk = process_map(jackknife_for_processmap, range(chunk_size),
                                max_workers=config.num_processes,
                                chunksize=10,
                                desc=f'Running LDSC for {chunk_size} cells in chunk-{chunk_index}.'
                                )
        out_chunk = pd.DataFrame.from_records(out_chunk, columns=['beta', 'se', ],
                                              index=ref_ld_spatial_cnames)
        out_chunk['z'] = out_chunk.beta / out_chunk.se
        out_chunk['p'] = norm.sf(out_chunk['z'])
        # Concat results
        out_all = pd.concat([out_all, out_chunk], axis=0)
    # Save the results
    out_dir = Path(config.ldsc_save_dir)
    out_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
    out_file_name = out_dir / f'{data_name}_{trait_name}.new.csv.gz'
    out_all['spot'] = out_all.index
    out_all = out_all[['spot', 'beta', 'se', 'z', 'p']]
    out_all.to_csv(out_file_name, compression='gzip', index=False)
    print(f'Done! Results saved to {out_file_name}')


if __name__ == '__main__':
    # Main function of analysis
    parser = argparse.ArgumentParser(
        description="Run Spatial LD Score Regression (LDSC) analysis for GWAS and spatial transcriptomic data."
    )
    parser = add_spatial_ldsc_args(parser)
    TEST = True
    if TEST:
        gwas_root = "/storage/yangjianLab/songliyang/GWAS_trait/LDSC"
        gwas_trait = "/storage/yangjianLab/songliyang/GWAS_trait/GWAS_Public_Use_MaxPower.csv"
        root = "/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/processed/h5ad"

        name = 'Cortex_151507'
        spe_name = name
        # ld_pth = f"/storage/yangjianLab/songliyang/SpatialData/Data/Brain/Human/Nature_Neuroscience_2021/annotation/{spe_name}/snp_annotation"
        ld_pth = f"/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/snake_workdir/{name}/generate_ldscore"
        out_pth = f"/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/snake_workdir/{name}/ldsc"
        gwas_file = "ADULT1_ADULT2_ONSET_ASTHMA"
        # Prepare the arguments list using f-strings
        args_list = [
            "--h2", f"{gwas_root}/{gwas_file}.sumstats.gz",
            "--w_file", "/storage/yangjianLab/sharedata/LDSC_resource/LDSC_SEG_ldscores/weights_hm3_no_hla/weights.",
            "--sample_name", spe_name,
            "--num_processes", '4',
            "--ldscore_input_dir", ld_pth,
            "--ldsc_save_dir", out_pth,
            '--trait_name', 'adult1_adult2_onset_asthma'
        ]
        args = parser.parse_args(args_list)
    else:
        args = parser.parse_args()
    config = SpatialLDSCConfig(**vars(args))
    run_spatial_ldsc(config)
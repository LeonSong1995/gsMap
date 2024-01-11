import argparse
import logging
import multiprocessing
from collections import defaultdict
from pathlib import Path
import os
from scipy.stats import norm
from tqdm.contrib.concurrent import process_map

import GPS.jackknife as jk
from GPS.config import add_spatial_ldsc_args, SpatialLDSCConfig
from GPS.regression_read import _read_sumstats, _read_w_ld, _read_ref_ld_v2, _read_M_v2

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
    est_ = jknife.est[0, 0] / Nbar
    se_ = jknife.jknife_se[0, 0] / Nbar
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
def get_weight(sumstats, ref_ld_baseline, spatial_annotation, w_ld, M_annot, intercept=1):
    M_tot = np.sum(M_annot)
    x_tot = ref_ld_baseline.sum(axis=1).values + spatial_annotation.sum(axis=1).values
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
    x_focal = np.concatenate((np.reshape(spatial_annotation[:, spot_id], (-1, 1)),
                              baseline_annotation), axis=1)
    # random shuffle x and y
    jknife = jk.LstsqJackknifeFast(x_focal, y, n_blocks)
    return _coef_new(jknife)


# Updated function
def get_weight_optimized(sumstats, x_tot_precomputed, M_tot, w_ld, intercept=1):
    tot_agg = aggregate(sumstats.chisq, x_tot_precomputed, sumstats.N, M_tot, intercept)
    initial_w = weights(x_tot_precomputed, w_ld.LD_weights.values, sumstats.N.values, M_tot, tot_agg, intercept)
    initial_w = np.sqrt(initial_w)
    return initial_w


def _preprocess_sumstats(trait_name, sumstat_file_path, baseline_and_w_ld_common_snp: pd.Index, chisq_max=None):
    # Load the gwas summary statistics
    sumstats = _read_sumstats(fh=sumstat_file_path, alleles=False, dropna=False)
    sumstats.set_index('SNP', inplace=True)
    sumstats = sumstats.astype(np.float32)
    sumstats = filter_sumstats_by_chisq(sumstats, chisq_max)

    # NB: The intersection order is essential for keeping the same order of SNPs by its BP location
    common_snp = baseline_and_w_ld_common_snp.intersection(sumstats.index)
    if len(common_snp) < 200000:
        logger.warning(f'WARNING: number of SNPs less than 200k; for {trait_name} this is almost always bad.')

    sumstats = sumstats.loc[common_snp]
    return sumstats


def _get_sumstats_from_sumstats_dict(sumstats_config_dict: dict, baseline_and_w_ld_common_snp: pd.Index,
                                     chisq_max=None):
    # first validate if all sumstats file exists
    logger.info('Validating sumstats files...')
    for trait_name, sumstat_file_path in sumstats_config_dict.items():
        if not os.path.exists(sumstat_file_path):
            raise FileNotFoundError(f'{sumstat_file_path} not found')
    # then load all sumstats
    sumstats_cleaned_dict = {}
    for trait_name, sumstat_file_path in sumstats_config_dict.items():
        sumstats_cleaned_dict[trait_name] = _preprocess_sumstats(trait_name, sumstat_file_path,
                                                                 baseline_and_w_ld_common_snp, chisq_max)
    logger.info('cleaned sumstats loaded')
    return sumstats_cleaned_dict


def run_spatial_ldsc(config: SpatialLDSCConfig):
    global spatial_annotation, baseline_annotation, y, n_blocks, Nbar
    # config
    n_blocks = config.n_blocks
    sample_name = config.sample_name
    num_cpus = min(multiprocessing.cpu_count(), config.num_processes)

    print(f'------Running Spatial LDSC for {sample_name}...')
    # Load the regression weights
    w_ld = _read_w_ld(config.w_file)
    w_ld_cname = w_ld.columns[1]
    w_ld.set_index('SNP', inplace=True)

    # Load the baseline annotations
    ld_file_baseline = f'{config.ldscore_input_dir}/baseline/baseline.'
    ref_ld_baseline = _read_ref_ld_v2(ld_file_baseline)
    n_annot_baseline = len(ref_ld_baseline.columns)
    M_annot_baseline = _read_M_v2(ld_file_baseline, n_annot_baseline, config.not_M_5_50)

    # common snp between baseline and w_ld
    baseline_and_w_ld_common_snp = ref_ld_baseline.index.intersection(w_ld.index)
    if len(baseline_and_w_ld_common_snp) < 200000:
        logger.warning(f'WARNING: number of SNPs less than 200k; for {sample_name} this is almost always bad.')
    ref_ld_baseline = ref_ld_baseline.loc[baseline_and_w_ld_common_snp]
    w_ld = w_ld.loc[baseline_and_w_ld_common_snp]

    # Clean the sumstats
    sumstats_cleaned_dict = _get_sumstats_from_sumstats_dict(config.sumstats_config_dict, baseline_and_w_ld_common_snp,
                                                             chisq_max=config.chisq_max)

    # Detect avalable chunk files
    all_file = os.listdir(config.ldscore_input_dir)
    if config.all_chunk is None:
        all_chunk = sum('chunk' in name for name in all_file)
        print(f'\t')
        print(f'Find {all_chunk} chunked files')
    else:
        all_chunk = config.all_chunk
        print(f'using {all_chunk} chunked files by provided argument')
        print(f'\t')
        print(f'Input {all_chunk} chunked files')

    # Process each chunk
    output_dict = defaultdict(list)
    for chunk_index in range(1, all_chunk + 1):
        print(f'------Processing chunk-{chunk_index}')

        # Load the spatial annotations for this chunk
        ld_file_spatial = f'{config.ldscore_input_dir}/{sample_name}_chunk{chunk_index}/{sample_name}.'
        ref_ld_spatial = _read_ref_ld_v2(ld_file_spatial)
        ref_ld_spatial = ref_ld_spatial.loc[baseline_and_w_ld_common_snp]
        ref_ld_spatial = ref_ld_spatial.astype(np.float32, copy=False)

        # check the variance of the design matrix
        # ref_ld_spatial, spatial_annotation = _check_variance_v2(ref_ld_spatial, ref_ld_baseline)

        # load M file of spatial annotation
        n_annot_spatial = len(ref_ld_spatial.columns)
        n_annot = n_annot_baseline + n_annot_spatial
        M_annot_spatial = _read_M_v2(ld_file_spatial, n_annot_spatial, config.not_M_5_50)
        M_annot = np.concatenate((M_annot_baseline, M_annot_spatial), axis=1)

        # precompute x_tot and M_tot
        x_tot_precomputed = ref_ld_baseline.sum(axis=1) + ref_ld_spatial.sum(axis=1)
        M_tot = np.sum(M_annot)

        for trait_name, sumstats in sumstats_cleaned_dict.items():
            logger.info(f'Processing {trait_name}...')

            # filter ldscore by common snp
            common_snp = sumstats.index
            spatial_annotation = ref_ld_spatial.loc[common_snp].astype(np.float32, copy=False)
            spatial_annotation_cnames = spatial_annotation.columns
            baseline_annotation = ref_ld_baseline.loc[common_snp].astype(np.float32, copy=False)
            w_ld_common_snp = w_ld.loc[common_snp].astype(np.float32, copy=False)
            x_tot_precomputed_common_snp = x_tot_precomputed.loc[common_snp].values

            initial_w = (
                get_weight_optimized(sumstats, x_tot_precomputed_common_snp, M_tot, w_ld_common_snp, intercept=1)
                .astype(np.float32)
                .reshape((-1, 1)))
            assert np.any(initial_w > 0), 'Weights must be > 0'

            baseline_annotation = baseline_annotation * sumstats.N.values.reshape((-1, 1)) / sumstats.N.mean()
            # append intercept
            baseline_annotation = append_intercept(baseline_annotation)

            # apply weight
            initial_w_scaled = initial_w / np.sum(initial_w)
            baseline_annotation *= initial_w_scaled
            spatial_annotation *= initial_w_scaled
            spatial_annotation = spatial_annotation.to_numpy(dtype=np.float32)
            CHISQ = sumstats.chisq.to_numpy(dtype=np.float32).reshape((-1, 1)).copy()
            y = CHISQ * initial_w_scaled
            Nbar = sumstats.N.mean()

            # Run the jackknife
            chunk_size = spatial_annotation.shape[1]
            out_chunk = process_map(jackknife_for_processmap, range(chunk_size),
                                    max_workers=config.num_processes,
                                    chunksize=10,
                                    desc=f'LDSC chunk-{chunk_index}: {trait_name}')

            # cache the results
            out_chunk = pd.DataFrame.from_records(out_chunk,
                                                  columns=['beta', 'se', ],
                                                  index=spatial_annotation_cnames)
            out_chunk['z'] = out_chunk.beta / out_chunk.se
            out_chunk['p'] = norm.sf(out_chunk['z'])
            output_dict[trait_name].append(out_chunk)

    # Save the results
    out_dir = Path(config.ldsc_save_dir)
    out_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
    for trait_name, out_chunk_list in output_dict.items():
        out_all = pd.concat(out_chunk_list, axis=0)
        out_file_name = out_dir / f'{sample_name}_{trait_name}.csv.gz'
        out_all['spot'] = out_all.index
        out_all = out_all[['spot', 'beta', 'se', 'z', 'p']]
        out_all.to_csv(out_file_name, compression='gzip', index=False)
        logger.info(f'Output saved to {out_file_name} for {trait_name}')
    logger.info(f'------Spatial LDSC for {sample_name} finished!')


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
    import os

    os.chdir('/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/macaque/representative_slices2')
    config = SpatialLDSCConfig(**{'all_chunk': 2,
                                  'chisq_max': None,
                                  # 'sumstats_file': '/storage/yangjianLab/songliyang/GWAS_trait/LDSC/ADULT1_ADULT2_ONSET_ASTHMA.sumstats.gz',
                                  'ldsc_save_dir': 'T135_macaque1/spatial_ldsc',
                                  'ldscore_input_dir': 'T135_macaque1/generate_ldscore',
                                  'n_blocks': 200,
                                  'not_M_5_50': False,
                                  'num_processes': 4,
                                  'sample_name': 'T135_macaque1',
                                  # 'trait_name': 'ADULT1_ADULT2_ONSET_ASTHMA',
                                  'sumstats_config_file': '/storage/yangjianLab/chenwenhao/projects/202312_GPS/src/GPS/example/sumstats_config_sub.yaml',
                                  'w_file': '/storage/yangjianLab/sharedata/LDSC_resource/LDSC_SEG_ldscores/weights_hm3_no_hla/weights.'
                                  })
    # config = SpatialLDSCConfig(**vars(args))
    run_spatial_ldsc(config)

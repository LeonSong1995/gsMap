import argparse
import gc
import multiprocessing
import os
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

import GPS.jackknife as jk
from GPS.regression_read import _read_sumstats, _read_w_ld, _read_ref_ld, _read_M, _merge_and_log, _check_variance, \
    _read_ref_ld_v2, _read_M_v2, _check_variance_v2
from GPS.config import add_spatial_ldsc_args, SpatialLDSCConfig


# %load_ext autoreload
# %autoreload 2


# Set regression weight
class Regression_weight:
    """
    Weight the LDSC regression
    """

    def __init__(self, sumstats, ref_ld_cnames, w_ld_cname, M_annot, n_baseline, chisq_max=None):

        self.sumstats = sumstats
        self.ref_ld_cnames = ref_ld_cnames
        self.w_ld_cname = w_ld_cname
        self.M_annot = M_annot
        self.n_baseline = n_baseline
        self.intercept = 1

        self.n_snp = len(self.sumstats)
        self.n_annot = len(self.ref_ld_cnames)
        s = lambda x: np.array(x).reshape((self.n_snp, 1))

        # Convert the data
        self.ref_ld = np.array(self.sumstats[self.ref_ld_cnames])

        # Check the input data
        # self._check_ld_condnum()
        self._warn_length()

        # Remove SNPs with high chisq
        if chisq_max is None:
            self.chisq_max = max(0.001 * sumstats.N.max(), 80)
        else:
            self.chisq_max = chisq_max

        self.chisq = s(self.sumstats.Z ** 2)
        ii = np.ravel(self.chisq < self.chisq_max)
        self.sumstats = self.sumstats.iloc[ii, :]
        print('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
            C=self.chisq_max, N=np.sum(ii), M=self.n_snp - np.sum(ii)))
        self.n_snp = np.sum(ii)
        self.ref_ld = np.array(self.sumstats[self.ref_ld_cnames])
        self.chisq = self.chisq[ii].reshape((self.n_snp, 1))

        # Update Inpute data
        self.w = s(self.sumstats[self.w_ld_cname])
        self.N = s(self.sumstats.N)
        self.n_snp, self.n_annot = self.ref_ld.shape
        self.M_tot = float(np.sum(self.M_annot))

        self.initial_w = self.get_weight()

    # Get the weight
    def get_weight(self):
        self.x_tot = np.sum(self.ref_ld, axis=1).reshape((self.n_snp, 1))

        self.tot_agg = self.aggregate(self.chisq, self.x_tot, self.N, self.M_tot, self.intercept)
        initial_w = self.weights(self.x_tot, self.w, self.N, self.M_tot, self.tot_agg, self.intercept)
        initial_w = np.sqrt(initial_w)

        return initial_w

    def weight_yx(self):
        basic_annotation = self.ref_ld[:, 0:self.n_baseline]
        cell_annotation = self.ref_ld[:, self.n_baseline:]

        Nbar = np.mean(self.N)
        basic_annotation = np.multiply(self.N, basic_annotation) / Nbar
        cell_annotation = np.multiply(self.N, cell_annotation) / Nbar

        # Append intercept
        basic_annotation = self.append_intercept(basic_annotation)

        # Apply weight
        basic_annotation = self.apply_weights(basic_annotation, self.initial_w)
        cell_annotation = self.apply_weights(cell_annotation, self.initial_w)
        y = self.apply_weights(self.chisq, self.initial_w)

        return y, basic_annotation, cell_annotation, Nbar

    def _check_ld_condnum(self):
        '''Check condition number'''
        if len(self.ref_ld.shape) >= 2:
            cond_num = int(np.linalg.cond(self.ref_ld))
            if cond_num > 100000:
                if self.args.invert_anyway:
                    warn = "WARNING: LD Score matrix condition number is {C}. "
                    warn += "Inverting anyway because the --invert-anyway flag is set."
                    print(warn.format(C=cond_num))
                else:
                    warn = "WARNING: LD Score matrix condition number is {C}. "
                    warn += "Remove collinear LD Scores. "
                    raise ValueError(warn.format(C=cond_num))

    def _warn_length(self):
        if len(self.sumstats) < 200000:
            print('WARNING: number of SNPs less than 200k; this is almost always bad.')

    def aggregate(self, y, x, N, M, intercept=1):
        num = M * (np.mean(y) - intercept)
        denom = np.mean(np.multiply(x, N))
        return num / denom

    def append_intercept(self, x):
        n_row = x.shape[0]
        intercept = np.ones((n_row, 1))
        x_new = np.concatenate((x, intercept), axis=1)

        return x_new

    def weights(self, ld, w_ld, N, M, hsq, intercept=1):
        M = float(M)
        hsq = np.clip(hsq, 0.0, 1.0)
        ld = np.maximum(ld, 1.0)
        w_ld = np.maximum(w_ld, 1.0)
        c = hsq * N / M
        het_w = 1.0 / (2 * np.square(intercept + np.multiply(c, ld)))
        oc_w = 1.0 / w_ld
        w = np.multiply(het_w, oc_w)
        return w

    def apply_weights(self, x, w):
        if np.any(w <= 0):
            raise ValueError('Weights must be > 0')
        (n, p) = x.shape
        if w.shape != (n, 1):
            raise ValueError(
                'w has shape {S}. w must have shape (n, 1).'.format(S=w.shape))
        # -
        w = w / float(np.sum(w))
        x_new = np.multiply(x, w)
        return x_new

    # Fun for running LDSC


def _coef(jknife, Nbar):
    '''Get coefficient estimates + cov from the jackknife.'''
    coef = jknife.est[0, 0:n_annot] / Nbar
    coef_cov = jknife.jknife_cov[0:n_annot, 0:n_annot] / Nbar ** 2
    coef_se = np.sqrt(np.diag(coef_cov))
    z = coef / coef_se
    return coef, coef_cov, coef_se, z


def process_columns_cpu_worker(args):
    """LDSC for one spot"""
    i, n_blocks, Nbar, n_snp = args

    # Block jknife for current design matrix
    x_focal = np.concatenate((np.reshape(spatial_annotation[:, i], (n_snp, 1)), baseline_annotation), axis=1)
    jknife = jk.LstsqJackknifeFast(x_focal, y, n_blocks)
    coef, coef_cov, coef_se, z = _coef(jknife, Nbar)

    # Return the regression results
    return [coef[0], coef_se[0], z[0]]


def process_columns_cpu(n_blocks, Nbar, n_snp, chunk_index, num_processes=2):
    """LDSC for all spots"""
    chunk_size = spatial_annotation.shape[1]
    output = []

    # Process cells in parallel
    print(f'Running LDSC for {chunk_size} cells in chunk-{chunk_index}.')
    with Pool(num_processes) as p:
        # with tqdm(total=chunk_size, desc=f"Running ldsc of chunk-{chunk_index} ") as progress_bar:
        for result in p.imap(process_columns_cpu_worker,
                             [(i, n_blocks, Nbar, n_snp) for i in range(chunk_size)]
                             ):
            # progress_bar.update()
            output.append(result)

    return output


def run_spatial_ldsc(config: SpatialLDSCConfig):
    global data_name, name, all_chunk, n_annot, y, baseline_annotation, spatial_annotation, out_chunk
    trait_name = config.trait_name
    print(f'------Running Spatial LDSC for {trait_name}...')
    data_name = config.sample_name
    num_cpus = min(multiprocessing.cpu_count(), config.num_processes)
    # Load the gwas summary statistics
    sumstats = _read_sumstats(fh=config.h2, alleles=False, dropna=False)
    sumstats.set_index('SNP', inplace=True)
    # Load the regression weights
    w_ld = _read_w_ld(config.w_file)
    w_ld_cname = w_ld.columns[1]
    w_ld.set_index('SNP', inplace=True)
    # Load the baseline annotations
    ld_file_baseline = f'{config.ldscore_input_dir}/baseline/baseline.'
    ref_ld_baseline = _read_ref_ld_v2(ld_file_baseline)
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
        ref_ld_spatial_cnames = ref_ld_spatial.columns

        n_annot_spatial = len(ref_ld_spatial.columns)
        M_annot_spatial = _read_M_v2(ld_file_spatial, n_annot_spatial, config.not_M_5_50)

        # Merge the spatial annotations and baseline annotations
        # ref_ld = pd.concat([ref_ld_baseline, ref_ld_spatial], axis=1)
        n_annot = n_annot_baseline + n_annot_spatial
        M_annot = np.concatenate((M_annot_baseline, M_annot_spatial), axis=1)
        ref_ld_cnames = ref_ld_baseline.columns.to_list() + ref_ld_spatial.columns.to_list()

        # # Check the variance of the design matrix
        M_annot, ref_ld_spatial = _check_variance_v2(M_annot, ref_ld_spatial)

        # Merge gwas summary statistics with annotations, and regression weights
        sumstats_chunk = pd.concat([sumstats, ref_ld_baseline, ref_ld_spatial, w_ld], axis=1, join='inner')

        # Weight gwas summary statistics
        re = Regression_weight(sumstats_chunk, ref_ld_cnames, w_ld_cname, M_annot, n_annot_baseline)
        y, baseline_annotation, spatial_annotation, Nbar = re.weight_yx()
        del sumstats_chunk, w_ld, M_annot, ld_file_spatial
        gc.collect()

        # Run LDSC
        out_chunk = process_columns_cpu(config.n_blocks, Nbar, re.n_snp, chunk_index, num_cpus)
        out_chunk = pd.DataFrame(out_chunk)
        out_chunk.index = ref_ld_spatial_cnames
        out_chunk.columns = ['beta', 'se', 'z']
        out_chunk['p'] = norm.sf(out_chunk['z'])

        # Concat results
        out_all = pd.concat([out_all, out_chunk], axis=0)
    # Save the results
    print(f'------Saving the results...')
    out_dir = Path(config.ldsc_save_dir)
    out_dir.mkdir(parents=True, exist_ok=True, mode=0o777)
    out_file_name = out_dir / f'{data_name}_{trait_name}.csv.gz'
    out_all['spot'] = out_all.index
    out_all = out_all[['spot', 'beta', 'se', 'z', 'p']]
    out_all.to_csv(out_file_name, compression='gzip', index=False)


# Main function of analysis
if __name__ == '__main__':
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
            '--trait_name','adult1_adult2_onset_asthma'
        ]
        args = parser.parse_args(args_list)
    else:
        args = parser.parse_args()
    config = SpatialLDSCConfig(**vars(args))
    run_spatial_ldsc(config)

import pprint
from dataclasses import dataclass
from pathlib import Path

import submitit

import gsMap.diagnosis
from gsMap.cauchy_combination_test import run_Cauchy_combination
from gsMap.config import GenerateLDScoreConfig, SpatialLDSCConfig, VisualizeConfig, \
    LatentToGeneConfig, FindLatentRepresentationsConfig, CauchyCombinationConfig
from gsMap.find_latent_representation import run_find_latent_representation
from gsMap.generate_ldscore import run_generate_ldscore
from gsMap.latent_to_gene import run_latent_to_gene
from gsMap.spatial_ldsc_multiple_sumstats import run_spatial_ldsc
from gsMap.visualize import run_Visualize

# %% ========================================
# 1.1 Configuration
# %% ========================================
SPATIAL_LDSC_PROCESS_NUM = 10

@dataclass
class Config_Mouse:
    # Constants and configuration
    WORKDIR: str = './example/Mouse_Embryo'
    SAMPLE_NAME: str = "E16.5_E1S1.MOSTA"

    # Input data
    # HDF5_PATH: str = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/example_data/ST/E16.5_E1S1.MOSTA.h5ad"
    HDF5_PATH: str = "/storage/yangjianLab/songliyang/SpatialData/Data/Embryo/Mice/Cell_MOSTA/h5ad/E16.5_E1S1.MOSTA.h5ad"
    ANNOTATION: str = "annotation"
    DATA_TYPE: str = 'count'

    # === ORTHOLOGOUS PARAMETERS ===
    species = 'MOUSE_GENE_SYM'
    gs_species = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/GPS_resource/homologs/mouse_human_homologs.txt'

    ## === FIXED PARAMETERS ===
    # Running Dependencies and Resources
    GTFFILE: str = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/gsMap_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
    ALL_ENHANCER_FILE: str = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/gsMap_resource/genome_annotation/enhancer/by_tissue/ALL/ABC_roadmap_merged.bed"
    BFILE_ROOT: str = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/gsMap_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
    KEEP_SNP_ROOT: str = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/gsMap_resource/LDSC_resource/hapmap3_snps/hm"
    W_FILE: str = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/gsMap_resource/LDSC_resource/weights_hm3_no_hla/weights."


@dataclass
class Mouse_Without_Denoise(Config_Mouse):
    WORKDIR: str = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/20240817_vanilla_pipeline_mouse_embryo_v4'


# Example of how to use the dataclasses
config = Mouse_Without_Denoise()

find_latent_config = FindLatentRepresentationsConfig(
    input_hdf5_path=config.HDF5_PATH,
    output_hdf5_path=f"{config.WORKDIR}/{config.SAMPLE_NAME}/find_latent_representations/{config.SAMPLE_NAME}_add_latent.h5ad",
    sample_name=config.SAMPLE_NAME,
    annotation=config.ANNOTATION,
    type=config.DATA_TYPE
)
latent_to_gene_config = LatentToGeneConfig(
        input_hdf5_with_latent_path=find_latent_config.output_hdf5_path,
        sample_name=config.SAMPLE_NAME,
        output_feather_path=f"{config.WORKDIR}/{config.SAMPLE_NAME}/latent_to_gene/{config.SAMPLE_NAME}_gene_marker_score.feather",
        annotation=config.ANNOTATION,
        type=config.DATA_TYPE,
        latent_representation='latent_GVAE',
        num_neighbour=51,
        num_neighbour_spatial=201,
        species=config.species,
        gs_species=config.gs_species,

    )

ldscore_configs = [
    GenerateLDScoreConfig(
        sample_name=config.SAMPLE_NAME,
        chrom=chrom,
        ldscore_save_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/generate_ldscore",
        mkscore_feather_file=latent_to_gene_config.output_feather_path,
        bfile_root=config.BFILE_ROOT,
        keep_snp_root=config.KEEP_SNP_ROOT,
        gtf_annotation_file=config.GTFFILE,
        spots_per_chunk=5_000,
        ldscore_save_format="feather",
    ) for chrom in range(1, 23)
]

spatial_ldsc_config = SpatialLDSCConfig(
    # sumstats_file="/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/docs_test/example_data/GWAS/IQ_NG_2018.sumstats.gz",
    # trait_name="IQ",
    w_file=config.W_FILE,
    sample_name=config.SAMPLE_NAME,
    ldscore_input_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/generate_ldscore",
    ldsc_save_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/spatial_ldsc",
    num_processes=10,
    sumstats_config_file="/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GWAS_data_Collection/sumstats_configs/sumstats_config_scPaGWAS.yaml",
    ldscore_save_format="feather",
)

# log_folder = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/20240605_without_denoise/mouse_embryo/submitit/log%j"
LOG_FOLDER = "/storage/yangjianLab/chenwenhao/projects/202312_GPS/test/20240815_without_celltype_annotation/submitit/log%j"
Path(LOG_FOLDER).parent.mkdir(exist_ok=True, parents=True)
executor = submitit.AutoExecutor(folder=LOG_FOLDER)
# executor = submitit.DebugExecutor(folder=log_folder)
# executor = submitit.LocalExecutor(folder=log_folder)

executor.update_parameters(slurm_partition='intel-sc3,amd-ep2,amd-ep2-short',
                           slurm_cpus_per_task=4,
                           slurm_qos='huge',
                           slurm_time='1-00:00:00',
                           slurm_exclude='amdepyc06',
                           )
# %%
print("Start running Pipeline")
# print("1. Find latent representations")
executor.update_parameters(
    slurm_mem='60G',
    slurm_cpus_per_task=4,
    job_name=f'1_find_latent_{config.SAMPLE_NAME}', )

if Path(find_latent_config.output_hdf5_path).exists():
    print(f"Skip {config.SAMPLE_NAME} find latent")
else:
    print(f"Run {config.SAMPLE_NAME} find latent")
    find_latent_job = executor.submit(run_find_latent_representation, find_latent_config)
    find_latent_job.wait()

    if find_latent_job.state == 'OUT_OF_MEMORY':
        print(f"Re-run {config.SAMPLE_NAME} find latent")
        executor_local = submitit.DebugExecutor(folder=LOG_FOLDER,
                                                )
        executor_local.update_parameters(
            timeout_min=60 * 24, )

        find_latent_job = executor_local.submit(run_find_latent_representation, find_latent_config)
        find_latent_job.result()
    else:
        find_latent_job.result()

print("2. Latent to gene")
# pprint.pprint(latent_to_gene_config)
executor.update_parameters(
    slurm_mem='80G',
    job_name=f'2_latent_to_gene_{config.SAMPLE_NAME}', )
# with memray.Tracker(Path(f"{config.WORKDIR}/{config.SAMPLE_NAME}/latent_to_gene")/'memray.log'):
if Path(latent_to_gene_config.output_feather_path).exists():
    print(f"Skip {config.SAMPLE_NAME} latent to gene")
else:
    print(f"Run {config.SAMPLE_NAME} latent to gene")
    latent_to_gene_job = executor.submit(run_latent_to_gene, latent_to_gene_config)
    latent_to_gene_job.wait()

    if latent_to_gene_job.state == 'OUT_OF_MEMORY':
        print(f"Re-run {config.SAMPLE_NAME} latent to gene")
        executor_local = submitit.DebugExecutor(folder=LOG_FOLDER,
                                                )
        executor_local.update_parameters(
            timeout_min=60 * 24, )

        latent_to_gene_job = executor_local.submit(run_latent_to_gene, latent_to_gene_config)
        latent_to_gene_job.result()
    else:
        latent_to_gene_job.result()

# %%
import yaml
with open(spatial_ldsc_config.sumstats_config_file, 'r') as f:
    sumstats_config = yaml.load(f, Loader=yaml.FullLoader)
# %% run generate ldscore
jobs = []
executor.update_parameters(
    slurm_job_name=f'3_generate_ldscore_{config.SAMPLE_NAME}',
    slurm_cpus_per_task=10,
    slurm_mem='35G',
)
with executor.batch():
    for ldscore_config in ldscore_configs:
        pprint.pprint(ldscore_config)
        job = executor.submit(run_generate_ldscore, ldscore_config)
        jobs.append(job)
# %%
for i, job in enumerate(jobs):
    job.result()
(Path(ldscore_configs[0].ldscore_save_dir) / 'generate_ldscore.done').touch()

if (Path(spatial_ldsc_config.ldsc_save_dir) / 'spatial_ldsc.done').exists():
    print(f"Skip {config.SAMPLE_NAME} spatial ldsc")
else:
    # submit each trait
    spatial_ldsc_jobs = []
    for trait_name in sumstats_config:
        print(f'run spatial ldsc on {trait_name}')
        spatial_ldsc_config = SpatialLDSCConfig(
            sumstats_file=sumstats_config[trait_name],
            trait_name=trait_name,
            w_file=config.W_FILE,
            sample_name=config.SAMPLE_NAME,
            ldscore_input_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/generate_ldscore",
            ldsc_save_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/spatial_ldsc",
            num_processes=SPATIAL_LDSC_PROCESS_NUM,
            ldscore_save_format="feather",
        )
        executor.update_parameters(slurm_mem='45G',
                                   slurm_job_name=f'4_spatial_ldsc_{config.SAMPLE_NAME}',
                                   slurm_cpus_per_task=SPATIAL_LDSC_PROCESS_NUM,
                                   # slurm_dependency=f'afterok:{",".join([str(job.job_id) for job in out_of_memory_jobs + jobs])}'
                                   # if len(out_of_memory_jobs) > 0 else None
                                   )
        spatial_ldsc_job = executor.submit(run_spatial_ldsc, spatial_ldsc_config)
        spatial_ldsc_jobs.append(spatial_ldsc_job)

    for spatial_ldsc_job in spatial_ldsc_jobs:
        spatial_ldsc_job.result()
    (Path(spatial_ldsc_config.ldsc_save_dir) / 'spatial_ldsc.done').touch()

# run_spatial_ldsc(spatial_ldsc_config)
# %% run visualize
print('start_visualization')
for trait_name in sumstats_config:
    visualize_config = VisualizeConfig(
        input_hdf5_path=config.HDF5_PATH,
        input_ldsc_dir=spatial_ldsc_config.ldsc_save_dir,
        output_figure_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/visualize",
        sample_name=config.SAMPLE_NAME,
        trait_name=trait_name,
        annotation=config.ANNOTATION,
        fig_title=trait_name,
        point_size=2,
    )
    run_Visualize(visualize_config)

# %% run cauchy combination
for trait_name in sumstats_config:
    cauchy_config = CauchyCombinationConfig(
        input_hdf5_path=config.HDF5_PATH,
        input_ldsc_dir=spatial_ldsc_config.ldsc_save_dir,
        sample_name=config.SAMPLE_NAME,
        annotation=config.ANNOTATION,
        output_cauchy_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/cauchy_combination",
        trait_name=trait_name,
    )
    if Path(Path(cauchy_config.output_cauchy_dir) / f'{config.SAMPLE_NAME}_{trait_name}.Cauchy.csv.gz').exists():
        print(f"Skip {config.SAMPLE_NAME} cauchy combination {trait_name}")
    else:
        run_Cauchy_combination(cauchy_config)

# %% run diagnosis
'''  config = DiagnosisConfig(
        sample_name='E16.5_E1S1.MOSTA',
        input_hdf5_path='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/example_data/ST/E16.5_E1S1.MOSTA.h5ad',
        annotation='annotation',
        mkscore_feather_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/latent_to_gene/E16.5_E1S1.MOSTA_gene_marker_score.feather',
        bfile_root='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/GPS_resource/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC',
        input_ldsc_dir='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/spatial_ldsc',
        trait_name='GIANT_EUR_Height_2022_Nature',
        gtf_annotation_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/GPS_resource/genome_annotation/gtf/gencode.v39lift37.annotation.gtf',
        gene_window_size=50000,
        top_corr_genes=10,
        selected_genes=['COL11A1', 'MECOM'],
        sumstats_file='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/example_data/GWAS/GIANT_EUR_Height_2022_Nature.sumstats.gz',
        diagnosis_save_dir='/mnt/e/0_Wenhao/7_Projects/20231213_GPS_Liyang/test/20240902_gsMap_Local_Test/E16.5_E1S1.MOSTA/diagnosis'
    )'''

for trait_name in sumstats_config:
    diagnosis_config = gsMap.diagnosis.DiagnosisConfig(
        sample_name=config.SAMPLE_NAME,
        input_hdf5_path=config.HDF5_PATH,
        annotation=config.ANNOTATION,
        mkscore_feather_file=latent_to_gene_config.output_feather_path,
        input_ldsc_dir=spatial_ldsc_config.ldsc_save_dir,
        trait_name=trait_name,
        plot_type='manhattan',
        ldscore_save_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/generate_ldscore",

        gene_window_size=50000,
        top_corr_genes=50,
        selected_genes=None,
        sumstats_file=sumstats_config[trait_name],
        diagnosis_save_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/diagnosis"

    )
    gsMap.diagnosis.generate_manhattan_plot(diagnosis_config)

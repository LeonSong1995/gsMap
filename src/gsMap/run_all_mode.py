import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional

from gsMap.diagnosis import DiagnosisConfig, generate_manhattan_plot, run_Diagnosis
from gsMap.cauchy_combination_test import run_Cauchy_combination
from gsMap.config import GenerateLDScoreConfig, SpatialLDSCConfig, VisualizeConfig, LatentToGeneConfig, \
    FindLatentRepresentationsConfig, CauchyCombinationConfig
from gsMap.find_latent_representation import run_find_latent_representation
from gsMap.generate_ldscore import run_generate_ldscore
from gsMap.latent_to_gene import run_latent_to_gene
from gsMap.report import run_Report
from gsMap.spatial_ldsc_multiple_sumstats import run_spatial_ldsc
from gsMap.visualize import run_Visualize
import yaml

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.FileHandler("pipeline.log"),
                        logging.StreamHandler()
                    ])
logger = logging.getLogger('gsMap Pipeline')



@dataclass
class RunAllModeConfig:
    workdir: str
    sample_name: str

    gsMap_resource_dir: str

    # == ST DATA PARAMETERS ==
    hdf5_path: str
    annotation: str
    data_layer: str = 'X'

    # ==GWAS DATA PARAMETERS==
    trait_name: Optional[str] = None
    sumstats_file: Optional[str] = None
    sumstats_config_file: Optional[str] = None

    # === homolog PARAMETERS ===
    homolog_file: Optional[str] = None

    num_processes: int = 5

    def __post_init__(self):
        logger.info(f"Running pipeline with configuration: {self}")

        self.gtffile = f"{self.gsMap_resource_dir}/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
        self.bfile_root = f"{self.gsMap_resource_dir}/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
        self.keep_snp_root = f"{self.gsMap_resource_dir}/LDSC_resource/hapmap3_snps/hm"
        self.w_file = f"{self.gsMap_resource_dir}/LDSC_resource/weights_hm3_no_hla/weights."

        if self.homolog_file is not None:
            logger.info(f"User provided homolog file to map gene names to human: {self.homolog_file}")
            # check the format of the homolog file
            with open(self.homolog_file, 'r') as f:
                first_line = f.readline().strip()
                _n_col = len(first_line.split())
                if _n_col != 2:
                    raise ValueError(f"Invalid homolog file format. Expected 2 columns, first column should be other species gene name, second column should be human gene name. "
                                     f"Got {_n_col} columns in the first line.")
                else:
                    first_col_name, second_col_name = first_line.split()
                    logger.info(f"Successfully parse homolog file with format and will map {first_col_name} to {second_col_name}")
        else:
            logger.info("No homolog file provided. Run in human mode.")

        # check the existence of the input files and resources files
        for file in [self.hdf5_path, self.gtffile]:
                raise FileNotFoundError(f"File {file} does not exist.")


def format_duration(seconds):
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    return f"{hours}h {minutes}m"


def run_pipeline(config: RunAllModeConfig):
    logger.info("Starting pipeline with configuration: %s", config)

    find_latent_config = FindLatentRepresentationsConfig(
        input_hdf5_path=config.hdf5_path,
        output_hdf5_path=f"{config.workdir}/{config.sample_name}/find_latent_representations/{config.sample_name}_add_latent.h5ad",
        sample_name=config.sample_name,
        annotation=config.annotation,
        type=config.data_layer
    )

    latent_to_gene_config = LatentToGeneConfig(
        input_hdf5_with_latent_path=find_latent_config.output_hdf5_path,
        sample_name=config.sample_name,
        output_feather_path=f"{config.workdir}/{config.sample_name}/latent_to_gene/{config.sample_name}_gene_marker_score.feather",
        annotation=config.annotation,
        type=config.data_layer,
        latent_representation='latent_GVAE',
        num_neighbour=51,
        num_neighbour_spatial=201,
        homolog_file=config.homolog_file
    )

    ldscore_config = GenerateLDScoreConfig(
        sample_name=config.sample_name,
        chrom='all',
        ldscore_save_dir=f"{config.workdir}/{config.sample_name}/generate_ldscore",
        mkscore_feather_file=latent_to_gene_config.output_feather_path,
        bfile_root=config.bfile_root,
        keep_snp_root=config.keep_snp_root,
        gtf_annotation_file=config.gtffile,
        spots_per_chunk=5_000,
        ldscore_save_format="feather",
    )

    spatial_ldsc_config = SpatialLDSCConfig(
        w_file=config.w_file,
        sample_name=config.sample_name,
        ldscore_input_dir=f"{config.workdir}/{config.sample_name}/generate_ldscore",
        ldsc_save_dir=f"{config.workdir}/{config.sample_name}/spatial_ldsc",
        num_processes=config.num_processes,
        trait_name=config.trait_name,
        sumstats_config_file=config.sumstats_config_file,
        sumstats_file=config.sumstats_file,
    )

    pipeline_start_time = time.time()

    # Step 1: Find latent representations
    start_time = time.time()
    logger.info("Step 1: Finding latent representations")
    if Path(find_latent_config.output_hdf5_path).exists():
        logger.info("Latent representations already found, skipping...")
    else:
        run_find_latent_representation(find_latent_config)
    end_time = time.time()
    logger.info(f"Step 1 completed in {format_duration(end_time - start_time)}.")

    # Step 2: Latent to gene
    start_time = time.time()
    logger.info("Step 2: Mapping latent representations to genes")
    if Path(latent_to_gene_config.output_feather_path).exists():
        logger.info("Latent to gene mapping already done, skipping...")
    else:
        run_latent_to_gene(latent_to_gene_config)
    end_time = time.time()
    logger.info(f"Step 2 completed in {format_duration(end_time - start_time)}.")

    # Step 3: Generate LDScores
    start_time = time.time()
    logger.info("Step 3: Generating LDScores")
    run_generate_ldscore(ldscore_config)
    end_time = time.time()
    logger.info(f"Step 3 completed in {format_duration(end_time - start_time)}.")

    # Step 4: Spatial LDSC
    start_time = time.time()
    logger.info("Step 4: Running spatial LDSC")
    with open(spatial_ldsc_config.sumstats_config_file, 'r') as f:
        sumstats_config = yaml.load(f, Loader=yaml.FullLoader)

    for trait_name in sumstats_config:
        logger.info("Running spatial LDSC for trait: %s", trait_name)
        spatial_ldsc_config_trait = SpatialLDSCConfig(
            sumstats_file=sumstats_config[trait_name],
            trait_name=trait_name,
            w_file=config.w_file,
            sample_name=config.sample_name,
            ldscore_input_dir=spatial_ldsc_config.ldscore_input_dir,
            ldsc_save_dir=spatial_ldsc_config.ldsc_save_dir,
            num_processes=spatial_ldsc_config.num_processes,
            ldscore_save_format=spatial_ldsc_config.ldscore_save_format,
        )
        run_spatial_ldsc(spatial_ldsc_config_trait)
    end_time = time.time()
    logger.info(f"Step 4 completed in {format_duration(end_time - start_time)}.")

    # Step 5: Visualization
    start_time = time.time()
    logger.info("Step 5: Visualization")
    for trait_name in sumstats_config:
        visualize_config = VisualizeConfig(
            input_hdf5_path=config.hdf5_path,
            input_ldsc_dir=spatial_ldsc_config.ldsc_save_dir,
            output_figure_dir=f"{config.workdir}/{config.sample_name}/visualize",
            sample_name=config.sample_name,
            trait_name=trait_name,
            annotation=config.annotation,
            fig_title=trait_name,
            point_size=2,
        )
        run_Visualize(visualize_config)
    end_time = time.time()
    logger.info(f"Step 5 completed in {format_duration(end_time - start_time)}.")

    # Step 6: Cauchy combination test
    start_time = time.time()
    logger.info("Step 6: Running Cauchy combination test")
    for trait_name in sumstats_config:
        cauchy_config = CauchyCombinationConfig(
            input_hdf5_path=config.hdf5_path,
            input_ldsc_dir=spatial_ldsc_config.ldsc_save_dir,
            sample_name=config.sample_name,
            annotation=config.annotation,
            output_cauchy_dir=f"{config.workdir}/{config.sample_name}/cauchy_combination",
            trait_name=trait_name,
        )
        run_Cauchy_combination(cauchy_config)
    end_time = time.time()
    logger.info(f"Step 6 completed in {format_duration(end_time - start_time)}.")

    # Step 7: Diagnosis
    start_time = time.time()
    logger.info("Step 7: Running diagnosis")
    for trait_name in sumstats_config:
        diagnosis_config = DiagnosisConfig(
            sample_name=config.sample_name,
            input_hdf5_path=config.hdf5_path,
            annotation=config.annotation,
            mkscore_feather_file=latent_to_gene_config.output_feather_path,
            input_ldsc_dir=spatial_ldsc_config.ldsc_save_dir,
            trait_name=trait_name,
            plot_type='manhattan',
            ldscore_save_dir=f"{config.workdir}/{config.sample_name}/generate_ldscore",
            top_corr_genes=50,
            selected_genes=None,
            sumstats_file=sumstats_config[trait_name],
            diagnosis_save_dir=f"{config.workdir}/{config.sample_name}/diagnosis"
        )
        run_Diagnosis(diagnosis_config)

        diagnosis_config.plot_type = 'GSS'
        diagnosis_config.top_corr_genes = 20
        run_Diagnosis(diagnosis_config)
    end_time = time.time()
    logger.info(f"Step 7 completed in {format_duration(end_time - start_time)}.")


    # Step 8: Report

    for trait_name in sumstats_config:
        logger.info("Running final report generation for trait: %s", trait_name)

        # Create the run parameters dictionary for each trait
        run_parameter_dict = {
            "sample_name": config.sample_name,
            "trait_name": trait_name,
            "sumstats_file": sumstats_config[trait_name],
            "hdf5_path": config.hdf5_path,
            "annotation": config.annotation,

            "num_processes": spatial_ldsc_config.num_processes,
            "ldscore_dir": ldscore_config.ldscore_save_dir,
            "w_file": config.w_file,
            "gtf_annotation_file": config.gtffile,
            "bfile_root": config.bfile_root,
            "keep_snp_root": config.keep_snp_root,
            "mkscore_feather_file": latent_to_gene_config.output_feather_path,
            "spatial_ldsc_save_dir": spatial_ldsc_config.ldsc_save_dir,
            "cauchy_dir": f"{config.workdir}/{config.sample_name}/cauchy_combination",
            'visualize_dir': f"{config.workdir}/{config.sample_name}/visualize",
            "diagnosis_dir": f"{config.workdir}/{config.sample_name}/diagnosis",

            "Spending_time": format_duration(time.time() - pipeline_start_time),
            "Finish_time": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

        }

        # Pass the run parameter dictionary to the report generation function
        run_Report(
            result_dir=f"{config.workdir}/{config.sample_name}",
            sample_name=config.sample_name,
            trait_name=trait_name,
            run_parameters=run_parameter_dict
        )


    logger.info("Pipeline completed successfully.")

if __name__ == '__main__':
    # Example usage:
    config = Mouse_Without_Denoise()
    run_pipeline(config)

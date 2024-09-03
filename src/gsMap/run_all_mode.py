import logging
from dataclasses import dataclass
from pathlib import Path

from gsMap.diagnosis import DiagnosisConfig, generate_manhattan_plot
from gsMap.cauchy_combination_test import run_Cauchy_combination
from gsMap.config import GenerateLDScoreConfig, SpatialLDSCConfig, VisualizeConfig, LatentToGeneConfig, FindLatentRepresentationsConfig, CauchyCombinationConfig
from gsMap.find_latent_representation import run_find_latent_representation
from gsMap.generate_ldscore import run_generate_ldscore
from gsMap.latent_to_gene import run_latent_to_gene
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
logger = logging.getLogger(__name__)

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


def run_pipeline(config: Config_Mouse):
    logger.info("Starting pipeline with configuration: %s", config)

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

    ldscore_config = GenerateLDScoreConfig(
            sample_name=config.SAMPLE_NAME,
            chrom='all',
            ldscore_save_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/generate_ldscore",
            mkscore_feather_file=latent_to_gene_config.output_feather_path,
            bfile_root=config.BFILE_ROOT,
            keep_snp_root=config.KEEP_SNP_ROOT,
            gtf_annotation_file=config.GTFFILE,
            spots_per_chunk=5_000,
            ldscore_save_format="feather",
        )

    spatial_ldsc_config = SpatialLDSCConfig(
        w_file=config.W_FILE,
        sample_name=config.SAMPLE_NAME,
        ldscore_input_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/generate_ldscore",
        ldsc_save_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/spatial_ldsc",
        num_processes=10,
        sumstats_config_file="/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GWAS_data_Collection/sumstats_configs/sumstats_config_scPaGWAS.yaml",
        ldscore_save_format="feather",
    )

    # Step 1: Find latent representations
    logger.info("Step 1: Finding latent representations")
    if Path(find_latent_config.output_hdf5_path).exists():
        logger.info("Latent representations already found, skipping...")
    else:
        run_find_latent_representation(find_latent_config)

    # Step 2: Latent to gene
    logger.info("Step 2: Mapping latent representations to genes")
    if Path(latent_to_gene_config.output_feather_path).exists():
        logger.info("Latent to gene mapping already done, skipping...")
    else:
        run_latent_to_gene(latent_to_gene_config)

    # Step 3: Generate LDScores
    logger.info("Step 3: Generating LDScores")
    run_generate_ldscore(ldscore_config)

    # Step 4: Spatial LDSC
    logger.info("Step 4: Running spatial LDSC")
    with open(spatial_ldsc_config.sumstats_config_file, 'r') as f:
        sumstats_config = yaml.load(f, Loader=yaml.FullLoader)

    for trait_name in sumstats_config:
        logger.info("Running spatial LDSC for trait: %s", trait_name)
        spatial_ldsc_config_trait = SpatialLDSCConfig(
            sumstats_file=sumstats_config[trait_name],
            trait_name=trait_name,
            w_file=config.W_FILE,
            sample_name=config.SAMPLE_NAME,
            ldscore_input_dir=spatial_ldsc_config.ldscore_input_dir,
            ldsc_save_dir=spatial_ldsc_config.ldsc_save_dir,
            num_processes=spatial_ldsc_config.num_processes,
            ldscore_save_format=spatial_ldsc_config.ldscore_save_format,
        )
        run_spatial_ldsc(spatial_ldsc_config_trait)

    # Step 5: Visualization
    logger.info("Step 5: Visualization")
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

    # Step 6: Cauchy combination test
    logger.info("Step 6: Running Cauchy combination test")
    for trait_name in sumstats_config:
        cauchy_config = CauchyCombinationConfig(
            input_hdf5_path=config.HDF5_PATH,
            input_ldsc_dir=spatial_ldsc_config.ldsc_save_dir,
            sample_name=config.SAMPLE_NAME,
            annotation=config.ANNOTATION,
            output_cauchy_dir=f"{config.WORKDIR}/{config.SAMPLE_NAME}/cauchy_combination",
            trait_name=trait_name,
        )
        run_Cauchy_combination(cauchy_config)

    # Step 7: Diagnosis
    logger.info("Step 7: Running diagnosis")
    for trait_name in sumstats_config:
        diagnosis_config = DiagnosisConfig(
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
        generate_manhattan_plot(diagnosis_config)

    logger.info("Pipeline completed successfully.")


# Example usage:
config = Mouse_Without_Denoise()
run_pipeline(config)

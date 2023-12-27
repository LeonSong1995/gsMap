import argparse
import logging
from typing import Optional, Literal

from omegaconf import OmegaConf, II
from dataclasses import dataclass, field
from omegaconf import MISSING


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)

def find_latent_representations(spe_path, spe_name, annotation, type, spe_out):
    """
    Find Latent Representations
    :param spe_path: Path to the SPE
    :param spe_name: SPE name
    :param annotation: Annotation type
    :param type: Type of analysis
    :param spe_out: Output path for SPE
    :return:
    """
    pass  # Replace with actual function code

def latent_to_gene(latent_representation, spe_path, spe_name, num_processes, type, annotation, num_neighbour, spe_out):
    """
    Latent to Gene V2
    :param latent_representation: Latent representation
    :param spe_path: Path to the SPE
    :param spe_name: SPE name
    :param num_processes: Number of processes
    :param type: Type of analysis
    :param annotation: Annotation type
    :param num_neighbour: Number of neighbours
    :param spe_out: Output path
    :return:
    """
    pass  # Replace with actual function code

def make_snp_annotations(mk_score_file, gtf_file, bfile_root, annot_root, keep_snp, annot_name, const_max_size, chr, ld_wind_cm):
    """
    Make Annotations
    :param mk_score_file: Path to the MK score file
    :param gtf_file: Path to the GTF file
    :param bfile_root: Bfile root
    :param annot_root: Annotation root
    :param keep_snp: Keep SNP file
    :param annot_name: Annotation name
    :param const_max_size: Constant max size
    :param chr: Chromosome ID
    :param ld_wind_cm: LD window size in centiMorgans
    :return:
    """
    pass  # Replace with actual function code

def spatial_ldsc(h2, w_file, data_name, num_processes, ld_file, out_file):
    """
    Spatial LDSC Analysis
    :param h2: Path to the GWAS sumstats file
    :param w_file: Weight file
    :param data_name: Data name
    :param num_processes: Number of processes
    :param ld_file: LD file path
    :param out_file: Output file path
    :return:
    """
    pass  # Replace with actual function code

def cauchy_combination(ldsc_path, ldsc_name, spe_path, spe_name, annotation):
    """
    Cauchy Combination Analysis
    :param ldsc_path: LDSC path
    :param ldsc_name: LDSC name
    :param spe_path: Path to SPE
    :param spe_name: SPE name
    :param annotation: Annotation type
    :return:
    """
    pass  # Replace with actual function code


def main():
    # Main parser

    parser = argparse.ArgumentParser(description="Spatial Data Analysis Tool")

    # Subparsers
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # Find_Latent_Representations subparser
    flr_parser = subparsers.add_parser("find_latent_representations", help="Find Latent Representations")
    flr_parser.add_argument("--spe_path", required=True, help="Path to SPE")
    flr_parser.add_argument("--spe_name", required=True, help="SPE Name")
    flr_parser.add_argument("--annotation", default="layer_guess", help="Annotation type")
    flr_parser.add_argument("--type", default="count", help="Type of analysis")
    flr_parser.add_argument("--spe_out", required=True, help="Output path for SPE")
    flr_parser.set_defaults(func=find_latent_representations)

    # Latent_to_Gene_V2 subparser
    ltg_parser = subparsers.add_parser("latent_to_gene", help="Latent to Gene V2")
    ltg_parser.add_argument("--latent_representation", default="latent_GVAE", help="Latent representation")
    ltg_parser.add_argument("--spe_path", required=True, help="Path to SPE")
    ltg_parser.add_argument("--spe_name", required=True, help="SPE Name")
    ltg_parser.add_argument("--num_processes", type=int, default=4, help="Number of processes")
    ltg_parser.add_argument("--type", default="count", help="Type of analysis")
    ltg_parser.add_argument("--annotation", default="layer_guess", help="Annotation type")
    ltg_parser.add_argument("--num_neighbour", type=int, default=51, help="Number of neighbours")
    ltg_parser.add_argument("--spe_out", required=True, help="Output path")
    ltg_parser.set_defaults(func=latent_to_gene)

    # Make_Annotations_V2 subparser
    ma_parser = subparsers.add_parser("make_annotations_v2", help="Make Annotations V2")
    ma_parser.add_argument("--mk_score_file", required=True, help="Path to the MK score file")
    ma_parser.add_argument("--gtf_file", required=True, help="Path to the GTF file")
    ma_parser.add_argument("--bfile_root", required=True, help="Bfile root")
    ma_parser.add_argument("--annot_root", required=True, help="Annotation root")
    ma_parser.add_argument("--keep_snp", required=True, help="Keep SNP file")
    ma_parser.add_argument("--annot_name", required=True, help="Annotation name")
    ma_parser.add_argument("--const_max_size", type=int, default=500, help="Constant max size")
    ma_parser.add_argument("--chr", required=True, help="Chromosome ID")
    ma_parser.add_argument("--ld_wind_cm", type=float, default=1, help="LD window size in centiMorgans")
    ma_parser.set_defaults(func=make_snp_annotations)

    # Spatial_LDSC subparser
    sldsc_parser = subparsers.add_parser("spatial_ldsc", help="Spatial LDSC Analysis")
    sldsc_parser.add_argument("--h2", required=True, help="Path to the GWAS sumstats file")
    sldsc_parser.add_argument("--w_file", required=True, help="Weight file")
    sldsc_parser.add_argument("--data_name", required=True, help="Data name")
    sldsc_parser.add_argument("--num_processes", type=int, default=3, help="Number of processes")
    sldsc_parser.add_argument("--ld_file", required=True, help="LD file path")
    sldsc_parser.add_argument("--out_file", required=True, help="Output file path")
    sldsc_parser.set_defaults(func=spatial_ldsc)

    # Cauchy_Combination subparser
    cc_parser = subparsers.add_parser("cauchy_combination", help="Cauchy Combination Analysis")
    cc_parser.add_argument("--ldsc_path", required=True, help="LDSC path")
    cc_parser.add_argument("--ldsc_name", required=True, help="LDSC name")
    cc_parser.add_argument("--spe_path", required=True, help="Path to SPE")
    cc_parser.add_argument("--spe_name", required=True, help="SPE name")
    cc_parser.add_argument("--annotation", default="layer_guess", help="Annotation type")
    cc_parser.set_defaults(func=cauchy_combination)

    # Parse arguments
    args = parser.parse_args()

    # Call the appropriate function based on the subcommand
    if args.command == "find_latent_representations":
        # Call the function for find_latent_representations
        pass
    elif args.command == "latent_to_gene":
        # Call the function for latent_to_gene
        pass  # Replace with actual function call
    elif args.command == "make_annotations_v2":
        # Call the function for make_annotations_v2
        pass
    elif args.command == "spatial_ldsc":
        # Call the function for spatial_ldsc
        pass
    elif args.command == "cauchy_combination":
        # Call the function for cauchy_combination
        pass
    else:
        print('Use "python GPS_main.py <subcommand> -h" for more information on a specific subcommand')
        parser.print_help()
        exit(1)

@dataclass
class ST_SAMPLE_INFO:
    sample_hdf5:str = MISSING
    sample_name:str = MISSING
    annotation_layer_name: Optional[str] = None
    is_count: bool = True


@dataclass
class FIND_LATENT_REPRESENTATIONS:
    sample_info: ST_SAMPLE_INFO = II("sample_info")
    # epochs: int = 300
    # feat_hidden1: int = 256
    # feat_hidden2: int = 128
    # feat_cell: int = 3000
    # n_comps: int = 300
    # n_neighbors: int = 11
    # p_drop: float = 0.1
    # convergence_threshold: float = 1e-4


@dataclass
class LATENT_TO_GENE:
    sample_info: ST_SAMPLE_INFO = II("sample_info")
    method: str = 'rank'
    num_neighbour: int = 21
    num_neighbour_spatial: int = 101
    num_processes: int = 4
    fold: float = 1.0
    pst: float = 0.2
    species: Optional[str] = None
    gs_species: Optional[str] = None
    gM_slices: Optional[str] = None

@dataclass
class MAKE_ANNOTATION:
    gtf_file: str = MISSING
    bfile_root: str = MISSING
    baseline_annotation: Optional[str] = None
    keep_snp_root: Optional[str] = None
    chr: Optional[int] = None
    window_size: int = 50_000
    chunk_size: int = 100
    ld_wind: int = 1
    ld_wind_unit: str = 'cm'
    r2_cache_dir: Optional[str] = None

@dataclass
class GPS:
    sample_info: ST_SAMPLE_INFO = field(default_factory=ST_SAMPLE_INFO)
    find_latent_representations: FIND_LATENT_REPRESENTATIONS = field(default_factory=FIND_LATENT_REPRESENTATIONS)
    latent_to_gene: LATENT_TO_GENE = field(default_factory=LATENT_TO_GENE)
    make_annotation: MAKE_ANNOTATION = field(default_factory=MAKE_ANNOTATION)
    output_dir: str = MISSING

cf: GPS = OmegaConf.structured(GPS)
config_yaml_path='/storage/yangjianLab/chenwenhao/projects/202312_GPS/src/GPS/test/config.yaml'
OmegaConf.to_yaml(cf)
OmegaConf.save(config=cf, f=config_yaml_path)
if __name__ == "__main__":
    main()

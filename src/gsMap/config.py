import argparse
import logging
from collections import OrderedDict, namedtuple
from dataclasses import dataclass
from pathlib import Path
from pprint import pprint
from typing import Callable
from typing import Union, Literal, Tuple, Optional, List
from functools import wraps
import pyfiglet

from gsMap.__init__ import __version__

# Global registry to hold functions
cli_function_registry = OrderedDict()
subcommand = namedtuple('subcommand', ['name', 'func', 'add_args_function', 'description'])
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:6s} {message}', style='{'))
logger.addHandler(handler)


# Decorator to register functions for cli parsing
def register_cli(name: str, description: str, add_args_function: Callable) -> Callable:
    def decorator(func: Callable) -> Callable:
        def wrapper(*args, **kwargs):
            name.replace('_', ' ')
            gsMap_main_logo = pyfiglet.figlet_format("gsMap", font='doom', width=80, justify='center', ).rstrip()
            print(gsMap_main_logo, )
            version_number = 'Version: ' + __version__
            print(version_number.center(80))
            print('=' * 80)
            logger.info(f"Running {name}...")
            func(*args, **kwargs)
            logger.info(f"Finished running {name}.")

        cli_function_registry[name] = subcommand(name=name, func=wrapper, add_args_function=add_args_function,
                                                 description=description)
        return wrapper

    return decorator


def add_find_latent_representations_args(parser):
    parser.add_argument('--input_hdf5_path', required=True, type=str, help='Path to the input hdf5 file.')
    parser.add_argument('--output_hdf5_path', required=True, type=str, help='Path to the output hdf5 file.')
    parser.add_argument('--sample_name', required=True, type=str, help='Name of the sample.')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    parser.add_argument('--type', default=None, type=str,
                        help="The gene expression layer of the input data (e.g., 'count', 'counts', 'log1p').")
    parser.add_argument('--epochs', default=300, type=int,
                        help="Number of training epochs for the GNN-VAE model. Default is 300.")

    parser.add_argument('--feat_hidden1', default=256, type=int,
                        help="Number of neurons in the first hidden layer of the feature extraction network. Default is 256.")
    parser.add_argument('--feat_hidden2', default=128, type=int,
                        help="Number of neurons in the second hidden layer of the feature extraction network. Default is 128.")
    parser.add_argument('--feat_cell', default=3000, type=int,
                        help="Number of top variable genes to select. Default is 3000.")
    parser.add_argument('--gcn_hidden1', default=64, type=int,
                        help="Number of units in the first hidden layer of the GCN. Default is 64.")
    parser.add_argument('--gcn_hidden2', default=30, type=int,
                        help="Number of units in the second hidden layer of the GCN. Default is 30.")
    parser.add_argument('--p_drop', default=0.1, type=float,
                        help="Dropout rate used in the GNN-VAE model. Default is 0.1.")
    parser.add_argument('--gcn_lr', default=0.001, type=float,
                        help="Learning rate for the GCN network. Default is 0.001.")
    parser.add_argument('--gcn_decay', default=0.01, type=float,
                        help="Weight decay (L2 penalty) for the GCN network. Default is 0.01.")
    parser.add_argument('--n_neighbors', default=11, type=int,
                        help="Number of neighbors to consider for graph construction in GCN. Default is 11.")
    parser.add_argument('--label_w', default=1, type=float,
                        help="Weight of the label loss in the loss function. Default is 1.")
    parser.add_argument('--rec_w', default=1, type=float,
                        help="Weight of the reconstruction loss in the loss function. Default is 1.")
    parser.add_argument('--n_comps', default=300, type=int,
                        help="Number of principal components to keep if PCA is performed. Default is 300.")
    parser.add_argument('--weighted_adj', action='store_true',
                        help="Use a weighted adjacency matrix in GCN. Default is False.")
    parser.add_argument('--nheads', default=3, type=int,
                        help="Number of heads in the attention mechanism of the GNN. Default is 3.")
    parser.add_argument('--var', action='store_true',
                        help="Enable var. Use --var to enable. Default is False.")
    parser.add_argument('--convergence_threshold', default=1e-4, type=float,
                        help="Threshold for convergence during training. Training stops if the loss change is below this threshold. Default is 1e-4.")
    parser.add_argument('--hierarchically', action='store_true',
                        help="Find latent representations hierarchically. Use --hierarchically to enable. Default is False.")


def chrom_choice(value):
    if value.isdigit():
        ivalue = int(value)
        if 1 <= ivalue <= 22:
            return ivalue
    elif value.lower() == 'all':
        return value
    else:
        raise argparse.ArgumentTypeError(f"'{value}' is an invalid chromosome choice. Choose from 1-22 or 'all'.")


def filter_args_for_dataclass(args_dict, data_class: dataclass):
    return {k: v for k, v in args_dict.items() if k in data_class.__dataclass_fields__}


def get_dataclass_from_parser(args: argparse.Namespace, data_class: dataclass):
    remain_kwargs = filter_args_for_dataclass(vars(args), data_class)
    print(f'Using the following arguments for {data_class.__name__}:')
    pprint(remain_kwargs)
    return data_class(**remain_kwargs)


def add_generate_ldscore_args(parser):
    parser.add_argument('--sample_name', type=str, required=True, help='Sample name')
    parser.add_argument('--chrom', type=str, required=True, help='Chromosome number (1-22) or "all"')
    parser.add_argument('--ldscore_save_dir', type=str, required=True, help='Directory to save ld score files')
    parser.add_argument('--mkscore_feather_file', type=str, required=True, help='Mkscore feather file path')

    # additional baseline annotation
    parser.add_argument('--additional_baseline_annotation_dir_path', type=str, default=None, )

    # reference panel
    parser.add_argument('--bfile_root', type=str, required=True, help='Bfile root path')
    parser.add_argument('--keep_snp_root', type=str, required=True, help='Keep SNP root path')

    # Annotation by gene distance
    parser.add_argument('--gtf_annotation_file', type=str, required=True, help='GTF file path')
    parser.add_argument('--gene_window_size', type=int, default=50000, help='Gene window size')

    # Enhancer annotation
    parser.add_argument('--enhancer_annotation_file', type=str, default=None,
                        help='Enhancer annotation bed file path, optional.')
    parser.add_argument('--snp_multiple_enhancer_strategy', type=str, default='max_mkscore',
                        choices=['max_mkscore', 'nearest_TSS'], help='Strategy for multiple enhancers per SNP')
    parser.add_argument('--gene_window_enhancer_priority', type=str, default=None,
                        choices=['gene_window_first', 'enhancer_first', 'enhancer_only'],
                        help='Priority between gene window and enhancer')

    # Arguments for calculating ld score
    parser.add_argument('--spots_per_chunk', type=int, default=5_000, help='Number of spots per chunk')
    parser.add_argument('--ld_wind', type=int, default=1, help='LD window size')
    parser.add_argument('--ld_unit', type=str, default='CM', help='LD window unit (SNP/KB/CM)',
                        choices=['SNP', 'KB', 'CM'])


def add_latent_to_gene_args(parser):
    parser.add_argument('--input_hdf5_with_latent_path', type=str, required=True,
                        help='Path to the input HDF5 file which contains latent representations.')
    parser.add_argument('--sample_name', type=str, required=True, help='Name of the sample.')
    parser.add_argument('--output_feather_path', type=str, required=True,
                        help='Path to save output gene marker score feather file.')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts').")

    # no_expression_fraction
    parser.add_argument('--no_expression_fraction', action='store_true', default=False,
                        help='Flag to not use expression fraction as filter when calculate the maker score. Default is False.')

    parser.add_argument('--latent_representation', type=str, default='latent_GVAE',
                        choices=['latent_GVAE', 'latent_PCA'],
                        help='Type of latent representation. Default is "latent_GVAE".')
    parser.add_argument('--num_neighbour', type=int, default=21,
                        help='Number of neighbours to consider. Default is 21.')
    parser.add_argument('--num_neighbour_spatial', type=int, default=101,
                        help='Number of spatial neighbours to consider. Default is 101.')
    parser.add_argument('--species', type=str, default=None,
                        help='Species name, which is the column name in the homolog gene conversion file.')
    parser.add_argument('--gs_species', type=str, default=None,
                        help='Homologous gene conversion file path, if applicable.')
    parser.add_argument('--gM_slices', type=str, default=None,
                        help='The mean Gene marker scores path across multiple slices.')


def add_spatial_ldsc_args(parser):
    # Group for GWAS input data
    parser.add_argument('--sample_name', required=True, help="Name of the spatial transcriptomic dataset.")

    parser.add_argument('--sumstats_file', default=None, help="Path to GWAS summary statistics file.")
    parser.add_argument('--sumstats_config_file', default=None, help="Path to GWAS summary statistics config file.")
    parser.add_argument('--w_file', required=True, help="Path to regression weight file.")
    parser.add_argument('--ldscore_input_dir', required=True, help="Input directory for LD Score files.")
    parser.add_argument('--ldsc_save_dir', required=True, help="Directory to save Spatial LDSC results.")
    parser.add_argument('--trait_name', default=None, help="Name of the trait.")
    parser.add_argument('--not_M_5_50', action='store_true', help="Flag to not use M 5 50 in calculations.")
    parser.add_argument('--n_blocks', type=int, default=200, help="Number of blocks for jackknife resampling.")
    parser.add_argument('--chisq_max', type=int, help="Maximum chi-square value for filtering SNPs.")
    parser.add_argument('--all_chunk', type=int, help="Number of chunks for processing spatial data.")

    # if use additional baseline annotation
    parser.add_argument('--disable_additional_baseline_annotation', action='store_true', default=False, )

    parser.add_argument('--num_processes', type=int, default=4, help="Number of processes for parallel computing.")

    return parser


def add_Cauchy_combination_args(parser):
    # Required arguments
    parser.add_argument('--input_hdf5_path', required=True, type=str, help='Path to the HDF5 file')
    parser.add_argument('--input_ldsc_dir', required=True, type=str, help='Directory containing LDSC results')
    parser.add_argument('--output_cauchy_dir', required=True, type=str,
                        help='Output directory for Cauchy combination results')
    parser.add_argument('--sample_name', required=True, type=str, help='Name of the sample')
    parser.add_argument('--trait_name', required=True, type=str, help='Name of the trait')
    parser.add_argument('--annotation', required=True, type=str, help='Annotation layer name')

    # Optional arguments
    parser.add_argument('--meta', default=None, type=str, )
    parser.add_argument('--slide', default=None, type=str, )


def add_Visualization_args(parser):
    # Required arguments
    parser.add_argument('--input_hdf5_path', required=True, type=str, help='Path to the HDF5 file')
    parser.add_argument('--input_ldsc_dir', required=True, type=str, help='Directory containing LDSC results')
    parser.add_argument('--output_figure_dir', required=True, type=str, help='Output directory for figures')
    parser.add_argument('--sample_name', required=True, type=str, help='Name of the sample')
    parser.add_argument('--trait_name', required=True, type=str, help='Name of the trait')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')

    # Arguments with defaults
    parser.add_argument('--fig_title', type=str, default=None, help='Title of figure')
    parser.add_argument('--fig_height', type=int, default=600, help='Height of figure in pixels')
    parser.add_argument('--fig_width', type=int, default=800, help='Width of figure in pixels')
    parser.add_argument('--point_size', type=int, default=None, help='Point size of figure')
    parser.add_argument('--fig_style', type=str, default='light', choices=['dark', 'light'],
                        help='Plot style of figure')


def add_all_mode_args(parser):
    parser.add_argument('--input_hdf5_path', required=True, type=str, help='Path to the input hdf5 file.')
    parser.add_argument('--save_dir', required=True, type=str, help='Path to the running results.')
    # output
    # parser.add_argument('--output_hdf5_path', required=True, type=str, help='Path to the output hdf5 file.')
    parser.add_argument('--sample_name', required=True, type=str, help='Name of the sample.')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    parser.add_argument('--type', default=None, type=str,
                        help="The gene expression layer of the input data (e.g., 'count', 'counts', 'log1p').")

    # latent_to_gene
    # input
    # parser.add_argument('--input_hdf5_path', type=str, required=True, help='Path to the input HDF5 file.')
    # parser.add_argument('--sample_name', type=str, required=True, help='Name of the sample.')
    # output
    # parser.add_argument('--output_feather_path', type=str, required=True,
    #                     help='Path to save output gene marker score feather file.')
    # parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    # parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts').")

    # no_expression_fraction
    # no_expression_fraction
    parser.add_argument('--no_expression_fraction', action='store_true', default=False,
                        help='Flag to not use expression fraction as filter when calculate the maker score. Default is False.')

    parser.add_argument('--latent_representation', type=str, default='latent_GVAE',
                        choices=['latent_GVAE', 'latent_PCA'],
                        help='Type of latent representation. Default is "latent_GVAE".')
    parser.add_argument('--num_neighbour', type=int, default=21,
                        help='Number of neighbours to consider. Default is 21.')
    parser.add_argument('--num_neighbour_spatial', type=int, default=101,
                        help='Number of spatial neighbours to consider. Default is 101.')
    parser.add_argument('--species', type=str, default=None,
                        help='Species name, which is the column name in the homolog gene conversion file.')
    parser.add_argument('--gs_species', type=str, default=None,
                        help='Homologous gene conversion file path, if applicable.')
    parser.add_argument('--gM_slices', type=str, default=None,
                        help='The mean Gene marker scores path across multiple slices.')

    # generate_ldscore
    # parser.add_argument('--sample_name', type=str, required=True, help='Sample name')
    # should be all
    # parser.add_argument('--chrom', type=chrom_choice, required=True, help='Chromosome number (1-22) or "all"')
    # output
    # parser.add_argument('--ldscore_save_dir', type=str, required=True, help='Directory to save ld score files')

    # reference panel
    parser.add_argument('--bfile_root', type=str, required=True, help='Bfile root path')
    parser.add_argument('--keep_snp_root', type=str, required=True, help='Keep SNP root path')

    # Annotation by gene distance
    parser.add_argument('--gtf_annotation_file', type=str, required=True, help='GTF file path')
    parser.add_argument('--gene_window_size', type=int, default=50000, help='Gene window size')

    # Enhancer annotation
    parser.add_argument('--enhancer_annotation_file', type=str, default=None,
                        help='Enhancer annotation bed file path, optional.')
    parser.add_argument('--snp_multiple_enhancer_strategy', type=str, default='max_mkscore',
                        choices=['max_mkscore', 'nearest_TSS'], help='Strategy for multiple enhancers per SNP')
    parser.add_argument('--gene_window_enhancer_priority', type=str, default=None,
                        choices=['gene_window_first', 'enhancer_first', 'enhancer_only'],
                        help='Priority between gene window and enhancer')

    # Arguments for calculating ld score
    parser.add_argument('--spots_per_chunk', type=int, default=5_000, help='Number of spots per chunk')
    parser.add_argument('--ld_wind', type=int, default=1, help='LD window size')
    parser.add_argument('--ld_unit', type=str, default='CM', help='LD window unit (SNP/KB/CM)',
                        choices=['SNP', 'KB', 'CM'])

    # spatial ldsc args:
    parser.add_argument('--sumstats_file', default=None, help="Path to GWAS summary statistics file.")
    parser.add_argument('--sumstats_config_file', default=None, help="Path to GWAS summary statistics config file.")
    parser.add_argument('--w_file', required=True, help="Path to regression weight file.")
    parser.add_argument('--ldscore_input_dir', required=True, help="Input directory for LD Score files.")
    parser.add_argument('--ldsc_save_dir', required=True, help="Directory to save Spatial LDSC results.")
    parser.add_argument('--trait_name', default=None, help="Name of the trait.")
    parser.add_argument('--not_M_5_50', action='store_true', help="Flag to not use M 5 50 in calculations.")
    parser.add_argument('--n_blocks', type=int, default=200, help="Number of blocks for jackknife resampling.")
    parser.add_argument('--chisq_max', type=int, help="Maximum chi-square value for filtering SNPs.")
    parser.add_argument('--all_chunk', type=int, help="Number of chunks for processing spatial data.")


def get_runall_mode_config(args: argparse.Namespace):
    # output
    args.output_hdf5_path = f'{args.save_dir}/{args.sample_name}/find_latent_representations/{args.sample_name}_add_latent.h5ad'
    args.output_feather_path = f'{args.save_dir}/{args.sample_name}/latent_to_gene/{args.sample_name}_gene_marker_score.feather'
    args.ldscore_save_dir = f'{args.save_dir}/{args.sample_name}/generate_ldscore'
    args.ldsc_save_dir = f'{args.save_dir}/{args.sample_name}/spatial_ldsc'
    args.output_cauchy_dir = f'{args.save_dir}/{args.sample_name}/cauchy_combination/'

    # input
    args.input_hdf5_with_latent_path = args.output_hdf5_path
    args.mkscore_feather_file = args.output_feather_path
    args.ldscore_input_dir = args.ldscore_save_dir
    args.chrom = 'all'
    args.input_ldsc_dir = args.ldsc_save_dir
    args.input_spatial_ldsc = f'{args.save_dir}/{args.sample_name}/spatial_ldsc/{args.sample_name}_{args.trait_name}.gz'
    # find_latent_representations
    flr_config = get_dataclass_from_parser(args, FindLatentRepresentationsConfig)
    # latent_to_gene
    ltg_config = get_dataclass_from_parser(args, LatentToGeneConfig)
    # generate_ldscore
    gls_config = get_dataclass_from_parser(args, GenerateLDScoreConfig)
    # spatial ldsc
    ldsc_config = get_dataclass_from_parser(args, SpatialLDSCConfig)
    # cauchy combination
    cauchy_config = get_dataclass_from_parser(args, CauchyCombinationConfig)
    return RunAllModeConfig(flr_config=flr_config, ltg_config=ltg_config, gls_config=gls_config,
                            ldsc_config=ldsc_config, cauchy_config=cauchy_config)


def add_format_sumstats_args(parser):
    # Required arguments
    parser.add_argument('--sumstats', required=True, type=str,
                        help='Path to gwas summary data')
    parser.add_argument('--out', required=True, type=str,
                        help='Path to save the formatted gwas data')

    # Arguments for specify column name
    parser.add_argument('--snp', default=None, type=str,
                        help="Name of snp column (if not a name that gsMap understands)")
    parser.add_argument('--a1', default=None, type=str,
                        help="Name of effect allele column (if not a name that gsMap understands)")
    parser.add_argument('--a2', default=None, type=str,
                        help="Name of none-effect allele column (if not a name that gsMap understands)")
    parser.add_argument('--info', default=None, type=str,
                        help="Name of info column (if not a name that gsMap understands)")
    parser.add_argument('--beta', default=None, type=str,
                        help="Name of gwas beta column (if not a name that gsMap understands).")
    parser.add_argument('--se', default=None, type=str,
                        help="Name of gwas standar error of beta column (if not a name that gsMap understands)")
    parser.add_argument('--p', default=None, type=str,
                        help="Name of p-value column (if not a name that gsMap understands)")
    parser.add_argument('--frq', default=None, type=str,
                        help="Name of A1 ferquency column (if not a name that gsMap understands)")
    parser.add_argument('--n', default=None, type=str,
                        help="Name of sample size column (if not a name that gsMap understands)")
    parser.add_argument('--z', default=None, type=str,
                        help="Name of gwas Z-statistics column (if not a name that gsMap understands)")
    parser.add_argument('--OR', default=None, type=str,
                        help="Name of gwas OR column (if not a name that gsMap understands)")
    parser.add_argument('--se_OR', default=None, type=str,
                        help="Name of standar error of OR column (if not a name that gsMap understands)")

    # Arguments for convert SNP (chr, pos) to rsid
    parser.add_argument('--chr', default="Chr", type=str,
                        help="Name of SNP chromosome column (if not a name that gsMap understands)")
    parser.add_argument('--pos', default="Pos", type=str,
                        help="Name of SNP positions column (if not a name that gsMap understands)")
    parser.add_argument('--dbsnp', default=None, type=str,
                        help='Path to reference dnsnp file')
    parser.add_argument('--chunksize', default=1e+6, type=int,
                        help='Chunk size for loading dbsnp file')

    # Arguments for output format and quality
    parser.add_argument('--format', default='gsMap', type=str,
                        help='Format of output data', choices=['gsMap', 'COJO'])
    parser.add_argument('--info_min', default=0.9, type=float,
                        help='Minimum INFO score.')
    parser.add_argument('--maf_min', default=0.01, type=float,
                        help='Minimum MAF.')
    parser.add_argument('--keep_chr_pos', action='store_true', default=False,
                        help='Keep SNP chromosome and position columns in the output data')


def ensure_path_exists(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        if isinstance(result, Path):
            if result.suffix:
                result.parent.mkdir(parents=True, exist_ok=True, mode=0o755)
            else:  # It's a directory path
                result.mkdir(parents=True, exist_ok=True, mode=0o755)
        return result

    return wrapper


@dataclass
class ConfigWithAutoPaths:
    workdir: str
    sample_name: str

    def __post_init__(self):
        if self.workdir is None:
            raise ValueError('workdir must be provided.')

    @property
    @ensure_path_exists
    def hdf5_with_latent_path(self) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/find_latent_representations/{self.sample_name}_add_latent.h5ad')

    @property
    @ensure_path_exists
    def mkscore_feather_path(self) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/latent_to_gene/{self.sample_name}_gene_marker_score.feather')

    @property
    @ensure_path_exists
    def ldscore_save_dir(self) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/generate_ldscore')

    @property
    @ensure_path_exists
    def ldsc_save_dir(self) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/spatial_ldsc')

    @property
    @ensure_path_exists
    def cauchy_save_dir(self) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/cauchy_combination')

    @ensure_path_exists
    def get_report_dir(self, trait_name: str) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/report/{trait_name}')

    def get_gsMap_report_file(self, trait_name: str) -> Path:
        return self.get_report_dir(trait_name) / f'{self.sample_name}_{trait_name}_gsMap_Report.html'

    @ensure_path_exists
    def get_manhattan_html_plot_path(self, trait_name: str) -> Path:
        return Path(
            f'{self.workdir}/{self.sample_name}/report/{trait_name}/manhattan_plot/{self.sample_name}_{trait_name}_Diagnostic_Manhattan_Plot.html')

    @ensure_path_exists
    def get_GSS_plot_dir(self, trait_name: str) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/report/{trait_name}/GSS_plot')

    def get_GSS_plot_select_gene_file(self, trait_name: str) -> Path:
        return self.get_GSS_plot_dir(trait_name) / 'plot_genes.csv'

    @ensure_path_exists
    def get_ldsc_result_file(self, trait_name: str) -> Path:
        return Path(f'{self.ldsc_save_dir}/{self.sample_name}_{trait_name}.csv.gz')

    @ensure_path_exists
    def get_cauchy_result_file(self, trait_name: str) -> Path:
        return Path(f'{self.cauchy_save_dir}/{self.sample_name}_{trait_name}.Cauchy.csv.gz')

    @ensure_path_exists
    def get_gene_diagnostic_info_save_path(self, trait_name: str) -> Path:
        return Path(
            f'{self.workdir}/{self.sample_name}/report/{trait_name}/{self.sample_name}_{trait_name}_Gene_Diagnostic_Info.csv')

    @ensure_path_exists
    def get_gsMap_plot_save_dir(self, trait_name: str) -> Path:
        return Path(f'{self.workdir}/{self.sample_name}/report/{trait_name}/gsMap_plot')

    def get_gsMap_html_plot_save_path(self, trait_name: str) -> Path:
        return self.get_gsMap_plot_save_dir(trait_name) / f'{self.sample_name}_{trait_name}_gsMap_plot.html'


@dataclass
class FindLatentRepresentationsConfig(ConfigWithAutoPaths):
    input_hdf5_path: str
    # output_hdf5_path: str
    annotation: str = None
    data_layer: str = None

    epochs: int = 300
    feat_hidden1: int = 256
    feat_hidden2: int = 128
    feat_cell: int = 3000
    gcn_hidden1: int = 64
    gcn_hidden2: int = 30
    p_drop: float = 0.1
    gcn_lr: float = 0.001
    gcn_decay: float = 0.01
    n_neighbors: int = 11
    label_w: float = 1
    rec_w: float = 1
    input_pca: bool = True
    n_comps: int = 300
    weighted_adj: bool = False
    nheads: int = 3
    var: bool = False
    convergence_threshold: float = 1e-4
    hierarchically: bool = False

    def __post_init__(self):
        # self.output_hdf5_path = self.hdf5_with_latent_path
        if self.hierarchically:
            if self.annotation is None:
                raise ValueError('annotation must be provided if hierarchically is True.')
            logger.info(
                f'------Hierarchical mode is enabled. This will find the latent representations within each annotation.')

        # remind for not providing annotation
        if self.annotation is None:
            logger.warning(
                'annotation is not provided. This will find the latent representations for the whole dataset.')
        else:
            logger.info(f'------Find latent representations for {self.annotation}...')


@dataclass
class LatentToGeneConfig(ConfigWithAutoPaths):
    # input_hdf5_with_latent_path: str
    # output_feather_path: str
    no_expression_fraction: bool = False
    latent_representation: str = 'latent_GVAE'
    num_neighbour: int = 21
    num_neighbour_spatial: int = 101
    species: str = None
    homolog_file: str = None
    gM_slices: str = None
    annotation: str = None
    type: str = None

    def __post_init__(self):
        if self.homolog_file is not None:
            logger.info(f"User provided homolog file to map gene names to human: {self.homolog_file}")
            # check the format of the homolog file
            with open(self.homolog_file, 'r') as f:
                first_line = f.readline().strip()
                _n_col = len(first_line.split())
                if _n_col != 2:
                    raise ValueError(
                        f"Invalid homolog file format. Expected 2 columns, first column should be other species gene name, second column should be human gene name. "
                        f"Got {_n_col} columns in the first line.")
                else:
                    first_col_name, second_col_name = first_line.split()
                    self.species = first_col_name
                    logger.info(
                        f"Homolog file provided and will map gene name from column1:{first_col_name} to column2:{second_col_name}")
        else:
            logger.info("No homolog file provided. Run in human mode.")


@dataclass
class GenerateLDScoreConfig(ConfigWithAutoPaths):
    chrom: Union[int, str]

    bfile_root: str
    keep_snp_root: Optional[str]

    # annotation by gene distance
    gtf_annotation_file: str
    gene_window_size: int = 50000

    # annotation by enhancer
    enhancer_annotation_file: str = None
    snp_multiple_enhancer_strategy: Literal['max_mkscore', 'nearest_TSS'] = 'max_mkscore'
    gene_window_enhancer_priority: Optional[Literal['gene_window_first', 'enhancer_first', 'enhancer_only',]] = None

    # for calculating ld score
    additional_baseline_annotation_dir_path: str = None
    spots_per_chunk: int = 5_000
    ld_wind: int = 1
    ld_unit: str = 'CM'

    # zarr config
    ldscore_save_format: Literal['feather', 'zarr', 'quick_mode'] = 'feather'

    zarr_chunk_size: Tuple[int, int] = None

    # for pre calculating the SNP Gene ldscore Weight
    save_pre_calculate_snp_gene_weight_matrix: bool = False

    baseline_annotation_dir: Optional[str] = None
    SNP_gene_pair_dir: Optional[str] = None
    def __post_init__(self):
        # if self.mkscore_feather_file is None:
        #     self.mkscore_feather_file = self._get_mkscore_feather_path()

        if self.enhancer_annotation_file is not None and self.gene_window_enhancer_priority is None:
            logger.warning("enhancer_annotation_file is provided but gene_window_enhancer_priority is not provided. "
                           "by default, gene_window_enhancer_priority is set to 'enhancer_only', when enhancer_annotation_file is provided.")
            self.gene_window_enhancer_priority = 'enhancer_only'
        if self.enhancer_annotation_file is None and self.gene_window_enhancer_priority is not None:
            logger.warning("gene_window_enhancer_priority is provided but enhancer_annotation_file is not provided. "
                           "by default, gene_window_enhancer_priority is set to None, when enhancer_annotation_file is not provided.")
            self.gene_window_enhancer_priority = None
        assert self.gene_window_enhancer_priority in [None, 'gene_window_first', 'enhancer_first', 'enhancer_only', ], \
            f"gene_window_enhancer_priority must be one of None, 'gene_window_first', 'enhancer_first', 'enhancer_only', but got {self.gene_window_enhancer_priority}."
        if self.gene_window_enhancer_priority in ['gene_window_first', 'enhancer_first']:
            logger.info(f'Both gene_window and enhancer annotation will be used to calculate LD score. ')
            logger.info(
                f'SNP within +-{self.gene_window_size} bp of gene body will be used and enhancer annotation will be used to calculate LD score. If a snp maps to multiple enhancers, the strategy to choose by your select strategy: {self.snp_multiple_enhancer_strategy}.')
        elif self.gene_window_enhancer_priority == 'enhancer_only':
            logger.info(f'Only enhancer annotation will be used to calculate LD score. ')
        else:
            logger.info(
                f'Only gene window annotation will be used to calculate LD score. SNP within +-{self.gene_window_size} bp of gene body will be used. ')

        # remind for baseline annotation
        if self.additional_baseline_annotation_dir_path is None:
            logger.info(f'------Baseline annotation is not provided. Default baseline annotation will be used.')
        else:
            logger.info(
                f'------Baseline annotation is provided. Additional baseline annotation will be used with the default baseline annotation.')
            logger.info(f'------Baseline annotation directory: {self.additional_baseline_annotation_dir_path}')
            # check the existence of baseline annotation
            if self.chrom == 'all':
                for chrom in range(1, 23):
                    chrom = str(chrom)
                    baseline_annotation_path = Path(
                        self.additional_baseline_annotation_dir_path) / f'baseline.{chrom}.annot.gz'
                    if not baseline_annotation_path.exists():
                        raise FileNotFoundError(
                            f'baseline.{chrom}.annot.gz is not found in {self.additional_baseline_annotation_dir_path}.')
            else:
                baseline_annotation_path = Path(
                    self.additional_baseline_annotation_dir_path) / f'baseline.{self.chrom}.annot.gz'
                if not baseline_annotation_path.exists():
                    raise FileNotFoundError(
                        f'baseline.{self.chrom}.annot.gz is not found in {self.additional_baseline_annotation_dir_path}.')

        # set the default zarr chunk size
        if self.ldscore_save_format == 'zarr' and self.zarr_chunk_size is None:
            self.zarr_chunk_size = (10_000, self.spots_per_chunk)


@dataclass
class SpatialLDSCConfig(ConfigWithAutoPaths):
    w_file: str
    # ldscore_save_dir: str
    disable_additional_baseline_annotation: bool = False
    trait_name: Optional[str] = None
    sumstats_file: Optional[str] = None
    sumstats_config_file: Optional[str] = None
    num_processes: int = 4
    not_M_5_50: bool = False
    n_blocks: int = 200
    chisq_max: Optional[int] = None
    all_chunk: Optional[int] = None
    chunk_range: Optional[Tuple[int, int]] = None

    ldscore_save_format: Literal['feather', 'zarr', 'quick_mode'] = 'feather'

    spots_per_chunk_quick_mode: int = 1_000
    snp_gene_weight_adata_path: Optional[str] = None

    def __post_init__(self):
        super().__post_init__()
        if self.sumstats_file is None and self.sumstats_config_file is None:
            raise ValueError('One of sumstats_file and sumstats_config_file must be provided.')
        if self.sumstats_file is not None and self.sumstats_config_file is not None:
            raise ValueError('Only one of sumstats_file and sumstats_config_file must be provided.')
        if self.sumstats_file is not None and self.trait_name is None:
            raise ValueError('trait_name must be provided if sumstats_file is provided.')
        if self.sumstats_config_file is not None and self.trait_name is not None:
            raise ValueError('trait_name must not be provided if sumstats_config_file is provided.')
        self.sumstats_config_dict = {}
        # load the sumstats config file
        if self.sumstats_config_file is not None:
            import yaml
            with open(self.sumstats_config_file) as f:
                config = yaml.load(f, Loader=yaml.FullLoader)
            for trait_name, sumstats_file in config.items():
                assert Path(sumstats_file).exists(), f'{sumstats_file} does not exist.'
        # load the sumstats file
        elif self.sumstats_file is not None:
            self.sumstats_config_dict[self.trait_name] = self.sumstats_file
        else:
            raise ValueError('One of sumstats_file and sumstats_config_file must be provided.')

        for sumstats_file in self.sumstats_config_dict.values():
            assert Path(sumstats_file).exists(), f'{sumstats_file} does not exist.'

        # check if additional baseline annotation is exist
        self.use_additional_baseline_annotation = False
        self.process_additional_baseline_annotation()

    def process_additional_baseline_annotation(self):
        additional_baseline_annotation_dir_path = Path(self.ldscore_save_dir) / 'additional_baseline'
        dir_exists = additional_baseline_annotation_dir_path.exists()

        if not dir_exists:
            if self.use_additional_baseline_annotation:
                logger.warning(f"additional_baseline directory is not found in {self.ldscore_save_dir}.")
                print('''\
                    if you want to use additional baseline annotation, 
                    please provide additional baseline annotation when calculating ld score.
                    ''')
                raise FileNotFoundError(
                    f'additional_baseline directory is not found. You should disable use_additional_baseline_annotation')
            return

        self.use_additional_baseline_annotation = self.use_additional_baseline_annotation or True

        if self.disable_additional_baseline_annotation:
            logger.warning(
                f"additional_baseline directory is found in {self.ldscore_save_dir}, but use_additional_baseline_annotation is disabled.")
            print('''\
                if you want to use additional baseline annotation,
                please enable by not adding --disable_additional_baseline_annotation.
                ''')
            self.use_additional_baseline_annotation = False
        else:
            logger.info(
                f'------Additional baseline annotation is provided. It will be used with the default baseline annotation.')
            logger.info(f'------Additional baseline annotation directory: {additional_baseline_annotation_dir_path}')

            chrom_list = range(1, 23)
            for chrom in chrom_list:
                baseline_annotation_path = additional_baseline_annotation_dir_path / f'baseline.{chrom}.l2.ldscore.feather'
                if not baseline_annotation_path.exists():
                    raise FileNotFoundError(
                        f'baseline.{chrom}.annot.gz is not found in {additional_baseline_annotation_dir_path}.')


@dataclass
class CauchyCombinationConfig(ConfigWithAutoPaths):
    # input_hdf5_path: str
    # input_ldsc_dir: str
    # output_cauchy_dir: str
    # sample_name: str
    trait_name: str
    annotation: str
    meta: str = None
    slide: str = None


@dataclass
class VisualizeConfig(ConfigWithAutoPaths):
    # input_hdf5_path: str
    # input_ldsc_dir: str
    # output_figure_dir: str
    # sample_name: str
    trait_name: str

    annotation: str = None
    fig_title: str = None
    fig_height: int = 600
    fig_width: int = 800
    point_size: int = None
    fig_style: Literal['dark', 'light'] = 'light'


@dataclass
class DiagnosisConfig(ConfigWithAutoPaths):
    annotation: str
    # mkscore_feather_file: str

    trait_name: str
    sumstats_file: str
    plot_type: Literal['manhattan', 'GSS', 'gsMap', 'all'] = 'all'
    top_corr_genes: int = 50
    selected_genes: Optional[List[str]] = None


@dataclass
class ReportConfig(DiagnosisConfig):
    pass


@dataclass
class RunAllModeConfig(ConfigWithAutoPaths):
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

    max_processes: int = 10

    def __post_init__(self):
        super().__post_init__()
        self.gtffile = f"{self.gsMap_resource_dir}/genome_annotation/gtf/gencode.v39lift37.annotation.gtf"
        self.bfile_root = f"{self.gsMap_resource_dir}/LD_Reference_Panel/1000G_EUR_Phase3_plink/1000G.EUR.QC"
        self.keep_snp_root = f"{self.gsMap_resource_dir}/LDSC_resource/hapmap3_snps/hm"
        self.w_file = f"{self.gsMap_resource_dir}/LDSC_resource/weights_hm3_no_hla/weights."
        self.snp_gene_weight_adata_path = f"{self.gsMap_resource_dir}/LDSC_resource/quick_mode/snp_gene_weight_matrix.h5ad"
        self.baseline_annotation_dir = f"{self.gsMap_resource_dir}/LDSC_resource/quick_mode/baseline"
        self.SNP_gene_pair_dir = f"{self.gsMap_resource_dir}/LDSC_resource/quick_mode/SNP_gene_pair"
        # check the existence of the input files and resources files
        for file in [self.hdf5_path, self.gtffile]:
            if not Path(file).exists():
                raise FileNotFoundError(f"File {file} does not exist.")

        if self.sumstats_file is None and self.sumstats_config_file is None:
            raise ValueError('One of sumstats_file and sumstats_config_file must be provided.')
        if self.sumstats_file is not None and self.sumstats_config_file is not None:
            raise ValueError('Only one of sumstats_file and sumstats_config_file must be provided.')
        if self.sumstats_file is not None and self.trait_name is None:
            raise ValueError('trait_name must be provided if sumstats_file is provided.')
        if self.sumstats_config_file is not None and self.trait_name is not None:
            raise ValueError('trait_name must not be provided if sumstats_config_file is provided.')
        self.sumstats_config_dict = {}
        # load the sumstats config file
        if self.sumstats_config_file is not None:
            import yaml
            with open(self.sumstats_config_file) as f:
                config = yaml.load(f, Loader=yaml.FullLoader)
            for trait_name, sumstats_file in config.items():
                assert Path(sumstats_file).exists(), f'{sumstats_file} does not exist.'
        # load the sumstats file
        elif self.sumstats_file is not None and self.trait_name is not None:
            self.sumstats_config_dict[self.trait_name] = self.sumstats_file
        else:
            raise ValueError('One of sumstats_file and sumstats_config_file must be provided.')

        for sumstats_file in self.sumstats_config_dict.values():
            assert Path(sumstats_file).exists(), f'{sumstats_file} does not exist.'


@dataclass
class FormatSumstatsConfig:
    sumstats: str
    out: str
    dbsnp: str
    snp: str = None
    a1: str = None
    a2: str = None
    info: str = None
    beta: str = None
    se: str = None
    p: str = None
    frq: str = None
    n: str = None
    z: str = None
    OR: str = None
    se_OR: str = None
    format: str = None
    chr: str = None
    pos: str = None
    chunksize: int = 1e+7
    info_min: float = 0.9
    maf_min: float = 0.01
    keep_chr_pos: bool = False


@register_cli(name='run_find_latent_representations',
              description='Run Find_latent_representations \nFind the latent representations of each spot by running GNN-VAE',
              add_args_function=add_find_latent_representations_args)
def run_find_latent_representation_from_cli(args: argparse.Namespace):
    from gsMap.find_latent_representation import run_find_latent_representation
    config = get_dataclass_from_parser(args, FindLatentRepresentationsConfig)
    run_find_latent_representation(config)


@register_cli(name='run_latent_to_gene',
              description='Run Latent_to_gene \nEstimate gene marker gene scores for each spot by using latent representations from nearby spots',
              add_args_function=add_latent_to_gene_args)
def run_latent_to_gene_from_cli(args: argparse.Namespace):
    from gsMap.latent_to_gene import run_latent_to_gene
    config = get_dataclass_from_parser(args, LatentToGeneConfig)
    run_latent_to_gene(config)


@register_cli(name='run_generate_ldscore',
              description='Run Generate_ldscore \nGenerate LD scores for each spot',
              add_args_function=add_generate_ldscore_args)
def run_generate_ldscore_from_cli(args: argparse.Namespace):
    from gsMap.generate_ldscore import run_generate_ldscore
    config = get_dataclass_from_parser(args, GenerateLDScoreConfig)
    run_generate_ldscore(config)


@register_cli(name='run_spatial_ldsc',
              description='Run Spatial_ldsc \nRun spatial LDSC for each spot',
              add_args_function=add_spatial_ldsc_args)
def run_spatial_ldsc_from_cli(args: argparse.Namespace):
    from gsMap.spatial_ldsc_multiple_sumstats import run_spatial_ldsc
    config = get_dataclass_from_parser(args, SpatialLDSCConfig)
    run_spatial_ldsc(config)


@register_cli(name='run_cauchy_combination',
              description='Run Cauchy_combination for each annotation',
              add_args_function=add_Cauchy_combination_args)
def run_Cauchy_combination_from_cli(args: argparse.Namespace):
    from gsMap.cauchy_combination_test import run_Cauchy_combination
    config = get_dataclass_from_parser(args, CauchyCombinationConfig)
    run_Cauchy_combination(config)


@register_cli(name='run_visualize',
              description='Visualize the gsMap results',
              add_args_function=add_Visualization_args)
def run_Visualize_from_cli(args: argparse.Namespace):
    from gsMap.visualize import run_Visualize
    config = get_dataclass_from_parser(args, VisualizeConfig)
    run_Visualize(config)


@register_cli(name='run_all_mode',
              description='Run gsMap method (the full process)',
              add_args_function=add_all_mode_args)
def run_all_mode_from_cli(args: argparse.Namespace):
    from gsMap.find_latent_representation import run_find_latent_representation
    from gsMap.latent_to_gene import run_latent_to_gene
    from gsMap.generate_ldscore import run_generate_ldscore
    from gsMap.spatial_ldsc_multiple_sumstats import run_spatial_ldsc
    from gsMap.cauchy_combination_test import run_Cauchy_combination
    config = get_runall_mode_config(args)
    run_find_latent_representation(config.flr_config)
    run_latent_to_gene(config.ltg_config)
    run_generate_ldscore(config.gls_config)
    run_spatial_ldsc(config.ldsc_config)
    if args.annotation is not None:
        config.cauchy_config.annotation = args.annotation
        run_Cauchy_combination(config.cauchy_config)


@register_cli(name='format_sumstats',
              description='Format gwas summary statistics',
              add_args_function=add_format_sumstats_args)
def gwas_format_from_cli(args: argparse.Namespace):
    from gsMap.format_sumstats import gwas_format
    config = get_dataclass_from_parser(args, FormatSumstatsConfig)
    gwas_format(config)

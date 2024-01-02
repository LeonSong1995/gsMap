import argparse
from dataclasses import dataclass
from typing import Union

from collections import OrderedDict, namedtuple
from typing import Callable
import pyfiglet
# Global registry to hold functions
cli_function_registry = OrderedDict()
subcommand = namedtuple('subcommand', ['name', 'func', 'add_args_function', 'description'])


# Decorator to register functions for cli parsing
def register_cli(name: str, description: str, add_args_function: Callable) -> Callable:
    def decorator(func: Callable) -> Callable:
        cli_function_registry[name] = subcommand(name=name, func=func, add_args_function=add_args_function,
                                                 description=description)
        def wrapper(*args, **kwargs):
            print(pyfiglet.figlet_format('GPS'))
            print(pyfiglet.figlet_format('Genetics-informed pathogenic spatial mapping',
                                         font='contessa'))

            print(pyfiglet.figlet_format(name))

            print(f'Running {name}...')
            return func(*args, **kwargs)
        return wrapper

    return decorator


def add_find_latent_representations_args(parser):
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
    parser.add_argument('--input_pca', default=True, type=bool,
                        help="Whether to perform PCA on input features. Default is True.")
    parser.add_argument('--n_comps', default=300, type=int,
                        help="Number of principal components to keep if PCA is performed. Default is 300.")
    parser.add_argument('--weighted_adj', default=False, type=bool,
                        help="Whether to use a weighted adjacency matrix in GCN. Default is False.")
    parser.add_argument('--nheads', default=3, type=int,
                        help="Number of heads in the attention mechanism of the GNN. Default is 3.")
    parser.add_argument('--var', default=False, type=bool)
    parser.add_argument('--convergence_threshold', default=1e-4, type=float,
                        help="Threshold for convergence during training. Training stops if the loss change is below this threshold. Default is 1e-4.")
    parser.add_argument('--hierarchically', default=False, type=bool,
                        help="Whether to find latent representations hierarchically. Default is False.")

    parser.add_argument('--input_hdf5_path', required=True, type=str, help='Path to the input hdf5 file.')
    parser.add_argument('--output_hdf5_path', required=True, type=str, help='Path to the output hdf5 file.')
    parser.add_argument('--sample_name', required=True, type=str, help='Name of the sample.')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts').")


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
def get_dataclass_from_parser(parser, data_class: dataclass):
    return data_class(**filter_args_for_dataclass(vars(parser), data_class))

def add_generate_ldscore_args(parser):
    parser.add_argument('--sample_name', type=str, required=True, help='Sample name')
    parser.add_argument('--chrom', type=chrom_choice, required=True, help='Chromosome number (1-22) or "all"')
    parser.add_argument('--ldscore_save_dir', type=str, required=True, help='Directory to save ld score files')
    parser.add_argument('--gtf_file', type=str, required=True, help='GTF file path')
    parser.add_argument('--mkscore_feather_file', type=str, required=True, help='Mkscore feather file path')
    parser.add_argument('--bfile_root', type=str, required=True, help='Bfile root path')
    parser.add_argument('--keep_snp_root', type=str, required=True, help='Keep SNP root path')

    # Arguments with defaults
    parser.add_argument('--window_size', type=int, default=50000, help='Annotation window size for each gene')
    parser.add_argument('--spots_per_chunk', type=int, default=10000, help='Number of spots per chunk')
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

    parser.add_argument('--method', type=str, default='rank', choices=['rank', 'other_method'],
                        help='Method to be used. Default is "rank".')
    parser.add_argument('--latent_representation', type=str, default='latent_GVAE',
                        choices=['latent_GVAE', 'latent_PCA'],
                        help='Type of latent representation. Default is "latent_GVAE".')
    parser.add_argument('--num_neighbour', type=int, default=21,
                        help='Number of neighbours to consider. Default is 21.')
    parser.add_argument('--num_neighbour_spatial', type=int, default=101,
                        help='Number of spatial neighbours to consider. Default is 101.')
    parser.add_argument('--num_processes', type=int, default=4, help='Number of processes to use. Default is 4.')
    parser.add_argument('--fold', type=float, default=1.0, help='Fold change threshold. Default is 1.0.')
    parser.add_argument('--pst', type=float, default=0.2, help='PST value. Default is 0.2.')
    parser.add_argument('--species', type=str, default=None, help='Species name, if applicable.')
    parser.add_argument('--gs_species', type=str, default=None, help='Gene species file path, if applicable.')
    parser.add_argument('--gM_slices', type=str, default=None, help='Path to gene model slices file, if applicable.')



def add_all_mode_args(parser):
    parser.add_argument('--input_hdf5_path', required=True, type=str, help='Path to the input hdf5 file.')
    parser.add_argument('--save_dir', required=True, type=str, help='Path to the running results.')
    # output
    # parser.add_argument('--output_hdf5_path', required=True, type=str, help='Path to the output hdf5 file.')
    parser.add_argument('--sample_name', required=True, type=str, help='Name of the sample.')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts').")

    # find_latent_representations
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
    parser.add_argument('--input_pca', default=True, type=bool,
                        help="Whether to perform PCA on input features. Default is True.")
    parser.add_argument('--n_comps', default=300, type=int,
                        help="Number of principal components to keep if PCA is performed. Default is 300.")
    parser.add_argument('--weighted_adj', default=False, type=bool,
                        help="Whether to use a weighted adjacency matrix in GCN. Default is False.")
    parser.add_argument('--nheads', default=3, type=int,
                        help="Number of heads in the attention mechanism of the GNN. Default is 3.")
    parser.add_argument('--var', default=False, type=bool)
    parser.add_argument('--convergence_threshold', default=1e-4, type=float,
                        help="Threshold for convergence during training. Training stops if the loss change is below this threshold. Default is 1e-4.")
    parser.add_argument('--hierarchically', default=False, type=bool,
                        help="Whether to find latent representations hierarchically. Default is False.")

    # latent_to_gene
    # input
    # parser.add_argument('--input_hdf5_path', type=str, required=True, help='Path to the input HDF5 file.')
    # parser.add_argument('--sample_name', type=str, required=True, help='Name of the sample.')
    # output
    # parser.add_argument('--output_feather_path', type=str, required=True,
    #                     help='Path to save output gene marker score feather file.')
    # parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    # parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts').")

    parser.add_argument('--method', type=str, default='rank', choices=['rank', 'other_method'],
                        help='Method to be used. Default is "rank".')
    parser.add_argument('--latent_representation', type=str, default='latent_GVAE',
                        choices=['latent_GVAE', 'latent_PCA'],
                        help='Type of latent representation. Default is "latent_GVAE".')
    parser.add_argument('--num_neighbour', type=int, default=21,
                        help='Number of neighbours to consider. Default is 21.')
    parser.add_argument('--num_neighbour_spatial', type=int, default=101,
                        help='Number of spatial neighbours to consider. Default is 101.')
    parser.add_argument('--num_processes', type=int, default=4, help='Number of processes to use. Default is 4.')
    parser.add_argument('--fold', type=float, default=1.0, help='Fold change threshold. Default is 1.0.')
    parser.add_argument('--pst', type=float, default=0.2, help='PST value. Default is 0.2.')
    parser.add_argument('--species', type=str, default=None, help='Species name, if applicable.')
    parser.add_argument('--gs_species', type=str, default=None, help='Gene species file path, if applicable.')
    parser.add_argument('--gM_slices', type=str, default=None, help='Path to gene model slices file, if applicable.')

    # generate_ldscore
    # parser.add_argument('--sample_name', type=str, required=True, help='Sample name')
    # should be all
    # parser.add_argument('--chrom', type=chrom_choice, required=True, help='Chromosome number (1-22) or "all"')
    # output
    # parser.add_argument('--ldscore_save_dir', type=str, required=True, help='Directory to save ld score files')
    parser.add_argument('--gtf_file', type=str, required=True, help='GTF file path')
    # input
    # parser.add_argument('--mkscore_feather_file', type=str, required=True, help='Mkscore feather file path')
    parser.add_argument('--bfile_root', type=str, required=True, help='Bfile root path')
    parser.add_argument('--keep_snp_root', type=str, required=True, help='Keep SNP root path')

    # Arguments with defaults
    parser.add_argument('--window_size', type=int, default=50000, help='Annotation window size for each gene')
    parser.add_argument('--spots_per_chunk', type=int, default=10000, help='Number of spots per chunk')
    parser.add_argument('--ld_wind', type=int, default=1, help='LD window size')
    parser.add_argument('--ld_unit', type=str, default='CM', help='LD window unit (SNP/KB/CM)',
                        choices=['SNP', 'KB', 'CM'])


def get_runall_mode_config(args: argparse.ArgumentParser):
    # output
    args.output_hdf5_path = f'{args.save_dir}/{args.sample_name}/find_latent_representations/{args.sample_name}_add_latent.h5ad'
    args.output_feather_path = f'{args.save_dir}/{args.sample_name}/latent_to_gene/{args.sample_name}_gene_marker_score.feather'
    args.ldscore_save_dir = f'{args.save_dir}/{args.sample_name}/generate_ldscore'

    # input
    args.input_hdf5_with_latent_path = args.output_hdf5_path
    args.mkscore_feather_file = args.output_feather_path
    args.chrom = 'all'

    # find_latent_representations
    flr_config = get_dataclass_from_parser(args, FindLatentRepresentationsConfig)
    # latent_to_gene
    ltg_config = get_dataclass_from_parser(args, LatentToGeneConfig)
    # generate_ldscore
    gls_config = get_dataclass_from_parser(args, GenerateLDScoreConfig)

    return RunAllModeConfig(flr_config=flr_config, ltg_config=ltg_config, gls_config=gls_config)


@dataclass
class GenerateLDScoreConfig:
    sample_name: str
    chrom: Union[int, str]
    ldscore_save_dir: str
    gtf_file: str
    mkscore_feather_file: str
    bfile_root: str
    keep_snp_root: str
    window_size: int = 50000
    spots_per_chunk: int = 10_000
    ld_wind: int = 1
    ld_unit: str = 'CM'


@dataclass
class FindLatentRepresentationsConfig:
    input_hdf5_path: str
    output_hdf5_path: str
    sample_name: str
    annotation: str = None
    type: str = None

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


@dataclass
class LatentToGeneConfig:
    input_hdf5_with_latent_path: str
    sample_name: str
    output_feather_path: str

    method: str = 'rank'
    latent_representation: str = 'latent_GVAE'
    num_neighbour: int = 21
    num_neighbour_spatial: int = 101
    num_processes: int = 4
    fold: float = 1.0
    pst: float = 0.2
    species: str = None
    gs_species: str = None
    gM_slices: str = None
    annotation: str = None
    type: str = None


@dataclass
class RunAllModeConfig:
    flr_config: FindLatentRepresentationsConfig
    ltg_config: LatentToGeneConfig
    gls_config: GenerateLDScoreConfig


@register_cli(name='run_find_latent_representations',
              description='Run Find_latent_representations \nFind the latent representations of each spot by running GNN-VAE',
              add_args_function=add_find_latent_representations_args)
def run_find_latent_representation_from_cli(args: argparse.ArgumentParser):
    from GPS.find_latent_representation import run_find_latent_representation
    config = get_dataclass_from_parser(args, FindLatentRepresentationsConfig)
    run_find_latent_representation(config)


@register_cli(name='run_latent_to_gene',
              description='Run Latent_to_gene \nFind gene marker gene scores for each spot by using latent representations from nearby spots',
              add_args_function=add_latent_to_gene_args)
def run_latent_to_gene_from_cli(args: argparse.ArgumentParser):
    from GPS.latent_to_gene import run_latent_to_gene
    config = get_dataclass_from_parser(args, LatentToGeneConfig)
    run_latent_to_gene(config)

@register_cli(name='run_generate_ldscore',
                description='Run Generate_ldscore \nGenerate LD scores for each spot',
                add_args_function=add_generate_ldscore_args)
def run_generate_ldscore_from_cli(args: argparse.ArgumentParser):
    from GPS.generate_ldscore import run_generate_ldscore
    config = get_dataclass_from_parser(args, GenerateLDScoreConfig)
    run_generate_ldscore(config)

@register_cli(name='run_all_mode',
                description='Run GPS Pipeline \nGSP Pipeline (Run Find_latent_representations, Latent_to_gene, and Generate_ldscore) in order',
                add_args_function=add_all_mode_args)
def run_all_mode_from_cli(args: argparse.ArgumentParser):
    from GPS.find_latent_representation import run_find_latent_representation
    from GPS.latent_to_gene import run_latent_to_gene
    from GPS.generate_ldscore import run_generate_ldscore
    config = get_runall_mode_config(args)
    run_find_latent_representation(config.flr_config)
    run_latent_to_gene(config.ltg_config)
    run_generate_ldscore(config.gls_config)

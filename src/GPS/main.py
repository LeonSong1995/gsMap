import logging
from collections import OrderedDict, namedtuple
from GPS import (__version__)
from GPS.config import *
from GPS.find_latent_representation import run_find_latent_representation
from GPS.generate_ldscore import run_generate_ldscore
from GPS.latent_to_gene import run_latent_to_gene

subcommand_config = namedtuple('subcommand_config',
                               ['name',
                                'run_config',
                                'run_function',
                                'add_args_function',
                                'description'])

find_latent_subcommand_config = subcommand_config(
    name='find_latent_representations',
    run_config=FindLatentRepresentationsConfig,
    run_function=run_find_latent_representation,
    add_args_function=add_find_latent_representations_args,
    description='Find Latent Representations'
)
latent_to_gene_subcommand_config = subcommand_config(
    name='latent_to_gene',
    run_config=LatentToGeneConfig,
    run_function=run_latent_to_gene,
    add_args_function=add_latent_to_gene_args,
    description='Latent to Gene'
)

generate_ldscore_subcommand_config = subcommand_config(
    name='generate_ldscore',
    run_config=GenerateLDScoreConfig,
    run_function=run_generate_ldscore,
    add_args_function=add_generate_ldscore_args,
    description='Generate LD Sc ore'
)

subcommand_list = [find_latent_subcommand_config, latent_to_gene_subcommand_config, generate_ldscore_subcommand_config]
subcommand_ordered_dict = OrderedDict([(subcommand.name, subcommand) for subcommand in subcommand_list])

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def main():
    parser = argparse.ArgumentParser(description=" GPS: Genetics-informed pathogenic spatial mapping")
    parser.add_argument('--version', '-v', action='version', version=f'GPS version {__version__}')
    subparsers = parser.add_subparsers(dest="subcommand", help="Subcommands")

    for subcommand in subcommand_list:
        subcommand_parser = subparsers.add_parser(subcommand.name, help=subcommand.description)
        subcommand.add_args_function(subcommand_parser)
        subcommand_parser.set_defaults(subcommand=subcommand)
    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    if args.subcommand is None:
        parser.print_help()
        exit(1)
    args.func(
        subcommand_ordered_dict[args.subcommand].run_config(
            **vars(args)
        ))

from GPS import (__version__)
from GPS.config import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def main():
    parser = create_parser()
    args = parser.parse_args()
    if args.subcommand is None:
        parser.print_help()
        exit(1)
    args.func(
        args
    )

def create_parser():
    parser = argparse.ArgumentParser(description=" GPS: Genetics-informed pathogenic spatial mapping",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     prog='GPS'
                                     )
    parser.add_argument('--version', '-v', action='version', version=f'GPS version {__version__}')
    subparsers = parser.add_subparsers(dest="subcommand", help="Subcommands", title="Available subcommands")
    for subcommand in cli_function_registry.values():
        subcommand_parser = subparsers.add_parser(subcommand.name, help=subcommand.description,
                                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                                  )
        subcommand.add_args_function(subcommand_parser)
        subcommand_parser.set_defaults(func=subcommand.func)
    return parser


if __name__ == "__main__":
    main()

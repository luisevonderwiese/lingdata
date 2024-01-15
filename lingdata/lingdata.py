import argparse
import sys

import lingdata.database as database
import lingdata.converter as converter


def print_header():
    print(  "Lingdata\n"
            "Developed by: Luise Häuser\n"
            "Latest version: https://github.com/luisevonderwiese/lingdata\n"
            "Questions/problems/suggestions? Please open an issue on GitHub.\n",
    )


def main():
    print_header()

    parser = argparse.ArgumentParser(
        description="Parser for lingdata command line options."
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "--generate",
        required=False,
        action="store_true",
        help="Generate MSA from CLDF",
    )

    group.add_argument(
        "--download",
        required=False,
        action="store_true",
        help="Download native linguistic data from sources provided in config",
    )

    group.add_argument(
        "--compile",
        required=False,
        action="store_true",
        help="Generate data and setup database. Requires native data to be stored in the native_dir speicfied in the provided config. Run download-native first.",
    )

    parser.add_argument(
        "-c",
        "--config",
        required=('--compile' in sys.argv or '--download' in sys.argv),
        help="Path to json config",
    )

    parser.add_argument(
        "-i",
        "--input",
        required=('--generate' in sys.argv),
        help="Path to directory with input cldf files",
    )

    parser.add_argument(
        "-l",
        "--ling_type",
        required=('--generate' in sys.argv),
        help="Linguistic datatype of input data",
        choices=["cognate", "structural", "correspondence"]
    )

    parser.add_argument(
        "-o",
        "--output",
        required=('--generate' in sys.argv),
        help="Path to output MSA file",
    )

    parser.add_argument(
        "-m",
        "--msa_type",
        required=('--generate' in sys.argv),
        help="Type of MSA to be generated",
        choices=["bin", "multi", "catg_bin", "catg_multi", "ambig"]
    )

    args = parser.parse_args()
    if args.download:
        database.read_config(args.config)
        database.update_native()
    elif args.compile:
        database.read_config(args.config)
        database.generate_data()
    elif args.generate:
        converter.cldf_to_msa(args.input, args.ling_type, args.output, args.msa_type)

"""Maintainer entry point: python -m proteolyzer.unimod ...

build   download UniMod and (re)build the cached database
export  write the CSVs that ship in proteolyzer/resources
"""

import argparse
import sys

from . import database, database_path, refresh
from .export import UniModProcessor


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command")

    build = sub.add_parser("build", help="download UniMod and rebuild the cache")
    build.add_argument("--xml-source", help="override the UniMod XML URL or path")
    build.add_argument("--xsd-source", help="override the UniMod XSD URL or path")

    export = sub.add_parser("export", help="write the bundled reference CSVs")
    export.add_argument(
        "--db-file",
        default=None,
        help=f"database to export from (default: {database_path()})",
    )
    export.add_argument("--mods-output", help="output CSV for the modifications")
    export.add_argument("--aa-output", help="output CSV for the amino acids")

    args = parser.parse_args()

    if args.command == "build":
        print(f"Built {refresh(args.xml_source, args.xsd_source)}")
    elif args.command == "export":
        db_file = args.db_file or str(database())
        UniModProcessor(
            db_file=db_file,
            mods_output=args.mods_output,
            aa_output=args.aa_output,
        ).process_and_save()
    else:
        parser.print_help()
        sys.exit(2)


if __name__ == "__main__":  # pragma: no cover
    main()

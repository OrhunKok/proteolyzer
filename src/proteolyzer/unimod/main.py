# main.py
import argparse
from loader import UnimodDBLoader
from processor import UniModProcessor


def main() -> None:
    """Handles command-line arguments and orchestrates the execution of data loading or processing."""
    parser = argparse.ArgumentParser(
        description="Load UNIMOD XML data into an SQLite database or process and save modifications and amino acid data."
    )

    subparsers = parser.add_subparsers(dest="command")

    load_parser = subparsers.add_parser('load', help="Load data from UNIMOD XML into the SQLite database.")
    load_parser.add_argument('--xsd-source', type=str, default=UnimodDBLoader.XSD_SOURCE, help="URL for the XSD schema.")
    load_parser.add_argument('--xml-source', type=str, default=UnimodDBLoader.XML_SOURCE, help="URL for the XML file.")
    load_parser.add_argument('--db-conn', type=str, default=UnimodDBLoader.DB_CONN, help="Connection string for the SQLite database.")

    process_parser = subparsers.add_parser('process', help="Process and save data from the SQLite database.")
    process_parser.add_argument('--db-file', type=str, default=UniModProcessor.DB_FILE, help="SQLite database file path.")
    process_parser.add_argument('--mods-output', type=str, default=UniModProcessor.MODS_OUTPUT, help="Output CSV file for modifications.")
    process_parser.add_argument('--aa-output', type=str, default=UniModProcessor.AA_OUTPUT, help="Output CSV file for amino acids.")

    args = parser.parse_args()

    if args.command == "load":
        print("Initializing UnimodDBLoader...")
        loader = UnimodDBLoader(xsd_source=args.xsd_source, xml_source=args.xml_source, db_conn=args.db_conn)
        print("Loading and cleaning data...")
        loader.load_and_clean()

    elif args.command == "process":
        print("Initializing UniModProcessor...")
        processor = UniModProcessor(db_file=args.db_file, mods_output=args.mods_output, aa_output=args.aa_output)
        print("Processing and saving data...")
        processor.process_and_save()


if __name__ == "__main__":
    main()

import argparse
import sys
from loader import UnimodDBLoader
from processor import UniModProcessor

def main() -> None:
    """Handles command-line arguments and orchestrates the execution of data loading or processing."""
    parser = argparse.ArgumentParser(
        description="Load UNIMOD XML data into an SQLite database or process and save modifications and amino acid data."
    )

    subparsers = parser.add_subparsers(dest="command")

    # Load Command
    load_parser = subparsers.add_parser('load', help="Load data from UNIMOD XML into the SQLite database.")
    load_parser.add_argument(
        '--db-output', 
        required=True, 
        type=str, 
        help="Path to the SQLite database file. (Required)"
    )
    load_parser.add_argument(
        '--xsd-source', 
        type=str, 
        default=UnimodDBLoader.XSD_SOURCE, 
        help=f"Source URL or path for the XSD schema file. Default: {UnimodDBLoader.XSD_SOURCE}"
    )
    load_parser.add_argument(
        '--xml-source', 
        type=str, 
        default=UnimodDBLoader.XML_SOURCE, 
        help=f"Source URL or path for the XML data file. Default: {UnimodDBLoader.XML_SOURCE}"
    )

    # Process Command
    process_parser = subparsers.add_parser('process', help="Process and save data from the SQLite database.")
    process_parser.add_argument(
        '--db-file', 
        required=True, 
        type=str, 
        help="SQLite database file path. (Required)"
    )
    process_parser.add_argument(
        '--mods-output', 
        type=str, 
        help="Output CSV file for modifications. If not specified, defaults to the same location as the database file."
    )
    process_parser.add_argument(
        '--aa-output', 
        type=str, 
        help="Output CSV file for amino acids. If not specified, defaults to the same location as the database file."
    )

    args = parser.parse_args()

    if args.command == "load":
        print("Initializing UnimodDBLoader...")
        loader = UnimodDBLoader(
            db_output=args.db_output,
            xsd_source=args.xsd_source,
            xml_source=args.xml_source
        )
        print("Loading and cleaning data...")
        try:
            loader.load_and_clean()
            print("Data loaded and cleaned successfully.")
        except Exception as e:
            print(f"Error while loading data: {e}")
            sys.exit(1)

    elif args.command == "process":
        print("Initializing UniModProcessor...")
        processor = UniModProcessor(
            db_file=args.db_file,
            mods_output=args.mods_output,
            aa_output=args.aa_output
        )
        print("Processing and saving data...")
        try:
            processor.process_and_save()
            print("Data processed and saved successfully.")
        except Exception as e:
            print(f"Error while processing data: {e}")
            sys.exit(1)

if __name__ == "__main__":
    main()

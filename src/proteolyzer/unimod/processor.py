import sqlite3
import pandas as pd
import numpy as np
import os
import sys
import argparse
from typing import Optional


class UniModProcessor:
    """
    A class to extract, process, and save modification and amino acid data
    from a UniMod SQLite database.
    """
    MODS_OUTPUT: str = 'unimod_modifications.csv'
    AA_OUTPUT: str = 'unimod_amino_acids.csv'

    def __init__(self, db_file: str, mods_output: Optional[str] = None, aa_output: Optional[str] = None) -> None:
        """
        Initializes the UniModDataProcessor with a database file path.
        """
        
        self.db_file = db_file
        db_dir = os.path.dirname(os.path.abspath(self.db_file))
        
        self.mods_output = mods_output if mods_output else os.path.join(db_dir, self.MODS_OUTPUT)
        self.aa_output = aa_output if aa_output else os.path.join(db_dir, self.AA_OUTPUT)
        
        if not os.path.exists(self.db_file):
            print(f"Error: Database file not found at '{self.db_file}'.")
            sys.exit(1)
        
        self.connection: Optional[sqlite3.Connection] = None

    def connect_to_db(self) -> None:
        """Establishes a connection to the SQLite database."""
        try:
            self.connection = sqlite3.connect(self.db_file)
            print(f"Connected to database: {self.db_file}")
        except sqlite3.Error as e:
            print(f"Database connection error: {e}")
            sys.exit(1)

    def close_connection(self) -> None:
        """Closes the database connection."""
        if self.connection:
            self.connection.close()
            print("Database connection closed.")

    # --- Data Access Methods ---
    def get_modifications(self) -> pd.DataFrame:
        """Fetches approved modification data."""
        mod_query = """
        SELECT
            s.one_letter,
            p.position,
            m.record_id AS unimod_id,
            m.mono_mass,
            m.full_name,
            m.code_name,
            m.composition,
            c.classification
        FROM specificity AS s
        JOIN modifications AS m ON s.mod_key = m.record_id
        JOIN positions AS p ON s.position_key = p.record_id
        JOIN classifications AS c ON s.classifications_key = c.record_id
        WHERE
            m.username_of_poster = 'unimod' OR m.approved = 1
        """
        return pd.read_sql_query(mod_query, self.connection)

    def get_amino_acids(self) -> pd.DataFrame:
        """Fetches standard amino acid composition data."""
        aa_query = """
        SELECT
            a.one_letter,
            a.three_letter,
            a.full_name,
            a.num_H,
            a.num_O,
            a.num_C,
            a.num_N,
            a.num_S,
            a.num_Se
        FROM amino_acids AS a
        WHERE
            a.one_letter != '-'
        """
        return pd.read_sql_query(aa_query, self.connection)

    def get_elements(self) -> pd.DataFrame:
        """Fetches element mono-mass data."""
        elements_query = """
        SELECT
            e.element,
            e.full_name,
            e.mono_mass
        FROM elements AS e
        """
        return pd.read_sql_query(elements_query, self.connection)

    # --- Data Processing Methods ---
    def calculate_aa_masses(self, amino_acids_df: pd.DataFrame, elements_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculates mono-mass for amino acids and aggregates Isoleucine/Leucine into Xle ('J').
        """
        element_vector: np.ndarray = elements_df.set_index('element').loc[['H', 'O', 'C', 'N', 'S', 'Se'], 'mono_mass'].values
        composition_cols: list[str] = ['num_H', 'num_O', 'num_C', 'num_N', 'num_S', 'num_Se']
        
        amino_acids_df['mono_mass'] = amino_acids_df.set_index(
            ['one_letter', 'three_letter', 'full_name']
        )[composition_cols].dot(element_vector).values
        
        # Add Xle ('J')
        xle_rows = amino_acids_df[amino_acids_df['one_letter'] == 'L'].copy()
        xle_rows['one_letter'] = 'J'
        xle_rows['three_letter'] = 'Xle'
        xle_rows['full_name'] = 'Isoleucine/Leucine'
        processed_aas = pd.concat([amino_acids_df, xle_rows], ignore_index=True)

        return processed_aas

    # --- Main Execution Method ---
    def process_and_save(self) -> None:
        """Executes the entire data extraction, processing, and saving process."""
        print("Starting data extraction and processing...")

        self.connect_to_db()

        mods_df = self.get_modifications()
        amino_acids_df = self.get_amino_acids()
        elements_df = self.get_elements()

        print("Calculating amino acid masses...")
        processed_aa_df = self.calculate_aa_masses(amino_acids_df, elements_df)

        self.close_connection()
        
        print(f"Saving results to CSV files: {self.mods_output}, {self.aa_output}")
        mods_df.to_csv(self.mods_output, index=False)
        processed_aa_df.to_csv(self.aa_output, index=False)

        print(f"Results saved:\n- {self.mods_output}\n- {self.aa_output}")


def main() -> None:
    """Handles command-line arguments and orchestrates the execution of data processing."""
    parser = argparse.ArgumentParser(
        description="Extracts and processes modification and amino acid data from a UniMod SQLite database."
    )

    parser.add_argument(
        '--db-file',
        required=True,
        type=str,
        help=f"Path to the UniMod SQLite database file. (Required)"
    )

    parser.add_argument(
        '--mods-output',
        type=str,
        help=f"Output file for the modifications data. If not specified, defaults to the same location as the database file."
    )

    parser.add_argument(
        '--aa-output',
        type=str,
        help=f"Output file for the amino acids data. If not specified, defaults to the same location as the database file."
    )
    
    args = parser.parse_args()
    processor = UniModProcessor(db_file=args.db_file, mods_output=args.mods_output, aa_output=args.aa_output)
    processor.process_and_save()


if __name__ == '__main__':
    main()

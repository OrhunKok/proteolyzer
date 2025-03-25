import pandas as pd
from functools import partial
import csv
import pyarrow.parquet as pq
from proteoboost.utils.models import Data
from proteoboost.utils.logging import MetaLogging

class DataLoader(metaclass=MetaLogging):
    """Loads data from various file formats."""

    __slots__ = ('file', 'data', 'INPUT_TYPE', 'cols_rename_mapping', 'logger')

    def __init__(self, file: Data):
        """Initializes the DataLoader."""
        self.file = file
        self.INPUT_TYPE = file.input_type
        self.cols_rename_mapping = file.cols_rename_mapping

        self.data = self._auto_load()

        if self.cols_rename_mapping:
            self.data = self._rename_cols()

        self._memory_check(self.data)

    def _auto_load(self) -> pd.DataFrame:
        """Automatically loads data based on file extension."""
        load_methods = {
            '.csv': partial(self._load_csv, delimiter=','),
            '.tsv': partial(self._load_csv, delimiter='\t'),
            '.txt': self._load_csv,
            '.xls': self._load_excel,
            '.xlsx': self._load_excel,
            '.parquet': self._load_parquet,
        }

        loader = load_methods.get(self.file.file_extension)
        if not loader:
            self.logger.error(f"Unsupported file format: {self.file.file_extension} for file: {self.file.file_name}")
            raise ValueError(f"Unsupported file format: {self.file.file_extension}")
        return loader()

    def _rename_cols(self) -> pd.DataFrame:
        """Renames columns based on the mapping."""
        self.logger.info(f'Renaming columns in {self.INPUT_TYPE} input to match mapping')
        return self.data.rename(columns=self.cols_rename_mapping)

    def _get_delimiter(self, default_delimiter='\t', min_sample_size=524288, sample_percent=0.01) -> str:
        """Detects the delimiter of a CSV-like file."""
        sample_size = max(min_sample_size, int(self.file.file_stats['Size (Bytes)'] * sample_percent))

        try:
            with open(self.file.file_path, 'r', newline='', encoding='utf-8') as csvfile:
                sample = csvfile.read(sample_size)
            delimiter = csv.Sniffer().sniff(sample).delimiter
        except Exception as e:
            delimiter = default_delimiter
            self.logger.error(f"{e} for: {self.file.file_path}. Falling back to default delimiter {repr(delimiter)}", stacklevel=2)

        return delimiter

    def _cols_to_load(self, all_cols: set) -> list:
        """Determines columns to load based on settings."""
        if not self.file.load_all_columns and self.file.cols_subset is not None:
            return list(all_cols & self.file.cols_subset)
        return list(all_cols)

    def _load_csv(self, delimiter: str = None) -> pd.DataFrame:
        """Loads a CSV or TSV file."""
        if delimiter is None:
            delimiter = self._get_delimiter()
        self.logger.info(f"Loading {self.file.file_path} with delimiter '{delimiter}'")

        try:
            df = pd.read_csv(self.file.file_path, delimiter=delimiter, nrows=0)
            cols_to_load = self._cols_to_load(set(df.columns))
            return pd.read_csv(self.file.file_path, delimiter=delimiter, usecols=cols_to_load)
        except FileNotFoundError:
            self.logger.error(f"File not found: {self.file.file_path}")
            raise
        except pd.errors.ParserError:
            self.logger.error(f"Error parsing CSV file: {self.file.file_path}")
            raise
        except Exception as e:
            self.logger.error(f"An unexpected error occured when loading {self.file.file_path}: {e}")
            raise

    def _load_excel(self) -> pd.DataFrame:
        """Loads an Excel file."""
        try:
            df = pd.read_excel(self.file.file_path, nrows=0)
            cols_to_load = self._cols_to_load(set(df.columns))
            return pd.read_excel(self.file.file_path, usecols=cols_to_load)
        except FileNotFoundError:
            self.logger.error(f"File not found: {self.file.file_path}")
            raise
        except Exception as e:
            self.logger.error(f"An unexpected error occured when loading {self.file.file_path}: {e}")
            raise

    def _load_parquet(self) -> pd.DataFrame:
        """Loads a Parquet file."""
        try:
            schema = pq.ParquetFile(self.file.file_path).schema
            cols_to_load = self._cols_to_load(set(schema.names))
            return pd.read_parquet(self.file.file_path, columns=cols_to_load)
        except FileNotFoundError:
            self.logger.error(f"File not found: {self.file.file_path}")
            raise
        except Exception as e:
            self.logger.error(f"An unexpected error occured when loading {self.file.file_path}: {e}")
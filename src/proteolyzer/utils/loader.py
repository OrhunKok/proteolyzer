"""Data loading helpers for common proteolyzer inputs.

This module implements robust, documented helpers to read CSV/TSV, Excel,
and other domain-specific export formats. Helpers validate required
columns and return pandas DataFrame objects.
"""

import pandas as pd
from functools import partial
import csv
import magic
import chardet
import pyarrow.parquet as pq
from proteolyzer.utils.models import Data
from proteolyzer.utils.logging import MetaLogging


class DataLoader(metaclass=MetaLogging):
    """Loads data from various file formats."""

    __slots__ = ("file", "data", "INPUT_TYPE", "cols_rename_mapping", "logger")

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
            ".csv": partial(self._load_csv, delimiter=","),
            ".tsv": partial(self._load_csv, delimiter="\t"),
            ".txt": self._load_txt,
            ".xls": self._load_excel,
            ".xlsx": self._load_excel,
            ".parquet": self._load_parquet,
        }

        loader = load_methods.get(self.file.file_extension, self._load_txt)
        if not loader:
            self.logger.error(
                f"Unsupported file format: {self.file.file_extension} for file: {self.file.file_name}"
            )
            raise ValueError(f"Unsupported file format: {self.file.file_extension}")
        return loader()

    def _rename_cols(self) -> pd.DataFrame:
        """Renames columns based on the mapping."""
        self.logger.info(
            f"Renaming columns in {self.INPUT_TYPE} input to match mapping"
        )
        return self.data.rename(columns=self.cols_rename_mapping)

    def _get_delimiter(
        self,
        default_delimiter: str = "\t",
        encoding: str = "utf-8",
        min_sample_size: int = 524288,
        sample_percent: float = 0.01,
    ) -> str:
        """Detects the delimiter of a CSV-like file."""
        sample_size = max(
            min_sample_size, int(self.file.file_stats["Size (Bytes)"] * sample_percent)
        )

        try:
            with open(
                self.file.file_path, "r", newline="", encoding=encoding
            ) as csvfile:
                sample = csvfile.read(sample_size)
            delimiter = csv.Sniffer().sniff(sample).delimiter
        except Exception as e:
            delimiter = default_delimiter
            self.logger.error(
                f"{e} for: {self.file.file_path}. Falling back to default delimiter {repr(delimiter)}",
                stacklevel=2,
            )

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
            return pd.read_csv(
                self.file.file_path, delimiter=delimiter, usecols=cols_to_load
            )
        except FileNotFoundError:
            self.logger.error(f"File not found: {self.file.file_path}")
            raise
        except pd.errors.ParserError:
            self.logger.error(f"Error parsing CSV file: {self.file.file_path}")
            raise
        except Exception as e:
            self.logger.error(
                f"An unexpected error occured when loading {self.file.file_path}: {e}"
            )
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
            self.logger.error(
                f"An unexpected error occured when loading {self.file.file_path}: {e}"
            )
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
            self.logger.error(
                f"An unexpected error occured when loading {self.file.file_path}: {e}"
            )

    def _load_txt(self) -> pd.DataFrame:
        """Loads a plaintext file with detected MIME type and encoding."""

        try:
            mime_type, encoding = self._detect_file_type_and_encoding()
        except Exception as e:
            self.logger.error(f"Failed to detect MIME type or encoding: {e}")
            raise

        if mime_type != "text/plain":
            self.logger.error(
                f"Unsupported MIME type {mime_type} for .txt file. Only 'text/plain' is supported."
            )
            raise

        if self.INPUT_TYPE == "MaxQuant":
            self.logger.info(
                "Detected MaxQuant format, treating .txt as structured CSV."
            )
            return self._load_csv()

        try:
            with open(self.file.file_path, encoding=encoding) as file:
                lines = [line for line in file]
            df = pd.DataFrame({"line": lines})
            self.logger.info(
                f"Loaded plaintext file as DataFrame with {len(df)} lines."
            )
            return df
        except FileNotFoundError:
            self.logger.error(f"File not found: {self.file.file_path}")
            raise
        except Exception as e:
            self.logger.error(
                f"An unexpected error occurred while loading {self.file.file_path}: {e}"
            )
            raise

    def _detect_file_type_and_encoding(self):
        """Strictly detects MIME type and text encoding for the file. Raises if undetectable."""

        file_path = str(self.file.file_path)

        mime_type = magic.from_file(file_path, mime=True)

        if mime_type:
            with open(file_path, "rb") as f:
                raw_data = f.read(100_000)
                encoding_result = chardet.detect(raw_data)
                encoding = encoding_result.get("encoding")
                confidence = encoding_result.get("confidence", 0.0)

            if encoding:
                self.logger.info(
                    f"Detected MIME type: {mime_type}, Encoding: {encoding}, Encoding detection reliability {confidence:.2f}, for file: {file_path}"
                )
            else:
                self.logger.error(
                    f"Detected MIME type: {mime_type}, but unable to detect encoding for file: {file_path}"
                )
                raise

        else:
            self.logger.error(
                f"Failed to detect MIME type and encoding for file: {file_path}"
            )
            raise

        return mime_type, encoding

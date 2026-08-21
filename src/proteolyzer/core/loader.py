"""Data loading helpers for common proteolyzer inputs.

This module implements robust, documented helpers to read CSV/TSV, Excel,
and other domain-specific export formats. Helpers validate required
columns and return pandas DataFrame objects.
"""

import csv
from collections.abc import Callable
from functools import partial
from pathlib import Path

import chardet
import pandas as pd
import pyarrow.parquet as pq

from .formats import Config
from .logging import Logged
from .models import Data

CONFIG = Config()


class DataLoader(Logged):
    """Loads data from various file formats with robust handling."""

    __slots__ = ("data", "file")

    def __init__(self, file: Data, verbose: bool = False):
        """
        Initializes DataLoader.

        Parameters
        ----------
        file : Data
            The source descriptor. `file.source` may be a Path, string path,
            or file-like object. It is held rather than copied field by field,
            so the two cannot drift apart.
        verbose : bool, default False
            If True, run optional diagnostic steps after loading (currently:
            log the in-memory size of the loaded data). Off by default because
            memory_usage(deep=True) walks every cell of every object column,
            which is expensive on large frames.
        """
        self.file = file
        self.data = self._auto_load()

        if self.cols_rename_mapping:
            self.data = self._rename_cols()

        if verbose:
            self._memory_check(self.data)

    @property
    def source(self):
        return self.file.source

    @property
    def INPUT_TYPE(self) -> str:
        return self.file.input_type

    @property
    def cols_rename_mapping(self) -> dict:
        return self.file.cols_rename_mapping

    @property
    def cols_subset(self):
        return self.file.cols_subset

    @property
    def file_name(self) -> str:
        return self.file.file_name

    @property
    def file_extension(self) -> str:
        """Lower-cased suffix; taken from Data so named streams dispatch too."""
        return self.file.file_extension.lower()

    @property
    def is_path(self) -> bool:
        return isinstance(self.source, (Path, str))

    def _auto_load(self) -> pd.DataFrame:
        """Automatically dispatch to the correct loader based on file extension."""
        load_methods: dict[str, Callable[[], pd.DataFrame]] = {
            ".csv": partial(self._load_csv, delimiter=","),
            ".tsv": partial(self._load_csv, delimiter="\t"),
            ".txt": self._load_txt,
            ".xls": self._load_excel,
            ".xlsx": self._load_excel,
            ".parquet": self._load_parquet,
        }

        loader = load_methods.get(self.file_extension, self._load_txt)
        return loader()

    def _rewind(self) -> None:
        """Return a file-like source to the start so it can be read again.

        Peeking at the header and then reading the body means the source is read
        twice; paths reopen transparently but streams have to be rewound.
        """
        if self.is_path:
            return
        seek = getattr(self.source, "seek", None)
        if callable(seek):
            seek(0)

    def _rename_cols(self, df: pd.DataFrame | None = None) -> pd.DataFrame:
        df = self.data if df is None else df
        if not self.cols_rename_mapping:
            return df
        self.logger.info(f"Renaming columns in {self.INPUT_TYPE} input")
        return df.rename(columns=self.cols_rename_mapping)

    def _cols_to_load(self, all_cols) -> list:
        """Intersect the file's columns with the configured subset.

        Keeps the order the columns appear in the file, so repeated loads of the
        same file always produce the same column order.
        """
        all_cols = list(all_cols)
        if not self.cols_subset:
            return all_cols
        wanted = set(self.cols_subset)
        return [col for col in all_cols if col in wanted]

    # ---------------- CSV / TSV ----------------

    def _load_csv(self, delimiter: str | None = None) -> pd.DataFrame:
        if delimiter is None:
            delimiter = self._get_delimiter()

        self.logger.info(f"Loading {self.source} with delimiter '{delimiter}'")
        try:
            # Peek at columns
            df = pd.read_csv(self.source, delimiter=delimiter, nrows=0)
            cols_to_load = self._cols_to_load(df.columns)
            self._rewind()
            return self._read_delimited(delimiter, cols_to_load)
        except Exception as e:
            self.logger.error(f"Error loading CSV: {self.source}, {e}")
            raise

    def _read_delimited(self, delimiter: str, cols_to_load: list) -> pd.DataFrame:
        """Read the body, preferring pyarrow's multithreaded CSV reader.

        A search report is the largest thing proteolyzer reads, and on a
        million-row report pyarrow is ~4x faster than the stock parser whole,
        ~11x when reading a subset of the columns, for identical dtypes and
        values. It is stricter, though -- it rejects ragged rows the stock
        parser pads -- so that one stays as a fallback.

        The two parsers agree on column order only because
        :meth:`_cols_to_load` hands over the columns in the order the file has
        them: pyarrow returns them in the order asked for, the stock parser
        always in file order.
        """
        try:
            df = pd.read_csv(
                self.source,
                delimiter=delimiter,
                usecols=cols_to_load,
                engine="pyarrow",
            )
        except Exception as exc:
            reason = f"pyarrow could not parse it ({exc})"
        else:
            # A column pyarrow cannot decode as UTF-8 comes back as raw bytes
            # rather than raising, which would leave bytes where the rest of
            # proteolyzer expects text. Decoded text arrives as str dtype, so
            # an object column means something was left undecoded.
            undecoded = [col for col in df.columns if df[col].dtype == object]
            if not undecoded:
                return df
            reason = f"pyarrow left {undecoded} undecoded"

        self.logger.info(f"{reason}; re-reading with the default parser.")
        self._rewind()
        return pd.read_csv(self.source, delimiter=delimiter, usecols=cols_to_load)

    def _get_delimiter(
        self, default_delimiter="\t", sample_size=524288, sample_percent=0.01
    ) -> str:
        """Detect delimiter for CSV-like files."""
        if not self.is_path:
            return default_delimiter  # Cannot sniff a stream
        file_path = str(self.source)
        try:
            size = max(
                sample_size, int(Path(file_path).stat().st_size * sample_percent)
            )
            with open(file_path, encoding="utf-8") as f:
                sample = f.read(size)
            return csv.Sniffer().sniff(sample).delimiter
        except Exception as e:
            self.logger.warning(
                f"Could not detect delimiter, using default '{default_delimiter}': {e}"
            )
            return default_delimiter

    # ---------------- Excel ----------------

    def _load_excel(self) -> pd.DataFrame:
        try:
            df = pd.read_excel(self.source, nrows=0)
            cols_to_load = self._cols_to_load(df.columns)
            self._rewind()
            return pd.read_excel(self.source, usecols=cols_to_load)
        except Exception as e:
            self.logger.error(f"Error loading Excel: {self.source}, {e}")
            raise

    # ---------------- Parquet ----------------

    def _load_parquet(self) -> pd.DataFrame:
        try:
            schema = pq.ParquetFile(self.source).schema
            cols_to_load = self._cols_to_load(schema.names)
            self._rewind()
            return pd.read_parquet(self.source, columns=cols_to_load)
        except Exception as e:
            self.logger.error(f"Error loading Parquet: {self.source}, {e}")
            raise

    # ---------------- TXT ----------------

    def _load_txt(self) -> pd.DataFrame:
        """Load a text file. Detect encoding and MIME type if Path."""
        if self.INPUT_TYPE == "MaxQuant":
            return self._load_csv()

        encoding = self._detect_encoding() if self.is_path else "utf-8"

        try:
            if self.is_path:
                with open(self.source, encoding=encoding) as f:
                    lines = list(f)
            else:
                lines = list(self.source)
            df = pd.DataFrame({"line": lines})
            self.logger.info(f"Loaded {len(df)} lines from plaintext file.")
            return df
        except Exception as e:
            self.logger.error(f"Error loading TXT: {self.source}, {e}")
            raise

    # ---------------- Utilities ----------------

    def _detect_encoding(self, sample_size: int = 100_000) -> str:
        """Detect the character encoding of a file path."""
        if not self.is_path:
            raise ValueError("Cannot detect encoding for non-path source")
        file_path = str(self.source)
        with open(file_path, "rb") as f:
            enc_info = chardet.detect(f.read(sample_size))
        encoding = enc_info.get("encoding")
        if not encoding:
            raise ValueError(f"Cannot detect encoding for {file_path}")
        self.logger.info(f"Detected encoding: {encoding}")
        return encoding

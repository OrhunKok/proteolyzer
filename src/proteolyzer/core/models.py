"""Typed models describing proteolyzer inputs and outputs.

:class:`Data` describes *where* data comes from and *what* should be read from
it; :class:`LoadedData` and :class:`ProcessedData` are thin ``DataFrame``
subclasses that carry the metadata needed by the next pipeline stage.
"""

import datetime
import logging
from functools import cached_property
from pathlib import Path
from typing import IO

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, computed_field, field_validator

from .formats import Config

CONFIG = Config()
logger = logging.getLogger(__name__)

SourceType = Path | IO[str] | IO[bytes]


class Data(BaseModel):
    """A description of a single input file (or file-like object)."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    source: object = Field(..., description="Path or file-like object.")
    load_all_columns: bool = Field(
        False, description="Whether to load all columns from the file."
    )
    extra_cols_to_load: set[str] | None = Field(
        None, description="Additional columns to load."
    )
    INPUT_TYPE: str | None = Field(
        None,
        description="Manually set the input data type (e.g., 'DIANN', 'MaxQuant').",
    )

    @property
    def is_path(self) -> bool:
        return isinstance(self.source, Path)

    @property
    def is_file_like(self) -> bool:
        return callable(getattr(self.source, "read", None))

    @property
    def _path(self) -> Path:
        """The source as a Path. Only meaningful when :attr:`is_path`."""
        if not isinstance(self.source, Path):
            raise TypeError(f"source is not a path: {type(self.source).__name__}")
        return self.source

    @field_validator("source", mode="after")
    @classmethod
    def _validate_source(cls, value: SourceType) -> SourceType:
        # Convert string paths to Path automatically
        if isinstance(value, str):
            value = Path(value)

        # If it's a Path, validate existence
        if isinstance(value, Path):
            if not value.exists():
                raise ValueError(f"Path does not exist: {value}")

        # If it's file-like, ensure it has a read method
        elif not hasattr(value, "read"):
            raise TypeError(
                "source must be a Path, string path, or file-like object with a "
                ".read() method"
            )

        return value

    @field_validator("extra_cols_to_load", mode="before")
    @classmethod
    def _validate_extra_cols_to_load(cls, value: object) -> set[str] | None:
        if value is None:
            return None
        if isinstance(value, str):
            return {value}
        if isinstance(value, (list, set, tuple, frozenset)) and all(
            isinstance(item, str) for item in value
        ):
            return set(value)
        raise TypeError(f"Invalid input type: {type(value)}")

    @computed_field
    @cached_property
    def file_name(self) -> str:
        """Stem of the source file, used to recognize known input formats."""
        if self.is_path:
            return self._path.stem
        name = getattr(self.source, "name", None)
        if not isinstance(name, str) or not name:
            return "in_memory"
        return Path(name).stem

    @computed_field
    @cached_property
    def file_extension(self) -> str:
        if self.is_path:
            return self._path.suffix
        name = getattr(self.source, "name", "")
        return Path(name).suffix if isinstance(name, str) and name else ""

    @computed_field
    @cached_property
    def file_stats(self) -> dict | None:
        if not self.is_path:
            return None

        stat = self._path.stat()

        def _utc(timestamp: float) -> str:
            return datetime.datetime.fromtimestamp(timestamp, tz=datetime.UTC).strftime(
                "%Y-%m-%d %H:%M:%S"
            )

        return {
            "Size (Bytes)": stat.st_size,
            "Created": _utc(stat.st_ctime),
            "Last Modified": _utc(stat.st_mtime),
            "Last Accessed": _utc(stat.st_atime),
        }

    @computed_field
    @cached_property
    def input_type(self) -> str:
        """The search engine that produced this file: DIANN, MaxQuant or Unknown."""
        user_override = self.INPUT_TYPE

        is_diann = (
            self.file_name in CONFIG.DIANN.FILES
            and self.file_extension in CONFIG.DIANN.FILE_EXTENSIONS
        )
        is_maxquant = (
            self.file_name in CONFIG.MaxQuant.FILES
            and self.file_extension in CONFIG.MaxQuant.FILE_EXTENSIONS
        )

        if is_diann and is_maxquant:
            raise ValueError(
                f"File {self.file_name} with extension {self.file_extension} "
                "matches multiple categories."
            )

        if is_diann:
            auto_type = "DIANN"
        elif is_maxquant:
            auto_type = "MaxQuant"
        else:
            auto_type = "Unknown"

        if user_override not in (None, "Unknown"):
            if auto_type != "Unknown" and user_override != auto_type:
                logger.warning(
                    f"User input '{user_override}' conflicts with file type "
                    f"'{auto_type}'. Recommend using auto-detected type."
                )
            logger.debug(f"Using manually set input type: {user_override}")
            return user_override

        if auto_type == "Unknown":
            logger.warning(
                f"{self.file_name} source program could not be determined, certain "
                "optimizations will not be performed."
            )
        else:
            logger.debug(f"{self.file_name} determined to be {auto_type} output")

        return auto_type

    @computed_field
    @cached_property
    def cols_subset(self) -> set[str] | None:
        """Columns to read, or ``None`` to read everything."""
        if self.input_type == "Unknown" or self.load_all_columns:
            return None

        cols = getattr(getattr(CONFIG, self.input_type, None), "LOAD_COLS", {}).get(
            self.file_name
        )

        if cols is None:
            return None

        if self.extra_cols_to_load:
            return set(cols) | set(self.extra_cols_to_load)

        return set(cols)

    @computed_field
    @cached_property
    def cols_rename_mapping(self) -> dict:
        config_block = getattr(CONFIG, self.input_type, None)
        return getattr(config_block, "COLS_RENAME_MAPPING", {})

    def load(self) -> LoadedData:
        """Read the source into memory."""
        from .loader import DataLoader

        return LoadedData(DataLoader(self))


class ProcessedData(pd.DataFrame):
    """
    A specialized DataFrame to hold processed data and its metadata.
    Inherits all pandas DataFrame methods and attributes.
    """

    _metadata = [
        "ID_COL",
        "LABEL_FREE",
        "LABELS_COMPLETE",
        "LABEL_GROUP_CAPTURE",
        "PROTEASE",
    ]

    @property
    def _constructor(self):
        return ProcessedData

    def __init__(
        self,
        data=None,
        ID_COL=None,
        LABEL_FREE=None,
        LABELS_COMPLETE=None,
        LABEL_GROUP_CAPTURE=None,
        PROTEASE=None,
        **kwargs,
    ):
        super().__init__(data, **kwargs)

        self.ID_COL = ID_COL
        self.LABEL_FREE = LABEL_FREE
        #: False when channel information could not be derived in full.
        self.LABELS_COMPLETE = LABELS_COMPLETE
        self.LABEL_GROUP_CAPTURE = LABEL_GROUP_CAPTURE
        self.PROTEASE = PROTEASE

    @property
    def unique_runs(self) -> set:
        if "Run" not in self.columns:
            return set()
        return set(self["Run"].unique())

    @property
    def unique_ids(self) -> int:
        if self.ID_COL is None or self.ID_COL not in self.columns:
            return 0
        return self[self.ID_COL].nunique()


class LoadedData(pd.DataFrame):
    """Raw data as read from disk, plus the loader that produced it."""

    # Deliberately no ``_constructor`` override: slicing a LoadedData yields a
    # plain DataFrame, because this class' constructor takes a loader, not data.
    _metadata = ["loader"]

    def __init__(self, loader):
        super().__init__(loader.data.copy())
        self.loader = loader

    def process(self, **kwargs) -> ProcessedData:
        """
        Initiates processing.
        Any kwargs passed here are forwarded
        to the DataProcessor constructor.
        """
        from .processor import DataProcessor

        processor = DataProcessor(self.loader, **kwargs)
        return processor.process()

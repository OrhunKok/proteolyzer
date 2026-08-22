"""Typed models describing proteolyzer inputs and outputs.

:class:`Data` describes *where* data comes from and *what* should be read from
it. :class:`Report` is what comes back: a frame, the source it was read from,
and -- once :meth:`Report.process` has run -- a :class:`Processing` record of
what was done to it.
"""

import datetime
import logging
from dataclasses import dataclass, fields, replace
from functools import cached_property
from pathlib import Path
from typing import IO

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, computed_field, field_validator

from .formats import Config

CONFIG = Config()

#: The search engines the config describes, in the order it lists them. A format
#: is a block on Config with FILES and FILE_EXTENSIONS; nothing else has to know
#: its name.
_ENGINES: tuple[str, ...] = tuple(
    field.name
    for field in fields(CONFIG)
    if hasattr(getattr(CONFIG, field.name), "FILES")
)
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
    cols_to_load: set[str] | None = Field(
        None,
        description=(
            "Columns to load, replacing the configured subset rather than "
            "adding to it. For a caller that wants less than LOAD_COLS names."
        ),
    )
    rename: bool = Field(
        True,
        description=(
            "Whether to rename the columns onto proteolyzer's names. False "
            "keeps the file's own, for a caller written against them."
        ),
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
        """The search engine that produced this file, or Unknown.

        Whichever of the engines on :class:`~proteolyzer.core.formats.Config`
        claims the file name and extension.
        """
        user_override = self.INPUT_TYPE

        # Over the engines the config carries, rather than naming two of them:
        # adding a format is then a block in formats.py and nothing else.
        matched = [
            name
            for name in _ENGINES
            if self.file_name in getattr(CONFIG, name).FILES
            and self.file_extension in getattr(CONFIG, name).FILE_EXTENSIONS
        ]

        if len(matched) > 1:
            raise ValueError(
                f"File {self.file_name} with extension {self.file_extension} "
                f"matches multiple categories: {matched}."
            )

        auto_type = matched[0] if matched else "Unknown"

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
        """Columns to read, or ``None`` to read everything.

        ``load_all_columns`` wins, then ``cols_to_load``, which replaces the
        configured subset; otherwise the configured subset with
        ``extra_cols_to_load`` added. Whichever it is, the loader intersects it
        with the columns the file actually has.
        """
        if self.load_all_columns:
            return None

        if self.cols_to_load:
            return set(self.cols_to_load)

        if self.input_type == "Unknown":
            # Nothing describes this file, so a subset would be guesswork --
            # unless the caller asked for one above.
            return set(self.extra_cols_to_load) if self.extra_cols_to_load else None

        cols = getattr(getattr(CONFIG, self.input_type, None), "LOAD_COLS", {}).get(
            self.file_name
        )

        if cols is None:
            # No subset configured for this file. Asked-for columns are still
            # asked for: reading everything instead is how a caller after four
            # columns of an allPeptides.txt ends up reading sixty.
            return set(self.extra_cols_to_load) if self.extra_cols_to_load else None

        if self.extra_cols_to_load:
            return set(cols) | set(self.extra_cols_to_load)

        return set(cols)

    @computed_field
    @cached_property
    def cols_rename_mapping(self) -> dict:
        if not self.rename:
            return {}

        config_block = getattr(CONFIG, self.input_type, None)
        return getattr(config_block, "COLS_RENAME_MAPPING", {})

    def load(self) -> Report:
        """Read the source into memory."""
        from .loader import DataLoader

        return Report(frame=DataLoader(self).data, source=self)


@dataclass(frozen=True)
class Processing:
    """How a frame was processed, and what that revealed about it."""

    #: Column identifying a precursor, used for the identification count.
    id_col: str
    #: No labelling groups were found in the identifiers.
    label_free: bool
    #: Regex capturing a labelling group out of an identifier.
    label_group_capture: str
    #: Protease the missed-cleavage flag was computed for.
    protease: str
    #: False when channel information could not be derived in full.
    labels_complete: bool = True
    #: Large-magnitude float columns were rounded to integers.
    rounded_large_floats: bool = False
    #: float64 columns were narrowed to float32 where their values allowed it.
    narrowed_floats: bool = True


@dataclass(frozen=True)
class Report:
    """A frame of proteomics data, with where it came from and what was done.

    Composition rather than a ``DataFrame`` subclass: pandas returns plain
    frames from most operations, so a subclass silently loses its metadata the
    first time anyone slices it, and the pandas internals it has to hook into
    are not a stable API. Use :attr:`frame` for anything pandas; the few
    pass-throughs below are for interactive work.
    """

    #: The data. Anything pandas goes through here.
    frame: pd.DataFrame
    #: What the frame was read from.
    source: Data
    #: Set once :meth:`process` has run.
    processing: Processing | None = None

    @property
    def is_processed(self) -> bool:
        return self.processing is not None

    @property
    def columns(self) -> pd.Index:
        return self.frame.columns

    def __len__(self) -> int:
        return len(self.frame)

    def __getitem__(self, key):
        """Column access, for interactive use. Returns a plain pandas object."""
        return self.frame[key]

    def _repr_html_(self) -> str:  # pragma: no cover - notebook display
        state = "processed" if self.is_processed else "raw"
        return (
            f"<p><code>Report</code> ({state}), {len(self)} rows "
            f"from <code>{self.source.file_name}</code></p>"
            # Delegated so the frame is truncated as pandas would truncate it;
            # pandas-stubs does not declare this one.
            f"{self.frame._repr_html_()}"  # type: ignore[operator]
        )

    def process(self, **kwargs) -> Report:
        """Normalize the frame: dtypes, derived columns, labelling information.

        Keyword arguments are passed to :class:`~proteolyzer.core.processor
        .DataProcessor`. Returns a new Report; this one is unchanged.
        """
        from .processor import DataProcessor

        return DataProcessor(self, **kwargs).process()

    def matrix(self, values: str, index: list[str], columns: list[str]):
        """Pivot to a quantitative matrix. Returns a MatrixBuilder to chain on."""
        from .matrix import MatrixBuilder

        return MatrixBuilder(self).matrix_generation(values, index, columns)

    @property
    def runs(self) -> set:
        """The distinct runs present, empty if the frame has no Run column."""
        if "Run" not in self.frame.columns:
            return set()
        return set(self.frame["Run"].unique())

    def memory(self) -> pd.DataFrame:
        """What each column costs in memory, largest first.

        Counted deeply, so the strings behind a text column are measured
        rather than the pointers to them. This is what :meth:`process` works
        to bring down -- narrowing numeric columns and replacing repeated
        strings with categories -- and where to look when a report is larger
        than expected.
        """
        usage = self.frame.memory_usage(deep=True).drop("Index", errors="ignore")
        total = usage.sum()
        breakdown = pd.DataFrame(
            {
                "Dtype": [str(self.frame[col].dtype) for col in usage.index],
                "Bytes": usage.astype("int64"),
                "Share": usage / total if total else 0.0,
            }
        )
        breakdown.index.name = "Column"
        return breakdown.sort_values("Bytes", ascending=False)

    def summary(self) -> pd.DataFrame:
        """Identification counts per run.

        The first table to look at for a new report: what each run
        contributed, and how evenly. Counts are of *distinct* values, so a
        precursor seen in two channels of one run counts once. Levels the
        frame does not carry are left out, and a frame with no ``Run`` column
        is summarized as a single group.
        """
        proteins = (
            "Leading.Razor.Protein"
            if "Leading.Razor.Protein" in self.frame.columns
            else "Protein.Group"
        )
        levels = {
            "Precursors": self.processing.id_col if self.processing else "Precursor.Id",
            "Peptides": "Stripped.Sequence",
            "Proteins": proteins,
        }

        run = (
            self.frame["Run"]
            if "Run" in self.frame.columns
            else pd.Series("all", index=self.frame.index, name="Run")
        )
        grouped = self.frame.groupby(run, observed=True)

        summary = pd.DataFrame({"Rows": grouped.size()})
        for name, column in levels.items():
            if column in self.frame.columns:
                summary[name] = grouped[column].nunique()
        return summary

    @property
    def n_identifications(self) -> int:
        """Distinct precursors, or 0 before processing."""
        if self.processing is None:
            return 0
        id_col = self.processing.id_col
        if id_col not in self.frame.columns:
            return 0
        return int(self.frame[id_col].nunique())

    def with_frame(self, frame: pd.DataFrame) -> Report:
        """The same report around a different frame, keeping its metadata."""
        return replace(self, frame=frame)

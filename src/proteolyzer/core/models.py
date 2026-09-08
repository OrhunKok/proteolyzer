"""Typed models describing proteolyzer inputs and outputs.

:class:`Data` describes *where* data comes from and *what* should be read from
it. :class:`Report` is what comes back: a frame, the source it was read from,
and -- once :meth:`Report.process` has run -- a :class:`Processing` record of
what was done to it.
"""

import datetime
import logging
import re
from collections.abc import Collection
from dataclasses import dataclass, fields, replace
from functools import cached_property
from pathlib import Path
from typing import IO, Any, cast

import pandas as pd
import pyarrow.parquet as pq
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


#: How many of a format's ``COLUMN_SIGNATURE`` columns a file has to carry
#: before its contents identify it. Two rather than one: a single distinctive
#: name turns up in frames people derive from a report and write back out.
SIGNATURE_THRESHOLD: int = 2

#: Extensions whose first line is a header, and what separates it.
_DELIMITED: dict[str, str] = {".tsv": "\t", ".csv": ",", ".txt": "\t"}


def _shaped(block: Any, columns: Collection[str]) -> bool:
    """Whether `columns` are the whole shape of one of `block`'s narrow tables.

    Every declared column has to be there and the table has to be about as
    narrow as declared. See :class:`~proteolyzer.core.formats.Narrow`.
    """
    held = set(columns)
    return any(
        narrow.columns <= held and len(held) <= narrow.width
        for narrow in getattr(block, "NARROW_SIGNATURES", ())
    )


def _signed(block: Any, columns: Collection[str]) -> bool:
    """Whether `columns` carry enough of `block`'s signature to identify it."""
    signature: frozenset[str] = getattr(block, "COLUMN_SIGNATURE", frozenset())
    if not signature:
        return False
    return len(signature & set(columns)) >= SIGNATURE_THRESHOLD


def _claims(block: Any, file_name: str, extension: str) -> bool:
    """Whether a format block recognizes a file by this name and extension.

    By exact name where an engine names its own output, and by pattern where it
    does not. A pattern has to match the stem *in full*: matching from the start
    would take a ``..._Report.setup`` beside a report for the report itself.
    """
    if extension not in block.FILE_EXTENSIONS:
        return False

    if file_name in block.FILES:
        return True

    return any(
        re.fullmatch(pattern, file_name)
        for pattern in getattr(block, "FILE_PATTERNS", ())
    )


class Data(BaseModel):
    """A description of a single input file (or file-like object)."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    source: object = Field(..., description="Path or file-like object.")
    load_all_columns: bool = Field(
        False,
        description=(
            "Load every column, overriding any the caller named. Redundant "
            "unless one was named: reading the file whole is the default."
        ),
    )
    extra_cols_to_load: set[str] | None = Field(
        None,
        description=(
            "Columns to add to cols_to_load. On its own it widens a request "
            "that is already every column, so it changes nothing."
        ),
    )
    cols_to_load: set[str] | None = Field(
        None,
        description=(
            "The columns to read, for a project that wants fewer than the "
            "file has. Which those are is the project's to say -- the core "
            "keeps no list of its own. Names the file lacks are ignored."
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

    def _rewind(self) -> None:
        """Put a file-like source back to its start after peeking at it."""
        seek = getattr(self.source, "seek", None)
        if callable(seek):
            seek(0)

    def peek_columns(self) -> tuple[str, ...]:
        """The source's column names, read as cheaply as the format allows.

        A parquet file carries its schema in the footer and a delimited file
        gives its header up in one line, so this is cheap -- and it is only
        reached when the name settled nothing.

        Empty for anything that cannot be looked at without reading it, and
        empty on any failure: choosing a reader is not the place to raise about
        a file. The reader that follows will, and will say what it was doing.
        """
        extension = self.file_extension.lower()
        try:
            if extension == ".parquet":
                names = tuple(pq.ParquetFile(self.source).schema_arrow.names)
            elif extension in _DELIMITED:
                # `source` is declared `object` on the model; the readers
                # take what they are given, as the loader's do.
                names = tuple(
                    pd.read_csv(
                        cast(Any, self.source),
                        delimiter=_DELIMITED[extension],
                        nrows=0,
                    ).columns
                )
            else:
                return ()
        except Exception:
            return ()
        finally:
            # Peeking consumes a stream, and the loader is about to read the
            # same source from its start.
            if not self.is_path:
                self._rewind()

        return names

    @computed_field
    @cached_property
    def input_type(self) -> str:
        """The search engine that produced this file, or Unknown.

        Whichever of the engines on :class:`~proteolyzer.core.formats.Config`
        claims the file: by name where the engine names its own output, and by
        the columns inside where it does not. See
        ``docs/notes/recognising-a-format.md``.
        """
        user_override = self.INPUT_TYPE

        # Over the engines the config carries, rather than naming two of them:
        # adding a format is then a block in formats.py and nothing else.
        matched = [
            name
            for name in _ENGINES
            if _claims(getattr(CONFIG, name), self.file_name, self.file_extension)
        ]

        if len(matched) > 1:
            raise ValueError(
                f"File {self.file_name} with extension {self.file_extension} "
                f"matches multiple categories: {matched}."
            )

        # Nothing claimed the name, which for a format whose output the analyst
        # names is the ordinary case rather than a failure. Ask the file what it
        # is. Only formats that offered a signature can answer, and only a file
        # whose columns can be read without reading the file is asked.
        if not matched:
            columns = self.peek_columns()
            if columns:
                # Whole shape first, being the stronger claim, then an overlap
                # of names.
                matched = [
                    name for name in _ENGINES if _shaped(getattr(CONFIG, name), columns)
                ] or [
                    name for name in _ENGINES if _signed(getattr(CONFIG, name), columns)
                ]
                if len(matched) > 1:
                    # Not the error the name clash above is: that one is the
                    # config contradicting itself, this one is a fact about
                    # somebody's file, and Unknown is what that means.
                    logger.warning(
                        f"The columns of {self.file_name} match {matched}, so "
                        "which engine wrote it cannot be told from them. Reading "
                        "it as an unknown format; pass INPUT_TYPE= to say."
                    )
                    matched = []
                elif matched:
                    logger.debug(
                        f"{self.file_name} identified as {matched[0]} by its "
                        "columns; its name matched no format."
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

        ``None`` -- every column -- unless the caller named some. Which columns
        matter is the caller's to say and nothing here knows it: the core no
        longer carries a subset per file, because one list cannot be right for a
        dashboard and a pipeline at once. ``cols_to_load`` is the list;
        ``extra_cols_to_load`` adds to it, and on its own only widens a request
        that is already everything.

        Whichever it is, the loader intersects it with the columns the file
        actually has, so naming one that is not there is not an error.
        """
        if self.load_all_columns:
            return None

        wanted = set(self.cols_to_load or ())
        if self.extra_cols_to_load:
            # Extras on their own mean nothing now that the base is everything,
            # and a caller that only passes them still gets everything.
            wanted = (wanted | set(self.extra_cols_to_load)) if wanted else set()

        return wanted or None

    @computed_field
    @cached_property
    def cols_rename_mapping(self) -> dict:
        if not self.rename:
            return {}

        config_block = getattr(CONFIG, self.input_type, None)
        return getattr(config_block, "COLS_RENAME_MAPPING", {})

    @computed_field
    @cached_property
    def numeric_cols(self) -> frozenset[str]:
        """Columns this format writes as text that hold numbers.

        Named as the *file* names them, and not conditioned on ``rename``:
        keeping the file's own column names is a question about names, and a
        q-value that cannot be compared to a float is no more use under one
        than another.
        """
        return getattr(
            getattr(CONFIG, self.input_type, None), "NUMERIC_COLS", frozenset()
        )

    @computed_field
    @cached_property
    def boolean_cols(self) -> frozenset[str]:
        """Columns this format writes as text that hold ``True``/``False``."""
        return getattr(
            getattr(CONFIG, self.input_type, None), "BOOLEAN_COLS", frozenset()
        )

    @computed_field
    @cached_property
    def built_cols(self) -> dict:
        """Canonical columns to build after the rename, and what out of.

        For a format that writes no column of its own for something the rest of
        the package keys on. Empty when the file keeps its own column names:
        these are stated in the core's vocabulary, which such a caller has
        asked not to be given.
        """
        if not self.rename:
            return {}

        config_block = getattr(CONFIG, self.input_type, None)
        return getattr(config_block, "BUILT_COLS", {})

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

"""Higher-level processing pipelines and orchestrators.

This module composes loader, transformer and model utilities into
reusable processing functions and light-weight pipelines.
"""

import logging
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import cast

import numpy as np
import pandas as pd

from proteolyzer import reference

from .formats import Config
from .logging import Logged
from .models import Processing, Report
from .operations import per_distinct

CONFIG = Config()


@dataclass
class Narrower:
    """The dtype narrowing, over a frame and nothing else.

    Separate from :class:`DataProcessor` because the narrowing is the part of
    processing a caller wants on its own: a dashboard that reads its own
    columns, derives its own, and only wants the frame to fit in memory has no
    use for the renames, the derived columns, the labelling or the miscleavage
    counts, and had to build a Report from a source it never read to reach
    these four steps. Asked for by streamlit-DO-MS; see #24.

    The steps are the same ones :class:`DataProcessor` runs, and it runs them
    through here.
    """

    input_type: str = "Unknown"
    narrow_floats: bool = True
    round_large_floats: bool = False
    logger: logging.Logger = field(
        default_factory=lambda: logging.getLogger("proteolyzer.Narrower")
    )

    def narrow(self, df: pd.DataFrame) -> pd.DataFrame:
        """Every step, in the order they have to run in.

        Integers before floats, because a column of whole numbers should end up
        an integer rather than a narrowed float; categoricals last, since they
        are decided by measuring what the numeric steps left behind.
        """
        df = (
            df.pipe(self.convert_float_columns_to_int)
            .pipe(self.narrow_float_columns)
            .pipe(self.narrow_integer_columns)
        )

        return self.convert_columns_to_categorical(df)

    def convert_float_columns_to_int(self, df: pd.DataFrame) -> pd.DataFrame:
        """Narrows float columns that already hold whole numbers to integers.

        Integrality is tested exactly. ``np.allclose`` is not usable here: its
        default relative tolerance of 1e-5 accepts a deviation of ~1000 at the
        1e8 magnitudes typical of intensity columns, which silently rounded
        quantitative values.

        Columns that are *not* integral are left alone unless
        ``ROUND_LARGE_FLOATS`` is set; see :meth:`_round_large_floats`.
        """
        converted = {}
        for col in df.select_dtypes(include=["float"]).columns:
            values = df[col].to_numpy(dtype="float64", na_value=np.nan)
            present = ~np.isnan(values)

            if not np.array_equal(values[present], np.round(values[present])):
                continue

            df[col] = self._to_narrow_int(df[col], values, bool(present.all()))
            converted[col] = str(df[col].dtype)

        if converted:
            self.logger.info(f"Narrowed columns to integer dtype: {converted}.")

        return self._round_large_floats(df)

    @staticmethod
    def _to_narrow_int(column: pd.Series, values: np.ndarray, complete: bool):
        """Return `column` as the smallest integer dtype that holds it exactly.

        A plain numpy integer dtype is used when nothing is missing; otherwise a
        nullable one, which needs an extra byte per value for its mask. Picking
        the width matters: unconditionally using Int64 makes the float32 columns
        DIA-NN writes *larger*, defeating the point of narrowing them.
        """
        rounded = np.round(values)
        if complete:
            return pd.to_numeric(rounded, downcast="integer")

        if np.isnan(rounded).all():
            # Nothing to narrow, and reducing over an empty slice would warn.
            return column

        nullable = pd.array(rounded, dtype="Int64")
        low, high = np.nanmin(rounded), np.nanmax(rounded)
        for dtype in ("Int8", "Int16", "Int32"):
            info = np.iinfo(dtype.lower())
            if info.min <= low and high <= info.max:
                return nullable.astype(dtype)
        return pd.Series(nullable, index=column.index)

    def _round_large_floats(self, df: pd.DataFrame) -> pd.DataFrame:
        """Optionally round large-magnitude float columns to integers.

        Opt-in (``ROUND_LARGE_FLOATS``) because it discards real precision at
        the low end of a quantitative column while saving nothing: the loss is
        below float32 resolution only above ~1.7e7, and the resulting nullable
        integer dtype is wider than the float32 it replaces.
        """
        if not self.round_large_floats:
            return df

        rounded = []
        for col in df.select_dtypes(include=["float"]).columns:
            values = df[col].to_numpy(dtype="float64", na_value=np.nan)
            if np.nanmedian(values) > CONFIG.COL_MEDIAN_THRESHOLD:
                present = ~np.isnan(values)
                df[col] = self._to_narrow_int(df[col], values, bool(present.all()))
                rounded.append(col)

        if rounded:
            self.logger.warning(
                f"Rounded columns to integer dtype, discarding their fractional "
                f"part: {rounded}."
            )

        return df

    def narrow_float_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Narrows float64 columns to float32 where the values allow it.

        DIA-NN stores these columns as float32 in its own parquet output, so a
        report read from parquet is already single precision while the same
        report read from the TSV is not -- the identical data costs 44% more
        memory depending on which file it came from. This makes the text path
        match, at a worst-case relative error of 6e-8, float32's epsilon,
        measured across every float column of a real report.

        Runs after :meth:`convert_float_columns_to_int`, so columns holding
        whole numbers are already integers and keep their exact values;
        ``Ms1.Total.Signal.*`` reach 1e10, far beyond the 2**24 that float32
        can represent exactly. It runs after :meth:`extra_info` too, so
        derived columns are computed before anything is narrowed.

        A column is left alone if any value falls outside float32's normal
        range, where narrowing would turn it into an infinity or a zero rather
        than round it.

        The bound is on each value, not on what is later computed from it:
        subtracting two nearly equal narrowed values amplifies their rounding
        by the ratio between them and their difference. Pass
        ``narrow_floats=False`` when doing that kind of arithmetic on the
        frame directly.
        """
        if not self.narrow_floats:
            return df

        narrowed = {}
        for col in df.select_dtypes(include=["float64"]).columns:
            if not self._fits_in_float32(df[col]):
                continue
            if str(df[col].dtype) == "Float64":
                # Nullable narrows to nullable: numpy float32 would turn its
                # pd.NA into nan and drop the mask.
                df[col] = df[col].astype("Float32")
            else:
                df[col] = df[col].astype("float32")
            narrowed[col] = str(df[col].dtype)

        if narrowed:
            self.logger.info(f"Narrowed columns to single precision: {narrowed}.")

        return df

    def narrow_integer_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Narrows integer columns to the smallest dtype that holds them.

        Unlike the float narrowing this is exact: an integer either fits or it
        does not, and the width is chosen per column. A charge state, a run
        index or a peptide length arrives as int64 and needs one byte, while
        ``Ms1.Total.Signal`` reaches 1e10 and keeps all eight.

        Runs after :meth:`extra_info` so the integer columns derived there are
        narrowed too. Nullable dtypes stay nullable.

        The values are unchanged, but the headroom is not: adding a scalar
        that takes a narrow column past its range wraps around rather than
        raising, as it did before for the columns
        :meth:`convert_float_columns_to_int` narrows. Reductions such as
        ``sum`` accumulate in int64 and are unaffected.
        """
        narrowed = {}
        for col in df.select_dtypes(include=["integer"]).columns:
            narrower = pd.to_numeric(df[col], downcast="integer")
            if narrower.dtype != df[col].dtype:
                df[col] = narrower
                narrowed[col] = str(narrower.dtype)

        if narrowed:
            self.logger.info(f"Narrowed integer columns: {narrowed}.")

        return df

    @staticmethod
    def _fits_in_float32(column: pd.Series) -> bool:
        """Whether every value of `column` is within float32's normal range."""
        values = column.to_numpy(dtype="float64", na_value=np.nan)
        magnitudes = np.abs(values[np.isfinite(values)])
        magnitudes = magnitudes[magnitudes != 0]
        if not magnitudes.size:
            return True

        limits = np.finfo("float32")
        return bool(magnitudes.max() <= limits.max and magnitudes.min() >= limits.tiny)

    def convert_columns_to_categorical(self, df: pd.DataFrame) -> pd.DataFrame:
        """Converts columns to categorical where that actually saves memory.

        Decided by measuring both representations rather than by a cardinality
        ratio, which is a poor proxy for the saving: on a real report the
        protein, gene and name columns hold 0.3-0.4 distinct values per row
        and still halve in size, while identifier columns approach 1.0, save
        nothing, and can come out *larger* than the strings they replace. A
        ratio tight enough to exclude the second group excludes most of the
        first as well.

        Numeric columns are never converted, however few distinct values they
        hold: a categorical of numbers no longer supports arithmetic, so
        q-value or intensity columns would stop working downstream.
        """
        try:
            exclude_from_conversion = getattr(
                CONFIG, self.input_type
            ).EXCLUDE_CAT_CONVERSION
        except AttributeError:
            exclude_from_conversion = set()

        candidates = df.select_dtypes(exclude=["number", "bool"])
        converted = {}
        for col in candidates.columns:
            if col in exclude_from_conversion:
                continue
            saving, as_category = self._as_categorical(df[col])
            if saving >= CONFIG.MIN_CATEGORICAL_SAVING:
                df[col] = as_category
                converted[col] = saving

        if converted:
            self.logger.info(
                "Converted columns to categorical dtype, saving: "
                + ", ".join(f"{col} {pct:.0%}" for col, pct in converted.items())
            )

        return df

    @staticmethod
    def _as_categorical(column: pd.Series) -> tuple[float, pd.Series]:
        """`column` as a categorical, and the fraction of memory that saves.

        Measured by building it rather than estimated: estimating means
        guessing the width of the codes and the size of the dictionary, and
        the factorize it costs is what counting distinct values would cost
        anyway. The conversion is handed back so a caller that decides to keep
        it does not build it a second time. The saving is negative when the
        categorical is the larger of the two.
        """
        as_category = column.astype("category")
        current = column.memory_usage(deep=True, index=False)
        if not current:
            return 0.0, as_category
        saving = 1 - as_category.memory_usage(deep=True, index=False) / current
        return saving, as_category


class DataProcessor(Logged):
    """Processes raw data into a structured DataFrame."""

    __slots__ = (
        "report",
        "data",
        "input_type",
        "id_col",
        "label_group_capture",
        "label_free",
        "labels_complete",
        "protease",
        "round_large_floats",
        "narrow_floats",
        "_narrower",
    )

    def __init__(
        self,
        report: Report,
        id_col: str = "Precursor.Id",
        label_group_capture: str = r"\(((?:mTRAQ|SILAC|TMT)[^()]*)\)",
        protease: str = "Trypsin",
        round_large_floats: bool = False,
        narrow_floats: bool = True,
    ):
        """Initializes the DataProcessor.

        Parameters
        ----------
        report
            The report to process. Its frame is copied, so the input is
            untouched.
        round_large_floats : bool, default False
            Round float columns whose median exceeds
            ``Config.COL_MEDIAN_THRESHOLD`` to integers, discarding their
            fractional part. Off by default because it loses real precision in
            the low range of quantitative columns.
        narrow_floats : bool, default True
            Narrow float64 columns to float32. See
            :meth:`narrow_float_columns`; pass False to keep double precision.
        """
        self.report = report
        self.data = report.frame.copy()
        self.input_type = report.source.input_type
        self.id_col = id_col
        self.label_group_capture = label_group_capture
        self.protease = protease
        self.round_large_floats = round_large_floats
        self.narrow_floats = narrow_floats

        # The narrowing steps live on Narrower now, which a caller can use on
        # its own. The four methods below stay because they are public, and
        # because process() pipes through them.
        self._narrower = Narrower(
            input_type=self.input_type,
            narrow_floats=narrow_floats,
            round_large_floats=round_large_floats,
            logger=self.logger,
        )
        # True unless a label matrix has to be skipped; see _LabelGenerator.
        self.labels_complete = True
        self._check_labelfree()

    def process(self, verbose: bool = False) -> Report:
        """Processes the data and returns a new :class:`Report`.

        Parameters
        ----------
        verbose : bool, default False
            If True, run optional diagnostic steps after processing (currently:
            log the in-memory size of the processed data).
        """
        self.data = (
            self.data.pipe(self.drop_identical_cols)
            .pipe(self.convert_float_columns_to_int)
            .pipe(self.rename_columns)
            .pipe(self.extra_info)
            # After extra_info, so the columns derived there are computed at
            # full precision: RT.Width is a difference of two nearly equal
            # retention times, and narrowing its inputs first would amplify
            # their rounding by the ratio between the times and the width --
            # 6e-8 becomes 3e-5.
            .pipe(self.narrow_float_columns)
            .pipe(self.narrow_integer_columns)
        )

        self.data = self.convert_columns_to_categorical(self.data)

        if not self.label_free and self.input_type == "DIANN":
            labels = _LabelGenerator(self)
            self.data = labels.data
            self.labels_complete = labels.labels_complete
            if not self.labels_complete:
                self.logger.warning(
                    "Labelling information is incomplete; some channel columns "
                    "are missing. Check Report.processing.labels_complete "
                    "before using channel-level results."
                )

        self.data = self.miscleavages(self.data, protease=self.protease)

        if verbose:
            self._memory_check(self.data)

        return Report(
            frame=self.data,
            source=self.report.source,
            processing=Processing(
                id_col=self.id_col,
                label_free=self.label_free,
                label_group_capture=self.label_group_capture,
                protease=self.protease,
                labels_complete=self.labels_complete,
                rounded_large_floats=self.round_large_floats,
                narrowed_floats=self.narrow_floats,
            ),
        )

    def convert_float_columns_to_int(self, df: pd.DataFrame) -> pd.DataFrame:
        """See :meth:`Narrower.convert_float_columns_to_int`."""
        return self._narrower.convert_float_columns_to_int(df)

    def narrow_float_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """See :meth:`Narrower.narrow_float_columns`."""
        return self._narrower.narrow_float_columns(df)

    def narrow_integer_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """See :meth:`Narrower.narrow_integer_columns`."""
        return self._narrower.narrow_integer_columns(df)

    def convert_columns_to_categorical(self, df: pd.DataFrame) -> pd.DataFrame:
        """See :meth:`Narrower.convert_columns_to_categorical`."""
        return self._narrower.convert_columns_to_categorical(df)

    def _check_labelfree(self) -> None:
        """Checks if the data is label-free."""
        if self.id_col not in self.data.columns:
            raise ValueError(
                f"id_col '{self.id_col}' is not a column of the loaded data. "
                f"Available columns: {sorted(self.data.columns)}"
            )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            extracted_matches = (
                self.data[self.id_col]
                .str.contains(self.label_group_capture, regex=True)
                .any()
            )

        if extracted_matches:
            self.label_free = False
            self.logger.info(
                "Data appears to be labelled. Proceeding with generating labelling information..."
            )
        else:
            self.label_free = True
            self.logger.info(
                "No labelling groups found in data. Data appears to be label-free or labels are in a custom format that is not recognized."
            )

    def drop_identical_cols(self, df: pd.DataFrame) -> pd.DataFrame:
        """Drops columns holding the same value in every row.

        Missing counts as a value, which reverses two cases that ``nunique``
        alone gets backwards: it ignores NA, so a column that is entirely
        empty counts 0 distinct values and survives, while one holding a
        single value in a handful of rows counts 1 and is dropped -- losing
        which rows those were.
        """
        cols_to_drop = [col for col in df.columns if df[col].nunique(dropna=False) == 1]
        if cols_to_drop:
            self.logger.info(
                f"Columns dropped for having identical values in all rows: {cols_to_drop}."
            )
        return df.drop(columns=cols_to_drop, errors="ignore")

    def rename_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Renames columns based on given alias mapping."""
        try:
            config_block = getattr(CONFIG, self.input_type)
            rename_mapping = config_block.COLS_RENAME_MAPPING
        except AttributeError:
            rename_mapping = {}

        if rename_mapping:
            return df.rename(columns=rename_mapping)
        else:
            return df

    def extra_info(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds extra information columns."""
        transformations: dict[str, tuple[list[str], Callable[..., pd.Series]]] = {
            "Leading.Razor.Protein": (
                ["Protein.Group"],
                # A regex split, which pandas runs per element, on a column
                # that repeats. Of the rest, RT.Width and Peptide.Length are
                # vectorized already, and Label.Free spans two columns, where
                # factorizing the combination costs about what it saves.
                per_distinct(lambda x: x.str.split(";").str[0]),
            ),
            "Peptide.Length": (["Stripped.Sequence"], lambda x: x.str.len()),
            "Label.Free": (
                ["Stripped.Sequence", "Precursor.Charge"],
                lambda x, y: x.astype(str) + y.astype(str),
            ),
            "RT.Width": (["RT.Stop", "RT.Start"], lambda x, y: x - y),
        }

        for new_col, (cols, func) in transformations.items():
            if all(col in set(df.columns) for col in cols):
                df[new_col] = func(*(df[col] for col in cols))

        return df

    def miscleavages(
        self,
        df: pd.DataFrame,
        seq_col: str = "Stripped.Sequence",
        protease: str = "Trypsin",
    ) -> pd.DataFrame:
        """Flags peptides whose residue counts are inconsistent with full cleavage."""
        rules = reference.protease(protease).allowed_counts

        def flag(seqs: pd.Series) -> pd.Series:
            terminal_aa = seqs.str[-1]
            fully_cleaved = pd.Series(False, index=seqs.index)
            for aa, count in rules.items():
                fully_cleaved |= (terminal_aa == aa) & (seqs.str.count(aa) == count)
            return ~fully_cleaved

        # Counting residues is per-element work in pandas and a peptide recurs
        # once per run per charge per channel, so this runs on the distinct
        # sequences: an order of magnitude less work on a multi-run report.
        df[f"{protease}.Miscleavages"] = per_distinct(flag)(df[seq_col])

        return df


class _LabelGenerator(Logged):
    """Generates label information for DIA-NN data.

    Every column derived here is a function of the precursor identifier alone,
    so the regex extraction and the pivots that follow it run over the
    *distinct* identifiers and the result is gathered back out. A report holds
    one row per identifier per run per channel, so on a 40-run experiment the
    regex extraction -- which pandas applies element by element in Python, on
    categorical columns too -- does a fortieth of the work.
    """

    __slots__ = (
        "data",
        "id_col",
        "id_codes",
        "label_group_capture",
        "extracted_matches",
        "sorted_matches",
        "UNIQUE_LABELS",
        "UNIQUE_CHANNELS",
        "labels_complete",
        "logger",
    )

    def __init__(self, processed_data: DataProcessor):
        """Initializes the LabelGenerator."""
        self.data = processed_data.data
        self.id_col = processed_data.id_col
        self.label_group_capture = processed_data.label_group_capture
        self.labels_complete = True

        if "Run" not in self.data.columns:
            raise ValueError(
                "Labelled data needs a 'Run' column to build channel identifiers; "
                f"available columns: {sorted(self.data.columns)}"
            )

        distinct_ids, self.id_codes = self._distinct_ids()
        self.extracted_matches = distinct_ids.str.extractall(self.label_group_capture)
        self.UNIQUE_LABELS = sorted(
            self.extracted_matches[0].str.split("-").str[0].unique()
        )

        sorted_matches = self._generate_sorted_matches()
        labelled_data = self._add_label_info(self.data, sorted_matches)
        self.data = self._generate_run_channels(labelled_data)
        self.logger.info("Data overwritten to include labelling information.")

    def _distinct_ids(self) -> tuple[pd.Series, np.ndarray]:
        """The distinct identifiers, and which one each row holds.

        Keyed by position rather than by the identifier itself: the frames
        derived from these get a plain integer index, so the pivots below sort
        a few thousand integers instead of a few hundred thousand long
        strings, and attaching the result is a positional gather rather than a
        join on those strings.
        """
        codes, uniques = pd.factorize(self.data[self.id_col])
        return pd.Series(pd.Index(uniques).astype(str)), codes

    def _validate_matrix_shape(
        self, matrix: pd.DataFrame | None, categorical: bool = True
    ) -> pd.DataFrame | None:
        """Validates the shape of the label matrix.

        `categorical` is off for the count matrix: it holds numbers, and a
        categorical of numbers no longer supports arithmetic.
        """
        if matrix is None or not matrix.index.is_unique:
            self.labels_complete = False
            self.logger.error(
                "Label matrix is not the expected shape, requires custom manipulation! Skipping labelling info generation..."
            )
            return None
        return matrix.astype("category") if categorical else matrix

    def _label_matrix(self, sorted_matches: pd.DataFrame) -> pd.DataFrame | None:
        """Generates the label matrix."""
        output_matrix = []
        for label in self.UNIQUE_LABELS:
            label_subset = sorted_matches[sorted_matches["Label"] == label]
            label_matrix = (
                label_subset.pivot_table(
                    index=["Index"], columns=["Location", "Offset"], aggfunc="size"
                )
                .fillna(0)
                .astype(int)
            )
            label_matrix.columns = label_matrix.columns.map("".join)
            name_matrix = np.char.multiply(
                np.array(label_matrix.columns, dtype=str), label_matrix.values
            )
            label_ids = ["".join(row) for row in name_matrix.tolist()]
            output_matrix.append(
                pd.DataFrame(
                    label_ids,
                    columns=[label + ".Label"],
                    index=label_matrix.index,
                    dtype="category",
                )
            )

        combined = pd.concat(output_matrix, axis=1)

        return self._validate_matrix_shape(combined)

    def _label_counts(self, sorted_matches: pd.DataFrame) -> pd.DataFrame | None:
        """Generates the label counts matrix."""
        label_count = sorted_matches.groupby(["Index", "Label"]).size().reset_index()
        label_count.columns = ["Index", "Label", "Size"]
        label_count["Count"] = label_count["Label"] + ".Count"
        label_count = pd.pivot(
            label_count, index="Index", columns="Count", values="Size"
        ).astype("Int8")

        return self._validate_matrix_shape(label_count, categorical=False)

    def _label_offset(
        self, sorted_matches: pd.DataFrame, label_count: pd.DataFrame | None
    ) -> pd.DataFrame | None:
        """Generates the label offset matrix."""
        if label_count is None:
            # The counts are needed to average the offsets; without them the
            # channel cannot be assigned. _label_counts has already reported it.
            return None

        sorted_matches["Offset"] = sorted_matches["Offset"].astype(int)
        offset_sum = sorted_matches.pivot_table(
            index="Index", columns="Label", values="Offset", aggfunc="sum"
        )
        offset_sum.columns = [col + ".Channel" for col in offset_sum.columns]

        if np.all(
            offset_sum.columns.str.split(".").str[0]
            == label_count.columns.str.split(".").str[0]
        ):
            # Dividing a frame by an array, and np.trunc of a frame, both keep
            # the frame at runtime; the numpy stubs say ndarray.
            offsets: pd.DataFrame = cast(
                pd.DataFrame, offset_sum / label_count.astype(float).values
            )
            rounded_offsets = cast(pd.DataFrame, np.trunc(offsets))

            if np.allclose(offsets, rounded_offsets, equal_nan=True):
                offsets = rounded_offsets.astype("Int8").astype(str)
            else:
                self.labels_complete = False
                self.logger.error(
                    "Channel offsets not uniform on peptide, cannot assign channel!"
                )
                return None
        else:
            self.labels_complete = False
            self.logger.error(
                "Label count and Label offset dfs do not have the same column order, cannot assign channel!"
            )
            return None

        return self._validate_matrix_shape(offsets)

    def _generate_sorted_matches(self) -> pd.DataFrame:
        """Generates sorted matches DataFrame."""
        sorted_matches = (
            self.extracted_matches[0].str.split("-", expand=True).reset_index()
        )
        sorted_matches.columns = ["Index", "Match", "Label", "Location", "Offset"]
        return sorted_matches

    def _add_label_info(
        self, df: pd.DataFrame, sorted_matches: pd.DataFrame
    ) -> pd.DataFrame:
        """Attaches the per-identifier label information to the data.

        The offsets are derived last: they need the counts to average over,
        and they convert a column of `sorted_matches` in place.
        """
        matrix = self._label_matrix(sorted_matches)
        counts = self._label_counts(sorted_matches)
        offsets = self._label_offset(sorted_matches, counts)

        blocks = [block for block in (matrix, counts, offsets) if block is not None]
        if not blocks:
            return df

        # One row per distinct identifier, gathered back out to one row per
        # data row. Identifiers with no group for a given label are absent
        # from the block and come through as NA, as they did when the blocks
        # were built per row and aligned on the index.
        info = pd.concat(blocks, axis=1).reindex(self.id_codes)
        info.index = df.index
        return pd.concat([df, info], axis=1)

    def _generate_run_channels(self, df: pd.DataFrame) -> pd.DataFrame:
        """Generates run channel information."""
        label_columns = [
            f"{label}.Channel"
            for label in self.UNIQUE_LABELS
            if f"{label}.Channel" in df.columns
        ]

        for label_column in label_columns:
            run_channel = df["Run"].astype(str) + "-" + df[label_column].astype(str)
            df[f"Run.{label_column}"] = run_channel.astype("category")

        # Concatenating column by column (rather than joining a row at a time)
        # keeps precursors with no channel for a label as NA instead of raising,
        # matching the per-label columns above.
        full_channel = df["Run"].astype(str)
        for label_column in label_columns:
            full_channel = full_channel + "-" + df[label_column].astype(str)
        df["Run.Full.Channel"] = full_channel.astype("category")

        return df


def narrow(
    frame: pd.DataFrame,
    *,
    input_type: str = "Unknown",
    floats: bool = True,
    round_large_floats: bool = False,
) -> pd.DataFrame:
    """Narrow a frame's dtypes, without the pipeline around them.

    Parameters
    ----------
    frame
        The frame to narrow. Modified in place and returned, as the steps are.
    input_type
        The engine the frame came from, which decides only which columns are
        kept out of the categorical conversion.
    floats
        Whether to narrow float64 to float32. Pass False when the caller will
        subtract nearly equal values from the frame afterwards; see
        :meth:`Narrower.narrow_float_columns`.
    round_large_floats
        Round large-magnitude float columns to integers, discarding their
        fractional part.
    """
    return Narrower(
        input_type=input_type,
        narrow_floats=floats,
        round_large_floats=round_large_floats,
    ).narrow(frame)

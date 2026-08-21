"""Higher-level processing pipelines and orchestrators.

This module composes loader, transformer and model utilities into
reusable processing functions and light-weight pipelines.
"""

import warnings
from collections.abc import Callable
from typing import cast

import numpy as np
import pandas as pd

from proteolyzer import reference

from .formats import Config
from .logging import Logged
from .models import Processing, Report
from .operations import per_distinct

CONFIG = Config()


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
    )

    def __init__(
        self,
        report: Report,
        id_col: str = "Precursor.Id",
        label_group_capture: str = r"\(((?:mTRAQ|SILAC|TMT)[^()]*)\)",
        protease: str = "Trypsin",
        round_large_floats: bool = False,
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
        """
        self.report = report
        self.data = report.frame.copy()
        self.input_type = report.source.input_type
        self.id_col = id_col
        self.label_group_capture = label_group_capture
        self.protease = protease
        self.round_large_floats = round_large_floats
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
            ),
        )

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

    def convert_columns_to_categorical(self, df: pd.DataFrame) -> pd.DataFrame:
        """Converts eligible low-cardinality columns to categorical type.

        Numeric columns are never converted, however few distinct values they
        hold: a categorical of numbers no longer supports arithmetic, so
        q-value or intensity columns would stop working downstream.
        """
        candidates = df.select_dtypes(exclude=["number", "bool"])
        cardinality = candidates.nunique() / len(df)
        cat_cols = list(cardinality[cardinality < CONFIG.CARDINALITY_THRESHOLD].index)

        try:
            config_block = getattr(CONFIG, self.input_type)
            exclude_from_conversion = config_block.EXCLUDE_CAT_CONVERSION
        except AttributeError:
            exclude_from_conversion = set()

        cat_cols = [col for col in cat_cols if col not in exclude_from_conversion]
        for col in cat_cols:
            df[col] = df[col].astype("category")

        self.logger.info(f"Converted columns: {cat_cols} to categorical dtype.")

        return df

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

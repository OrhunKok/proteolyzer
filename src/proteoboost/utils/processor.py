import pandas as pd
import numpy as np
import warnings
from ..utils.logging import MetaLogging
from proteoboost.utils.models import ProcessedData
import proteoboost.utils.constants as Constants
from proteoboost.utils.loader import DataLoader


class DataProcessor(metaclass=MetaLogging):
    """Processes raw data into a structured DataFrame."""

    __slots__ = (
        "data_loader",
        "data",
        "INPUT_TYPE",
        "logger",
        "ID_COL",
        "LABEL_GROUP_CAPTURE",
        "LABEL_FREE",
    )

    def __init__(
        self,
        data_loader: DataLoader,
        ID_COL: str = "Precursor.Id",
        LABEL_GROUP_CAPTURE: str = r"\((?!.*\bUniMod\b)(.*?)\)",
    ):
        """Initializes the DataProcessor."""
        self.data = data_loader.data
        self.INPUT_TYPE = data_loader.INPUT_TYPE
        self.ID_COL = ID_COL
        self.LABEL_GROUP_CAPTURE = LABEL_GROUP_CAPTURE
        self._check_labelfree()

    def process(self) -> ProcessedData:
        """Processes the data and returns a ProcessedData object."""
        self.data = (
            self.data.pipe(self.drop_identical_cols)
            .pipe(self.convert_float_columns_to_int)
            .pipe(self.rename_columns)
            .pipe(self.extra_info)
        )

        self.data = self.convert_columns_to_categorical(self.data)

        if not self.LABEL_FREE and self.INPUT_TYPE == "DIANN":
            self.data = _LabelGenerator(self).data

        self.data = self.miscleavages(self.data)

        self._memory_check(self.data)

        return ProcessedData(self)

    def _check_labelfree(self) -> None:
        """Checks if the data is label-free."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            extracted_matches = (
                self.data[self.ID_COL]
                .str.contains(self.LABEL_GROUP_CAPTURE, regex=True)
                .any()
            )

        if extracted_matches:
            self.LABEL_FREE = False
            self.logger.info(
                "Data appears to be labelled. Proceeding with generating labelling information..."
            )
        else:
            self.LABEL_FREE = True
            self.logger.info(
                "No labelling groups found in data. Data appears to be label-free."
            )

    def drop_identical_cols(self, df: pd.DataFrame) -> pd.DataFrame:
        """Drops columns with identical values."""
        cols_to_drop = [col for col in df.columns if df[col].nunique() == 1]
        if cols_to_drop:
            self.logger.info(
                f"Columns dropped for having identical values in all rows: {cols_to_drop}."
            )
        return df.drop(columns=cols_to_drop, errors="ignore")

    def convert_float_columns_to_int(self, df: pd.DataFrame) -> pd.DataFrame:
        """Converts eligible float columns to integer type."""
        float_cols = df.select_dtypes(include=["float"]).columns

        converted_cols = []
        for col in float_cols:
            columns_vals = df[col].values
            col_median = np.nanmedian(columns_vals)

            rounded_column_vals = np.round(
                pd.to_numeric(df[col], errors="coerce")
            ).astype("Int64")
            if (
                np.allclose(columns_vals, rounded_column_vals, equal_nan=True)
                or col_median > Constants.COL_MEDIAN_THRESHOLD
            ):
                converted_cols += [col]
                df[col] = np.round(pd.to_numeric(df[col], errors="coerce")).astype(
                    "Int64"
                )

        if converted_cols:
            self.logger.info(f"Converted columns: {converted_cols} to integer dtype.")

        return df

    def convert_columns_to_categorical(self, df: pd.DataFrame) -> pd.DataFrame:
        """Converts eligible columns to categorical type."""
        cardinality = df.apply(lambda x: x.nunique(), axis=0) / len(df)
        cat_cols = list(
            cardinality[cardinality < Constants.CARDINALITY_THRESHOLD].index
        )

        for col in cat_cols:
            df[col] = df[col].astype("category")

        self.logger.info(f"Converted columns: {cat_cols} to categorical dtype.")

        return df

    def rename_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Renames columns based on given alias mapping."""
        rename_mapping = Constants.COLS_RENAME_MAPPING.get(self.INPUT_TYPE)
        if rename_mapping:
            return df.rename(columns=rename_mapping)
        else:
            return df

    def extra_info(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds extra information columns."""
        transformations = {
            "Leading.Razor.Protein": (
                ["Protein.Group"],
                lambda x: x.str.split(";").str[0],
            ),
            "Peptide.Length": (["Stripped.Sequence"], lambda x: x.str.len()),
            "Label.Free": (
                ["Stripped.Sequence", "Precursor.Charge"],
                lambda x, y: x + y.astype(str),
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
        """Calculates miscleavages and adds a boolean column."""
        if protease not in Constants.PROTEASE_RULES:
            raise ValueError(
                f"Invalid protease: '{protease}'. Must be one of: {list(Constants.PROTEASE_RULES.keys())}"
            )

        rules = Constants.PROTEASE_RULES[protease]

        seqs = np.array(df[seq_col].values, dtype=str)
        terminal_aa = np.array([x[-1] for x in seqs])

        bool_arrays = np.ones(len(seqs), dtype=bool)
        for aa, count in rules.items():
            bool_arrays ^= (terminal_aa == aa) & (np.char.count(seqs, aa) == count)

        self.logger.info(
            f"{round(sum(bool_arrays) / len(bool_arrays) * 100, 1)}% Miscleavage rate according to {protease} rules"
        )
        df[f"{protease}.Miscleavages"] = bool_arrays

        return df


class _LabelGenerator(metaclass=MetaLogging):
    """Generates label information for DIA-NN data."""

    __slots__ = (
        "data",
        "ID_COL",
        "LABEL_GROUP_CAPTURE",
        "extracted_matches",
        "sorted_matches",
        "UNIQUE_LABELS",
        "UNIQUE_CHANNELS",
        "logger",
    )

    def __init__(self, processed_data: ProcessedData):
        """Initializes the LabelGenerator."""
        self.data = processed_data.data
        self.ID_COL = processed_data.ID_COL
        self.LABEL_GROUP_CAPTURE = processed_data.LABEL_GROUP_CAPTURE
        self.extracted_matches = self.data[self.ID_COL].str.extractall(
            self.LABEL_GROUP_CAPTURE
        )
        self.UNIQUE_LABELS = sorted(
            self.extracted_matches[0].str.split("-").str[0].unique()
        )

        sorted_matches = self._generate_sorted_matches()
        labelled_data = self._add_label_info(self.data, sorted_matches)
        self.data = self._generate_run_channels(labelled_data)
        self.logger.info("Data overwritten to include labelling information.")

    @property
    def label_checked_data(self):
        """Returns the label-checked data."""
        return self._label_checked_data

    def _validate_matrix_shape(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """Validates the shape of the label matrix."""
        if matrix is None or not matrix.index.is_unique:
            self.logger.error(
                "Label matrix is not the expected shape, requires custom manipulation! Skipping labelling info generation..."
            )
            return None
        return matrix.astype("category")

    def _label_matrix(self, sorted_matches: pd.DataFrame) -> pd.DataFrame:
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

        output_matrix = pd.concat(output_matrix, axis=1)

        return self._validate_matrix_shape(output_matrix)

    def _label_counts(self, sorted_matches: pd.DataFrame) -> pd.DataFrame:
        """Generates the label counts matrix."""
        label_count = sorted_matches.groupby(["Index", "Label"]).size().reset_index()
        label_count.columns = ["Index", "Label", "Size"]
        label_count["Count"] = label_count["Label"] + ".Count"
        label_count = pd.pivot(
            label_count, index="Index", columns="Count", values="Size"
        ).astype("Int8")

        return self._validate_matrix_shape(label_count)

    def _label_offset(
        self, sorted_matches: pd.DataFrame, label_count: pd.DataFrame
    ) -> pd.DataFrame:
        """Generates the label offset matrix."""
        sorted_matches["Offset"] = sorted_matches["Offset"].astype(int)
        offset_sum = sorted_matches.pivot_table(
            index="Index", columns="Label", values="Offset", aggfunc="sum"
        )
        offset_sum.columns = [col + ".Channel" for col in offset_sum.columns]

        if np.all(
            offset_sum.columns.str.split(".").str[0]
            == label_count.columns.str.split(".").str[0]
        ):
            offsets = offset_sum / label_count.astype(float).values
            rounded_offsets = np.trunc(offsets)

            if np.allclose(offsets, rounded_offsets, equal_nan=True):
                offsets = rounded_offsets.astype("Int8").astype(str)
            else:
                self.logger.error(
                    "Channel offsets not uniform on peptide, cannot assign channel!"
                )
                return None
        else:
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
        """Adds label information to the DataFrame."""
        label_matrix = self._label_matrix(sorted_matches)
        if label_matrix is not None:
            df = pd.concat([df, label_matrix], axis=1)

        label_counts = self._label_counts(sorted_matches)
        if label_counts is not None:
            df = pd.concat([df, label_counts], axis=1)

        label_offsets = self._label_offset(sorted_matches, label_counts)
        if label_offsets is not None:
            df = pd.concat([df, label_offsets], axis=1)

        return df

    def _generate_run_channels(self, df: pd.DataFrame) -> pd.DataFrame:
        """Generates run channel information."""
        for label in self.UNIQUE_LABELS:
            label_column = f"{label}.Channel"
            if label_column in df.columns:
                df[f"Run.{label}.Channel"] = (
                    df["Run"].astype(str) + "-" + df[label_column].astype(str)
                )
                df[f"Run.{label}.Channel"] = df[f"Run.{label}.Channel"].astype(
                    "category"
                )

        label_columns = [
            f"{label}.Channel"
            for label in self.UNIQUE_LABELS
            if f"{label}.Channel" in df.columns
        ]
        df["Run.Full.Channel"] = (
            df["Run"].astype(str)
            + "-"
            + df[label_columns].astype(str).agg("-".join, axis=1)
        )
        df["Run.Full.Channel"] = df["Run.Full.Channel"].astype("category")

        return df

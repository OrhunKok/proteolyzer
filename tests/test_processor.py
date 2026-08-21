import warnings

import numpy as np
import pandas as pd
import pytest

from proteolyzer.core.loader import DataLoader
from proteolyzer.core.models import Data
from proteolyzer.core.processor import DataProcessor


def _process(frame: pd.DataFrame, tmp_path, **kwargs):
    path = tmp_path / "report.parquet"
    frame.to_parquet(path)
    return Data(source=path, load_all_columns=True).load().process(**kwargs)


def test_label_free_data_is_detected(label_free_report, tmp_path):
    processed = _process(label_free_report, tmp_path)
    assert processed.LABEL_FREE is True
    assert not [c for c in processed.columns if c.endswith(".Channel")]


def test_derived_columns_are_added(label_free_report, tmp_path):
    processed = _process(label_free_report, tmp_path)
    assert "Peptide.Length" in processed.columns
    assert "RT.Width" in processed.columns
    assert "Leading.Razor.Protein" in processed.columns
    assert processed["RT.Width"].gt(0).all()


def test_labelled_data_gets_channel_columns(labelled_report, tmp_path):
    processed = _process(labelled_report, tmp_path)
    assert processed.LABEL_FREE is False
    assert "mTRAQ.Channel" in processed.columns
    assert "Run.mTRAQ.Channel" in processed.columns
    assert "Run.Full.Channel" in processed.columns


def test_unlabelled_precursor_does_not_break_channel_generation(
    labelled_report, tmp_path
):
    """Regression: joining a row containing NA used to raise TypeError."""
    processed = _process(labelled_report, tmp_path)
    full = processed.set_index("Precursor.Id")["Run.Full.Channel"]

    assert full["MYSEQK(mTRAQ-K-8)2"] == "run2-8"
    assert pd.isna(full["UNLABELLEDPEPTIDEK2"])


def test_missing_id_column_is_reported(label_free_report, tmp_path):
    with pytest.raises(ValueError, match="ID_COL 'Nope' is not a column"):
        _process(label_free_report, tmp_path, ID_COL="Nope")


@pytest.mark.parametrize(
    ("sequence", "expected_miscleavage"),
    [
        ("AAAK", False),  # one K, at the C-terminus
        ("AAAR", False),  # one R, at the C-terminus
        ("AKAK", True),  # a second K, so one is internal
        ("AAAG", True),  # does not end in K or R
        # Documents current behaviour: each cleavage residue is counted on its
        # own, so an internal K is not flagged on a peptide ending in R.
        ("AKAAR", False),
    ],
)
def test_miscleavage_flag(sequence, expected_miscleavage, tmp_path, label_free_report):
    frame = label_free_report.copy()
    frame.loc[0, "Stripped.Sequence"] = sequence
    processed = _process(frame, tmp_path)
    assert bool(processed.loc[0, "Trypsin.Miscleavages"]) is expected_miscleavage


def test_unknown_protease_lists_the_valid_ones(label_free_report, tmp_path):
    with pytest.raises(
        ValueError, match=r"Must be one of: \['ArgC', 'LysC', 'Trypsin'\]"
    ):
        _process(label_free_report, tmp_path, PROTEASE="Chymotrypsin")


def test_constant_columns_are_dropped(label_free_report, tmp_path):
    frame = label_free_report.copy()
    frame["Constant"] = "same"
    processed = _process(frame, tmp_path)
    assert "Constant" not in processed.columns


def test_integer_like_floats_are_narrowed(label_free_report, tmp_path):
    frame = label_free_report.copy()
    frame["Counts"] = np.arange(len(frame), dtype="float64")
    processed = _process(frame, tmp_path)

    assert pd.api.types.is_integer_dtype(processed["Counts"])
    # Narrowed to the smallest dtype that holds the values, not always Int64:
    # a nullable Int64 is wider than the float32 columns DIA-NN writes.
    assert processed["Counts"].dtype == np.int8
    assert processed["Counts"].tolist() == list(range(len(frame)))


def test_integer_like_floats_with_gaps_stay_nullable(label_free_report, tmp_path):
    frame = label_free_report.copy()
    frame["Counts"] = np.arange(len(frame), dtype="float64")
    frame.loc[0, "Counts"] = np.nan
    processed = _process(frame, tmp_path)

    assert processed["Counts"].dtype == "Int8"
    assert pd.isna(processed.loc[0, "Counts"])


def test_fractional_columns_keep_their_precision(label_free_report, tmp_path):
    """Regression: large-median float columns were silently rounded to integers."""
    frame = label_free_report.copy()
    frame["Ms1.Area"] = np.linspace(10.5, 1e7, len(frame), dtype="float32")
    processed = _process(frame, tmp_path)

    assert processed["Ms1.Area"].dtype == "float32"
    assert processed.loc[0, "Ms1.Area"] == pytest.approx(10.5)


def test_rounding_large_floats_is_opt_in(label_free_report, tmp_path, caplog):
    frame = label_free_report.copy()
    frame["Ms1.Area"] = np.linspace(150.5, 1e7, len(frame), dtype="float64")
    processed = _process(frame, tmp_path, ROUND_LARGE_FLOATS=True)

    assert pd.api.types.is_integer_dtype(processed["Ms1.Area"])
    assert processed.loc[0, "Ms1.Area"] == 150
    assert "discarding their fractional part" in caplog.text


@pytest.mark.parametrize("column", ["PEP", "Precursor.Charge", "Proteotypic"])
def test_numeric_columns_are_never_categorical(label_free_report, tmp_path, column):
    """Regression: q-values and charges became categorical and lost arithmetic."""
    processed = _process(label_free_report, tmp_path)

    assert not isinstance(processed[column].dtype, pd.CategoricalDtype)
    processed[column].mean()  # would raise on a categorical


def test_low_cardinality_strings_are_categorical(report_parquet):
    processor = DataProcessor(DataLoader(Data(source=report_parquet)))
    frame = pd.DataFrame(
        {
            "Condition": ["a", "b"] * 25,
            "Intensity": np.arange(50, dtype="float64"),
        }
    )

    converted = processor.convert_columns_to_categorical(frame)

    assert isinstance(converted["Condition"].dtype, pd.CategoricalDtype)
    # Same cardinality, but numeric: left alone.
    assert not isinstance(converted["Intensity"].dtype, pd.CategoricalDtype)


def test_label_counts_stay_numeric(labelled_report, tmp_path):
    """Label counts are numbers, so they must remain summable."""
    processed = _process(labelled_report, tmp_path)
    counts = processed["mTRAQ.Count"]
    assert not isinstance(counts.dtype, pd.CategoricalDtype)
    assert counts.sum() == 6  # five labelled precursors, one carrying two labels


def test_processed_metadata_is_propagated(label_free_report, tmp_path):
    processed = _process(label_free_report, tmp_path)
    assert processed.ID_COL == "Precursor.Id"
    assert processed.PROTEASE == "Trypsin"
    assert processed.unique_runs == {"run1", "run2"}


def test_processor_accepts_a_loader_directly(report_parquet):
    loader = DataLoader(Data(source=report_parquet))
    processed = DataProcessor(loader).process()
    assert len(processed) == 6


def test_all_nan_float_columns_are_left_alone(label_free_report, tmp_path):
    """Regression: narrowing an empty slice warned and produced a pointless Int64."""
    frame = label_free_report.copy()
    frame["Empty"] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        processed = _process(frame, tmp_path)

    assert processed["Empty"].dtype == "float64"
    assert processed["Empty"].isna().all()

import warnings

import numpy as np
import pandas as pd
import pytest

from proteolyzer.core.models import Data
from proteolyzer.core.processor import DataProcessor


def _process(frame: pd.DataFrame, tmp_path, **kwargs):
    """Process `frame` as if it had been read as a DIA-NN report."""
    path = tmp_path / "report.parquet"
    frame.to_parquet(path)
    report = Data(source=path, load_all_columns=True).load().process(**kwargs)
    return report.frame


def _report(frame: pd.DataFrame, tmp_path, **kwargs):
    """As `_process`, but keeping the Report so its metadata can be checked."""
    path = tmp_path / "report.parquet"
    frame.to_parquet(path)
    return Data(source=path, load_all_columns=True).load().process(**kwargs)


def test_label_free_data_is_detected(label_free_report, tmp_path):
    processed = _report(label_free_report, tmp_path)
    assert processed.processing.label_free is True
    assert not [c for c in processed.columns if c.endswith(".Channel")]


def test_derived_columns_are_added(label_free_report, tmp_path):
    processed = _process(label_free_report, tmp_path)
    assert "Peptide.Length" in processed.columns
    assert "RT.Width" in processed.columns
    assert "Leading.Razor.Protein" in processed.columns
    assert processed["RT.Width"].gt(0).all()


def test_labelled_data_gets_channel_columns(labelled_report, tmp_path):
    processed = _report(labelled_report, tmp_path)
    assert processed.processing.label_free is False
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


def test_a_precursor_measured_in_many_runs_is_labelled_in_every_row(
    tmp_path, report_frame
):
    """Label columns are derived per distinct identifier and gathered back out.

    Every row of a repeated precursor has to come back with that precursor's
    channel, not just the first one.
    """
    ids = ["MYSEQK(mTRAQ-K-8)2", "GGGGR(mTRAQ-R-0)2"] * 4
    runs = ["run1", "run1", "run2", "run2", "run3", "run3", "run4", "run4"]
    frame = report_frame(ids, runs)

    processed = _process(frame, tmp_path)

    by_id = processed.groupby("Precursor.Id", observed=True)
    assert by_id["mTRAQ.Channel"].nunique().eq(1).all()
    channels = processed.set_index(["Precursor.Id", "Run"])["mTRAQ.Channel"]
    assert channels["MYSEQK(mTRAQ-K-8)2"].eq("8").all()
    assert channels["GGGGR(mTRAQ-R-0)2"].eq("0").all()
    # And the run-specific column still varies with the run.
    assert processed["Run.mTRAQ.Channel"].nunique() == 8


def test_derived_columns_survive_repeated_values(tmp_path, report_frame):
    """Columns derived per distinct value must line up row by row."""
    ids = ["AKAK2", "MYSEQK2", "AKAK2", "VLDATRK3"]
    frame = report_frame(ids, ["run1", "run1", "run2", "run2"])
    frame["Protein.Group"] = ["P1;P2", "Q9", "P1;P2", "R4;R5"]
    frame["Stripped.Sequence"] = ["AKAK", "MYSEQK", "AKAK", "VLDATRK"]

    processed = _process(frame, tmp_path)

    assert processed["Leading.Razor.Protein"].tolist() == ["P1", "Q9", "P1", "R4"]
    # AKAK has an internal K; the other two end on their only cleavage residue.
    assert processed["Trypsin.Miscleavages"].tolist() == [True, False, True, False]


def test_missing_id_column_is_reported(label_free_report, tmp_path):
    with pytest.raises(ValueError, match="id_col 'Nope' is not a column"):
        _process(label_free_report, tmp_path, id_col="Nope")


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
        _process(label_free_report, tmp_path, protease="Chymotrypsin")


def test_constant_columns_are_dropped(label_free_report, tmp_path):
    frame = label_free_report.copy()
    frame["Constant"] = "same"
    processed = _process(frame, tmp_path)
    assert "Constant" not in processed.columns


def test_an_entirely_empty_column_is_dropped(label_free_report, tmp_path):
    """nunique() ignores NA, so an all-NA column used to count 0 and survive."""
    frame = label_free_report.copy()
    frame["Empty"] = pd.NA
    assert "Empty" not in _process(frame, tmp_path).columns


def test_a_sparsely_populated_column_is_kept(label_free_report, tmp_path):
    """One value plus gaps is not an identical column: which rows have it is data."""
    frame = label_free_report.copy()
    frame["Flag"] = pd.NA
    frame.loc[0, "Flag"] = "seen"

    processed = _process(frame, tmp_path)

    assert "Flag" in processed.columns
    assert processed["Flag"].notna().sum() == 1


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
    processed = _process(frame, tmp_path, round_large_floats=True)

    assert pd.api.types.is_integer_dtype(processed["Ms1.Area"])
    assert processed.loc[0, "Ms1.Area"] == 150
    assert "discarding their fractional part" in caplog.text


# ------------------------------------------------- single-precision narrowing


def test_float64_columns_are_narrowed_to_float32(label_free_report, tmp_path, caplog):
    """DIA-NN stores these as float32 itself, so the text path should match."""
    frame = label_free_report.copy()
    frame["Ms1.Area"] = np.linspace(10.5, 1e7, len(frame), dtype="float64")

    processed = _process(frame, tmp_path)

    assert processed["Ms1.Area"].dtype == "float32"
    assert processed["Ms1.Area"].iloc[0] == pytest.approx(10.5)
    assert "single precision" in caplog.text


def test_narrowing_floats_can_be_turned_off(label_free_report, tmp_path):
    frame = label_free_report.copy()
    frame["Ms1.Area"] = np.linspace(10.5, 1e7, len(frame), dtype="float64")

    report = _report(frame, tmp_path, narrow_floats=False)

    assert report.frame["Ms1.Area"].dtype == "float64"
    assert report.processing.narrowed_floats is False


@pytest.mark.parametrize(
    ("value", "reason"),
    [
        (1e300, "overflows to inf"),
        (1e-300, "underflows to zero"),
    ],
)
def test_values_outside_single_precision_keep_double(
    label_free_report, tmp_path, value, reason
):
    """Narrowing must round values, never replace them with inf or 0."""
    frame = label_free_report.copy()
    frame["Extreme"] = np.linspace(0.5, 1.5, len(frame))
    frame.loc[0, "Extreme"] = value

    processed = _process(frame, tmp_path)

    assert processed["Extreme"].dtype == "float64", reason
    assert processed.loc[0, "Extreme"] == value


def test_integer_columns_are_narrowed_to_the_width_they_need(
    label_free_report, tmp_path, caplog
):
    """Exact, so nothing is traded away: charges need one byte, not eight."""
    frame = label_free_report.copy()
    frame["Run.Index"] = np.arange(len(frame), dtype="int64")
    frame["Huge"] = np.full(len(frame), 10**10, dtype="int64")
    frame.loc[0, "Huge"] = 1

    processed = _process(frame, tmp_path)

    assert processed["Run.Index"].dtype == "int8"
    assert processed["Precursor.Charge"].dtype == "int8"
    # Derived in extra_info, so narrowing has to run after it.
    assert processed["Peptide.Length"].dtype == "int8"
    # Too large to narrow: left alone rather than wrapped.
    assert processed["Huge"].dtype == "int64"
    assert processed["Huge"].max() == 10**10
    assert "Narrowed integer columns" in caplog.text


def test_narrowed_integers_keep_their_values(label_free_report, tmp_path):
    frame = label_free_report.copy()
    frame["Counts"] = np.arange(-5, len(frame) - 5, dtype="int64")

    processed = _process(frame, tmp_path)

    assert processed["Counts"].tolist() == list(range(-5, len(frame) - 5))


def test_nullable_integers_stay_nullable_when_narrowed(report_parquet):
    processor = DataProcessor(Data(source=report_parquet).load())
    frame = pd.DataFrame({"Count": pd.Series([1, None, 3], dtype="Int64")})

    narrowed = processor.narrow_integer_columns(frame)

    assert narrowed["Count"].dtype == "Int8"
    assert narrowed["Count"].iloc[1] is pd.NA


def test_booleans_are_not_treated_as_integers(report_parquet):
    processor = DataProcessor(Data(source=report_parquet).load())
    frame = pd.DataFrame({"Flag": [True, False, True]})

    assert processor.narrow_integer_columns(frame)["Flag"].dtype == "bool"


def test_nullable_floats_keep_their_mask_when_narrowed(report_parquet):
    """Float64 narrows to Float32; numpy float32 would turn pd.NA into nan."""
    processor = DataProcessor(Data(source=report_parquet).load())
    frame = pd.DataFrame({"Quant": pd.Series([1.5, None, 2.5], dtype="Float64")})

    narrowed = processor.narrow_float_columns(frame)

    assert narrowed["Quant"].dtype == "Float32"
    assert narrowed["Quant"].iloc[1] is pd.NA


def test_derived_columns_are_computed_before_narrowing(label_free_report, tmp_path):
    """RT.Width is a difference of two close values, so order matters.

    Narrowing the retention times first would amplify their rounding by the
    ratio between a retention time and a peak width -- about 500x here.
    """
    frame = label_free_report.copy()
    frame["RT.Start"] = np.linspace(20.0, 30.0, len(frame))
    frame["RT.Stop"] = frame["RT.Start"] + 0.05

    processed = _process(frame, tmp_path)

    exact = 0.05
    error = abs(float(processed["RT.Width"].iloc[0]) - exact) / exact
    assert error < 1e-6, "RT.Width was computed from already-narrowed inputs"


@pytest.mark.parametrize("column", ["PEP", "Precursor.Charge", "Proteotypic"])
def test_numeric_columns_are_never_categorical(label_free_report, tmp_path, column):
    """Regression: q-values and charges became categorical and lost arithmetic."""
    processed = _process(label_free_report, tmp_path)

    assert not isinstance(processed[column].dtype, pd.CategoricalDtype)
    processed[column].mean()  # would raise on a categorical


def test_low_cardinality_strings_are_categorical(report_parquet):
    processor = DataProcessor(Data(source=report_parquet).load())
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


def test_conversion_is_decided_by_the_saving_not_the_cardinality(report_parquet):
    """A column can repeat too little for a ratio rule and still halve in size.

    Two thirds of the values here are distinct -- far above any workable
    cardinality threshold -- but the strings are long enough that replacing
    them with codes still saves most of the column.
    """
    processor = DataProcessor(Data(source=report_parquet).load())
    long_names = [f"sp|P{i:05d}|PROTEIN_{i}_HUMAN" for i in range(200)]
    frame = pd.DataFrame({"Protein": long_names * 3 + long_names[:100]})
    assert frame["Protein"].nunique() / len(frame) == pytest.approx(0.29, abs=0.01)

    converted = processor.convert_columns_to_categorical(frame)

    assert isinstance(converted["Protein"].dtype, pd.CategoricalDtype)


def test_columns_that_would_not_shrink_are_left_alone(report_parquet):
    """A column of distinct values gains nothing and can come out larger."""
    processor = DataProcessor(Data(source=report_parquet).load())
    frame = pd.DataFrame({"Precursor.Id": [f"PEPTIDEK{i}2" for i in range(500)]})

    converted = processor.convert_columns_to_categorical(frame)

    assert not isinstance(converted["Precursor.Id"].dtype, pd.CategoricalDtype)


def test_the_measured_saving_is_reported(report_parquet, caplog):
    processor = DataProcessor(Data(source=report_parquet).load())
    processor.convert_columns_to_categorical(pd.DataFrame({"Run": ["a", "b"] * 50}))
    assert "Run" in caplog.text and "%" in caplog.text


def test_label_counts_stay_numeric(labelled_report, tmp_path):
    """Label counts are numbers, so they must remain summable."""
    processed = _process(labelled_report, tmp_path)
    counts = processed["mTRAQ.Count"]
    assert not isinstance(counts.dtype, pd.CategoricalDtype)
    assert counts.sum() == 6  # five labelled precursors, one carrying two labels


def test_processing_metadata_is_recorded(label_free_report, tmp_path):
    report = _report(label_free_report, tmp_path)
    assert report.processing.id_col == "Precursor.Id"
    assert report.processing.protease == "Trypsin"
    assert report.processing.labels_complete is True
    assert report.runs == {"run1", "run2"}


def test_the_processor_can_be_driven_directly(report_parquet):
    report = Data(source=report_parquet).load()
    processed = DataProcessor(report).process()
    assert len(processed) == 6
    # The input report is untouched.
    assert not report.is_processed


def test_all_nan_float_columns_are_left_alone(label_free_report, tmp_path):
    """Regression: narrowing an empty slice warned and produced a pointless Int64.

    Driven through the method rather than the pipeline, which now drops an
    entirely empty column before the narrowing ever sees it.
    """
    frame = label_free_report.copy()
    frame["Empty"] = np.nan
    processor = DataProcessor(_report(frame, tmp_path))

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        narrowed = processor.convert_float_columns_to_int(frame)

    assert narrowed["Empty"].dtype == "float64"
    assert narrowed["Empty"].isna().all()

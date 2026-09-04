import io
import warnings

import numpy as np
import pandas as pd
import pytest

import proteolyzer as pz
from proteolyzer.core.models import Data, Report
from proteolyzer.core.processor import DataProcessor, _LabelGenerator


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


def test_a_spectronaut_report_processes_end_to_end(spectronaut_report, tmp_path):
    """The identifier the whole pipeline keys on is built as the file is read,
    because DataProcessor asks for it before any of its own steps run."""
    path = tmp_path / "20260901_164751_GluC_30min_Report.tsv"
    spectronaut_report.to_csv(path, sep="\t", index=False)

    processed = Data(source=path).load().process()

    assert processed.source.input_type == "Spectronaut"
    assert processed.processing.label_free is True
    assert processed.n_identifications == len(spectronaut_report)
    assert {"Peptide.Length", "Trypsin.Miscleavages"} <= set(processed.columns)
    assert processed.summary()["Precursors"].sum() == len(spectronaut_report)


def test_a_spectronaut_quantity_keeps_its_low_end(spectronaut_report, tmp_path):
    """2.54 and 400,000 share the column, so what a DIA-NN area could afford
    this cannot: rounding the low end throws away a fifth of it. Nothing has to
    be configured for that -- round_large_floats is off for every format -- and
    this pins the consequence for the one it would cost the most."""
    path = tmp_path / "20260901_164751_GluC_30min_Report.tsv"
    spectronaut_report.to_csv(path, sep="\t", index=False)

    quantity = Data(source=path).load().process().frame["Precursor.Quantity"]

    assert quantity.min() == pytest.approx(2.54)
    assert quantity.max() == pytest.approx(400000.0)


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


def test_narrowing_a_frame_needs_nothing_but_the_frame():
    """The four steps, without a Report or an id column; see #24."""
    frame = pd.DataFrame(
        {
            "Run": ["run1", "run2"] * 50,
            "Ms1.Area": np.linspace(1e5, 1e7, 100),
            "Precursor.Charge": np.full(100, 2.0),
            "PEP": np.linspace(0.001, 0.05, 100),
        }
    )
    narrowed = pz.narrow(frame, input_type="DIANN", floats=False)

    assert narrowed["Precursor.Charge"].dtype == np.int8
    assert isinstance(narrowed["Run"].dtype, pd.CategoricalDtype)
    # A number is never made categorical, however few values it takes, and a
    # fractional column is left as it is when floats are not being narrowed.
    assert narrowed["PEP"].dtype == np.float64


def test_narrowing_floats_is_optional():
    frame = pd.DataFrame({"PEP": np.linspace(0.001, 0.05, 10)})
    assert pz.narrow(frame.copy(), floats=True)["PEP"].dtype == np.float32
    assert pz.narrow(frame.copy(), floats=False)["PEP"].dtype == np.float64


def test_the_processor_runs_the_same_steps_as_the_narrower():
    """The pipeline goes through Narrower, so the two cannot drift apart."""
    frame = pd.DataFrame(
        {
            "Run": ["run1"] * 10,
            "Ms1.Area": np.linspace(1e5, 1e6, 10),
            "Precursor.Charge": np.full(10, 3.0),
        }
    )
    through_processor = DataProcessor(
        Report(
            frame=frame.copy(), source=Data(source=io.BytesIO(), INPUT_TYPE="DIANN")
        ),
        id_col="Run",
    )
    assert through_processor.narrow_integer_columns is not None
    assert (
        through_processor.convert_float_columns_to_int(frame.copy())[
            "Precursor.Charge"
        ].dtype
        == pz.narrow(frame.copy(), input_type="DIANN")["Precursor.Charge"].dtype
    )


def _label_generator(labelled_report, tmp_path) -> _LabelGenerator:
    """A label generator over labelled data, for the defensive paths below.

    These branches are what the labelling gives up on rather than guesses at, and
    valid input never reaches them -- so they are reached by handing the method
    the shape it refuses, which is the only way to see that it refuses.
    """
    path = tmp_path / "report.parquet"
    labelled_report.to_parquet(path)
    processor = DataProcessor(Data(source=path, load_all_columns=True).load())
    processor._check_labelfree()
    return _LabelGenerator(processor)


def test_a_matrix_with_a_repeated_index_is_refused(labelled_report, tmp_path, caplog):
    """One identifier cannot hold two rows of labelling: which one is the label?"""
    generator = _label_generator(labelled_report, tmp_path)
    repeated = pd.DataFrame({"mTRAQ.Channel": ["0", "8"]}, index=["A", "A"])

    with caplog.at_level("ERROR"):
        assert generator._validate_matrix_shape(repeated) is None

    assert generator.labels_complete is False
    assert "not the expected shape" in caplog.text


def test_no_matrix_at_all_is_refused(labelled_report, tmp_path):
    generator = _label_generator(labelled_report, tmp_path)

    assert generator._validate_matrix_shape(None) is None
    assert generator.labels_complete is False


def test_offsets_that_disagree_on_one_peptide_are_refused(
    labelled_report, tmp_path, caplog
):
    """Two rows of one identifier offset differently average to something that is
    not an offset -- 4 and 5 give 4.5, and no channel sits there."""
    generator = _label_generator(labelled_report, tmp_path)
    matches = pd.DataFrame(
        {"Index": [0, 0], "Label": ["mTRAQ", "mTRAQ"], "Offset": [4, 5]}
    )
    counts = pd.DataFrame({"mTRAQ.Count": [2]}, index=pd.Index([0], name="Index"))

    with caplog.at_level("ERROR"):
        assert generator._label_offset(matches, counts) is None

    assert generator.labels_complete is False
    assert "not uniform" in caplog.text


def test_counts_and_offsets_about_different_labels_are_refused(
    labelled_report, tmp_path, caplog
):
    """The offsets are divided by the counts positionally, so counting one label
    and offsetting another would divide the wrong pair and say nothing."""
    generator = _label_generator(labelled_report, tmp_path)
    matches = pd.DataFrame({"Index": [0], "Label": ["mTRAQ"], "Offset": [8]})
    counts = pd.DataFrame({"TMT.Count": [1]}, index=pd.Index([0], name="Index"))

    with caplog.at_level("ERROR"):
        assert generator._label_offset(matches, counts) is None

    assert generator.labels_complete is False
    assert "same column order" in caplog.text


def test_offsets_without_counts_give_up_quietly(labelled_report, tmp_path):
    """_label_counts has already said why; saying it twice is noise."""
    generator = _label_generator(labelled_report, tmp_path)
    matches = pd.DataFrame({"Index": [0], "Label": ["mTRAQ"], "Offset": [8]})

    assert generator._label_offset(matches, None) is None

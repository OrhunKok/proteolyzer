import dataclasses
import io
from pathlib import Path

import pandas as pd
import pytest

from proteolyzer.core.models import Data, Processing, Report


def test_string_source_is_coerced_to_path(report_parquet):
    data = Data(source=str(report_parquet))
    assert isinstance(data.source, Path)
    assert data.is_path
    assert not data.is_file_like


def test_missing_path_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="Path does not exist"):
        Data(source=tmp_path / "nope.parquet")


def test_non_readable_source_is_rejected():
    with pytest.raises(Exception, match="file-like"):
        Data(source=42)


def test_file_like_source_is_accepted():
    data = Data(source=io.StringIO("a\tb\n1\t2\n"))
    assert data.is_file_like
    assert data.file_name == "in_memory"
    assert data.file_extension == ""


def test_input_type_detected_from_name_and_extension(report_parquet):
    data = Data(source=report_parquet)
    assert data.file_name == "report"
    assert data.file_extension == ".parquet"
    assert data.input_type == "DIANN"


def test_unknown_input_type_loads_every_column(tmp_path, label_free_report):
    path = tmp_path / "something_else.parquet"
    label_free_report.to_parquet(path)
    data = Data(source=path)
    assert data.input_type == "Unknown"
    assert data.cols_subset is None
    assert data.cols_rename_mapping == {}


def test_cols_subset_comes_from_config(report_parquet):
    subset = Data(source=report_parquet).cols_subset
    assert "Precursor.Id" in subset
    assert "Q.Value" not in subset


def test_extra_cols_are_added_to_the_subset(report_parquet):
    data = Data(source=report_parquet, extra_cols_to_load="Q.Value")
    assert "Q.Value" in data.cols_subset
    assert "Precursor.Id" in data.cols_subset


@pytest.mark.parametrize("extra", [["Q.Value"], {"Q.Value"}, ("Q.Value",)])
def test_extra_cols_accept_any_string_collection(report_parquet, extra):
    assert Data(source=report_parquet, extra_cols_to_load=extra).cols_subset >= {
        "Q.Value"
    }


def test_extra_cols_reject_non_strings(report_parquet):
    with pytest.raises(Exception, match="Invalid input type"):
        Data(source=report_parquet, extra_cols_to_load=[1, 2])


def test_load_all_columns_overrides_the_subset(report_parquet):
    assert Data(source=report_parquet, load_all_columns=True).cols_subset is None


def test_file_stats_only_exist_for_paths(report_parquet):
    assert set(Data(source=report_parquet).file_stats) == {
        "Size (Bytes)",
        "Created",
        "Last Modified",
        "Last Accessed",
    }
    assert Data(source=io.StringIO("x")).file_stats is None


def test_manual_input_type_wins(report_parquet, caplog):
    data = Data(source=report_parquet, INPUT_TYPE="MaxQuant")
    assert data.input_type == "MaxQuant"
    assert "conflicts with file type" in caplog.text


def _report(frame, source, processed=True) -> Report:
    processing = (
        Processing(
            id_col="Precursor.Id",
            label_free=True,
            label_group_capture="",
            protease="Trypsin",
        )
        if processed
        else None
    )
    return Report(frame=frame, source=source, processing=processing)


def test_a_report_summarises_its_frame(report_parquet):
    frame = pd.DataFrame({"Run": ["a", "a", "b"], "Precursor.Id": ["x", "x", "y"]})
    report = _report(frame, Data(source=report_parquet))

    assert report.runs == {"a", "b"}
    assert report.n_identifications == 2
    assert report.processing.protease == "Trypsin"
    assert len(report) == 3
    assert list(report.columns) == ["Run", "Precursor.Id"]
    assert report["Run"].tolist() == ["a", "a", "b"]


def test_summaries_tolerate_missing_columns(report_parquet):
    report = _report(pd.DataFrame({"other": [1, 2]}), Data(source=report_parquet))
    assert report.runs == set()
    assert report.n_identifications == 0


def test_summary_counts_identifications_per_run(report_parquet):
    frame = pd.DataFrame(
        {
            "Run": ["a", "a", "a", "b"],
            "Precursor.Id": ["p1", "p1", "p2", "p1"],
            "Stripped.Sequence": ["PEPK", "PEPK", "SEQR", "PEPK"],
            "Leading.Razor.Protein": ["P1", "P1", "P1", "P1"],
        }
    )
    summary = _report(frame, Data(source=report_parquet)).summary()

    assert summary.index.tolist() == ["a", "b"]
    assert summary.loc["a"].to_dict() == {
        "Rows": 3,
        # p1 appears twice in run a -- distinct values, so it counts once.
        "Precursors": 2,
        "Peptides": 2,
        "Proteins": 1,
    }
    assert summary.loc["b", "Rows"] == 1


def test_summary_leaves_out_levels_the_frame_lacks(report_parquet):
    frame = pd.DataFrame({"Run": ["a", "b"], "Precursor.Id": ["p1", "p2"]})
    summary = _report(frame, Data(source=report_parquet)).summary()
    assert summary.columns.tolist() == ["Rows", "Precursors"]


def test_summary_falls_back_to_the_protein_group(report_parquet):
    """An unprocessed report has no razor protein yet."""
    frame = pd.DataFrame({"Run": ["a"], "Protein.Group": ["P1;P2"]})
    assert _report(frame, Data(source=report_parquet)).summary()[
        "Proteins"
    ].tolist() == [1]


def test_summary_of_a_frame_without_runs_is_a_single_group(report_parquet):
    frame = pd.DataFrame({"Precursor.Id": ["p1", "p2"]})
    summary = _report(frame, Data(source=report_parquet)).summary()
    assert summary.index.tolist() == ["all"]
    assert summary.loc["all", "Precursors"] == 2


def test_an_unprocessed_report_can_still_be_summarised(report_parquet):
    """Before processing, the identifier column falls back to the DIA-NN name."""
    summary = Data(source=report_parquet).load().summary()
    assert summary["Precursors"].sum() == 6


def test_an_unprocessed_report_has_no_identification_count(report_parquet):
    report = Data(source=report_parquet).load()
    assert not report.is_processed
    assert report.processing is None
    assert report.n_identifications == 0


def test_a_report_is_immutable_and_replaceable(report_parquet):
    """The metadata cannot drift away from the frame it describes."""
    report = Data(source=report_parquet).load()

    with pytest.raises(dataclasses.FrozenInstanceError):
        report.frame = pd.DataFrame()

    swapped = report.with_frame(report.frame.head(2))
    assert len(swapped) == 2
    assert len(report) == 6
    assert swapped.source is report.source


def test_load_returns_a_report(report_parquet):
    report = Data(source=report_parquet).load()
    assert isinstance(report, Report)
    assert isinstance(report.frame, pd.DataFrame)
    assert type(report.frame) is pd.DataFrame
    assert report.source.input_type == "DIANN"


def test_read_is_the_entry_point(report_parquet):
    import proteolyzer as pz

    assert isinstance(pz.read(report_parquet), Report)
    assert len(pz.read(report_parquet, load_all_columns=True).columns) > len(
        pz.read(report_parquet).columns
    )

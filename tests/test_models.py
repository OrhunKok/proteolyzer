import io
from pathlib import Path

import pandas as pd
import pytest

from proteolyzer.core.models import Data, LoadedData, ProcessedData


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


def test_processed_data_keeps_metadata_and_summarises():
    frame = pd.DataFrame(
        {"Run": ["a", "a", "b"], "Precursor.Id": ["x", "x", "y"]},
    )
    processed = ProcessedData(frame, ID_COL="Precursor.Id", PROTEASE="Trypsin")

    assert processed.unique_runs == {"a", "b"}
    assert processed.unique_ids == 2
    assert processed.PROTEASE == "Trypsin"
    # Metadata survives operations that go through _constructor.
    assert processed.head(2).ID_COL == "Precursor.Id"


def test_processed_data_summaries_tolerate_missing_columns():
    processed = ProcessedData(pd.DataFrame({"other": [1, 2]}), ID_COL="Precursor.Id")
    assert processed.unique_runs == set()
    assert processed.unique_ids == 0


def test_loaded_data_exposes_its_loader(report_parquet):
    loaded = Data(source=report_parquet).load()
    assert isinstance(loaded, LoadedData)
    assert loaded.loader.INPUT_TYPE == "DIANN"

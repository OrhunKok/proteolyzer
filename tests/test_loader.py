import io

import pandas as pd
import pyarrow
import pytest

from proteolyzer.utils.loader import DataLoader
from proteolyzer.utils.models import Data


def test_parquet_load_honours_the_column_subset(report_parquet):
    """Regression: the configured subset never reached the loader."""
    loaded = Data(source=report_parquet).load()
    assert "Precursor.Id" in loaded.columns
    assert "Q.Value" not in loaded.columns


def test_parquet_load_keeps_file_column_order(report_parquet, label_free_report):
    loaded = Data(source=report_parquet).load()
    expected = [c for c in label_free_report.columns if c in set(loaded.columns)]
    assert list(loaded.columns) == expected


def test_load_all_columns(report_parquet, label_free_report):
    loaded = Data(source=report_parquet, load_all_columns=True).load()
    assert set(loaded.columns) == set(label_free_report.columns)


def test_tsv_load_and_subset(tmp_path, label_free_report):
    path = tmp_path / "report.tsv"
    label_free_report.to_csv(path, sep="\t", index=False)
    loaded = Data(source=path).load()
    assert "Q.Value" not in loaded.columns
    assert len(loaded) == len(label_free_report)


def test_csv_delimiter_is_sniffed(tmp_path, label_free_report):
    path = tmp_path / "unknown_export.csv"
    label_free_report.to_csv(path, index=False)
    loaded = Data(source=path).load()
    # Unknown file name, so every column is loaded.
    assert set(loaded.columns) == set(label_free_report.columns)


def test_excel_load(tmp_path, label_free_report):
    pytest.importorskip("openpyxl")
    path = tmp_path / "sheet.xlsx"
    label_free_report.to_excel(path, index=False)
    loaded = Data(source=path).load()
    assert len(loaded) == len(label_free_report)


def test_plaintext_falls_back_to_one_row_per_line(tmp_path):
    path = tmp_path / "run.log"
    path.write_text("first\nsecond\nthird\n")
    loaded = Data(source=path).load()
    assert list(loaded.columns) == ["line"]
    assert len(loaded) == 3


def test_maxquant_txt_is_parsed_as_a_table(tmp_path):
    path = tmp_path / "evidence.txt"
    pd.DataFrame(
        {
            "Experiment": ["e1", "e2"],
            "Sequence": ["AAAK", "BBBR"],
            "Charge": [2, 3],
            "Intensity": [10, 20],
        }
    ).to_csv(path, sep="\t", index=False)

    data = Data(source=path)
    assert data.input_type == "MaxQuant"

    loaded = data.load()
    # Columns are renamed to the canonical proteolyzer names.
    assert "Run" in loaded.columns
    assert "Stripped.Sequence" in loaded.columns
    assert "Experiment" not in loaded.columns


def test_unnamed_stream_falls_back_to_plaintext(label_free_report):
    buffer = io.StringIO(label_free_report.to_csv(sep="\t", index=False))
    loaded = Data(source=buffer).load()
    # No name, so no extension to dispatch on: one row per line, header included.
    assert list(loaded.columns) == ["line"]
    assert len(loaded) == len(label_free_report) + 1


def test_named_stream_is_read_as_a_table(label_free_report):
    """Peeking at the header consumes a stream, so it has to be rewound."""
    buffer = io.StringIO(label_free_report.to_csv(sep="\t", index=False))
    buffer.name = "report.tsv"

    loader = DataLoader(Data(source=buffer))

    assert len(loader.data) == len(label_free_report)
    assert "Q.Value" not in loader.data.columns


def test_loader_reports_unreadable_files(tmp_path, caplog):
    path = tmp_path / "report.parquet"
    path.write_bytes(b"not a parquet file")
    with pytest.raises(pyarrow.ArrowInvalid):
        Data(source=path).load()
    assert "Error loading Parquet" in caplog.text

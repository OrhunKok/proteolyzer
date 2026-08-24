import io
from dataclasses import replace

import pandas as pd
import pyarrow
import pytest

from proteolyzer.core import loader, models
from proteolyzer.core.loader import DataLoader
from proteolyzer.core.models import Data


def test_parquet_load_honours_the_column_subset(report_parquet):
    """Regression: the subset the caller asks for never reached the loader."""
    loaded = Data(source=report_parquet, cols_to_load={"Precursor.Id"}).load()
    assert "Precursor.Id" in loaded.columns
    assert "Q.Value" not in loaded.columns


def test_parquet_load_keeps_every_column_when_none_are_asked_for(
    report_parquet, label_free_report
):
    """The core holds no list of its own, so nothing is dropped by default."""
    loaded = Data(source=report_parquet).load()
    assert set(loaded.columns) == set(label_free_report.columns)


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
    loaded = Data(source=path, cols_to_load={"Run", "Precursor.Id"}).load()
    assert "Q.Value" not in loaded.columns
    assert len(loaded) == len(label_free_report)


def test_tsv_load_keeps_file_column_order(tmp_path, label_free_report):
    """The fast parser returns columns in the order asked for, not file order.

    They agree only because the subset is built in file order; this pins that.
    """
    path = tmp_path / "report.tsv"
    label_free_report.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path).load()

    expected = [c for c in label_free_report.columns if c in set(loaded.columns)]
    assert list(loaded.columns) == expected


def test_tsv_load_matches_the_stock_parser(tmp_path, label_free_report):
    """The fast parser is an optimization, so it must not change the frame."""
    path = tmp_path / "report.tsv"
    label_free_report.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path, load_all_columns=True).load()
    expected = pd.read_csv(path, delimiter="\t")

    pd.testing.assert_frame_equal(loaded.frame, expected)


def test_a_file_too_large_for_memory_uses_the_stock_parser(
    tmp_path, label_free_report, caplog, monkeypatch
):
    """The fast parser needs several times the file's size; the other does not."""
    path = tmp_path / "report.tsv"
    label_free_report.to_csv(path, sep="\t", index=False)
    monkeypatch.setattr(loader, "_available_memory", lambda: 1024)

    loaded = Data(source=path, load_all_columns=True).load()

    assert "reading it with the stock parser" in caplog.text
    pd.testing.assert_frame_equal(loaded.frame, pd.read_csv(path, delimiter="\t"))


def test_a_file_that_fits_uses_the_fast_parser(
    tmp_path, label_free_report, caplog, monkeypatch
):
    path = tmp_path / "report.tsv"
    label_free_report.to_csv(path, sep="\t", index=False)
    monkeypatch.setattr(loader, "_available_memory", lambda: 1024**4)

    Data(source=path, load_all_columns=True).load()

    assert "stock parser" not in caplog.text


def test_a_stream_is_already_in_memory(label_free_report):
    """There is no file size to weigh, so the choice does not arise."""
    data = Data(source=io.StringIO("Run\tPrecursor.Id\nrun1\tp1\n"))
    assert DataLoader(data)._fast_read_fits() is True


def test_available_memory_falls_back_where_there_is_no_sysconf(monkeypatch):
    """Windows has no os.sysconf at all, which is what the guard is for."""
    monkeypatch.delattr(loader.os, "sysconf", raising=False)
    assert loader._available_memory() == loader.ASSUMED_AVAILABLE_MEMORY


def test_available_memory_falls_back_when_sysconf_refuses(monkeypatch):
    """Present but without the keys asked for, as on some BSDs."""

    def unavailable(name):
        raise ValueError(name)

    # raising=False: on a platform with no sysconf there is nothing to replace.
    monkeypatch.setattr(loader.os, "sysconf", unavailable, raising=False)
    assert loader._available_memory() == loader.ASSUMED_AVAILABLE_MEMORY


def test_ragged_rows_fall_back_to_the_stock_parser(tmp_path, caplog):
    """pyarrow rejects a short row; the stock parser pads it, as before."""
    path = tmp_path / "report.tsv"
    path.write_text("Run\tPrecursor.Id\nrun1\tp1\nrun2\n")

    loaded = Data(source=path, load_all_columns=True).load()

    assert loaded["Run"].tolist() == ["run1", "run2"]
    assert pd.isna(loaded["Precursor.Id"].iloc[1])
    assert "re-reading with the default parser" in caplog.text


def test_undecodable_bytes_are_not_passed_off_as_data(tmp_path, caplog):
    """pyarrow reads what it cannot decode as bytes instead of failing.

    Falling back means the caller gets the stock parser's error rather than a
    column of bytes objects. The bad byte has to sit beyond the header peek,
    which is where it would be in a real file.
    """
    path = tmp_path / "report.tsv"
    padding = b"".join(b"run1\tPEPTIDEK%d\n" % i for i in range(60_000))
    path.write_bytes(
        b"Run\tPrecursor.Id\n" + padding + "run1\tcaf\xe9\n".encode("latin-1")
    )

    with pytest.raises(UnicodeDecodeError):
        Data(source=path, load_all_columns=True).load()

    assert "undecoded" in caplog.text


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

    loader = DataLoader(Data(source=buffer, cols_to_load={"Run", "Precursor.Id"}))

    assert len(loader.data) == len(label_free_report)
    assert "Q.Value" not in loader.data.columns


def test_loader_reports_unreadable_files(tmp_path, caplog):
    path = tmp_path / "report.parquet"
    path.write_bytes(b"not a parquet file")
    with pytest.raises(pyarrow.ArrowInvalid):
        Data(source=path).load()
    assert "Error loading Parquet" in caplog.text


def test_a_jmod_table_is_recognised_and_subset(tmp_path, jmod_ids):
    """A format is a block on the config; nothing else has to know its name."""
    path = tmp_path / "filtered_IDs.csv"
    jmod_ids.to_csv(path, index=False)

    data = Data(source=path, cols_to_load={"file_name", "stripped_seq", "z", "rt"})
    assert data.input_type == "JMod"

    loaded = data.load()
    assert "unused_column" not in loaded.columns
    # Onto the canonical names, as every other engine is.
    assert {"Run", "Stripped.Sequence", "Precursor.Charge", "RT"} <= set(loaded.columns)
    assert len(loaded) == len(jmod_ids)


def test_a_jmod_parquet_reads_the_same_way(tmp_path, jmod_ids):
    path = tmp_path / "filtered_IDs.parquet"
    jmod_ids.to_parquet(path, index=False)

    loaded = Data(source=path).load()
    assert loaded.source.input_type == "JMod"
    assert "unused_column" in loaded.columns


def test_a_fragpipe_psm_table_is_recognised_and_subset(tmp_path, fragpipe_psms):
    path = tmp_path / "psm.tsv"
    fragpipe_psms.to_csv(path, sep="\t", index=False)

    data = Data(source=path)
    assert data.input_type == "FragPipe"

    loaded = data.load()
    # Recognised, renamed, and nothing withheld: the caller says what it wants.
    assert "Unused Column" in loaded.columns
    assert {"Run", "Stripped.Sequence", "Precursor.Charge", "RT"} <= set(loaded.columns)


def test_a_maxquant_table_is_read_whole_unless_asked_otherwise(tmp_path):
    """allPeptides runs to several GB, and narrowing it is the caller's call --
    the core cannot know which of its columns this project plots."""
    frame = pd.DataFrame(
        {
            "Raw file": ["run1", "run2"],
            "Retention time": [10.0, 20.0],
            "MS/MS count": [1, 2],
            "Total ion current": [1e6, 2e6],
            "Ion injection time": [10.0, 12.0],
            **{f"Unused {n}": [1.0, 2.0] for n in range(20)},
        }
    )
    path = tmp_path / "msScans.txt"
    frame.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path).load()
    assert [column for column in loaded.columns if column.startswith("Unused")]

    asked = Data(source=path, cols_to_load={"Raw file", "Retention time"}).load()
    # Onto the canonical names: 'Raw file' is 'Run' and 'Retention time' is 'RT'.
    assert set(asked.columns) == {"Run", "RT"}
    assert "Run" in loaded.columns


def test_a_file_matching_two_engines_says_which(tmp_path, monkeypatch):
    """The two-engine clash is general now, not a hand-written pair."""
    frame = pd.DataFrame({"Raw file": ["run1"]})
    path = tmp_path / "psm.tsv"
    frame.to_csv(path, sep="\t", index=False)

    clashing = replace(
        models.CONFIG,
        DIANN=replace(
            models.CONFIG.DIANN,
            FILES=[*models.CONFIG.DIANN.FILES, "psm"],
            FILE_EXTENSIONS=[*models.CONFIG.DIANN.FILE_EXTENSIONS, ".tsv"],
        ),
    )
    monkeypatch.setattr(models, "CONFIG", clashing)

    with pytest.raises(ValueError, match="matches multiple categories"):
        _ = Data(source=path).input_type


def test_cols_to_load_replaces_the_configured_subset(report_parquet):
    """A caller naming the columns its own project reads."""
    loaded = Data(source=report_parquet, cols_to_load={"Run", "Precursor.Id"}).load()
    assert set(loaded.columns) == {"Run", "Precursor.Id"}


def test_extras_widen_cols_to_load_and_both_lose_to_load_all(report_parquet):
    loaded = Data(
        source=report_parquet,
        cols_to_load={"Run"},
        extra_cols_to_load={"Q.Value"},
    ).load()
    assert set(loaded.columns) == {"Run", "Q.Value"}

    everything = Data(
        source=report_parquet, cols_to_load={"Run"}, load_all_columns=True
    ).load()
    assert "Q.Value" in everything.columns


def test_the_file_can_keep_its_own_column_names(tmp_path, label_free_report):
    """For a caller written against the engine's names rather than ours."""
    frame = label_free_report.rename(columns={"Run": "Experiment"})
    path = tmp_path / "evidence.txt"
    frame.to_csv(path, sep="\t", index=False)

    renamed = Data(source=path).load()
    assert "Run" in renamed.columns and "Experiment" not in renamed.columns

    as_written = Data(source=path, rename=False).load()
    assert "Experiment" in as_written.columns and "Run" not in as_written.columns


def test_extras_on_their_own_still_read_the_file_whole(tmp_path):
    """There is no base subset for them to be extra to. Answering with the extras
    alone would drop every other column without the caller having said so."""
    frame = pd.DataFrame(
        {
            "Raw file": ["run1"],
            "Kept": [1.0],
            **{f"Unused {n}": [1.0] for n in range(20)},
        }
    )
    path = tmp_path / "peptides.txt"
    frame.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path, extra_cols_to_load={"Raw file", "Kept"}).load()
    assert {"Run", "Kept"} <= set(loaded.columns)
    assert [column for column in loaded.columns if column.startswith("Unused")]

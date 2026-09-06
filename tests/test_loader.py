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


#: What Spectronaut writes an export as: the date, the time, the name of the
#: analysis, and only then the table.
STAMPED_REPORT = "20260901_164751_2026-08-27_CF_PD_GluC_30min_Report.tsv"

#: The same export as parquet, which is what Spectronaut writes by default. The
#: stem is identical, which is the point -- the pattern does not care.
STAMPED_PARQUET = "20260901_164751_2026-08-27_CF_PD_GluC_30min_Report.parquet"


def test_a_spectronaut_report_is_recognised_by_its_ending(tmp_path, spectronaut_report):
    """No fixed name could match one: Spectronaut stamps the export with the
    moment it was written and the name of the analysis."""
    path = tmp_path / STAMPED_REPORT
    spectronaut_report.to_csv(path, sep="\t", index=False)

    data = Data(source=path)
    assert data.input_type == "Spectronaut"

    loaded = data.load()
    assert {"Run", "Stripped.Sequence", "Precursor.Charge", "RT"} <= set(loaded.columns)
    assert len(loaded) == len(spectronaut_report)


def test_a_spectronaut_report_that_was_never_stamped_is_recognised_too(
    tmp_path, spectronaut_report
):
    path = tmp_path / "Report.tsv"
    spectronaut_report.to_csv(path, sep="\t", index=False)
    assert Data(source=path).input_type == "Spectronaut"


def test_a_diann_report_is_not_taken_for_a_spectronaut_one(tmp_path, label_free_report):
    """They differ by one capital letter and share an extension, so the pattern
    is matched case-sensitively. Two claimants would be refused outright."""
    path = tmp_path / "report.tsv"
    label_free_report.to_csv(path, sep="\t", index=False)
    assert Data(source=path).input_type == "DIANN"


def test_the_file_beside_a_spectronaut_report_is_not_the_report(tmp_path):
    """Spectronaut writes a `..._Report.setup.txt` next to it. The pattern is
    matched against the whole stem, so a prefix of one is not a match."""
    path = tmp_path / "20260901_164751_GluC_30min_Report.setup.txt"
    path.write_text("some setup\n")
    assert Data(source=path).input_type == "Unknown"


def test_a_spectronaut_precursor_identifier_is_built(tmp_path, spectronaut_report):
    """The report has no one column for it -- EG.PrecursorId is not in every
    export -- so it comes from the modified sequence and the charge."""
    path = tmp_path / STAMPED_REPORT
    spectronaut_report.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path).load()

    assert loaded["Precursor.Id"].tolist() == [
        f"_PEPTIDEK{i}_{charge}"
        for i, charge in enumerate([2, 3] * (len(spectronaut_report) // 2))
    ]


def test_an_export_that_carries_its_own_identifier_keeps_it(
    tmp_path, spectronaut_report
):
    """A report can be configured to write EG.PrecursorId, and one that did is
    taken at its word rather than overwritten."""
    spectronaut_report["EG.PrecursorId"] = [
        f"_PEPTIDEK{i}_.2" for i in range(len(spectronaut_report))
    ]
    path = tmp_path / STAMPED_REPORT
    spectronaut_report.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path).load()

    # Nothing maps EG.PrecursorId, so it arrives under its own name and the
    # built column is the one the rest of the package keys on.
    assert loaded["Precursor.Id"].iloc[0] == "_PEPTIDEK0_2"
    assert loaded["EG.PrecursorId"].iloc[0] == "_PEPTIDEK0_.2"


def test_nothing_is_built_for_a_caller_keeping_the_file_s_own_names(
    tmp_path, spectronaut_report
):
    """`Precursor.Id` is a name in the core's vocabulary, and rename=False asks
    for the file's. The dashboard reads this format that way."""
    path = tmp_path / STAMPED_REPORT
    spectronaut_report.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path, rename=False).load()

    assert "Precursor.Id" not in loaded.columns
    assert {"R.FileName", "EG.ModifiedSequence", "FG.Charge"} <= set(loaded.columns)


def test_a_column_the_identifier_needs_can_simply_not_be_there(
    tmp_path, spectronaut_report, caplog
):
    """A subset is an intersection, and so is this: the caller is told what could
    not be built rather than handed a column built out of half of it."""
    path = tmp_path / STAMPED_REPORT
    spectronaut_report.to_csv(path, sep="\t", index=False)

    loaded = Data(
        source=path, cols_to_load={"R.FileName", "EG.ModifiedSequence"}
    ).load()

    assert "Precursor.Id" not in loaded.columns
    assert "Not building Precursor.Id" in caplog.text


def test_a_spectronaut_subset_may_name_a_column_the_export_lacks(
    tmp_path, spectronaut_report
):
    """A report is configurable column by column, so a list written against one
    lab's export names columns another's does not have."""
    path = tmp_path / STAMPED_REPORT
    spectronaut_report.to_csv(path, sep="\t", index=False)

    loaded = Data(
        source=path,
        cols_to_load={
            "R.FileName",
            "FG.Quantity",
            # Neither is in this export, and neither is an error.
            "FG.MS1Quantity",
            "EG.IonMobility",
        },
    ).load()

    assert set(loaded.columns) == {"Run", "Precursor.Quantity"}


def test_spectronaut_column_names_carrying_spaces_and_brackets_survive(
    tmp_path, spectronaut_report
):
    path = tmp_path / STAMPED_REPORT
    spectronaut_report.to_csv(path, sep="\t", index=False)

    loaded = Data(source=path).load()

    assert "PG.Cscore (Run-Wise)" in loaded.columns
    assert "EG.TotalQuantity (Settings)" in loaded.columns


def test_a_spectronaut_report_read_from_an_upload(spectronaut_report):
    """An upload is read twice -- once for the header, once for the body -- and
    is dispatched on the name it carries rather than on a path."""
    buffer = io.StringIO(spectronaut_report.to_csv(sep="\t", index=False))
    buffer.name = STAMPED_REPORT

    data = Data(
        source=buffer, cols_to_load={"R.FileName", "EG.ModifiedSequence", "FG.Charge"}
    )
    assert data.input_type == "Spectronaut"

    loaded = data.load()
    assert len(loaded) == len(spectronaut_report)
    assert loaded["Precursor.Id"].iloc[0] == "_PEPTIDEK0_2"


def test_a_spectronaut_parquet_is_recognised_and_read(tmp_path, spectronaut_report):
    """Parquet is what Spectronaut writes by default, so it is the shape most
    exports actually arrive in."""
    path = tmp_path / STAMPED_PARQUET
    spectronaut_report.to_parquet(path, index=False)

    data = Data(source=path)
    assert data.input_type == "Spectronaut"

    loaded = data.load()
    assert {"Run", "Stripped.Sequence", "Precursor.Charge", "RT"} <= set(loaded.columns)
    assert loaded["Precursor.Id"].iloc[0] == "_PEPTIDEK0_2"
    assert len(loaded) == len(spectronaut_report)


def test_a_diann_parquet_is_not_taken_for_a_spectronaut_one(
    tmp_path, label_free_report
):
    """Both formats claim `.parquet` now, so the only thing keeping DIA-NN's
    `report.parquet` out of Spectronaut's hands is the capital letter. A file
    two blocks claimed would be refused outright rather than guessed at, so this
    failing would take DIA-NN's most common input down, not Spectronaut's."""
    path = tmp_path / "report.parquet"
    label_free_report.to_parquet(path, index=False)

    assert Data(source=path).input_type == "DIANN"


def test_the_two_spectronaut_serializations_read_the_same(tmp_path, spectronaut_report):
    """One report in two containers, not two formats. Whichever a lab exports,
    the frame that comes back is the same one -- names, values and dtypes."""
    tsv = tmp_path / STAMPED_REPORT
    parquet = tmp_path / STAMPED_PARQUET
    spectronaut_report.to_csv(tsv, sep="\t", index=False)
    spectronaut_report.to_parquet(parquet, index=False)

    pd.testing.assert_frame_equal(
        Data(source=tsv).load().frame, Data(source=parquet).load().frame
    )


def test_a_spectronaut_parquet_subset_may_name_a_column_the_export_lacks(
    tmp_path, spectronaut_report
):
    """A report is configurable column by column whichever way it is written,
    so the intersection has to hold on this path too."""
    path = tmp_path / STAMPED_PARQUET
    spectronaut_report.to_parquet(path, index=False)

    loaded = Data(
        source=path,
        cols_to_load={"R.FileName", "FG.Quantity", "FG.MS1Quantity", "EG.IonMobility"},
    ).load()

    assert set(loaded.columns) == {"Run", "Precursor.Quantity"}


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


# --- Spectronaut, as its parquet export actually arrives -----------------------
#
# Spectronaut has no default output name: whoever runs the analysis names the
# export, so the file below is called what a real one was called, and nothing
# about the name says what it is.

ANALYST_NAMED = "GluC-30min.parquet"


def test_a_spectronaut_parquet_is_recognised_by_its_columns_not_its_name(
    tmp_path, spectronaut_parquet_report
):
    """The name settles nothing -- an analyst chose it -- so the columns do.
    The same call the cellenONE reader makes: which file is which is worked out
    from the file, because names are unreliable."""
    path = tmp_path / ANALYST_NAMED
    spectronaut_parquet_report.to_parquet(path, index=False)

    data = Data(source=path)
    assert data.input_type == "Spectronaut"

    loaded = data.load()
    assert {"Run", "Stripped.Sequence", "Precursor.Charge", "RT"} <= set(loaded.columns)
    assert loaded["Precursor.Id"].iloc[0] == "_PEPTIDEK0_2"


def test_the_parquet_spelling_of_the_names_is_mapped_too(
    tmp_path, spectronaut_parquet_report
):
    """`R_FileName`, not `R.FileName`. Every canonical column the tab-separated
    export gives has to arrive from the parquet one as well, or the mapping is
    right about a file nobody exports."""
    path = tmp_path / ANALYST_NAMED
    spectronaut_parquet_report.to_parquet(path, index=False)

    loaded = Data(source=path).load()

    assert {
        "Run",
        "Modified.Sequence",
        "Stripped.Sequence",
        "Precursor.Charge",
        "Precursor.Mz",
        "Precursor.Quantity",
        "RT",
        "Predicted.RT",
        "Q.Value",
        "PEP",
        "Protein.Group",
        "PG.Q.Value",
        "Missed.Cleavages",
        "Proteotypic",
        "Decoy",
    } <= set(loaded.columns)
    # None of the mapped names survive under the file's own spelling ...
    assert not {"R_FileName", "EG_ModifiedSequence", "FG_Charge"} & set(loaded.columns)
    # ... and a column the schema has no canonical name for keeps the file's,
    # rather than being mangled towards one or dropped for lacking one.
    assert "EG_TotalQuantity_(Settings)" in loaded.columns


def test_a_tab_separated_export_under_an_analyst_s_name_is_recognised_too(
    tmp_path, spectronaut_report
):
    """Nothing about this is parquet's: the text export has no default name
    either, and its columns say the same thing."""
    path = tmp_path / "GluC-30min.tsv"
    spectronaut_report.to_csv(path, sep="\t", index=False)

    assert Data(source=path).input_type == "Spectronaut"


def test_one_familiar_column_is_not_enough_to_claim_a_file(tmp_path):
    """A signature of one would claim any frame somebody derived from a report
    and wrote back out. Two together are what no other engine's output has."""
    frame = pd.DataFrame({"FG_Charge": [2, 3], "Something.Else": [1.0, 2.0]})
    path = tmp_path / "derived.parquet"
    frame.to_parquet(path, index=False)

    assert Data(source=path).input_type == "Unknown"


def test_a_file_that_is_not_a_report_is_left_unknown(tmp_path):
    """Looking inside must not claim anything that merely has columns. A frame
    belonging to no search engine is still answered with Unknown rather than
    with whichever signature came closest."""
    path = tmp_path / "something_else.parquet"
    pd.DataFrame({"height": [1.0], "colour": ["red"]}).to_parquet(path, index=False)

    assert Data(source=path).input_type == "Unknown"


def test_peeking_at_a_stream_leaves_it_where_the_loader_expects_it(
    spectronaut_parquet_report,
):
    """Reading the columns consumes an upload, and the loader is about to read
    the same source from its start. An upload is what the dashboard passes."""
    buffer = io.BytesIO()
    spectronaut_parquet_report.to_parquet(buffer, index=False)
    buffer.name = ANALYST_NAMED
    buffer.seek(0)

    data = Data(source=buffer)
    assert data.input_type == "Spectronaut"

    loaded = data.load()
    assert len(loaded) == len(spectronaut_parquet_report)
    assert loaded["Precursor.Id"].iloc[0] == "_PEPTIDEK0_2"


def test_a_file_that_cannot_be_peeked_at_is_not_an_error(tmp_path):
    """Deciding which reader to use is not the place to raise about a file. The
    reader that follows will, and will say what it was trying to do."""
    path = tmp_path / "truncated.parquet"
    path.write_bytes(b"not a parquet file at all")

    assert Data(source=path).input_type == "Unknown"


# --- Columns a format writes as text that are not text -------------------------


def test_a_q_value_written_as_text_comes_back_a_number(
    tmp_path, spectronaut_text_typed_report
):
    """Spectronaut writes EG.Qvalue as '1.99e-13' and PG.Qvalue beside it as a
    double. A string q-value is not a lesser q-value: it raises TypeError the
    first time anyone filters on it, which is the first thing anyone does."""
    path = tmp_path / ANALYST_NAMED
    spectronaut_text_typed_report.to_parquet(path, index=False)

    frame = Data(source=path).load().frame

    assert frame["Q.Value"].dtype == "float64"
    assert (frame["Q.Value"] < 0.01).all()
    assert frame["Missed.Cleavages"].dtype == "int64"
    assert frame["Proteotypic"].dtype == "bool"


def test_an_identifier_that_looks_like_a_number_is_left_alone(
    tmp_path, spectronaut_text_typed_report
):
    """Every value of FG.XICDBID parses as a number and it is a database key.
    This is why the columns are a list and not a rule -- a rule gets this one
    wrong, and turning an identifier into an integer is quiet damage."""
    path = tmp_path / ANALYST_NAMED
    spectronaut_text_typed_report.to_parquet(path, index=False)

    frame = Data(source=path).load().frame

    assert frame["FG_XICDBID"].iloc[0] == "42775050"
    assert not pd.api.types.is_numeric_dtype(frame["FG_XICDBID"])


def test_a_column_that_does_not_convert_whole_is_left_as_it_came(
    tmp_path, spectronaut_text_typed_report, caplog
):
    """All or nothing. `errors="coerce"` would hand back a numeric column
    shorter by however many values nobody was told about."""
    frame = spectronaut_text_typed_report.copy()
    frame.loc[frame.index[0], "EG_Qvalue"] = "Filtered"
    path = tmp_path / ANALYST_NAMED
    frame.to_parquet(path, index=False)

    loaded = Data(source=path).load().frame

    assert not pd.api.types.is_numeric_dtype(loaded["Q.Value"])
    assert loaded["Q.Value"].iloc[0] == "Filtered"
    assert "not every value in them converts whole" in caplog.text


def test_the_text_a_file_uses_for_a_gap_becomes_a_real_gap(
    tmp_path, spectronaut_text_typed_report
):
    """`pd.isna` says False of the string 'NaN', so a gap reads as data until
    somebody plots it. The nullable dtype is for the gap, not for its own sake."""
    frame = spectronaut_text_typed_report.copy()
    frame.loc[frame.index[0], "EG_Qvalue"] = "NaN"
    path = tmp_path / ANALYST_NAMED
    frame.to_parquet(path, index=False)

    loaded = Data(source=path).load().frame

    assert pd.isna(loaded["Q.Value"].iloc[0])
    assert loaded["Q.Value"].iloc[1:].notna().all()


def test_the_dtype_is_given_back_under_the_file_s_own_names_too(
    tmp_path, spectronaut_text_typed_report
):
    """rename=False asks about names. A q-value that cannot be compared to a
    float is no more use under one name than another, and the dashboard reads
    every format this way."""
    path = tmp_path / ANALYST_NAMED
    spectronaut_text_typed_report.to_parquet(path, index=False)

    frame = Data(source=path, rename=False).load().frame

    assert frame["EG_Qvalue"].dtype == "float64"
    assert frame["PEP_IsProteotypic"].dtype == "bool"
    assert "Q.Value" not in frame.columns


def test_an_export_that_stored_them_properly_is_not_undone(
    tmp_path, spectronaut_parquet_report
):
    """Another lab's export may write these as numbers already, and a reader
    that converts unconditionally would be round-tripping them for nothing."""
    frame = spectronaut_parquet_report.copy()
    frame["EG_Qvalue"] = 1.5e-13
    path = tmp_path / ANALYST_NAMED
    frame.to_parquet(path, index=False)

    loaded = Data(source=path).load().frame

    assert loaded["Q.Value"].dtype == "float64"
    assert loaded["Q.Value"].iloc[0] == 1.5e-13


def test_a_text_typed_parquet_still_reads_as_its_own_tab_separated_export(
    tmp_path, spectronaut_text_typed_report
):
    """Which container a lab exported must not change the frame it gets back,
    even though one of them stores these columns as text and the other does
    not. Over every column the format claims -- see below for the one it does
    not, which is the only place the two still part."""
    parquet = tmp_path / ANALYST_NAMED
    tsv = tmp_path / "GluC-30min.tsv"
    spectronaut_text_typed_report.to_parquet(parquet, index=False)
    spectronaut_text_typed_report.to_csv(tsv, sep="\t", index=False)

    from_parquet = Data(source=parquet).load().frame.drop(columns="FG_XICDBID")
    from_tsv = Data(source=tsv).load().frame.drop(columns="FG_XICDBID")

    pd.testing.assert_frame_equal(from_parquet, from_tsv)


def test_an_all_digit_identifier_is_a_number_from_text_and_a_string_from_parquet(
    tmp_path, spectronaut_text_typed_report
):
    """The one column the two serializations still disagree about, recorded
    rather than left to be discovered.

    `FG.XICDBID` is a database key of digits, so the text reader infers an
    integer from it and parquet hands back what it stored. Its *values* survive
    either way -- there are no leading zeros in it to lose -- so this is a dtype
    to know about rather than damage, and it is not one of the names the schema
    maps. Naming it in a list of columns-that-are-text-whatever-they-look-like
    would be a second list to keep true for this one column; if a second turns
    up, that is the point to write it.
    """
    parquet = tmp_path / ANALYST_NAMED
    tsv = tmp_path / "GluC-30min.tsv"
    spectronaut_text_typed_report.to_parquet(parquet, index=False)
    spectronaut_text_typed_report.to_csv(tsv, sep="\t", index=False)

    from_parquet = Data(source=parquet).load().frame["FG_XICDBID"]
    from_tsv = Data(source=tsv).load().frame["FG_XICDBID"]

    assert not pd.api.types.is_numeric_dtype(from_parquet)
    assert pd.api.types.is_integer_dtype(from_tsv)
    assert from_parquet.astype("int64").tolist() == from_tsv.tolist()


# --- A renamed file is still the file it was -----------------------------------


def test_a_renamed_diann_report_is_still_a_diann_report(tmp_path, label_free_report):
    """People rename what they download. `report.parquet` is what DIA-NN writes
    it as, not what it has to be called for the next six months."""
    path = tmp_path / "2026-08-27_experiment_three.parquet"
    label_free_report.to_parquet(path, index=False)

    assert Data(source=path).input_type == "DIANN"


def test_a_renamed_maxquant_table_is_still_maxquant(tmp_path):
    path = tmp_path / "run_five_evidence_backup.txt"
    pd.DataFrame(
        {
            "Raw file": ["run1"],
            "Modified sequence": ["_AAAK_"],
            "Sequence": ["AAAK"],
            "Charge": [2],
            "Intensity": [10.0],
        }
    ).to_csv(path, sep="\t", index=False)

    data = Data(source=path)
    assert data.input_type == "MaxQuant"
    # and it is read as MaxQuant, not merely labelled one
    assert "Run" in data.load().columns


def test_a_renamed_jmod_table_is_still_jmod(tmp_path, jmod_ids):
    path = tmp_path / "ids_after_the_rerun.csv"
    jmod_ids.to_csv(path, index=False)

    data = Data(source=path)
    assert data.input_type == "JMod"
    assert {"Run", "Stripped.Sequence"} <= set(data.load().columns)


def test_a_renamed_fragpipe_table_is_still_fragpipe(tmp_path, fragpipe_psms):
    path = tmp_path / "psms_2026_08.tsv"
    fragpipe_psms.to_csv(path, sep="\t", index=False)

    data = Data(source=path)
    assert data.input_type == "FragPipe"
    assert {"Run", "Stripped.Sequence"} <= set(data.load().columns)


def test_a_name_that_matches_is_still_taken_at_its_word(
    tmp_path, label_free_report, monkeypatch
):
    """The columns are the fallback, not the first question. A file named the
    way its engine names it is claimed without being opened -- which is cheaper,
    and is how every release before this one behaved.

    Checked by making the peek fail: if it were consulted, this would raise.
    """
    path = tmp_path / "report.parquet"
    label_free_report.to_parquet(path, index=False)

    def unreachable(self):
        raise AssertionError("the file was opened to identify a name that matched")

    monkeypatch.setattr(models.Data, "peek_columns", unreachable)

    assert Data(source=path).input_type == "DIANN"


def test_a_misleading_name_still_wins_over_the_columns(tmp_path, fragpipe_psms):
    """Worth pinning because it is a real limit rather than an oversight: a
    FragPipe table renamed to something DIA-NN writes is read as DIA-NN, because
    the name is asked first and answers. Nothing here can tell the difference
    between a rename and a report, and preferring the columns would mean opening
    every file to find out. What it costs is a wrong rename mapping on a file
    somebody deliberately misnamed."""
    path = tmp_path / "report.tsv"
    fragpipe_psms.to_csv(path, sep="\t", index=False)

    assert Data(source=path).input_type == "DIANN"


def test_a_renamed_xic_export_is_recognised_by_its_whole_shape(tmp_path):
    """The table no overlap of names could claim. `pr`, `feature`, `rt` and
    `value` each belong to nobody in particular and `rt` is JMod's as well, so
    what says which it is, is that all four are there and almost nothing else."""
    path = tmp_path / "traces_for_run_one.parquet"
    pd.DataFrame(
        {"pr": ["PEPTIDEK2"], "feature": ["ms1"], "rt": [10.2], "value": [1.0]}
    ).to_parquet(path, index=False)

    assert Data(source=path).input_type == "DIANN"


def test_a_jmod_table_is_not_taken_for_an_xic_export(tmp_path, jmod_ids):
    """They share `rt`, which is exactly why the shape has to be the whole
    schema and not an overlap. JMod's is thirty-odd columns wide."""
    path = tmp_path / "ids_after_the_rerun.csv"
    jmod_ids.to_csv(path, index=False)

    assert Data(source=path).input_type == "JMod"


def test_a_frame_carrying_two_engines_columns_is_unknown_not_an_error(tmp_path, caplog):
    """Somebody joins a DIA-NN report to a MaxQuant table for a figure, saves
    it, and reads it back. Two *names* colliding is this package contradicting
    itself and raises; two signatures matching is a fact about that file, and
    refusing to load it would be worse than saying it is nobody's."""
    path = tmp_path / "joined_for_the_figure.parquet"
    pd.DataFrame(
        {
            "Precursor.Id": ["PEPK2"],
            "Stripped.Sequence": ["PEPK"],
            "Ms1.Area": [1.0],
            "Raw file": ["run1"],
            "Modified sequence": ["_PEPK_"],
        }
    ).to_parquet(path, index=False)

    data = Data(source=path)

    assert data.input_type == "Unknown"
    assert "cannot be told from them" in caplog.text
    assert len(data.load()) == 1


def test_a_name_two_formats_claim_is_still_an_error(tmp_path, monkeypatch):
    """The other half of that: a file two blocks *name* means the config
    contradicts itself, and nothing but raising gets that noticed."""
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


# --- Text that means "no value" ------------------------------------------------


def test_text_standing_for_a_gap_is_read_as_a_gap(tmp_path):
    """A delimited reader does this itself; parquet stores the string it was
    given. So the two disagreed about which cells were empty, and parquet lied
    the worse way round -- `pd.isna` says False of the string "NaN", so the gap
    read as data all the way to whatever plotted it."""
    path = tmp_path / "GluC-30min.parquet"
    pd.DataFrame(
        {
            "R_FileName": ["run1", "run2"],
            "EG_ModifiedSequence": ["_PEPK_", "_SEQR_"],
            "FG_Charge": [2, 3],
            "EG_IsVerified": ["NaN", "NaN"],
            "EG_InSourceFragmentationParentID": ["", "EG_12"],
        }
    ).to_parquet(path, index=False)

    frame = Data(source=path, rename=False).load().frame

    assert frame["EG_IsVerified"].isna().all()
    assert pd.isna(frame["EG_InSourceFragmentationParentID"].iloc[0])
    assert frame["EG_InSourceFragmentationParentID"].iloc[1] == "EG_12"


def test_a_word_that_could_be_a_value_is_left_alone(tmp_path):
    """`NA` and `None` are words, and on the export this was written from
    `EG.InSourceFragmentationClass` is "None" in 168,532 rows of 170,795 with a
    real class in the rest -- a category, not an absence. A text reader nulls
    both by default, so the two serializations still disagree about these two;
    agreeing means deciding pandas' default is right, which is a wider change
    than this one and belongs to whoever knows the domain."""
    path = tmp_path / "GluC-30min.parquet"
    pd.DataFrame(
        {
            "R_FileName": ["run1", "run2"],
            "EG_ModifiedSequence": ["_PEPK_", "_SEQR_"],
            "FG_Charge": [2, 3],
            "R_Fraction": ["NA", "NA"],
            "EG_InSourceFragmentationClass": ["None", "SecondaryFragment"],
        }
    ).to_parquet(path, index=False)

    frame = Data(source=path, rename=False).load().frame

    assert frame["R_Fraction"].notna().all()
    assert frame["EG_InSourceFragmentationClass"].iloc[0] == "None"


def test_closing_a_gap_leaves_a_column_that_has_none_alone(tmp_path, jmod_ids):
    """It runs on every read, so it has to cost nothing on a file with no text
    sentinels in it -- which is every file any other engine writes."""
    path = tmp_path / "filtered_IDs.csv"
    jmod_ids.to_csv(path, index=False)

    loaded = Data(source=path, rename=False).load().frame
    expected = pd.read_csv(path)

    pd.testing.assert_frame_equal(loaded, expected)


def test_a_word_is_not_a_gap_on_the_text_path_either(tmp_path):
    """pandas nulls `NA`, `None`, `null` and `N/A` by default, and on a real
    Spectronaut export `EG.InSourceFragmentationClass` is a three-state
    classification whose third state is the word `None` -- 168,532 rows of
    170,795, against 1,117 `Likely Parent` and 1,146 `Likely Child`. Taking the
    default erased that state into "not recorded" for 99% of the report."""
    path = tmp_path / "GluC-30min.tsv"
    path.write_text(
        "R.FileName\tEG.ModifiedSequence\tFG.Charge\tEG.InSourceFragmentationClass\n"
        "run1\t_VVEAHVDQKNKVVTTPAFMCE_\t3\tLikely Parent\n"
        "run1\t_AHVDQKNKVVTTPAFMCE_\t3\tLikely Child\n"
        "run1\t_VERVLKE_\t2\tNone\n"
    )

    frame = Data(source=path, rename=False).load().frame

    assert frame["EG.InSourceFragmentationClass"].tolist() == [
        "Likely Parent",
        "Likely Child",
        "None",
    ]
    assert frame["EG.InSourceFragmentationClass"].notna().all()


def test_the_two_serializations_agree_about_a_word(tmp_path):
    """The invariant the change is for: one report, two containers, one answer
    about what is a value and what is an absence."""
    rows = {
        "R_FileName": ["run1", "run1", "run1"],
        "EG_ModifiedSequence": ["_A_", "_B_", "_C_"],
        "FG_Charge": [2, 3, 2],
        "EG_InSourceFragmentationClass": ["Likely Parent", "Likely Child", "None"],
        "R_Fraction": ["NA", "NA", "NA"],
    }
    tsv, parquet = tmp_path / "x.tsv", tmp_path / "x.parquet"
    pd.DataFrame(rows).to_csv(tsv, sep="\t", index=False)
    pd.DataFrame(rows).to_parquet(parquet, index=False)

    pd.testing.assert_frame_equal(
        Data(source=tsv, rename=False).load().frame,
        Data(source=parquet, rename=False).load().frame,
    )


def test_the_machine_artefacts_are_still_gaps_on_the_text_path(tmp_path):
    """Narrowing what counts as missing must not stop counting the things that
    are not words: `NaN` is `str(float('nan'))`, `#N/A` is a spreadsheet's."""
    path = tmp_path / "GluC-30min.tsv"
    path.write_text(
        "R.FileName\tEG.ModifiedSequence\tFG.Charge\tEG.IsVerified\tnote\n"
        "run1\t_A_\t2\tNaN\t#N/A\n"
        "run1\t_B_\t3\t\t<NA>\n"
    )

    frame = Data(source=path, rename=False).load().frame

    assert frame["EG.IsVerified"].isna().all()
    assert frame["note"].isna().all()


def test_a_number_column_with_empty_fields_still_reads_as_numbers(tmp_path):
    """The obvious way to get this wrong. An empty field is still a gap, so a
    quantitative column with holes in it stays quantitative rather than turning
    into text the moment `keep_default_na` is switched off."""
    path = tmp_path / "report.tsv"
    path.write_text("Run\tPrecursor.Id\tMs1.Area\nrun1\tp1\t1000.5\nrun1\tp2\t\n")

    frame = Data(source=path).load().frame

    assert pd.api.types.is_float_dtype(frame["Ms1.Area"])
    assert frame["Ms1.Area"].iloc[0] == 1000.5
    assert pd.isna(frame["Ms1.Area"].iloc[1])

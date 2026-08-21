"""CoordinatesMapping parses cellenOne exports and maps cells to wells.

The fixtures below reproduce the shape of a cellenOne run directory: one
isolation table per condition plus the dispense, label and pickup logs.
"""

import json

import pandas as pd
import pytest

from proteolyzer.cellenone import CoordinatesMapping
from proteolyzer.core.io import read_frame

LOG_HEADER = "01.01.25-10:00:00.000\tHumidity\t45.0\tTemperature\t20.0\t10.0\t19.5\t4.0"


def _droplet_line(timestamp, plate, platepos, nozzle, well, target, field, x, y):
    """One 'drops' line as written by the instrument."""
    return "\t".join(
        [
            timestamp,
            str(plate),
            str(platepos),
            str(nozzle),
            well,
            str(target),
            "1",
            str(field),
            "3 drops",
            str(x),
            str(y),
        ]
    )


#: Two cells per pickup well, one well per nozzle. Logged coordinates are one
#: less than the geoprops coordinates, which the parser accounts for.
CELLS = [(100, 100), (105, 105), (136, 100), (140, 105)]


@pytest.fixture
def run_dir(tmp_path):
    root = tmp_path / "cellenone_files"
    cells = root / "Cells" / "Condition_1.Run"
    cells.mkdir(parents=True)
    (root / "Label").mkdir()
    (root / "pickup").mkdir()

    # Isolated cell geometry, tab separated.
    pd.DataFrame(
        {
            "Date": ["01.01.25"] * len(CELLS),
            "Time": [f"10:00:0{i}" for i in range(len(CELLS))],
            "Plate": [1] * len(CELLS),
            "Well": ["A1"] * len(CELLS),
            "Target": [1] * len(CELLS),
            "Field": [1] * len(CELLS),
            "XPos": [x for x, _ in CELLS],
            "YPos": [y for _, y in CELLS],
            "Diameter": [15.0, 16.0, 17.0, 18.0],
            "ImageFile": ["a.bmp"] * len(CELLS),
            "Background": [0] * len(CELLS),
        }
    ).to_csv(cells / "Reordered_run_isolated.xls", sep="\t", index=False)

    (cells / "run_dispense.log").write_text(LOG_HEADER + "\n")

    # One label droplet per cell, alternating between the two mTRAQ channels.
    (root / "Label" / "label.log").write_text(
        "\n".join(
            [LOG_HEADER]
            + [
                _droplet_line(
                    f"01.01.25-10:00:0{i}.000", 1, 1, 1, f"G{i + 1}", 1, 1, x - 1, y - 1
                )
                for i, (x, y) in enumerate(CELLS)
            ]
        )
        + "\n"
    )

    # Only the first nozzle is logged; the parser derives the second.
    (root / "pickup" / "pickup_1_.log").write_text(
        "\n".join(
            [
                LOG_HEADER,
                _droplet_line("01.01.25-10:00:04.000", 1, 1, 1, "A1", 1, 1, 99, 99),
            ]
        )
        + "\n"
    )
    return root


def test_files_are_discovered_and_parsed(run_dir):
    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)
    assert set(mapper.parsed_data) == {"Geoprops", "Dispense", "Label", "Pickup"}
    assert len(mapper.parsed_data["Geoprops"]) == len(CELLS)


def test_geoprops_gets_a_timestamp_and_sample_name(run_dir):
    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)
    geo = mapper.parsed_data["Geoprops"]
    assert geo["SampleName"].unique().tolist() == ["Condition_1"]
    assert pd.api.types.is_datetime64_any_dtype(geo["Timestamp"])


def test_labels_are_mapped_to_plex_channels(run_dir):
    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)
    assert mapper.parsed_data["Label"]["Plex"].tolist() == ["0", "8", "0", "8"]


def test_pickup_gains_the_second_nozzle(run_dir):
    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)
    pickup = mapper.parsed_data["Pickup"]
    # One logged nozzle plus its mirrored partner.
    assert len(pickup) == 2
    assert pickup["Nozzle"].tolist() == ["1", 3]
    assert pickup["Well"].tolist() == ["A1", "C1"]
    assert pickup["XPos"].tolist() == [100, 136]


def test_map_data_assigns_pickup_wells_to_cells(run_dir):
    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)
    metadata = mapper.map_data()

    assert len(metadata) == len(CELLS)
    # Each nozzle picks up the two cells closest to it.
    assert metadata["Well.Pickup"].tolist() == ["A1", "A1", "C1", "C1"]
    assert set(metadata["Plate.Pickup"]) == {1}
    assert metadata["Plex.Label"].tolist() == ["0", "8", "0", "8"]
    assert not metadata["Label.Clash"].any()


def test_map_data_warns_about_well_clashes(run_dir, caplog):
    """Without a plex, two cells in one well is a configuration error."""
    mapper = CoordinatesMapping(str(run_dir))
    metadata = mapper.map_data()
    assert "Well.Clash" in metadata.columns
    assert metadata["Well.Clash"].all()
    assert "Well clashes detected" in caplog.text


def test_wells_are_assigned_by_distance_not_by_row_order(run_dir):
    """Regression: every cell used to take the last pickup row's well.

    Without a plex each pickup row keeps every candidate cell, so the claims
    overlap; resolving them by distance makes the result independent of
    iteration order and equal to the correctly configured run.
    """
    without_plex = CoordinatesMapping(str(run_dir)).map_data()
    with_plex = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2).map_data()

    assert without_plex["Well.Pickup"].tolist() == ["A1", "A1", "C1", "C1"]
    assert without_plex["Well.Pickup"].tolist() == with_plex["Well.Pickup"].tolist()


def test_map_stats_concatenates_environment_readings(run_dir):
    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)
    stats = mapper.map_stats()
    assert set(stats["Type"]) == {"Dispense", "Label", "Pickup"}
    assert {"Humidity", "Temperature", "Bath Temp"} <= set(stats.columns)


def test_missing_directory_is_reported(tmp_path, caplog):
    empty = tmp_path / "empty"
    empty.mkdir()
    mapper = CoordinatesMapping(str(empty))
    assert mapper.parsed_data == {}
    assert mapper.map_data() is None
    assert "Missing required files" in caplog.text


# ------------------------------------------------------- file classification


@pytest.mark.parametrize(
    ("directory", "filename", "expected"),
    [
        # Names taken from a real prep, where the labelling logs are called
        # labels_1, labeling_2, L9, l4, L1_final and L1_dispense.
        ("label_1_20260122_110200_121.Run", "labels_1_..._Logfile.log", "Label"),
        ("label_2_20260122_111053_767.Run", "labeling_2_..._Logfile.log", "Label"),
        ("label_9_20260122_135702_078.Run", "L9_..._Logfile.log", "Label"),
        ("label_11_20260122_135951_429.Run", "L11_..._Logfile.log", "Label"),
        ("label_4_20260122_135356_597.Run", "l4_..._Logfile.log", "Label"),
        ("label_1_final_20260122_140824_305.Run", "L1_final_..._Logfile.log", "Label"),
        ("label_1_dispense_2026.Run", "L1_dispense_..._Logfile.log", "Label"),
        ("pickup_1_20260122_163057_008.Run", "pickup_1_..._Logfile.log", "Pickup"),
        ("pickup_2_2026.Run", "P2_..._Logfile.log", "Pickup"),
        ("cells_01fbs.Run", "01_fbs_cells_..._Logfile.log", "Dispense"),
        ("quench_2_2026.Run", "HA2_..._Logfile.log", "Dispense"),
        ("digest_2_2026.Run", "digest_2_..._Logfile.log", "Dispense"),
        ("cells_01fbs.Run", "Reordered_run_isolated.xls", "Geoprops"),
        ("cells_01fbs.Run", "Tscatter.xls", None),
        ("cells_01fbs.Run", "run_isolated.xls", None),
        ("cells_01fbs.Run", "Tscatter.bmp", None),
    ],
)
def test_files_are_classified_by_directory_and_name(directory, filename, expected):
    """Regression: the file name alone misfiled 5 of 10 labelling logs."""
    assert CoordinatesMapping._classify_file(directory, filename) == expected


def test_each_log_lands_in_exactly_one_bucket():
    """A name mentioning both steps must not be counted twice."""
    assert CoordinatesMapping._classify_file("pickup_1.Run", "label_pickup.log") == (
        "Pickup"
    )


# ------------------------------------------------------- imaging channels


CHANNEL_MEASUREMENTS = ["X", "Y", "Diameter", "Elongation", "Circularity", "Intensity"]


def _geoprops_frame(channels=("Transmission", "Green")):
    """Geoprops as exported: one row per (cell, imaging channel).

    The fluorescence row repeats the cell's identity but measures zero, which
    is what used to double the cell count and halve every geometry average.
    """
    rows = []
    for drop, (x, y) in enumerate(CELLS, start=1):
        for channel in channels:
            real = channel == "Transmission"
            rows.append(
                {
                    "DropNo": drop,
                    "X": 177.4 if real else 0.0,
                    "Y": 169.3 if real else 0.0,
                    "Diameter": 20.0 if real else 0.0,
                    "Elongation": 1.8 if real else 0.0,
                    "Circularity": 1.1 if real else 0.0,
                    "Intensity": 89.7 if real else None,
                    "Teg": channel,
                    "Target": 1,
                    "Field": 1,
                    "XPos": x,
                    "YPos": y,
                    "Date        ": "01.01.25",
                    "Time    ": f"10:00:0{drop}",
                }
            )
    return pd.DataFrame(rows)


@pytest.fixture
def geoprops_file(tmp_path):
    def write(channels=("Transmission", "Green")):
        run = tmp_path / "cells_a.Run"
        run.mkdir(exist_ok=True)
        path = run / "Reordered_a_isolated.xls"
        _geoprops_frame(channels).to_csv(path, sep="\t", index=False)
        return path

    return write


def _mapper() -> CoordinatesMapping:
    """A CoordinatesMapping without running __init__ (no directory needed)."""
    return CoordinatesMapping.__new__(CoordinatesMapping)


def test_imaging_channels_become_columns_not_rows(geoprops_file):
    parsed = _mapper().xls_parse([str(geoprops_file())])

    assert len(parsed) == len(CELLS)
    assert parsed["Teg"].unique().tolist() == ["Transmission"]
    for column in CHANNEL_MEASUREMENTS:
        assert column in parsed.columns
        assert f"{column}.Green" in parsed.columns


def test_geometry_is_not_diluted_by_the_fluorescence_rows(geoprops_file):
    """Regression: averaging over both channels halved every measurement."""
    parsed = _mapper().xls_parse([str(geoprops_file())])

    assert parsed["Diameter"].mean() == pytest.approx(20.0)
    assert (parsed["Diameter.Green"] == 0).all()


def test_identity_columns_are_not_suffixed(geoprops_file):
    parsed = _mapper().xls_parse([str(geoprops_file())])

    # Shared between the channels, so they stay single columns.
    for column in ("DropNo", "Target", "Field", "XPos", "YPos", "Date", "Time"):
        assert column in parsed.columns
        assert f"{column}.Green" not in parsed.columns


def test_single_channel_exports_are_unchanged(geoprops_file):
    parsed = _mapper().xls_parse([str(geoprops_file(channels=("Transmission",)))])

    assert len(parsed) == len(CELLS)
    assert not [col for col in parsed.columns if col.endswith(".Green")]


def test_export_without_the_channel_column(tmp_path):
    frame = _geoprops_frame(channels=("Transmission",)).drop(columns=["Teg"])
    path = tmp_path / "Reordered_a_isolated.xls"
    frame.to_csv(path, sep="\t", index=False)

    parsed = _mapper().xls_parse([str(path)])
    assert len(parsed) == len(CELLS)


def test_a_cell_missing_from_the_second_channel_is_kept(geoprops_file):
    path = geoprops_file()
    frame = pd.read_csv(path, sep="\t")
    frame = frame.drop(frame.index[(frame["Teg"] == "Green")][:1])
    frame.to_csv(path, sep="\t", index=False)

    parsed = _mapper().xls_parse([str(path)])

    assert len(parsed) == len(CELLS)
    assert parsed["Diameter.Green"].isna().sum() == 1


def test_ambiguous_channel_rows_are_left_alone(tmp_path, caplog):
    """Without a usable cell key, merging could mangle the data; don't."""
    frame = _geoprops_frame()
    frame["DropNo"] = 1  # every row claims to be the same drop
    frame[["Target", "Field", "XPos", "YPos"]] = 1
    path = tmp_path / "Reordered_a_isolated.xls"
    frame.to_csv(path, sep="\t", index=False)

    parsed = _mapper().xls_parse([str(path)])

    assert len(parsed) == len(frame)
    assert "cannot tell the imaging channels of a cell apart" in caplog.text


def test_optional_columns_may_be_absent(tmp_path):
    """Regression: the Plate/Well/ImageFile/Background drop was unconditional."""
    frame = _geoprops_frame(channels=("Transmission",)).drop(columns=["Teg"])
    path = tmp_path / "Reordered_a_isolated.xls"
    frame.to_csv(path, sep="\t", index=False)

    assert len(_mapper().xls_parse([str(path)])) == len(CELLS)


# ------------------------------------------------------------ pickup plates


@pytest.mark.parametrize(
    ("directory", "filename", "expected"),
    [
        ("pickup_1_2026.Run", "pickup_1_2026_Logfile.log", 1),
        # Recognized by its directory, so the number must be looked for there.
        ("pickup_7_2026.Run", "P7_2026_Logfile.log", 7),
        ("PICKUP_12_2026.Run", "whatever.log", 12),
    ],
)
def test_pickup_plate_is_read_from_name_or_directory(directory, filename, expected):
    path = f"/data/{directory}/{filename}"
    assert CoordinatesMapping._pickup_plate(path) == expected


def test_unnumbered_pickup_log_is_reported():
    """Regression: this raised a bare IndexError from the regex result."""
    with pytest.raises(ValueError, match="Cannot read the pickup plate number"):
        CoordinatesMapping._pickup_plate("/data/pickup_unnumbered.Run/weird.log")


# ------------------------------------------------------- repeated dispenses


def test_single_dispense_reports_its_droplets(run_dir):
    metadata = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2).map_data()

    # "3 drops" per line in the fixture.
    assert metadata["Droplets.Label"].tolist() == [3, 3, 3, 3]
    assert metadata["Dispenses.Label"].tolist() == [1, 1, 1, 1]
    # The raw text is replaced by the numbers.
    assert "Drops.Label" not in metadata.columns


def test_a_repeated_labelling_run_does_not_duplicate_cells(run_dir):
    """Regression: a re-dispensed channel multiplied the merged metadata."""
    label_log = run_dir / "Label" / "label.log"
    x, y = CELLS[0]
    label_log.write_text(
        label_log.read_text()
        + _droplet_line("01.01.25-11:00:00.000", 1, 1, 1, "G1", 1, 1, x - 1, y - 1)
        + "\n"
    )

    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)
    assert len(mapper.parsed_data["Label"]) == len(CELLS) + 1

    metadata = mapper.map_data()
    assert len(metadata) == len(CELLS)

    repeated = metadata.set_index("Plex.Label").loc["0"]
    first = repeated.iloc[0] if hasattr(repeated, "iloc") else repeated
    assert first["Dispenses.Label"] == 2
    assert first["Droplets.Label"] == 6  # 3 + 3


def test_the_latest_dispense_represents_the_position(run_dir):
    label_log = run_dir / "Label" / "label.log"
    x, y = CELLS[0]
    label_log.write_text(
        label_log.read_text()
        + _droplet_line("01.01.25-11:00:00.000", 1, 1, 2, "G1", 1, 1, x - 1, y - 1)
        + "\n"
    )

    metadata = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2).map_data()

    # Nozzle 2 dispensed last, so its row is the one kept.
    assert metadata.loc[metadata["Dispenses.Label"] == 2, "Nozzle.Label"].tolist() == [
        "2"
    ]


def test_an_unknown_droplet_count_stays_unknown():
    """Summing missing counts must not report zero droplets delivered."""
    label = pd.DataFrame(
        {
            "XPos": [1, 1],
            "YPos": [2, 2],
            "Target": [1, 1],
            "Field": [1, 1],
            "Well": ["G1", "G1"],
            "Timestamp": pd.to_datetime(["2026-01-01 10:00", "2026-01-01 11:00"]),
        }
    )

    collapsed = _mapper()._collapse_label_dispenses(label)

    assert len(collapsed) == 1
    assert collapsed["Dispenses"].tolist() == [2]
    assert collapsed["Droplets"].isna().all()


def test_positions_without_a_well_are_left_alone():
    """Without the channel column there is nothing to collapse on."""
    label = pd.DataFrame({"XPos": [1], "YPos": [2], "Target": [1], "Field": [1]})
    assert len(_mapper()._collapse_label_dispenses(label)) == 1


# --------------------------------------------------------------------- saving


def test_save_writes_the_metadata_and_a_run_record(run_dir, tmp_path):
    out = tmp_path / "results"
    mapper = CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2)

    assert mapper.save(out) == out
    assert {p.name for p in out.iterdir()} == {
        "metadata.parquet",
        "instrument_stats.parquet",
        "provenance.jsonl",
    }
    assert len(read_frame(out / "metadata")) == len(CELLS)


def test_the_run_record_says_how_the_metadata_was_produced(run_dir, tmp_path):
    """The metadata frame alone does not say what was parsed or configured."""
    out = tmp_path / "results"
    CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=2).save(out)

    entry = json.loads((out / "provenance.jsonl").read_text().splitlines()[0])

    assert entry["step"] == "CoordinatesMapping"
    assert entry["proteolyzer_version"]
    assert entry["params"] == {
        "root_dir": str(run_dir),
        "label_type": "mTRAQ",
        "plex": 2,
    }
    assert entry["summary"]["cells"] == len(CELLS)
    assert entry["summary"]["parsed_rows"]["Geoprops"] == len(CELLS)
    # Which logs were picked up as what: the audit trail for a misfiled run.
    assert len(entry["inputs"]["Label"]) == 1
    assert len(entry["inputs"]["Pickup"]) == 1


def test_the_run_record_counts_the_clashes(run_dir, tmp_path):
    out = tmp_path / "results"
    CoordinatesMapping(str(run_dir)).save(out)  # no plex, so every well clashes

    entry = json.loads((out / "provenance.jsonl").read_text().splitlines()[0])
    assert entry["summary"]["clashes"]["Well.Clash"] == len(CELLS)


def test_saving_twice_appends_a_second_record(run_dir, tmp_path):
    out = tmp_path / "results"
    for plex in (2, 3):
        CoordinatesMapping(str(run_dir), label_type="mTRAQ", plex=plex).save(out)

    records = (out / "provenance.jsonl").read_text().splitlines()
    assert [json.loads(r)["params"]["plex"] for r in records] == [2, 3]


def test_saving_nothing_is_refused(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()

    with pytest.raises(ValueError, match="Nothing to save"):
        CoordinatesMapping(str(empty)).save(tmp_path / "results")

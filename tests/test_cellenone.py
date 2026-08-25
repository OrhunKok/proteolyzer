"""Reading a cellenONE run directory whose names say nothing about its steps.

The files here are written to the shape the instrument writes them in — a header
that mentions drops without being a droplet, a 'Run ID' line, chamber readings,
then the dispensing — and given names that say nothing about which step they came
from, which is the case this is for.

Instrument output is not in the repository, so these files are a reconstruction
of it. A real run directory can be checked against as well, which is where test
runs, second attempts and steps named after their reagent turn up:

    PROTEOLYZER_CELLENONE_RUN=path/to/run_directory pytest tests/test_cellenone.py

Ported from streamlit-DO-MS, along with the implementation.
"""

import io
import os
import pathlib

import pandas as pd
import pytest

from proteolyzer.cellenone import UNREAD, CoordinatesMapping


class Upload(io.BytesIO):
    """An uploaded file: bytes with a relative path, as a browser hands them over."""

    def __init__(self, name, text):
        super().__init__(text.encode("latin1"))
        self.name = name


def chamber(minute, humidity=52.0):
    """A chamber reading, the line the log parser pulls the statistics out of."""
    stamp = f"21.01.26-1{minute // 60}:{minute % 60:02d}:00.000"

    return "\t".join(
        [stamp, "Humidity", str(humidity), "Temperature", "15.2", "9.1", "15.0", "24.0"]
    )


def dispensing(minute, well, target, field, x, y, plate="1", nozzle="1"):
    """A dispensing line: eleven columns, one of which says how many drops."""
    stamp = f"21.01.26-1{minute // 60}:{minute % 60:02d}:00.000"

    return "\t".join(
        [
            stamp,
            plate,
            "A1",
            nozzle,
            well,
            str(target),
            "1",
            str(field),
            "3 drops",
            str(x),
            str(y),
        ]
    )


def log(started, lines, run_id="", method="Hum_60/DP3/ScanFieldsOnly"):
    """A log file, headed the way the instrument heads one.

    Three of these lines mention drops without being a droplet: the column header
    the table below it never uses, the summary of the method that ran, and a
    parameter. All three have to be left out of the dispensing, and none of them
    starts with a time, which is what tells them apart.
    """
    header = "\t".join(
        [
            "Plate",
            "PlatePos",
            "Nozzle",
            "Well",
            "Target",
            "Level",
            "Field",
            "No. of Drops",
            "XPos",
            "YPos",
        ]
    )

    return (
        "\n".join(
            [
                "cellenONE_v2.0_nPOP\t091803",
                "2.2.2.1113\t11:24 AM - 11/15/2022",
                "",
                f"Run ID: {run_id}",
                "",
                f"Task Names: {method}/PreSpot250drops/BeginLoop/EndLoop",
                "PreSpot250drops=",
                header,
                chamber(started),
                chamber(started + 1),
            ]
            + lines
        )
        + "\n"
    )


#: Positions in a log are counted from zero and the parser adds one to line them
#: up with the table, which counts from one. The fixtures follow the instrument.
SLIDE = [(x, y) for x in range(3) for y in range(2)]
DROPLETS = 2 * len(SLIDE)

#: The columns of the isolated-cells table, in the order the instrument exports.
CELLS_HEADER = [
    "Teg",
    "DropNo",
    "Date",
    "Time",
    "XPos",
    "YPos",
    "Field",
    "Target",
    "Diameter",
    "Circularity",
    "Elongation",
    "Intensity",
    "X",
    "Y",
    "Plate",
    "Well",
    "ImageFile",
    "Background",
]


def dmso_log(started=0, run_id=""):
    # DMSO comes off a reagent well, one of a handful, onto the whole slide.
    return log(
        started,
        [
            dispensing(started + 2, "D1", 1, field, x, y)
            for field in (1, 2)
            for x, y in SLIDE
        ],
        run_id=run_id,
    )


def label_log(started=200, run_id=""):
    # Labelling draws from the row of the reagent plate the label mapping covers.
    return log(
        started,
        [
            dispensing(started + 2, f"G{1 + (x + y) % 3}", 1, field, x, y)
            for field in (1, 2)
            for x, y in SLIDE
        ],
        run_id=run_id,
    )


def quench_log(started=400, run_id="", method="DewPointControl/Take_TEAB/mTRAQ_quench"):
    # The quench log names the label it is quenching in its method summary, which
    # is how a search of the whole text lands on the wrong step.
    return log(
        started,
        [
            dispensing(started + 2, "D3", 1, field, x, y)
            for field in (1, 2)
            for x, y in SLIDE
        ],
        run_id=run_id,
        method=method,
    )


def pickup_log(started=500, run_id=""):
    # One cell to one well, so the wells are as many as the lines and spread over
    # a plate.
    wells = [f"{letter}{number}" for letter in "ABEF" for number in range(1, 4)]

    return log(
        started,
        [
            dispensing(started + 2, well, 1, 1 + index // 6, *SLIDE[index % 6])
            for index, well in enumerate(wells)
        ],
        run_id=run_id,
    )


def pickup_field(field, wells, started=500, run_id=""):
    """One pickup file, all of it landing on a single field of the slide.

    Two files calling this with the same well names and different fields are
    two destination plates: the same well letter, placed by two different
    dispensing runs.
    """
    return log(
        started,
        [
            dispensing(started + 2, well, 1, field, x, y)
            for well, (x, y) in zip(wells, SLIDE, strict=True)
        ],
        run_id=run_id,
    )


def cells_log(started=100, run_id=""):
    """The imaging log: chamber readings and no dispensing at all."""
    return log(started, [], run_id=run_id)


def empty_log(run_id=""):
    """A run stopped before it started: no readings, no droplets, nothing."""
    return (
        "\n".join(
            [
                "cellenONE_v2.0_nPOP\t091803",
                f"Run ID: {run_id}",
                "Task Names: PreSpot250drops",
                "Aborted",
            ]
        )
        + "\n"
    )


def cells_table(channels=("Transmission",)):
    """The isolated, reordered table, by its columns rather than by its name.

    A row per cell per imaging channel, which is how the instrument writes it: a
    slide imaged in transmission and one fluorescence channel has two rows for
    every cell on it.
    """
    rows = [CELLS_HEADER]
    for field in (1, 2):
        for x in range(1, 4):
            for y in range(1, 3):
                for index, channel in enumerate(channels):
                    rows.append(
                        [
                            channel,
                            "1",
                            "21.01.2026",
                            "11:00:00",
                            str(x),
                            str(y),
                            str(field),
                            "1",
                            str(18.2 + index),
                            str(0.93),
                            str(1.1),
                            str(120 + index * 40),
                            "100",
                            "200",
                            "1",
                            "A1",
                            f"img_{channel}.png",
                            "10",
                        ]
                    )

    return "\n".join("\t".join(row) for row in rows) + "\n"


def droplet_table():
    """Every droplet the instrument saw, of which the placed ones are a handful.

    The same columns as the cells table, written alongside it, and the file that
    must not be read as the cells: a row per droplet detected, with no cell
    placed on the slide.
    """
    lines = cells_table().splitlines()
    header, placed = lines[0], lines[1:]
    # A droplet that was never placed stops *before* the columns that would say
    # where it went. It is a short row, not a full row carrying blanks, and a real
    # geoprops export is short rows almost to the last: one had 6797 of them in
    # the read window against 27 full ones.
    #
    # That distinction is the whole test. Writing these padded to full width is
    # what let a reader discard every short row as malformed and then take its
    # placed-row ratio over what was left -- which is the placed rows, giving 1
    # for this file and 1 for the cells, so a reordered geoprops export was read
    # as a second table of cells in every folder of a real run while this suite
    # stayed green. See DECISIONS.md.
    stub = "\t".join(["Transmission", "2", "21.01.2026", "11:00:00"])

    return "\n".join([header] + [stub] * (len(placed) * 200) + placed) + "\n"


def sorted_upload(files, stages=None):
    """What the sorter made of an upload, as {stage: [names]}."""
    mapping = CoordinatesMapping.__new__(CoordinatesMapping)
    mapping.root_dir = files
    mapping.label_type, mapping.plex = "mTRAQ", 3
    mapping.given = dict(stages or {})

    return (
        {
            stage: [file.name for file in group]
            for stage, group in mapping._output_file_paths().items()
        },
        mapping.stages,
    )


@pytest.fixture
def named_run():
    """A directory whose folder names say which step each file came from."""
    return [
        Upload("dmso.Run/dmso_1.log", dmso_log()),
        Upload("cells.Run/cells_1.log", cells_log()),
        Upload("cells.Run/sample_isolated_reordered.xls", cells_table()),
        Upload("label.Run/label_1.log", label_log()),
        Upload("quench.Run/quench_1.log", quench_log()),
        Upload("pickup.Run/pickup_1_a.log", pickup_log()),
    ]


@pytest.fixture
def blind_run():
    """The same run copied off the instrument and renamed: numbered, silent."""
    return [
        Upload("01/output.log", dmso_log()),
        Upload("02/output.log", cells_log()),
        Upload("02/table.xls", cells_table()),
        Upload("03/output.log", label_log()),
        Upload("04/output.log", quench_log()),
        Upload("05/output.log", pickup_log()),
    ]


def test_the_names_are_read_first_where_there_are_names(named_run):
    groups, _ = sorted_upload(named_run)

    assert {stage: len(files) for stage, files in groups.items()} == {
        "Geoprops": 2,
        "DMSO": 1,
        "Label": 1,
        "Quench": 1,
        "Pickup": 1,
    }


def test_every_step_is_found_with_nothing_in_the_names(blind_run):
    _, stages = sorted_upload(blind_run)

    # The cells table by its columns, and the log beside it comes with it.
    assert stages.get("02/table.xls") == "Geoprops"
    assert stages.get("02/output.log") == "Geoprops"
    # Pickup by its wells, one to a cell.
    assert stages.get("05/output.log") == "Pickup"
    # Labelling by the reagent row it draws from.
    assert stages.get("03/output.log") == "Label"
    # DMSO by going down before the label, quench by going down after it.
    assert stages.get("01/output.log") == "DMSO"
    assert stages.get("04/output.log") == "Quench"


def test_a_log_that_names_its_own_run_is_taken_at_its_word():
    """'Run ID:' is the folder name over again, so a renamed folder still says
    what it held. The method summary beside it is not: it lists every reagent the
    run touched, and the quench log names the label it is quenching."""
    spoken = [
        Upload("a.log", dmso_log(run_id="dmso_20260316_143133_606")),
        Upload("b.log", quench_log(run_id="teab_20260317_104444_710")),
        Upload("c.log", dmso_log(600, run_id="digest_argc_fields15to16")),
    ]
    _, stages = sorted_upload(spoken)

    assert stages["a.log"] == "DMSO"
    assert stages["b.log"] == "Quench"
    # A step with no plots anywhere stays a plain dispense record.
    assert stages["c.log"] == "Dispense"


def test_a_step_that_cannot_be_named_is_left_as_a_dispense_record():
    """Nothing in the name, nothing in the log, wells that say nothing, and DMSO
    and quench both already found: there is no answer left to give."""
    crowded = [
        Upload("01/dmso.log", dmso_log()),
        Upload("02/quench.log", quench_log()),
        Upload("03/output.log", dmso_log(700)),
    ]
    _, stages = sorted_upload(crowded)

    assert stages["03/output.log"] == "Dispense"


def test_the_files_that_are_not_the_cells_table_are_not_read_as_it():
    littered = [
        Upload("01/output.log", dmso_log()),
        Upload("02/table.xls", cells_table()),
        Upload("02/droplets.xls", droplet_table()),
        Upload(
            "02/Fscatter.xls",
            "Diameter (green) - Detected\tIntensity (green)\n8.35\t113.62\n",
        ),
        Upload("02/output.log", cells_log()),
        Upload("06/output.log", empty_log(run_id="pickup_test_20260317_124203")),
    ]
    _, stages = sorted_upload(littered)

    # The table of every droplet seen is not the table of cells placed, and the
    # scatter exports are not either.
    assert stages["02/droplets.xls"] == UNREAD
    assert stages["02/Fscatter.xls"] == UNREAD
    assert stages["02/table.xls"] == "Geoprops"
    # A run stopped before it started is read as nothing, name and all.
    assert stages["06/output.log"] == UNREAD


def test_a_step_can_be_corrected_since_a_directory_holds_test_runs(blind_run):
    corrected = {"05/output.log": UNREAD, "01/output.log": "Quench"}
    groups, stages = sorted_upload(blind_run, stages=corrected)

    assert "Pickup" not in groups
    assert stages["01/output.log"] == "Quench"
    # And a step nobody spoke for is still worked out.
    assert stages["03/output.log"] == "Label"


def test_the_run_reads_end_to_end_from_the_blind_upload(blind_run):
    mapping = CoordinatesMapping(blind_run, "mTRAQ", 3)
    droplets = mapping.parse_droplets()
    metadata = mapping.map_data()
    metastats = mapping.map_stats()

    assert set(droplets["Type"]) == {"DMSO", "Geoprops", "Label", "Quench"}
    # One row a cell: two fields of six.
    assert len(metadata) == 12
    assert {
        "Drops.DMSO",
        "Drops.Label",
        "Drops.Quench",
        "Diameter",
        "Well",
        "Plex",
    } <= set(metadata.columns)
    assert set(metadata["Plex"].dropna()) <= {"0", "4", "8"}
    # The chamber readings are stamped with the step they were taken during.
    assert set(metastats["Type"]) >= {"DMSO", "Label", "Quench"}
    assert pd.api.types.is_datetime64_any_dtype(metastats["Timestamp"])


def test_a_pickup_log_the_name_says_nothing_about_still_has_its_plate():
    plates = CoordinatesMapping(
        [
            Upload("05/output.log", pickup_log()),
            Upload("02/table.xls", cells_table()),
            Upload("02/output.log", cells_log()),
        ],
        "mTRAQ",
        3,
    )

    assert set(plates.parsed_data["Pickup"]["Plate"]) == {1}


def test_a_well_shared_by_two_plates_still_says_which_placed_it():
    """Two pickups, two destination plates, the same well name on each.

    `Plate` cannot tell them apart -- it is the plate's position in the run,
    which is 1 for both -- so the file that did the dispensing has to.
    """
    wells = ["A1", "A2", "A3", "B1", "B2", "B3"]
    mapping = CoordinatesMapping(
        [
            Upload("02/table.xls", cells_table()),
            Upload("02/output.log", cells_log()),
            Upload("05/pickup_plate1.log", pickup_field(1, wells)),
            Upload("06/pickup_plate2.log", pickup_field(2, wells)),
        ],
        "mTRAQ",
        3,
    )
    metadata = mapping.map_data()

    field1 = metadata[metadata["Field"] == 1]
    field2 = metadata[metadata["Field"] == 2]

    # Two different pickups placed these cells, and the well names alone say
    # nothing about which -- both fields draw from the same well letters.
    assert set(field1["Well"]) & set(wells)
    assert set(field2["Well"]) & set(wells)
    assert (field1["Pickup.Source"] == "05/pickup_plate1.log").all()
    assert (field2["Pickup.Source"] == "06/pickup_plate2.log").all()


def test_a_step_run_twice_is_a_step_redone():
    """The instrument was stopped part way — a label droplet missed, a cell would
    not pick up — and the step was run again over some of the slide. The second
    attempt is what became of those cells; the rest keep the first."""
    half = SLIDE[:3]
    first = log(
        0, [dispensing(2, "D1", 1, 1, x, y) for x, y in SLIDE], run_id="dmso_first"
    )
    redone = log(
        100, [dispensing(102, "D4", 1, 1, x, y) for x, y in half], run_id="dmso_again"
    )

    twice = CoordinatesMapping(
        [
            Upload("01/first.log", first),
            Upload("02/again.log", redone),
            Upload("03/table.xls", cells_table()),
            Upload("03/output.log", cells_log()),
        ],
        "mTRAQ",
        3,
    )
    dmso = twice.parsed_data["DMSO"]

    # Each cell recorded once, not once an attempt.
    assert len(dmso) == len(SLIDE)

    wells = dmso.set_index(["XPos", "YPos"])["Well"]
    assert all(wells[(x + 1, y + 1)] == "D4" for x, y in half)
    assert all(wells[(x + 1, y + 1)] == "D1" for x, y in SLIDE[3:])

    reached = dmso.set_index(["XPos", "YPos"])["Timestamp"]
    assert all(reached[(x + 1, y + 1)] == reached.max() for x, y in half)
    assert all(reached[(x + 1, y + 1)] == reached.min() for x, y in SLIDE[3:])

    # What each attempt is left accounting for.
    assert twice.contributed["01/first.log"] == len(SLIDE) - len(half)
    assert twice.contributed["02/again.log"] == len(half)


def test_the_paperwork_in_a_log_is_not_mistaken_for_droplets():
    """Every log opens with a column header, a summary of the method and a
    parameter, all three of which say 'drops' and none of which is a droplet."""
    counted = CoordinatesMapping(
        [
            Upload("01/output.log", dmso_log()),
            Upload("03/output.log", label_log()),
            Upload("02/table.xls", cells_table()),
            Upload("02/output.log", cells_log()),
        ],
        "mTRAQ",
        3,
    )

    assert len(counted.parsed_data["DMSO"]) == DROPLETS
    assert len(counted.parsed_data["Label"]) == DROPLETS


def test_a_fluorescence_channel_is_kept_beside_the_transmission():
    """A sorting experiment images the same cell twice, and the second image is
    the reason it was sorted. Counting rows counts each cell once a channel;
    keeping only the transmission row throws the sorting away."""
    two_channels = [
        Upload("02/table.xls", cells_table(channels=("Transmission", "Green"))),
        Upload("02/output.log", cells_log()),
        Upload("01/output.log", dmso_log()),
    ]
    mapping = CoordinatesMapping(two_channels, "mTRAQ", 3)
    cells = mapping.parsed_data["Geoprops"]

    # One row a cell, not one a channel.
    assert len(cells) == 12
    # The transmission measurements under their own names, the fluorescence
    # beside them under the channel's.
    assert {"Diameter", "Intensity"} <= set(cells.columns)
    assert {"Diameter.Green", "Intensity.Green"} <= set(cells.columns)
    assert cells["Diameter"].eq(18.2).all()
    assert cells["Diameter.Green"].eq(19.2).all()
    # And the column that named the channel has done its job.
    assert "Teg" not in cells.columns


def test_transmission_alone_comes_out_as_it_went_in():
    """The common case, and it must not grow columns it has no data for."""
    mapping = CoordinatesMapping(
        [
            Upload("02/table.xls", cells_table()),
            Upload("02/output.log", cells_log()),
            Upload("01/output.log", dmso_log()),
        ],
        "mTRAQ",
        3,
    )
    cells = mapping.parsed_data["Geoprops"]

    assert len(cells) == 12
    assert not [column for column in cells.columns if "." in column]


def test_fluorescence_alone_is_still_read():
    """A preparation imaged in one fluorescence channel and no transmission: the
    first channel is the cell rather than there being no cells."""
    mapping = CoordinatesMapping(
        [
            Upload("02/table.xls", cells_table(channels=("Green",))),
            Upload("02/output.log", cells_log()),
            Upload("01/output.log", dmso_log()),
        ],
        "mTRAQ",
        3,
    )
    cells = mapping.parsed_data["Geoprops"]

    assert len(cells) == 12
    assert cells["Diameter"].eq(18.2).all()


@pytest.mark.skipif(
    not os.environ.get("PROTEOLYZER_CELLENONE_RUN"),
    reason="set PROTEOLYZER_CELLENONE_RUN to a real run directory",
)
def test_a_real_run_directory_reads_end_to_end():
    """Against instrument output, where test runs and second attempts turn up."""
    root = pathlib.Path(os.environ["PROTEOLYZER_CELLENONE_RUN"])
    files = [
        Upload(str(path.relative_to(root)), path.read_bytes().decode("latin1"))
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.suffix.lower() in (".log", ".xls")
    ]

    mapping = CoordinatesMapping(files, "mTRAQ", 3)
    placed = {
        stage: [file.name for file in group]
        for stage, group in mapping._output_file_paths().items()
    }

    # Every file is either placed in a step or said to be read by none.
    assert all(name in mapping.stages for name in (file.name for file in files))
    # No folder yields two cells tables, however many spreadsheets it holds.
    folders = {name.rsplit("/", 1)[0] for name in placed.get("Geoprops", [])}
    for folder in folders:
        tables = [
            name
            for name in placed.get("Geoprops", [])
            if name.startswith(folder) and name.endswith(".xls")
        ]
        assert len(tables) <= 1, tables

    metadata = mapping.map_data()
    assert len(metadata) > 0
    assert not metadata.duplicated(["Target", "Field", "XPos", "YPos"]).any()
    assert metadata["Diameter"].notna().any()
    assert len(mapping.map_stats()) > 0


def test_the_same_upload_can_be_read_twice(blind_run):
    """The case the app is built on: read the run, correct a step, read it again.

    A parser that starts where the last one stopped reads an empty file the
    second time, which arrives as 'No columns to parse from file'.
    """
    first = CoordinatesMapping(blind_run, "mTRAQ", 3)
    again = CoordinatesMapping(blind_run, "mTRAQ", 3, stages={"05/output.log": UNREAD})

    assert len(again.map_data()) == len(first.map_data())
    assert again.stages["05/output.log"] == UNREAD


def test_a_preparation_that_skipped_steps_still_maps_its_cells():
    """Not every run does every step, and the missing ones are not an error.

    The run this was written from has its cells, one label and one pickup, and no
    DMSO, quench or dispense at all -- which is what the reader said of it, and
    what it should say. Every end-to-end case here uses `blind_run`, which has all
    six, so the absent side of each step was never read.
    """
    sparse = [
        Upload("02/output.log", cells_log()),
        Upload("02/table.xls", cells_table()),
        Upload("03/output.log", label_log()),
        Upload("05/output.log", pickup_log()),
    ]
    mapping = CoordinatesMapping(sparse, "mTRAQ", 3)
    metadata = mapping.map_data()

    assert "DMSO" not in mapping.parsed_data
    assert "Quench" not in mapping.parsed_data

    # The cells are all there, once each, with where they went.
    assert len(metadata) > 0
    assert not metadata.duplicated(["Target", "Field", "XPos", "YPos"]).any()
    assert metadata["Diameter"].notna().any()

    # And the steps that did run are still read, rather than dropped with them.
    assert "Drops.Label" in metadata.columns
    assert len(mapping.map_stats()) > 0

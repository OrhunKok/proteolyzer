"""Read a cellenONE run directory: what was dispensed where, and onto which cell.

The instrument writes a directory per run — a log per step, an isolated-cells
table beside it, and a great many images. What DO-MS and any other consumer want
out of it is one row per cell: where it landed, which label it got, how big it
was, and what the chamber was doing while that happened.

Moved here from streamlit-DO-MS, which had the more developed of two copies of
this: the same class name and method names existed in a second repository, both
had independently fixed the same imaging-channel bug, and neither knew about the
other. This is the surviving copy; see DECISIONS.md for what was kept and why.
"""

import os
import re
from functools import reduce
from typing import Literal

import numpy as np
import pandas as pd

from ..core.logging import Logged


def _column_to_datetime(column: pd.Series) -> pd.Series:
    return pd.to_datetime(column, format="%d.%m.%y-%H:%M:%S.%f", errors="coerce")


def _column_to_numeric(column: pd.Series) -> pd.Series:
    return pd.to_numeric(column, errors="coerce")


CELLEONE_MAPPING = {
    "mTRAQ": {
        2: {
            "G1": "0",
            "G2": "8",
            "G3": "0",
            "G4": "8",
            "G5": "0",
            "G6": "8",
            "G7": "0",
            "G8": "8",
        },
        3: {
            "G1": "0",
            "G2": "4",
            "G3": "8",
            "G4": "0",
            "G5": "4",
            "G6": "8",
            "G7": "0",
            "G8": "4",
            "G9": "8",
            "G10": "0",
            "G11": "4",
            "G12": "8",
        },
    },
    "TMT": {
        14: {
            "G1": "TMTpro-128C",
            "G2": "TMTpro-129N",
            "G3": "TMTpro-129C",
            "G4": "TMTpro-130N",
            "G5": "TMTpro-130C",
            "G6": "TMTpro-131N",
            "G7": "TMTpro-131C",
            "G8": "TMTpro-132N",
            "G9": "TMTpro-132C",
            "G10": "TMTpro-133N",
            "G11": "TMTpro-133C",
            "G12": "TMTpro-134N",
            "G13": "TMTpro-134C",
            "G14": "TMTpro-135N",
        }
    },
}

# What a file has to say about itself for DO-MS to place it in the preparation.
# Read against the upload path, and against the run's own name inside the log.
# 'digest' and 'dispense' are here so that a step DO-MS has no plots for stays a
# plain dispense record rather than being taken for one it does.
STAGE_WORDS = {
    "Pickup": ("pickup",),
    "Quench": ("quench", "teab"),
    "Label": ("label", "mtraq", "tmt"),
    "DMSO": ("dmso",),
    "Geoprops": ("isolated", "reordered", "cells", "geoprops"),
    "Dispense": ("digest", "dispense"),
}

# The steps that put something on the slide, and so have dispensing lines to be
# recognised by.
DISPENSING_STAGES = ("Pickup", "Quench", "Label", "DMSO")

# The steps a file can be filed under, and the answer for one DO-MS took nothing
# from: neither a step of the preparation nor a record of the chamber. Shown as
# itself rather than filed under a step it is not.
STAGES = ("Geoprops", "DMSO", "Dispense", "Label", "Quench", "Pickup")
UNREAD = "Not read"

# How a dispensing line starts. Every log also carries a header row and a summary
# of the method, both of which mention drops without being a droplet, so this is
# what separates the droplets from the paperwork.
DISPENSED_AT = re.compile(r"\d{2}\.\d{2}\.\d{2}-\d{2}:\d{2}:\d{2}\.\d+")

# How much of a spreadsheet has to be placed cells for it to be the table of
# isolated cells. cellenONE writes two tables with the same columns: one row per
# cell it printed, and one row per droplet it ever saw, of which the printed ones
# are a fraction of a percent.
PLACED_ENOUGH = 0.9

# What identifies a cell on the slide, and so what makes two records records of
# the same cell.
CELL_KEYS = ["Target", "Field", "XPos", "YPos"]

# The second nozzle's wells, which the instrument does not record: the same plate
# one letter along from the nozzle it does.
NOZZLE_WELLS = {
    "A": "C",
    "B": "D",
    "E": "G",
    "F": "H",
    "I": "K",
    "J": "L",
    "M": "O",
    "N": "P",
}


class CoordinatesMapping(Logged):
    def __init__(
        self,
        root_dir: str,
        label_type: Literal["mTRAQ", "TMT"] | None = None,
        plex: int | None = None,
        stages: dict | None = None,
    ):

        self.root_dir = root_dir
        self.label_type = label_type
        self.plex = plex
        # {file name: step}, from the user. A directory off the instrument holds
        # test runs and second attempts that look like the real thing and are not,
        # and no reading of the files can know which; what it can do is say what it
        # made of each one and take a correction. Anything not named here is worked
        # out as usual.
        self.given = {
            name: stage
            for name, stage in (stages or {}).items()
            if stage in STAGES + (UNREAD,)
        }
        # How many cells each file ends up accounting for, filled in as the steps
        # are processed. A run that was abandoned after two cells says so here.
        self.contributed: dict[str, int] = {}
        self.file_paths = self._output_file_paths()
        self.parsed_data, self.parsed_stats = self._files_parse(self.file_paths)
        self.parsed_data = self._data_process(self.parsed_data)
        self.parsed_stats = self._stats_process(self.parsed_stats)

    @staticmethod
    def _rewound(file):
        """`file`, positioned at its start.

        A parser should not depend on where the last one left the file. An upload
        is read more than once — the directory is peeked at to work out which step
        each file is, and read again to parse it, and the whole run is read again
        when a step is corrected and resubmitted — and the second read of a
        consumed buffer is an empty file, which arrives as 'No columns to parse'.
        """
        seek = getattr(file, "seek", None)
        if callable(seek):
            seek(0)

        return file

    @staticmethod
    def _peek(file, limit: int = 400_000) -> str:
        """The start of a file as text, without moving where it is being read from.

        The parsers below read each file themselves, from the beginning, so this
        cannot leave the position anywhere else. An upload keeps its bytes, so
        there is nothing to put back; anything else is read and wound back.
        """
        try:
            data = file.getvalue()
        except AttributeError:
            position = file.tell()
            data = file.read()
            file.seek(position)

        return data[:limit].decode("latin1", errors="ignore").lower()

    @staticmethod
    def _dispensing_table(text: str):
        """The dispensing lines of a log, in columns, or None if it holds none.

        Same rules log_parse works by: a dispensing line says 'drops' and starts
        with the time it was dispensed at, and the columns are the ones that are
        not empty on every such line. Eleven of them is the layout the parser
        names, so anything else is left alone rather than guessed at.
        """
        rows = [
            line.split("\t")
            for line in text.splitlines()
            if ("drops" in line) and DISPENSED_AT.match(line)
        ]
        if not rows:
            return None

        width = max(len(cells) for cells in rows)
        rows = [cells + [""] * (width - len(cells)) for cells in rows]
        filled = [
            index
            for index in range(width)
            if any(cells[index].strip() for cells in rows)
        ]

        if len(filled) != 11:
            return None

        return [[cells[index].strip() for index in filled] for cells in rows]

    @staticmethod
    def _run_id(text: str) -> str:
        """The name cellenONE gave the run, from the head of its log.

        The line reads 'Run ID: dmso_20260316_143133_606', which is the folder
        name over again — so a renamed folder still says what it held. It is read
        rather than the whole log because the log also lists the method it ran, and
        a method chain mentions every reagent it touches: the quench log names the
        label it is quenching, which is how a text search lands on the wrong step.
        """
        found = re.search(r"^run id:\s*(.+)$", text, flags=re.MULTILINE)

        return found.group(1).strip() if found else ""

    def _is_cells_table(self, file) -> bool:
        """Whether a spreadsheet is the table of isolated cells.

        cellenONE writes two tables with the same fifty-two columns: the cells it
        printed, and every droplet it detected on the way, the printed ones being
        a fraction of a percent of it. Both carry 'Teg', which is the column the
        parser filters by, so they are told apart by what they are mostly made of —
        rows that say which target, field and position a cell went to.
        """
        lines = self._peek(file).splitlines()
        if not lines:
            return False

        header = [cell.strip() for cell in lines[0].split("\t")]
        if "teg" not in header:
            return False

        try:
            where = [
                header.index(column)
                for column in ("teg", "target", "field", "xpos", "ypos")
            ]
        except ValueError:
            return False

        # The half line a capped read ends on is the last one, and only that one
        # is dropped. A short row is otherwise a droplet the instrument never
        # placed, which is what the table of all droplets is almost entirely made
        # of and exactly what the ratio below is measuring: discarding every short
        # row left the ratio computed over the placed rows alone, so it came to 1
        # for both tables and a run's geoprops table was read as a second table of
        # cells in every folder it had one.
        rows = [line.split("\t") for line in lines[1:]]
        if rows and len(rows[-1]) <= max(where):
            rows = rows[:-1]

        def value(row: list[str], column: int) -> str:
            """One cell of a row, empty where the row stops short of it."""
            return row[column].strip() if column < len(row) else ""

        # One channel's rows, because the instrument writes a row per cell per
        # imaging channel and the ratio below is over cells rather than over
        # images. Transmission where there is any, since that is the channel every
        # preparation has; otherwise whichever channel it imaged in, so a
        # fluorescence-only export is still recognised as the cells.
        channels = [
            value(row, where[0]).lower() for row in rows if value(row, where[0])
        ]
        if not channels:
            return False

        channel = (
            "transmission"
            if "transmission" in channels
            else max(set(channels), key=channels.count)
        )

        imaged = [row for row in rows if value(row, where[0]).lower() == channel]
        placed = [
            row for row in imaged if all(value(row, column) for column in where[1:])
        ]

        return bool(imaged) and ((len(placed) / len(imaged)) >= PLACED_ENOUGH)

    def _stage_of(self, file) -> str | None:
        """Which step of the preparation a file came from, or None if it will not say.

        The name is read first, the whole path of it: cellenONE writes the task
        into the file and folder names, and a name that says 'pickup' is not a
        guess. Where the name says nothing, the file itself is asked.
        """
        path = file.name.lower().replace("\\", "/")

        if path.endswith(".xls"):
            return "Geoprops" if self._is_cells_table(file) else UNREAD

        if not path.endswith(".log"):
            return UNREAD

        text = self._peek(file)
        table = self._dispensing_table(text)
        readings = any(
            ("humidity" in line) and ("temperature" in line)
            for line in text.splitlines()
        )

        # A run that was stopped before it started dispensed nothing and recorded
        # nothing, whatever it was called. There is no step to put it under.
        if (table is None) and not readings:
            return UNREAD

        named = [
            stage
            for stage, words in STAGE_WORDS.items()
            if any(word in path for word in words)
        ]
        if len(named) == 1:
            return named[0]

        spoken = [
            stage
            for stage, words in STAGE_WORDS.items()
            if any(word in self._run_id(text) for word in words)
        ]
        if len(spoken) == 1:
            return spoken[0]

        # Nothing was dispensed, so this is the imaging or a plain dispense log,
        # read for its chamber readings alone. Which of the two it is only decides
        # the name on a line in the chamber plots, and if the file will not say,
        # the folder it came in might.
        if table is None:
            return None

        return self._stage_from_dispensing(table)

    @staticmethod
    def _stage_from_dispensing(table: list) -> str | None:
        """The step a dispensing log came from, read off what went where.

        Pickup takes one cell to one well, so its wells are many and spread over a
        plate. Labelling draws from the row of the reagent plate that the label
        mapping is written for. What is left is a reagent going onto the slide,
        DMSO or quench, which the two cannot be told apart by here.
        """
        wells = {row[4] for row in table if row[4]}
        if not wells:
            return None

        if (len(wells) > 4) and all(
            re.fullmatch(r"[a-p]\d{1,2}", well) for well in wells
        ):
            return "Pickup"

        if all(re.fullmatch(r"g\d{1,2}", well) for well in wells):
            return "Label"

        return None

    @staticmethod
    def _first_dispense(table):
        """When a log first dispensed anything, for putting logs in order."""

        if not table:
            return pd.NaT

        return pd.to_datetime(
            table[0][0], format="%d.%m.%y-%H:%M:%S.%f", errors="coerce"
        )

    def _output_file_paths(self):
        """Sort an uploaded run directory into the steps of the preparation.

        A directory does not have to be named for its steps. Every file is asked
        what it is (see _stage_of), and the ones that will not say are settled
        against the rest of the upload: a file sitting in a folder whose other
        files all came from one step came from that step too, which is what pairs
        the cells log with the cells table; and a reagent dispense that named
        neither itself nor its wells is DMSO if it ran before the labelling and
        nothing else claims to be DMSO, or quench if it ran after.
        """
        # Each log read once. Its dispensing lines are wanted twice over — to
        # recognise the step and then to order the steps — and a log is megabytes.
        tables = {
            file: self._dispensing_table(self._peek(file))
            for file in self.root_dir
            if file.name.lower().endswith(".log")
        }

        stages, undecided = {}, []
        for file in self.root_dir:
            stage = self.given.get(file.name) or self._stage_of(file)
            stages[file] = stage

            # Only logs are settled by the passes below. A spreadsheet that is not
            # the cells table is read by nothing here, and putting it in a group on
            # the strength of the folder it came in would only hand it to a parser
            # that cannot read it.
            if (stage is None) and file.name.lower().endswith(".log"):
                undecided.append(file)

        # Files sitting together came off the same step. Only where every decided
        # file in the folder agrees: a flat upload of the whole run says nothing.
        for file in list(undecided):
            folder = os.path.dirname(file.name.replace("\\", "/"))
            alongside = {
                stage
                for other, stage in stages.items()
                if (stage is not None)
                and (stage != UNREAD)
                and (os.path.dirname(other.name.replace("\\", "/")) == folder)
            }

            if len(alongside) == 1:
                stages[file] = alongside.pop()
                undecided.remove(file)

        # DMSO goes onto the slide before the cells and quench after the label, so
        # the clock separates them. Only ever used to fill a step nothing else
        # claimed, never to move a file off one it did, and only for a log that
        # dispensed something: a log of chamber readings alone is neither.
        claimed = set(stages.values())
        labelled = min(
            (
                self._first_dispense(tables.get(file))
                for file, stage in stages.items()
                if stage == "Label"
            ),
            default=pd.NaT,
        )

        for file in sorted(
            undecided, key=lambda file: self._first_dispense(tables.get(file))
        ):
            started = self._first_dispense(tables.get(file))

            if tables.get(file) is None:
                stages[file] = "Dispense"
            elif ("DMSO" not in claimed) and (
                pd.isna(labelled) or pd.isna(started) or (started < labelled)
            ):
                stages[file] = "DMSO"
            elif (
                ("Quench" not in claimed)
                and (not pd.isna(started))
                and (not pd.isna(labelled))
                and (started > labelled)
            ):
                stages[file] = "Quench"
            else:
                stages[file] = "Dispense"

            claimed.add(stages[file])

        cellenone_files = {
            stage: [file for file, found in stages.items() if found == stage]
            for stage in STAGES
        }

        for key, value in cellenone_files.items():
            if len(value) == 0:
                print(f"No {key} files found in the upload.")
            else:
                print(f"Found {len(value)} {key} files.")

        # What each file was taken to be, for the app to show: a guess the user
        # cannot see is a guess they cannot correct. Files nothing could place are
        # in here too, saying so, rather than disappearing quietly.
        self.stages = {file.name: stage or "Not read" for file, stage in stages.items()}

        return {key: value for key, value in cellenone_files.items() if len(value) > 0}

    def _files_parse(self, file_paths: dict):

        parsed_data = {}
        parsed_stats = {}

        for key, value in file_paths.items():
            if key == "Geoprops":
                xls_files = [f for f in value if f.name.endswith(".xls")]
                log_files = [f for f in value if f.name.endswith(".log")]
                df, stats = self.xls_parse(xls_files), self.log_parse(key, log_files)[1]
            elif key == "DMSO":
                df, stats = self.log_parse(key, value)
            elif key == "Dispense":
                df, stats = self.log_parse(key, value)
            elif key == "Label":
                df, stats = self.log_parse(key, value)
            elif key == "Quench":
                df, stats = self.log_parse(key, value)
            elif key == "Pickup":
                df, stats = self.log_parse(key, value)

            parsed_data[key] = df
            parsed_stats[key] = stats

        return parsed_data, parsed_stats

    def _data_process(self, parsed_data: dict):

        for key, df in parsed_data.items():
            if key == "DMSO":
                df = df.drop(["Plate", "Nozzle"], axis=1)
                # The two rows dropped here were the header and the method summary,
                # which log_parse no longer hands over: it keeps the lines that
                # start with the time they were dispensed at, and those two do not.
                df["Drops"] = df["Drops"].str.strip(" drops").astype(int)
                df = df.reset_index(drop=True)

            # A cells log with no cells table beside it leaves nothing to process
            # here; its chamber readings are still worth keeping.
            if (key == "Geoprops") and ("Teg" in df.columns):
                df = self._merge_imaging_channels(df)
                df = df.drop(["DropNo"], axis=1)
                df.columns = [re.sub(r"\s+", "", col) for col in df.columns]
                df["Timestamp"] = pd.to_datetime(
                    df["Date"] + " " + df["Time"], format="mixed"
                )
                df = df.drop(["Date", "Time"], axis=1)

            if (key == "Label") and (self.label_type is not None):
                df["Drops"] = df["Drops"].str.strip(" drops").astype(int)
                df = self.label_well_plex(df)

            if (key == "Quench") and (self.label_type is not None):
                df = df.drop(["Plate", "PlatePos", "Nozzle", "Well"], axis=1)
                df["Drops"] = df["Drops"].str.strip(" drops").astype(int)
                df["Timestamp"] = pd.to_datetime(
                    df["Timestamp"], format="%d.%m.%y-%H:%M:%S.%f", errors="coerce"
                ).dt.floor("s")
                df = df.reset_index(drop=True)

            if key not in ["Geoprops", "Dispense"]:
                df["Timestamp"] = pd.to_datetime(
                    df["Timestamp"], format="%d.%m.%y-%H:%M:%S.%f", errors="coerce"
                ).dt.floor("s")
                df["Target"] = df["Target"].astype(int)
                df["Field"] = df["Field"].astype(int)
                df["XPos"] = df["XPos"].astype(int) + 1
                df["YPos"] = df["YPos"].astype(int) + 1

            # A step that ran twice ran again because the first attempt was
            # abandoned: a label droplet that missed, a cell that would not pick
            # up. Where both attempts recorded the same cell, the later one is what
            # became of it and the earlier record is history — and where the redo
            # only covered part of the slide, the cells it left alone keep theirs.
            if (key != "Dispense") and (
                set(CELL_KEYS + ["Timestamp"]) <= set(df.columns)
            ):
                df = df.sort_values("Timestamp").drop_duplicates(CELL_KEYS, keep="last")
                df = df.reset_index(drop=True)

            # What each file is left accounting for, before the column goes: a run
            # abandoned after two cells, or superseded entirely by a later one,
            # comes out at two cells or none.
            if "Source" in df.columns:
                self.contributed.update(
                    {
                        str(name): int(count)
                        for name, count in df["Source"].value_counts().items()
                    }
                )

                # A pickup is a destination plate, and the plate a run writes is
                # its position rather than its identity -- one destination plate
                # is mounted, so it reads 1 for every real pickup. The file that
                # did the dispensing is what tells two plates apart, so it has to
                # survive past here, onto the cell it placed.
                if key == "Pickup":
                    df["Pickup.Source"] = df["Source"]

                df = df.drop("Source", axis=1)

            if key == "Pickup":
                df["Drops"] = df["Drops"].str.strip(" drops").astype(int)

                # cellenONE only records the information for one nozzle; the
                # other put the same drops in the same wells one letter along,
                # and has to be written down here.
                other_nozzle = df.copy()
                other_nozzle["Nozzle"] = 3
                other_nozzle["XPos"] = other_nozzle["XPos"] + 36
                other_nozzle["Well"] = other_nozzle["Well"].apply(
                    lambda well: NOZZLE_WELLS.get(well[0]) + well[1:]
                )

                df = pd.concat([df, other_nozzle], axis=0)
                df = df.reset_index(drop=True)

            df["Type"] = key
            parsed_data[key] = df

        return parsed_data

    #: What a geoprops row measures, as against what it identifies. These are the
    #: columns that differ between imaging channels of the same cell, so these
    #: are what a fluorescence channel contributes.
    MEASURED = (
        "Diameter",
        "Circularity",
        "Elongation",
        "Intensity",
        "Background",
        "X",
        "Y",
        "ImageFile",
    )

    def _merge_imaging_channels(self, frame: pd.DataFrame) -> pd.DataFrame:
        """One row per cell, with each imaging channel's measurements beside it.

        cellenONE writes a geoprops row per cell *per imaging channel*, named in
        'Teg'. Counted as rows, a slide imaged in transmission and one
        fluorescence channel holds twice the cells it holds, and every geometry
        average is halved.

        Keeping only the transmission row fixes the count and throws the
        fluorescence away — which in a sorting experiment is the half of the
        record the experiment was for. So the transmission row is the cell, and
        every other channel joins onto it under its own name: Diameter.Green
        beside Diameter. A preparation imaged in transmission alone comes out
        exactly as it went in, minus the duplicate rows it never had.
        """
        channels = [
            channel
            for channel in frame["Teg"].dropna().unique()
            if str(channel).strip()
        ]
        keys = [key for key in CELL_KEYS if key in frame.columns]

        # Nothing to merge on, or nothing to merge: leave it as it is rather than
        # guess, and let the column go since it no longer says anything.
        if (not keys) or (len(channels) < 2):
            return frame.drop("Teg", axis=1).reset_index(drop=True)

        # Transmission is the geometry of record. A run that has none of it —
        # fluorescence-only imaging — takes its first channel as the base rather
        # than coming out empty.
        base_channel = "Transmission" if "Transmission" in channels else channels[0]
        base = frame.loc[frame["Teg"] == base_channel].drop("Teg", axis=1)

        measured = [column for column in self.MEASURED if column in frame.columns]
        merged = base
        for channel in channels:
            if channel == base_channel:
                continue

            beside = frame.loc[frame["Teg"] == channel, keys + measured]
            beside = beside.rename(
                columns={column: f"{column}.{channel}" for column in measured}
            )
            merged = merged.merge(beside, on=keys, how="left")

        self.logger.info(
            f"Imaging channels merged onto {base_channel}: "
            f"{[channel for channel in channels if channel != base_channel]}."
        )

        return merged.reset_index(drop=True)

    def _metadata_validate(self, metadata: pd.DataFrame):

        if all(col in metadata.columns for col in ["Plate.Pickup", "Well.Pickup"]) and (
            self.plex is None
        ):
            well_clash = metadata.groupby(["Plate.Pickup", "Well.Pickup"]).size() != 1
            metadata["Well.Clash"] = (
                metadata.set_index(["Plate.Pickup", "Well.Pickup"])
                .index.map(well_clash.to_dict())
                .fillna(False)
            )

            if well_clash.sum() > 0:
                self.logger.warning(
                    "Well clashes detected: more than one cell mapped to a well. "
                    "The Well.Clash column says which."
                )

        if all(
            col in metadata.columns
            for col in ["Plate.Pickup", "Well.Pickup", "Plex.Label"]
        ):
            # metadata.groupby(['Plate.Pickup', 'Well.Pickup'])['Plex.Label'].apply(lambda x: x.duplicated(keep = False))
            labels_clash = (
                metadata.groupby(["Plate.Pickup", "Well.Pickup"])[
                    "Plex.Label"
                ].nunique()
                != self.plex
            )
            metadata["Label.Clash"] = (
                metadata.set_index(["Plate.Pickup", "Well.Pickup"])
                .index.map(labels_clash.to_dict())
                .fillna(False)
            )

            if labels_clash.sum() > 0:
                self.logger.warning(
                    "Label clashes detected: a well carries more than one label. "
                    "The Label.Clash column says which."
                )

        return metadata

    def parse_droplets(self) -> pd.DataFrame:

        droplets = pd.DataFrame(columns=["XPos", "YPos", "Target", "Field"])

        if "DMSO" in self.parsed_data:
            dmso_df = self.parsed_data["DMSO"]
            droplets = pd.merge(
                dmso_df, droplets, on=["XPos", "YPos", "Target", "Field"], how="left"
            )

        if "Label" in self.parsed_data:
            label_df = self.parsed_data["Label"].copy()
            droplets = pd.concat([label_df, droplets], axis=0)

        if "Quench" in self.parsed_data:
            quench_df = self.parsed_data["Quench"].copy()
            droplets = pd.concat([quench_df, droplets], axis=0)

        if "Geoprops" in self.parsed_data:
            geo_df = self.parsed_data["Geoprops"]

            if "Pickup" in self.parsed_data:
                pickup_df = self.parsed_data["Pickup"].copy()

                pickup_df["MappedGeo"] = self._map_coords(geo_df, pickup_df)
                exploded = pickup_df.explode("MappedGeo")
                # strict: two columns of one frame, so the lengths are equal by
                # construction, and if they ever are not that is a bug worth
                # hearing about rather than a mapping quietly short of entries.
                well_mapping = dict(
                    zip(exploded["MappedGeo"], exploded["Well"], strict=True)
                )
                plate_mapping = dict(
                    zip(exploded["MappedGeo"], exploded["Plate"], strict=True)
                )
                # Which pickup placed the cell, so that a well shared by two
                # destination plates does not collapse into one key: see
                # DECISIONS.md and the `Pickup.Source` column it produces.
                source_mapping = dict(
                    zip(
                        exploded["MappedGeo"],
                        exploded["Pickup.Source"],
                        strict=True,
                    )
                )
                geo_df["Well"] = geo_df.index.map(well_mapping)
                geo_df["Plate"] = geo_df.index.map(plate_mapping)
                geo_df["Pickup.Source"] = geo_df.index.map(source_mapping)

            droplets = pd.concat([geo_df, droplets], axis=0)

        if "Dispense" in self.parsed_data:
            dispense_df = self.parsed_data["Dispense"].copy()
            droplets = pd.concat([dispense_df, droplets], axis=0)

        self.droplets = droplets.copy()

        return self.droplets

    def map_data(self) -> pd.DataFrame:

        # Parsed here if it has not been already: the droplets are what this maps,
        # and a caller that asks for the cells without asking for the droplets
        # first used to get an AttributeError about an attribute it has no reason
        # to know exists.
        if not hasattr(self, "droplets"):
            self.parse_droplets()

        droplets = self.droplets.copy()

        dfs = {
            "DMSO": droplets[droplets["Type"] == "DMSO"],
            "Geoprops": droplets[droplets["Type"] == "Geoprops"],
            "Label": droplets[droplets["Type"] == "Label"],
            "Quench": droplets[droplets["Type"] == "Quench"],
        }

        key_cols = ["XPos", "YPos", "Target", "Field"]

        processed = []
        for name, df in dfs.items():
            if df.empty:
                continue

            if name == "DMSO":
                df = df.drop(
                    ["Type", "Timestamp", "PlatePos", "Well", "Level"], axis=1
                ).dropna(how="all", axis=1)
                df = df.rename(columns={"Drops": "Drops.DMSO"})
            elif name == "Geoprops":
                df = df.drop(["Type", "Timestamp"], axis=1).dropna(how="all", axis=1)
            elif name == "Label":
                df = df.drop(
                    ["Type", "Timestamp", "Level", "Nozzle", "Well", "Plate"], axis=1
                ).dropna(how="all", axis=1)
                df = df.rename(columns={"Drops": "Drops.Label"})
            elif name == "Quench":
                df = df.drop(["Type", "Timestamp", "Level"], axis=1).dropna(
                    how="all", axis=1
                )
                df = df.rename(columns={"Drops": "Drops.Quench"})

            processed.append(df)

        metadata = (
            reduce(
                lambda left, right: left.merge(right, on=key_cols, how="outer"),
                processed,
            )
            if processed
            else pd.DataFrame()
        )

        if len(metadata) > 0:
            self._metadata_validate(metadata)
            return metadata
        else:
            self.logger.warning(
                "Missing the files metadata mapping needs; nothing was mapped. "
                f"Read: {sorted(self.parsed_data)}."
            )
            return metadata

    def _stats_process(self, parsed_stats: dict):

        for key, df in parsed_stats.items():
            df["Timestamp"] = df["Timestamp"].dt.floor("s")

            parsed_stats[key] = df

        return parsed_stats

    def map_stats(
        self,
    ):

        metastats = []
        for key, df in self.parsed_stats.items():
            df["Type"] = key

            metastats.append(df)

        metastats = pd.concat(metastats)

        metastats = metastats.sort_values("Timestamp").reset_index(drop=True)

        return metastats

    def xls_parse(self, file_paths: list):

        dispense_dfs = []
        for file in file_paths:
            df = pd.read_csv(self._rewound(file), sep="\t")
            df = df.loc[:, ~df.columns.str.contains("^Unnamed")]

            # The folder the table sat in names the sample. A flat upload has no
            # folder to read, so the file names itself.
            folder = os.path.basename(os.path.dirname(file.name.replace("\\", "/")))
            sample_name = (
                re.sub(r"\.Run$", "", folder)
                if folder
                else os.path.splitext(os.path.basename(file.name))[0]
            )
            df["SampleName"] = np.repeat(sample_name, len(df))
            df["Source"] = file.name
            df = df.drop(["Plate", "Well", "ImageFile", "Background"], axis=1)

            dispense_dfs.append(df)

        # Nothing to concatenate when a cells log turned up without its table.
        if not dispense_dfs:
            return pd.DataFrame()

        return pd.concat(dispense_dfs, axis=0).reset_index(drop=True)

    def log_parse(self, key: str, file_paths: list):

        dfs_list: list[pd.DataFrame] = []
        stats_dfs_list: list[pd.DataFrame] = []
        for file in file_paths:
            text = self._rewound(file).read().decode("latin1")

            # Which lines are wanted is decided on the lines themselves, before
            # anything is split into columns. Splitting a whole log into a frame of
            # every line by every column and then testing every cell of it was the
            # slowest thing in an upload by an order of magnitude: a log is tens of
            # thousands of lines, of which a few hundred are wanted.
            readings = [
                line.split("\t")
                for line in text.splitlines()
                if ("Humidity" in line) and ("Temperature" in line)
            ]

            temp_stats = (
                pd.DataFrame(readings).replace("", np.nan).dropna(how="all", axis=1)
            )
            time_col = temp_stats.apply(_column_to_datetime).dropna(how="all", axis=1)
            val_cols = temp_stats.apply(_column_to_numeric).dropna(how="all", axis=1)
            temp_stats = pd.concat([time_col, val_cols], axis=1).reset_index(drop=True)

            # Six columns is a chamber reading. A log of a run that was stopped
            # before the chamber settled has none, and naming six columns on a
            # frame that has none of them used to end the whole upload.
            if temp_stats.shape[1] == 6:
                temp_stats.columns = [
                    "Timestamp",
                    "Humidity",
                    "Temperature",
                    "Dew Point",
                    "Adj. Temp",
                    "Bath Temp",
                ]
                temp_stats["Timestamp"] = temp_stats["Timestamp"].dt.floor("s")
                stats_dfs_list.append(temp_stats)
            else:
                print(f"No chamber readings in {file.name}.")

            # The dispensing lines: they say 'drops' and they start with the time
            # they went down. Saying 'drops' is not enough on its own — the log
            # opens with a header row naming a 'No. of Drops' column and goes on to
            # list the method it ran, which has steps with names like
            # 'PreSpot250drops'. Neither is a droplet. The first two rows used to
            # be dropped by hand further down to get past them, which held only as
            # long as there were exactly two.
            #
            # The same reading of a log that decided which step it came from, so
            # the two cannot come to different answers about what a droplet is.
            table = self._dispensing_table(text)
            df = pd.DataFrame(table if table else [])

            if len(df) == 0:
                # Can exit from geoprops or dispense log here
                if key in ["Geoprops", "Dispense"]:
                    dfs_list.append(pd.DataFrame())
                else:
                    # In case log from aborted process is included in the output
                    print(
                        f"Skipping: {file.name}. Likely an aborted process during sample-prep."
                    )

            else:
                df.columns = [
                    "Timestamp",
                    "Plate",
                    "PlatePos",
                    "Nozzle",
                    "Well",
                    "Target",
                    "Level",
                    "Field",
                    "Drops",
                    "XPos",
                    "YPos",
                ]
                df = df.reset_index(drop=True)
                # Which run each droplet came from, for telling a redo from the
                # attempt it replaced. Dropped again once the steps are processed.
                df["Source"] = file.name

                if key == "Label":
                    df = df.drop(["PlatePos"], axis=1)

                elif key == "Pickup":
                    # The plate number where the file name carries it, which is
                    # where cellenONE puts it, and the log's own Plate column
                    # otherwise: the name is no longer required to say anything.
                    numbered = re.findall(
                        r"pickup[_\-. ]*(\d+)", os.path.basename(file.name).lower()
                    )
                    if numbered:
                        df["Plate"] = int(numbered[0])
                    else:
                        df["Plate"] = pd.to_numeric(df["Plate"], errors="coerce")

                dfs_list.append(df)

        dfs = pd.concat(dfs_list, axis=0) if dfs_list else pd.DataFrame()

        # An empty frame still has to carry the columns and the timestamp dtype
        # the chamber plots read, in case every log of a step was an aborted run.
        if stats_dfs_list:
            stats_dfs = pd.concat(stats_dfs_list, axis=0)
            stats_dfs = stats_dfs.groupby("Timestamp").mean().reset_index()
        else:
            stats_dfs = pd.DataFrame(
                {
                    "Timestamp": pd.Series(dtype="datetime64[ns]"),
                    **{
                        name: pd.Series(dtype=float)
                        for name in (
                            "Humidity",
                            "Temperature",
                            "Dew Point",
                            "Adj. Temp",
                            "Bath Temp",
                        )
                    },
                }
            )

        return dfs, stats_dfs

    def label_well_plex(self, label_df: pd.DataFrame):

        available_plex = (
            CELLEONE_MAPPING.get(self.label_type, {})
            if self.label_type is not None
            else {}
        )

        # The plex widget accepts any integer, so an unsupported combination has
        # to degrade to 'no label mapping' rather than raise a KeyError. The
        # plots that need it already guard on the Plex column being present.
        if self.plex not in available_plex:
            self.logger.warning(
                f"No {self.label_type} well layout for plex {self.plex}. "
                f"Supported: {sorted(available_plex)}. Skipping label mapping."
            )
            return label_df

        label_mapping = available_plex[self.plex]
        label_df["Plex"] = label_df["Well"].map(label_mapping)

        return label_df

    def _map_coords(self, geo_df, map_df, coord_cols=None, group_cols=None):

        # Not as defaults: a mutable default is shared between calls, and a tuple
        # handed to groupby is read as one key rather than as two.
        coord_cols = list(coord_cols or ("XPos", "YPos"))
        group_cols = list(group_cols or ("Target", "Field"))

        grouped_geo = geo_df.groupby(group_cols)
        grouped_map = map_df.groupby(group_cols)

        results = {}
        for group_key, geo_group in grouped_geo:
            # A field with no pickup against it was not picked up, which is a thing
            # that happens — fields left out of the run on purpose, most often —
            # and the pickup overlay is where it is shown. It was being written
            # onto the page as a message per field.
            if group_key not in grouped_map.groups:
                continue

            map_group = grouped_map.get_group(group_key)

            geo_coords = np.stack(geo_group[coord_cols].values)
            map_coords = np.stack(map_group[coord_cols].values)

            distances = np.linalg.norm(
                map_coords[:, np.newaxis, :] - geo_coords[np.newaxis, :, :], axis=2
            )
            sorted_indices = np.argsort(distances, axis=1)[:, : self.plex]
            closest_points = geo_group.index.values[sorted_indices]

            # Key on the source row so the result lines up with map_df when it is
            # assigned back as a column. Groups are visited in group_cols order,
            # not file order, so a positional result would silently mis-map wells
            # to cells (and shift everything after a skipped group).
            for i, map_idx in enumerate(map_group.index):
                results[map_idx] = closest_points[i]

        return pd.Series(results, dtype="object").reindex(map_df.index)

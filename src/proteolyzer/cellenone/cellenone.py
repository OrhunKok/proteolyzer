import os
import re
from typing import Literal

import numpy as np
import pandas as pd

from proteolyzer.utils.logging import Logged

from .config import (
    CELLEONE_MAPPING,
    CHANNEL_COL,
    DROPLET_COLS,
    MERGE_COLS,
    NOZZLE_WELL_MAPPING,
    PICKUP_NOZZLE_ID,
    PICKUP_NOZZLE_XPOS_OFFSET,
    PRIMARY_CHANNEL,
    TEMP_STATS_COLS,
)


class CoordinatesMapping(Logged):
    def __init__(
        self,
        root_dir: str,
        label_type: Literal["mTRAQ", "TMT"] | None = None,
        plex: int | None = None,
    ):
        self.root_dir = root_dir
        self.label_type = label_type
        self.plex = plex
        self.file_paths = self._output_file_paths()
        self.parsed_data, self.parsed_stats = self._files_parse(self.file_paths)
        self.parsed_data = self._data_process(self.parsed_data)
        self.parsed_stats = self._stats_process(self.parsed_stats)

    @staticmethod
    def _classify_file(dirpath: str, filename: str) -> str | None:
        """Which kind of output `filename` is, or None to ignore it.

        The run directory name counts as well as the file name. Operators name
        the logs inconsistently within one prep ("labels_1", "labeling_2",
        "L9", "l4"), while the ``.Run`` directory cellenONE creates keeps the
        step name, so classifying on the file name alone silently files
        labelling logs as dispense logs.
        """
        name = filename.lower()
        context = f"{os.path.basename(dirpath).lower()} {name}"

        if name.endswith(".xls"):
            if "isolated" in name and "reordered" in name:
                return "Geoprops"
            return None
        if not name.endswith(".log"):
            return None

        # Pickup first: a pickup log is never a labelling log.
        if "pickup" in context:
            return "Pickup"
        if "label" in context:
            return "Label"
        return "Dispense"

    def _output_file_paths(self):
        cellenone_files: dict[str, list[str]] = {
            "Geoprops": [],
            "Dispense": [],
            "Label": [],
            "Pickup": [],
        }

        for dirpath, _, filenames in os.walk(self.root_dir):
            for filename in sorted(filenames):
                kind = self._classify_file(dirpath, filename)
                if kind is None:
                    continue
                path = os.path.join(dirpath, filename)
                cellenone_files[kind].append(path)
                # Logged so a misfiled run can be spotted without guessing.
                self.logger.debug(f"{kind}: {path}")

        for key, value in cellenone_files.items():
            if len(value) == 0:
                self.logger.warning(
                    f"No {key} files found in the directory: {self.root_dir}"
                )
            else:
                self.logger.info(f"Found {len(value)} {key} files.")

        return {key: value for key, value in cellenone_files.items() if value}

    def _files_parse(self, file_paths: dict):
        parsed_data = {}
        parsed_stats = {}

        for key, value in file_paths.items():
            if key == "Geoprops":
                parsed_data[key] = self.xls_parse(value)
            else:
                parsed_data[key], parsed_stats[key] = self.log_parse(key, value)

        return parsed_data, parsed_stats

    def _data_process(self, parsed_data: dict):
        for key, df in parsed_data.items():
            if key == "Geoprops":
                df.columns = [re.sub(r"\s+", "", col) for col in df.columns]
                df["Timestamp"] = pd.to_datetime(
                    df["Date"] + " " + df["Time"], format="mixed"
                )

            if (key == "Label") and (self.label_type is not None):
                df = self.label_well_plex(df)

            if key not in ["Geoprops", "Dispense"]:
                df["Timestamp"] = pd.to_datetime(
                    df["Timestamp"], format="%d.%m.%y-%H:%M:%S.%f", errors="coerce"
                ).dt.floor("s")
                df["Target"] = df["Target"].astype(int)
                df["Field"] = df["Field"].astype(int)
                df["XPos"] = df["XPos"].astype(int) + 1
                df["YPos"] = df["YPos"].astype(int) + 1

            if key == "Pickup":
                other_nozzle = df.copy()
                other_nozzle["Nozzle"] = PICKUP_NOZZLE_ID
                other_nozzle["XPos"] = other_nozzle["XPos"] + PICKUP_NOZZLE_XPOS_OFFSET
                other_nozzle["Well"] = other_nozzle["Well"].apply(
                    lambda x: NOZZLE_WELL_MAPPING.get(x[0]) + x[1:]
                )

                df = pd.concat([df, other_nozzle], axis=0)
                df = df.reset_index(drop=True)

            parsed_data[key] = df

        return parsed_data

    def _metadata_validate(self, metadata: pd.DataFrame):
        if all(col in metadata.columns for col in ["Plate.Pickup", "Well.Pickup"]) and (
            self.plex is None
        ):
            well_clash = metadata.groupby(["Plate.Pickup", "Well.Pickup"]).size() != 1
            metadata["Well.Clash"] = (
                metadata.set_index(["Plate.Pickup", "Well.Pickup"])
                .index.map(well_clash)
                .fillna(False)
            )

            if well_clash.sum() > 0:
                self.logger.warning(
                    "Well clashes detected! Check Well.Clash column to view clashes."
                )

        if all(
            col in metadata.columns
            for col in ["Plate.Pickup", "Well.Pickup", "Plex.Label"]
        ):
            labels_clash = (
                metadata.groupby(["Plate.Pickup", "Well.Pickup"])[
                    "Plex.Label"
                ].nunique()
                != self.plex
            )
            metadata["Label.Clash"] = (
                metadata.set_index(["Plate.Pickup", "Well.Pickup"])
                .index.map(labels_clash)
                .fillna(False)
            )

            if labels_clash.sum() > 0:
                self.logger.warning(
                    "Label clashes detected! Check Label.Clash column to view clashes."
                )

        return metadata

    def map_data(self) -> pd.DataFrame:
        metadata = None

        if "Pickup" in self.parsed_data and "Geoprops" in self.parsed_data:
            pickup_df = self.parsed_data["Pickup"]
            geo_df = self.parsed_data["Geoprops"]

            mapped, distances = self._map_coords(geo_df, pickup_df)
            pickup_df["MappedGeo"] = mapped
            pickup_df["MappedGeoDistance"] = distances

            # A cell can be claimed by more than one pickup row (always so when
            # `plex` is unset, since then every candidate is kept). Resolve by
            # distance so the assignment does not depend on iteration order.
            claims = (
                pickup_df.explode(["MappedGeo", "MappedGeoDistance"])
                .dropna(subset=["MappedGeo"])
                .sort_values("MappedGeoDistance", kind="stable")
                .drop_duplicates(subset="MappedGeo", keep="first")
            )

            geo_df["Well.Pickup"] = geo_df.index.map(
                dict(zip(claims["MappedGeo"], claims["Well"], strict=True))
            )
            geo_df["Plate.Pickup"] = geo_df.index.map(
                dict(zip(claims["MappedGeo"], claims["Plate"], strict=True))
            )

            self.parsed_data["Geoprops"] = geo_df
            self.parsed_data["Pickup"] = pickup_df
            metadata = geo_df.copy()

        if "Geoprops" in self.parsed_data and "Label" in self.parsed_data:
            geoprop_df = self.parsed_data["Geoprops"].copy()
            label_df = self._collapse_label_dispenses(self.parsed_data["Label"].copy())

            geoprop_df = geoprop_df.rename(
                columns={
                    col: f"{col}.Geoprops"
                    for col in geoprop_df.columns
                    if col not in MERGE_COLS + ["Well.Pickup", "Plate.Pickup"]
                }
            )
            label_df = label_df.rename(
                columns={
                    col: f"{col}.Label"
                    for col in label_df.columns
                    if col not in MERGE_COLS
                }
            )

            metadata = pd.merge(
                geoprop_df,
                label_df,
                on=MERGE_COLS,
                how="inner",
            )

        if metadata is not None:
            self._metadata_validate(metadata)
            return metadata
        else:
            self.logger.warning(
                "Missing required files to do metadata mapping. Ensure all required inputs are present."
            )

    def _collapse_label_dispenses(self, label_df: pd.DataFrame) -> pd.DataFrame:
        """One row per labelled position, carrying the droplets delivered to it.

        A repeated labelling run dispenses the same channel onto a position
        again, so the log holds several rows for it. Keeping them all
        multiplies the merged metadata, so the latest dispense represents the
        position and the droplet counts are summed. The raw ``Drops`` text
        ("58 drops") is replaced by two numbers:

        Droplets
            Total droplets delivered to the position.
        Dispenses
            How many dispense events delivered them.
        """
        key = [*MERGE_COLS, "Well"]
        if not all(col in label_df.columns for col in key):
            return label_df

        counts = (
            label_df["Drops"].str.extract(r"(\d+)", expand=False).astype("Int64")
            if "Drops" in label_df.columns
            else pd.Series(pd.NA, index=label_df.index, dtype="Int64")
        )
        frame = label_df.assign(Droplets=counts).drop(columns="Drops", errors="ignore")

        if "Timestamp" in frame.columns:
            frame = frame.sort_values("Timestamp", kind="stable")

        grouped = frame.groupby(key, dropna=False)["Droplets"]
        totals = pd.DataFrame(
            {
                # min_count keeps an unknown total unknown: summing all-missing
                # counts would otherwise report zero droplets delivered.
                "Droplets": grouped.sum(min_count=1),
                "Dispenses": grouped.size(),
            }
        )
        collapsed = (
            frame.drop_duplicates(subset=key, keep="last")
            .drop(columns="Droplets")
            .merge(totals.reset_index(), on=key, how="left")
        )

        repeated = int((totals["Dispenses"] > 1).sum())
        if repeated:
            self.logger.info(
                f"Collapsed {len(frame)} label dispenses onto {len(collapsed)} "
                f"positions ({repeated} were dispensed to more than once)."
            )
        return collapsed

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

        return metastats

    def xls_parse(self, file_paths: list):
        dispense_dfs = []
        for file in file_paths:
            df = pd.read_csv(file, sep="\t")
            df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
            # Normalized here rather than downstream: the export pads some
            # headers ("Date        "), and the steps below match on names.
            df.columns = [re.sub(r"\s+", "", col) for col in df.columns]

            # Replace with something smarter later
            sample_name = os.path.basename(os.path.dirname(file))
            sample_name = re.sub(".Run", "", sample_name)
            df["SampleName"] = sample_name
            df = df.drop(
                ["Plate", "Well", "ImageFile", "Background"], axis=1, errors="ignore"
            )

            dispense_dfs.append(self._merge_imaging_channels(df, sample_name))

        return pd.concat(dispense_dfs, axis=0).reset_index(drop=True)

    def _merge_imaging_channels(
        self, df: pd.DataFrame, sample_name: str
    ) -> pd.DataFrame:
        """Collapse the per-imaging-channel rows of a cell onto one row.

        cellenONE writes one geoprops row per (cell, channel): the
        Transmission row carries the geometry, and each fluorescence channel
        adds a row whose measurements are zero where nothing was detected.
        Counting those as separate cells doubles the cell count and halves
        every geometry average, so the extra channels become extra columns
        (``Diameter.Green``) instead of extra rows.
        """
        if CHANNEL_COL not in df.columns:
            return df

        channels = list(pd.unique(df[CHANNEL_COL].dropna()))
        if len(channels) < 2:
            return df

        key = self._cell_key(df)
        if key is None:
            self.logger.error(
                f"{sample_name}: cannot tell the imaging channels of a cell apart, "
                f"so channels {channels} are left as separate rows. Cell counts and "
                "geometry averages will include every channel."
            )
            return df

        primary = PRIMARY_CHANNEL if PRIMARY_CHANNEL in channels else min(channels)
        per_channel = self._per_channel_columns(df, key)

        merged = df[df[CHANNEL_COL] == primary].copy()
        for channel in sorted(c for c in channels if c != primary):
            extra = df.loc[df[CHANNEL_COL] == channel, [*key, *per_channel]].rename(
                columns={col: f"{col}.{channel}" for col in per_channel}
            )
            merged = merged.merge(extra, on=key, how="left")

        self.logger.info(
            f"{sample_name}: merged imaging channels {channels} onto {len(merged)} "
            f"cells (primary channel '{primary}')."
        )
        return merged

    @staticmethod
    def _cell_key(df: pd.DataFrame) -> list[str] | None:
        """Columns identifying one cell, such that (cell, channel) is unique."""
        for key in (
            ["SampleName", "DropNo"],
            ["SampleName", "Target", "Field", "XPos", "YPos"],
        ):
            if (
                all(col in df.columns for col in key)
                and not df.duplicated(subset=[*key, CHANNEL_COL]).any()
            ):
                return key
        return None

    @staticmethod
    def _per_channel_columns(df: pd.DataFrame, key: list[str]) -> list[str]:
        """Columns measured per channel, i.e. those differing within a cell.

        Derived from the data rather than hardcoded, so an export with more or
        fewer measurement columns still splits correctly.
        """
        candidates = [col for col in df.columns if col not in {*key, CHANNEL_COL}]
        varies = df.groupby(key, dropna=False)[candidates].nunique(dropna=False) > 1
        return [col for col in candidates if varies[col].any()]

    def log_parse(self, key: str, file_paths: list):
        dfs = []
        stats_dfs = []
        for file_path in file_paths:
            with open(file_path, encoding="latin1") as handle:
                df = handle.readlines()
                df = [re.sub("\n", "", ele) for ele in df]
                df = pd.DataFrame(df)
                df = df[0].str.split("\t", expand=True)
                df = df.replace("", np.nan)
                df = df.dropna(how="all", axis=1)

                temp_stats = df[
                    (df.apply(lambda row: row.str.contains("Humidity")).any(axis=1))
                    & (
                        df.apply(lambda row: row.str.contains("Temperature")).any(
                            axis=1
                        )
                    )
                ]
                temp_stats = temp_stats.dropna(how="all", axis=1)
                time_col = temp_stats.apply(
                    lambda x: pd.to_datetime(
                        x, format="%d.%m.%y-%H:%M:%S.%f", errors="coerce"
                    )
                ).dropna(how="all", axis=1)
                val_cols = temp_stats.apply(
                    lambda x: pd.to_numeric(x, errors="coerce")
                ).dropna(how="all", axis=1)
                temp_stats = pd.concat([time_col, val_cols], axis=1).reset_index(
                    drop=True
                )
                temp_stats.columns = TEMP_STATS_COLS
                temp_stats["Timestamp"] = temp_stats["Timestamp"].dt.floor("s")
                stats_dfs.append(temp_stats)

                # Subset the log messages (relevant for label and pickup)
                df = df[df.apply(lambda row: row.str.contains("drops")).any(axis=1)]
                df = df.dropna(how="all", axis=1)

                if len(df) == 0:
                    # Can exit from dispense log here
                    if key == "Dispense":
                        dfs.append(pd.DataFrame())
                    else:
                        # In case log from aborted process is included in the output
                        self.logger.warning(
                            f"Skipping: {file_path}. Likely an aborted process "
                            "during sample-prep."
                        )

                else:
                    df.columns = DROPLET_COLS
                    df = df.reset_index(drop=True)

                    if key == "Label":
                        df = df.drop(["PlatePos"], axis=1)

                    # Quick fix for plate number for now
                    if key == "Pickup":
                        df["Plate"] = self._pickup_plate(file_path)

                    dfs.append(df)

        dfs = pd.concat(dfs, axis=0)
        stats_dfs = pd.concat(stats_dfs, axis=0)
        return dfs, stats_dfs

    @staticmethod
    def _pickup_plate(file_path: str) -> int:
        """Plate number for a pickup log, from its name or its run directory.

        A pickup log is recognized by either, so the number has to be looked
        for in both: `P7_..._Logfile.log` inside `pickup_7_....Run` carries it
        only in the directory.
        """
        for candidate in (
            os.path.basename(file_path),
            os.path.basename(os.path.dirname(file_path)),
        ):
            match = re.search(r"pickup_(\d+)", candidate, flags=re.IGNORECASE)
            if match:
                return int(match.group(1))

        raise ValueError(
            f"Cannot read the pickup plate number from {file_path}. Expected "
            "'pickup_<number>' in the file or run directory name."
        )

    def label_well_plex(self, label_df: pd.DataFrame):
        label_mapping = CELLEONE_MAPPING[self.label_type][self.plex]
        label_df["Plex"] = label_df["Well"].map(label_mapping)

        return label_df

    def _map_coords(
        self,
        geo_df,
        map_df,
        coord_cols=("XPos", "YPos"),
        group_cols=("Target", "Field"),
    ) -> tuple[pd.Series, pd.Series]:
        """Match each row of `map_df` to its nearest `plex` rows in `geo_df`.

        Returns the matched geoprops indices and their distances, both as
        Series of arrays indexed like `map_df` so they can be assigned straight
        back onto it. With `plex` unset every candidate is kept, and the caller
        resolves cells claimed more than once by distance.
        """
        coord_cols, group_cols = list(coord_cols), list(group_cols)
        grouped_geo = geo_df.groupby(group_cols)
        grouped_map = map_df.groupby(group_cols)

        closest, closest_distance = {}, {}
        for group_key, geo_group in grouped_geo:
            if group_key not in grouped_map.groups:
                self.logger.info(f"{group_key} not in pickup groups, skipping.")
                continue

            map_group = grouped_map.get_group(group_key)

            geo_coords = np.stack(geo_group[coord_cols].values)
            map_coords = np.stack(map_group[coord_cols].values)

            distances = np.linalg.norm(
                map_coords[:, np.newaxis, :] - geo_coords[np.newaxis, :, :], axis=2
            )
            sorted_indices = np.argsort(distances, axis=1)[:, : self.plex]
            closest_points = geo_group.index.values[sorted_indices]

            for i, map_idx in enumerate(map_group.index):
                closest[map_idx] = closest_points[i]
                closest_distance[map_idx] = distances[i][sorted_indices[i]]

        return (
            pd.Series(closest, dtype=object).reindex(map_df.index),
            pd.Series(closest_distance, dtype=object).reindex(map_df.index),
        )

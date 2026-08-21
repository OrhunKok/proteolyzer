"""Conversion of raw search-engine output into the parquet inputs of the pipeline.

Each search directory is filtered down to the columns and rows the downstream
stages need, then written next to the original text export as parquet.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from .base import Stage
from .config import Config
from .utils import column_mapping

CONFIG = Config()


class Preprocessor:
    """Namespace grouping the per-search-engine preprocessors."""

    class MaxQuant(Stage):
        FILES = CONFIG.MaxQuant.FILES
        FILES_NEEDED = CONFIG.MaxQuant.FILES_NEEDED
        LOAD_COLS = CONFIG.MaxQuant.LOAD_COLS
        RENAME_MAP = CONFIG.MaxQuant.COLS_RENAME_MAPPING

        def __init__(self, params, queue=None):
            super().__init__(params, queue)
            self.translated_frames = self.params["Translation"][
                "Translated Frames Folder"
            ]
            # Dependent peptides are filtered here rather than in the detection
            # stage, so this stage owns the configured threshold.
            self.detection_pep = float(self.params["Detection"]["Detection PEP"])

        def run(self):
            self.record_run()
            self.translated_frames.mkdir(parents=True, exist_ok=True)
            (self.output_dir / "ALT").mkdir(parents=True, exist_ok=True)
            (self.output_dir / "SAAP").mkdir(parents=True, exist_ok=True)

            converted = sum(
                self._convert_directory(subdir) for subdir in self._search_directories()
            )

            if not converted:
                expected = sorted(
                    {f for group in self.FILES_NEEDED.values() for f in group}
                )
                raise FileNotFoundError(
                    f"No MaxQuant search output found under {self.data_dir}. "
                    "Expected one directory per sample, or the search files "
                    f"directly in it, containing any of: {expected}."
                )

        def _search_directories(self):
            """Every directory that may hold search output, `data_dir` included.

            Yielding only the child names would skip a flat layout, where the
            text files sit directly in the data folder.
            """
            data_dir = Path(self.data_dir)
            if not data_dir.is_dir():
                raise FileNotFoundError(
                    f"Data Folder does not exist or is not a directory: {data_dir}"
                )

            yield data_dir
            for root, dirs, _ in os.walk(data_dir):
                for dirname in sorted(dirs):
                    yield Path(root) / dirname

        def _convert_directory(self, subdir: Path) -> int:
            """Converts the needed search files in `subdir`; returns how many."""
            files = {entry.name for entry in subdir.iterdir() if entry.is_file()}

            # Detect type of search
            search_type = (
                "Detection" if "dependentPeptides.txt" in files else "Validation"
            )
            needed = {f + ".txt" for f in self.FILES_NEEDED[search_type]}

            count = 0
            for filename in sorted(files & needed):
                txt_path = subdir / filename
                base = filename.removesuffix(".txt")
                parquet_path = subdir / f"{base}.parquet"

                cols2keep = self.LOAD_COLS[base]
                df = pd.read_csv(txt_path, delimiter="\t")

                if (
                    base == "evidence"
                    and self.label_setup == "TMT"
                    and "_val" in subdir.name
                ):
                    cols2keep = cols2keep + [
                        f"Reporter intensity corrected {i}"
                        for i in range(1, self.label_plex + 1)
                    ]

                df = column_mapping(df, cols2keep)
                df = self._filter(base, df)

                df = df.reset_index(drop=True)
                df.to_parquet(parquet_path, engine="fastparquet")
                self.queue.put(
                    ("stdout", f"{subdir.name} converted {txt_path} -> {parquet_path}")
                )
                count += 1

            return count

        def _filter(self, base: str, df: pd.DataFrame) -> pd.DataFrame:
            """Apply the per-file row filtering and derived columns."""
            filters = {
                "allPeptides": self._allpeptides,
                "evidence": self._evidence,
                "msms": self._msms,
                "peptides": self._peptides,
            }
            transform = filters.get(base)
            return df if transform is None else transform(df)

        def _evidence(self, df: pd.DataFrame) -> pd.DataFrame:
            return df[(df["Reverse"].isna()) & (df["Potential contaminant"].isna())]

        def _allpeptides(self, df: pd.DataFrame) -> pd.DataFrame:
            df = df.copy()
            df["Leading.Razor.DP.Protein"] = df["DP Proteins"].str.split(";").str[0]
            df["Con.DP.Protein"] = df["DP Proteins"].str.contains("CON__", na=False)

            df = df[(df["Reverse"].isna()) & (~df["Con.DP.Protein"])]

            df = df.loc[
                df["DP Mass Difference"].notna() & (df["DP PEP"] <= self.detection_pep)
            ]

            df = df[df["MSMS Scan Numbers"].notna()].copy()
            df["MSMS Scan Numbers"] = (
                df["MSMS Scan Numbers"]
                .str.split(";")
                .apply(lambda x: list(map(int, x)))
            )

            return df.rename(columns={"MSMS Scan Numbers": "MS/MS scan number"})

        def _msms(self, df: pd.DataFrame) -> pd.DataFrame:
            df = df[df["Reverse"].isna()].copy()
            df = df.rename(columns={"Scan number": "MS/MS scan number"})

            df["Matches"] = df["Matches"].str.split(";")
            df = df[df["Matches"].notna()].copy()

            df["Matches"] = df["Matches"].apply(
                lambda matches: [
                    frag
                    for frag in matches
                    if "a" not in frag
                    and "(" not in frag
                    and "NH3" not in frag
                    and "H2O" not in frag
                ]
            )

            df["Frag.Type"] = df["Matches"].apply(lambda m: [frag[0] for frag in m])
            df["Frag.Number"] = df["Matches"].apply(
                lambda m: [int(frag[1:]) for frag in m]
            )

            return df.drop(columns="Matches").explode(["Frag.Type", "Frag.Number"])

        def _peptides(self, df: pd.DataFrame) -> pd.DataFrame:
            df = df[(df["Reverse"].isna()) & (df["Potential contaminant"].isna())]
            df = df[
                (df["Start position"] == 1) | (df["Amino acid after"] == "-")
            ].copy()

            df["Terminus"] = np.where(df["Amino acid after"] == "-", "C", "N")
            return df

    class DIANN:
        FILES = CONFIG.DIANN.FILES
        FILE_EXT = CONFIG.DIANN.FILE_EXTENSIONS

        def __init__(self, data_folder: str):
            self.data_folder = Path(data_folder)

        def run(self):
            raise NotImplementedError(
                "No DIA-NN preprocessing pipeline is implemented yet."
            )

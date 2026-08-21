"""Read back what the AAS pipeline produced.

The stages write one frame per sample per step, named after the step that
wrote them and split across ``SAAP/`` and ``ALT/``. That is fine for the
pipeline and unhelpful for the person asking "what did we find?", who has to
know that the answer is ``SAAP/<sample>_SAAP_Quant.parquet`` and that stage 2
is validated while stage 1 is not.

:class:`Results` is the way in: it finds the samples, says which steps
completed for each, and combines a step's frames across samples.

    from proteolyzer import aas

    results = aas.Results("out/")
    results.samples                     # what is in there
    results.summary()                   # rows per sample per step
    results.combined("quantified")      # every sample, with a Sample column
    results.provenance()                # how it was produced
"""

import json
from pathlib import Path

import pandas as pd

from proteolyzer.core.io import frame_exists, read_frame
from proteolyzer.core.logging import Logged
from proteolyzer.core.pipeline import PROVENANCE_FILE

#: What the stages write, under a name describing the result rather than the
#: step: (subdirectory, file stem). The order is the order of the pipeline.
ARTEFACTS: dict[str, tuple[str, str]] = {
    "candidates": ("SAAP", "{sample}_SAAP"),
    "alt": ("ALT", "{sample}_ALT"),
    "filtered": ("SAAP", "{sample}_SAAP_Filtered_Stage_1"),
    "fasta_entries": ("SAAP", "{sample}_FASTA"),
    "validated": ("SAAP", "{sample}_SAAP_Filtered_Stage_2"),
    "evidence": ("SAAP", "{sample}_Val_Evidence_Filtered_Stage_2"),
    "quantified": ("SAAP", "{sample}_SAAP_Quant"),
}

#: Column the sample identifier is added under when frames are combined.
SAMPLE_COL = "Sample"


class Results(Logged):
    """The output folder of an AAS run."""

    def __init__(self, output_dir: str | Path):
        self.output_dir = Path(output_dir)

    @classmethod
    def from_params(cls, params) -> Results:
        """Open the results of the run described by a parameter file or dict."""
        from .params import load_params

        return cls(load_params(params)["Utils"]["Output Folder"])

    def path(self, artefact: str, sample: str) -> Path:
        """Where `artefact` for `sample` lives, whether or not it exists."""
        if artefact not in ARTEFACTS:
            raise ValueError(
                f"No artefact '{artefact}'. Available: {sorted(ARTEFACTS)}"
            )
        subdir, stem = ARTEFACTS[artefact]
        return self.output_dir / subdir / stem.format(sample=sample)

    def has(self, artefact: str, sample: str) -> bool:
        """Whether the step that writes `artefact` got as far as `sample`."""
        return frame_exists(self.path(artefact, sample))

    @property
    def samples(self) -> list[str]:
        """Every sample with at least one artefact, in a stable order."""
        found: set[str] = set()
        for subdir, stem in ARTEFACTS.values():
            suffix = stem.replace("{sample}", "")
            directory = self.output_dir / subdir
            if not directory.is_dir():
                continue
            for path in directory.iterdir():
                # Match the stem exactly: sample ids contain underscores, so
                # splitting on them would truncate them.
                name = path.name.removesuffix(path.suffix)
                if name.endswith(suffix) and path.suffix == ".parquet":
                    found.add(name.removesuffix(suffix))
        return sorted(found)

    def load(self, artefact: str, sample: str) -> pd.DataFrame:
        """One sample's frame for `artefact`."""
        return read_frame(self.path(artefact, sample))

    def combined(self, artefact: str, samples: list[str] | None = None):
        """Every sample's `artefact` in one frame, with a Sample column.

        Samples the step did not reach are skipped and reported, so a partial
        run still returns what it did produce.
        """
        wanted = self.samples if samples is None else samples
        frames, missing = [], []
        for sample in wanted:
            if not self.has(artefact, sample):
                missing.append(sample)
                continue
            frame = self.load(artefact, sample)
            frames.append(frame.assign(**{SAMPLE_COL: sample}))

        if missing:
            self.logger.warning(
                f"No '{artefact}' for {len(missing)} of {len(wanted)} samples: "
                f"{missing}"
            )
        if not frames:
            return pd.DataFrame(columns=[SAMPLE_COL])

        combined = pd.concat(frames, ignore_index=True)
        # Lead with the sample so the frame reads as one table of results.
        return combined[[SAMPLE_COL, *combined.columns.drop(SAMPLE_COL)]]

    def summary(self) -> pd.DataFrame:
        """Rows per sample per step; NA where the step did not run.

        Reading down a column shows where the pipeline stopped.
        """
        samples = self.samples
        rows = {
            sample: {
                artefact: (
                    len(self.load(artefact, sample))
                    if self.has(artefact, sample)
                    else pd.NA
                )
                for artefact in ARTEFACTS
            }
            for sample in samples
        }
        summary = pd.DataFrame.from_dict(rows, orient="index", dtype="Int64")
        summary.index.name = SAMPLE_COL
        return summary

    def provenance(self) -> pd.DataFrame:
        """What ran in this folder, one row per stage run, newest last."""
        log = self.output_dir / PROVENANCE_FILE
        if not log.exists():
            self.logger.warning(f"No {PROVENANCE_FILE} in {self.output_dir}")
            return pd.DataFrame()
        entries = [
            json.loads(line) for line in log.read_text().splitlines() if line.strip()
        ]
        return pd.DataFrame(entries)

"""The AAS pipeline's shared stage behaviour.

:class:`proteolyzer.core.pipeline.Stage` supplies what any pipeline needs —
parameters, a progress queue, provenance. This adds what is specific to this
one: parameters come from the AAS parameter file, and the samples to process
come from an experiment metadata spreadsheet.
"""

from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from proteolyzer.core.pipeline import PROVENANCE_FILE, NullQueue
from proteolyzer.core.pipeline import Stage as _CoreStage

from .config import Config
from .params import load_params

CONFIG = Config()

__all__ = ["PROVENANCE_FILE", "NullQueue", "Stage"]

#: Columns required from the experiment metadata spreadsheet.
METADATA_COLS: tuple[str, ...] = (
    "TMT plex",
    "TMT channel",
    "ParticipantID",
    "Group",
    "MQ",
    "sample_name",
    "sample_ID",
)


class Stage(_CoreStage):
    """Base class for the stages of the AAS pipeline.

    Parameters
    ----------
    params
        Path to a parameter YAML file, or an equivalent mapping.
    queue
        Optional queue receiving ``(stream, payload)`` progress messages.
    """

    def __init__(self, params, queue=None):
        # The unmerged form is kept so a stage can hand it to another stage
        # (see Detection._initialize_workflow) without re-merging merged params.
        self._raw_params = params
        super().__init__(load_params(params), queue)

        utils_params = self.params["Utils"]
        self.workflow: str = utils_params["Workflow"]
        self.data_dir: Path = utils_params["Data Folder"]
        self.output_dir: Path = utils_params["Output Folder"]
        self.metadata_file: Path = utils_params["Metadata File"]
        self.label_setup: str = utils_params["Labelling Setup"]
        self.label_plex: int = utils_params["Label Plex"]

        self.announce()

    def record_run(self) -> Path:
        """Log this run in the configured output folder."""
        return super().record_run(self.output_dir)

    @cached_property
    def metadata(self) -> pd.DataFrame:
        """Experiment metadata, with TMT channels mapped to MaxQuant indices."""
        metadata = pd.read_excel(
            self.metadata_file, engine="openpyxl", usecols=list(METADATA_COLS)
        )
        metadata["MQ"] = metadata["TMT channel"].map(CONFIG.TMT.MQ_TMT_MAP)
        return metadata

    @property
    def samples(self) -> np.ndarray:
        """The unique sample IDs to process, in a stable order."""
        return np.unique(self.metadata["sample_ID"])

    def run(self) -> None:
        """Process every sample in the metadata."""
        self.record_run()
        for sample in self.samples:
            self.process_sample(sample)

    def process_sample(self, sample: str) -> None:
        """Handle one sample; implemented by each stage."""
        raise NotImplementedError

    def _locate_sample_dir(self, sample: str, suffix: str = "") -> Path | None:
        """Find the directory holding `sample`'s search output, if any."""
        target = f"{sample}{suffix}"
        data_dir = Path(self.data_dir)
        if data_dir.name == target:
            return data_dir
        for candidate in data_dir.rglob(target):
            if candidate.is_dir():
                return candidate
        return None

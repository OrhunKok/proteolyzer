"""Shared plumbing for the AAS pipeline stages.

Every stage (preprocessing, translation, detection, validation,
quantification) is constructed the same way: from a parameter file or dict,
with an optional multiprocessing queue for progress reporting. :class:`Stage`
holds that plumbing so each stage module only contains its own analysis.
"""

import datetime
import json
from functools import cached_property
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Config
from .params import load_params

CONFIG = Config()

#: One JSON object per stage run, appended in the output folder.
PROVENANCE_FILE = "provenance.jsonl"

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


class NullQueue:
    """A 'do-nothing' queue to replace the real Queue when multiprocessing isn't used."""

    def put(self, item):
        pass


class Stage:
    """Base class for the stages of the AAS pipeline.

    Parameters
    ----------
    params
        Path to a parameter YAML file, or an equivalent mapping.
    queue
        Optional queue receiving ``(stream, payload)`` progress messages.
        Defaults to a :class:`NullQueue`, which discards them.
    """

    def __init__(self, params, queue=None):
        # The unmerged form is kept so a stage can hand it to another stage
        # (see Detection._initialize_workflow) without re-merging merged params.
        self._raw_params = params
        self.params = load_params(params)
        self.queue = queue if queue is not None else NullQueue()

        utils_params = self.params["Utils"]
        self.workflow: str = utils_params["Workflow"]
        self.data_dir: Path = utils_params["Data Folder"]
        self.output_dir: Path = utils_params["Output Folder"]
        self.metadata_file: Path = utils_params["Metadata File"]
        self.label_setup: str = utils_params["Labelling Setup"]
        self.label_plex: int = utils_params["Label Plex"]

        self.queue.put(("stdout", f"{type(self).__name__} initialized."))

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

    def record_run(self) -> Path:
        """Append this stage's parameters to the provenance log.

        Outputs are keyed only by sample, so a re-run with different thresholds
        overwrites the previous results. The log is what makes an output folder
        self-describing: which stage ran, when, from which version, and with
        which resolved parameters.
        """
        try:
            package_version = version("proteolyzer")
        except PackageNotFoundError:  # pragma: no cover - source checkout
            package_version = "unknown"

        entry = {
            "stage": type(self).__name__,
            "timestamp": datetime.datetime.now(datetime.UTC).isoformat(
                timespec="seconds"
            ),
            "proteolyzer_version": package_version,
            "params": self.params,
        }

        output_dir = Path(self.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        log = output_dir / PROVENANCE_FILE
        with open(log, "a") as f:
            # default=str renders the Path values the schema produces.
            f.write(json.dumps(entry, default=str, sort_keys=True) + "\n")
        return log

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

"""Shared plumbing for a pipeline stage.

A stage is one step of an analysis: it is constructed from resolved
parameters, reports progress on a queue, and records what it ran so a results
folder describes itself. None of that is specific to a particular pipeline, so
it lives here rather than in the subpackage that happened to need it first.

Subclasses add their own inputs and implement :meth:`run`.
"""

import datetime
import json
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

#: One JSON object per stage run, appended in the output folder.
PROVENANCE_FILE = "provenance.jsonl"


class NullQueue:
    """A 'do-nothing' queue to replace the real Queue when multiprocessing isn't used."""

    def put(self, item):
        pass


class Stage:
    """Base class for a pipeline stage.

    Parameters
    ----------
    params
        The resolved parameters for the whole pipeline. A subclass that reads
        them from a file resolves them before calling this.
    queue
        Optional queue receiving ``(stream, payload)`` progress messages.
        Defaults to a :class:`NullQueue`, which discards them.
    """

    def __init__(self, params: dict, queue=None):
        self.params = params
        self.queue = queue if queue is not None else NullQueue()

    def announce(self) -> None:
        """Report that the stage is ready, by name."""
        self.queue.put(("stdout", f"{type(self).__name__} initialized."))

    def run(self) -> None:  # pragma: no cover - implemented by subclasses
        raise NotImplementedError

    def record_run(self, output_dir: Path) -> Path:
        """Append this stage's parameters to the provenance log in `output_dir`.

        Outputs are usually keyed only by sample, so a re-run with different
        thresholds overwrites the previous results. This log is what makes an
        output folder self-describing: which stage ran, when, from which
        version, and with which resolved parameters.
        """
        entry = {
            "stage": type(self).__name__,
            "timestamp": datetime.datetime.now(datetime.UTC).isoformat(
                timespec="seconds"
            ),
            "proteolyzer_version": package_version(),
            "params": self.params,
        }

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        log = output_dir / PROVENANCE_FILE
        with open(log, "a") as f:
            # default=str renders the Path values a parameter schema produces.
            f.write(json.dumps(entry, default=str, sort_keys=True) + "\n")
        return log


def package_version() -> str:
    """The installed proteolyzer version, or "unknown" from a bare checkout."""
    try:
        return version("proteolyzer")
    except PackageNotFoundError:  # pragma: no cover - source checkout
        return "unknown"

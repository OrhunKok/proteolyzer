"""FrameTranslator writes the six-frame translations the detection stage reads.

quickdna ships as a platform-specific wheel that cannot be installed
everywhere, and the stage imports it lazily, so a stub stands in for it here.
That keeps everything except the third-party translation itself under test.
"""

import pickle
import sys
import types
from pathlib import Path

import pytest

pytest.importorskip("yaml")
pytest.importorskip("Bio")

from proteolyzer.aas.base import PROVENANCE_FILE  # noqa: E402
from proteolyzer.aas.translation import FrameTranslator  # noqa: E402

#: What the stub "translates" each record into, one string per frame.
STUB_FRAMES = ("MAIV", "AIVX", "IVXX", "VXXM", "XXMA", "XMAI")


@pytest.fixture
def stub_quickdna(monkeypatch):
    """Install a minimal quickdna whose translation is deterministic."""

    class DnaSequence:
        def __init__(self, seq):
            self.seq = seq

        def translate_all_frames(self):
            # Real quickdna returns Protein objects; only str() is used.
            return [f"{frame}{len(self.seq)}" for frame in STUB_FRAMES]

    module = types.ModuleType("quickdna")
    module.DnaSequence = DnaSequence
    monkeypatch.setitem(sys.modules, "quickdna", module)
    return module


@pytest.fixture
def genome(aas_params):
    """Two FASTA records, so the per-frame concatenation is observable."""
    path = Path(aas_params["Translation"]["Genome FASTA"])
    path.write_text(">chr1\nATGGCCATTGTAATG\n>chr2\nATGGCC\n")
    return path


def test_run_writes_both_variants_of_every_frame(aas_params, genome, stub_quickdna):
    translator = FrameTranslator(aas_params)
    translator.run()

    frames_dir = Path(aas_params["Translation"]["Translated Frames Folder"])
    written = {p.name for p in frames_dir.glob("*.p")}
    assert written == {
        name
        for frame in range(1, 7)
        for name in (f"frame_{frame}.p", f"frame_{frame}_il_ambigous.p")
    }


def test_frames_concatenate_every_record(aas_params, genome, stub_quickdna):
    FrameTranslator(aas_params).run()

    frames_dir = Path(aas_params["Translation"]["Translated Frames Folder"])
    with open(frames_dir / "frame_1.p", "rb") as f:
        frame_1 = pickle.load(f)

    # First frame of record 1 (length 15) followed by record 2 (length 6).
    assert frame_1 == "MAIV15MAIV6"


def test_isoleucine_is_collapsed_onto_leucine(aas_params, genome, stub_quickdna):
    FrameTranslator(aas_params).run()

    frames_dir = Path(aas_params["Translation"]["Translated Frames Folder"])
    with open(frames_dir / "frame_1.p", "rb") as f:
        exact = pickle.load(f)
    with open(frames_dir / "frame_1_il_ambigous.p", "rb") as f:
        ambiguous = pickle.load(f)

    assert "I" in exact
    assert "I" not in ambiguous
    assert ambiguous == exact.replace("I", "L")


def test_entry_count_drives_the_progress_reports(aas_params, genome, stub_quickdna):
    class RecordingQueue:
        def __init__(self):
            self.messages = []

        def put(self, item):
            self.messages.append(item)

    queue = RecordingQueue()
    FrameTranslator(aas_params, queue).run()

    progress = [payload for stream, payload in queue.messages if stream == "progress"]
    assert progress == [(1, 2), (2, 2)]


def test_run_is_recorded_in_the_provenance_log(aas_params, genome, stub_quickdna):
    FrameTranslator(aas_params).run()

    log = Path(aas_params["Utils"]["Output Folder"]) / PROVENANCE_FILE
    assert "FrameTranslator" in log.read_text()


def test_missing_quickdna_is_reported_only_when_translating(aas_params, genome):
    """Constructing the stage must not require the optional wheel."""
    translator = FrameTranslator(aas_params)  # no stub installed
    assert translator.frames_folder.is_dir()

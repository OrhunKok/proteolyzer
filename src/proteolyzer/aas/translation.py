"""Six-frame translation of a genome FASTA into searchable protein strings."""

import pickle

import numpy as np
from Bio import SeqIO

from .base import Stage


class FrameTranslator(Stage):
    """Translates every FASTA entry in all six frames.

    Two versions of each frame are written: the translation as-is, and one with
    isoleucine rewritten to leucine, since the two are isobaric and cannot be
    distinguished by mass.
    """

    #: str.translate table collapsing I onto L.
    ISOLEUCINE_LEUCINE = {73: 76}

    def __init__(self, params, queue=None):
        super().__init__(params, queue)

        self.path_to_fasta = self.params["Translation"]["Genome FASTA"]
        self.frames_folder = self.params["Translation"]["Translated Frames Folder"]
        self.frames_folder.mkdir(parents=True, exist_ok=True)

    def run(self):
        self.record_run()
        self.queue.put(("stdout", f"{type(self).__name__} commencing."))
        self._count_entries()
        self._translate_sequences()
        self._write_outputs()

    def _count_entries(self):
        with open(self.path_to_fasta) as handle:
            self.entry_counts = sum(1 for _ in SeqIO.parse(handle, "fasta"))

        self.translated_frames = np.zeros((self.entry_counts, 6), dtype=object)
        self.translated_frames_il_ambiguous = np.zeros(
            (self.entry_counts, 6), dtype=object
        )

    def _translate_sequences(self):
        # Imported here so the rest of the pipeline stays usable without the
        # platform-specific quickdna wheel installed.
        from quickdna import DnaSequence

        fasta = SeqIO.parse(self.path_to_fasta, "fasta")

        for i, record in enumerate(fasta):
            self.queue.put(("progress", (i + 1, self.entry_counts)))

            translation_frames = [
                str(frame)
                for frame in DnaSequence(seq=str(record.seq)).translate_all_frames()
            ]

            self.translated_frames[i, :] = translation_frames
            self.translated_frames_il_ambiguous[i, :] = [
                frame.translate(self.ISOLEUCINE_LEUCINE) for frame in translation_frames
            ]

    def _write_outputs(self):
        for frame in range(6):
            il_ambigous_frame = "".join(self.translated_frames_il_ambiguous[:, frame])
            translated_frame = "".join(self.translated_frames[:, frame])

            with open(
                self.frames_folder / f"frame_{frame + 1}_il_ambigous.p", "wb"
            ) as f:
                pickle.dump(il_ambigous_frame, f)

            with open(self.frames_folder / f"frame_{frame + 1}.p", "wb") as f:
                pickle.dump(translated_frame, f)

            self.queue.put(("stdout", f"Frame {frame + 1} written."))

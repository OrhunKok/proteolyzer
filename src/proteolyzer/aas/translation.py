from Bio import SeqIO
import numpy as np
import pickle
from typing import Dict, Union
from quickdna import DnaSequence
from multiprocessing import Queue
from .config import Config
from . import utils

CONFIG = Config()


class FrameTranslator:
    def __init__(self, params : Union[str, Dict], queue: Queue = None):
        if isinstance(params, str):
            self.params = utils._load_params(params)
        elif isinstance(params, dict):
            self.params = params
        else:
            raise ValueError("params must be a file path or dict")
        self.queue = queue if queue is not None else utils.NullQueue()

        self.path_to_fasta = self.params.get('Translation').get('Genome FASTA')
        self.frames_folder = self.params.get('Translation').get('Translated Frames Folder')
        self.frames_folder.mkdir(parents=True, exist_ok=True)
        self.isoleucine_leucine_ascii = {73 : 76}  # I → L

        self.queue.put(("stdout", f"{self.__class__.__name__} initialized."))

    def run(self):
        self.queue.put(("stdout", f"{self.__class__.__name__} commencing."))
        self._count_entries()
        self._translate_sequences()
        self._write_outputs()

    def _count_entries(self):
        
        self.entry_counts = 0
        with open(self.path_to_fasta, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            self.entry_counts = sum(1 for _ in fasta)

        self.translated_frames = np.zeros((self.entry_counts, 6), dtype=object)
        self.translated_frames_il_ambiguous = np.zeros((self.entry_counts, 6), dtype=object)

    def _translate_sequences(self):
        fasta = SeqIO.parse(self.path_to_fasta, "fasta")

        for i, record in enumerate(fasta):
            self.queue.put(("progress", (i + 1, self.entry_counts)))

            seq = str(record.seq)
            translation_frames = DnaSequence(seq=seq).translate_all_frames()
            translation_frames = [str(x) for x in translation_frames]

            self.translated_frames[i, :] = translation_frames
            il_ambiguous_translations = [x.translate(self.isoleucine_leucine_ascii) for x in translation_frames]
            self.translated_frames_il_ambiguous[i, :] = il_ambiguous_translations

    def _write_outputs(self):
        for frame in range(6):
            il_ambigous_frame = ''.join(self.translated_frames_il_ambiguous[:, frame])
            translated_frame = ''.join(self.translated_frames[:, frame])

            with open(self.frames_folder / f'frame_{frame + 1}_il_ambigous.p', 'wb') as f:
                pickle.dump(il_ambigous_frame, f)

            with open(self.frames_folder / f'frame_{frame + 1}.p', 'wb') as f:
                pickle.dump(translated_frame, f)

            self.queue.put(("stdout", f"Frame {frame + 1} written."))
from Bio import SeqIO
import numpy as np
import pickle
from quickdna import DnaSequence
from multiprocessing import Queue
from Utils import NullQueue


class FrameTranslator:
    def __init__(self, params : dict, queue: Queue = None):
        self.params = params
        self.queue = queue if queue is not None else NullQueue()

        self.path_to_fasta = params.get('Genome FASTA')
        self.output_dir = params.get('Output Folder')
        self.isoleucine_leucine_ascii = {73: 76}  # I → L
        (self.output_dir / 'TranslatedFrames').mkdir(parents=True, exist_ok=True)

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

            with open(self.output_dir / 'TranslatedFrames' / f'frame_{frame + 1}_il_ambigous.p', 'wb') as f:
                pickle.dump(il_ambigous_frame, f)

            with open(self.output_dir / 'TranslatedFrames' / f'frame_{frame + 1}.p', 'wb') as f:
                pickle.dump(translated_frame, f)

            self.queue.put(("stdout", f"Frame {frame + 1} written."))
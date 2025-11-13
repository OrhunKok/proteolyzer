from Bio import SeqIO
import numpy as np
import pickle
from quickdna import DnaSequence
from multiprocessing import Queue


class FrameTranslator:
    def __init__(self, params : dict, queue: Queue = None):
        self.params = params
        self.queue = queue

        self.path_to_fasta = params.get('Genome FASTA')
        self.output_dir = params.get('Output Folder')

        self.translation_arr_all_frames = None
        self.translated_arr_ambiguous = None
        self.entry_counts = 0
        self.isoleucine_leucine_ascii = {73: 76}  # I → L

        (self.output_dir / 'Translation_Frames').mkdir(parents=True, exist_ok=True)

        if self.queue:
            self.queue.put(("stdout", "FrameTranslator initialized."))
            self.queue.put(("stdout", "FrameTranslator commencing."))

    def run(self):
        self._count_entries()
        self._translate_sequences()
        self._write_outputs()

    def _count_entries(self):
        with open(self.path_to_fasta, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            self.entry_counts = sum(1 for _ in fasta)

        self.translation_arr_all_frames = np.zeros((self.entry_counts, 6), dtype=object)
        self.translated_arr_ambiguous = np.zeros((self.entry_counts, 6), dtype=object)

    def _translate_sequences(self):
        fasta = SeqIO.parse(self.path_to_fasta, "fasta")

        for i, record in enumerate(fasta):
            if self.queue:
                self.queue.put(("progress", (i + 1, self.entry_counts)))

            seq = str(record.seq)
            translation_frames = DnaSequence(seq=seq).translate_all_frames()
            translation_frames = [str(x) for x in translation_frames]

            self.translation_arr_all_frames[i, :] = translation_frames
            ambiguous_translations = [x.translate(self.isoleucine_leucine_ascii) for x in translation_frames]
            self.translated_arr_ambiguous[i, :] = ambiguous_translations

    def _write_outputs(self):
        for frame in range(6):
            W_aa = ''.join(self.translated_arr_ambiguous[:, frame])
            S_aa = ''.join(self.translation_arr_all_frames[:, frame])

            with open(self.output_dir / 'Translation_Frames' / f'W{frame + 1}_aa_ambig.p', 'wb') as f:
                pickle.dump(W_aa, f)

            with open(self.output_dir / 'Translation_Frames' / f's{frame + 1}a_ambig.p', 'wb') as f:
                pickle.dump(S_aa, f)

            if self.queue:
                self.queue.put(("stdout", f"Frame {frame + 1} written."))
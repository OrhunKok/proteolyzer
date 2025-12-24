import pandas as pd
import numpy as np
import os
import pickle
from queue import Queue
from typing import Dict, Union
from pathlib import Path
from . import utils
from .config import Config

CONFIG = Config()

class Quantification:
    def __init__(self, params: Union[str, Dict], queue: Queue = None):
        if isinstance(params, str):
            self.params = utils._load_params(params)
        elif isinstance(params, dict):
            self.params = params
        else:
            raise ValueError("params must be a file path or dict")
        self.queue = queue if queue is not None else utils.NullQueue()

        self.label_setup = self.params.get('Labelling Setup')
        self.metadata = pd.read_excel(self.params.get('Utils').get('Metadata File'), engine = 'openpyxl', usecols=['TMT plex', 'TMT channel', 'ParticipantID', 'Group', 'MQ', 'sample_name', 'sample_ID'])
        self.metadata['MQ'] = self.metadata['TMT channel'].map(CONFIG.TMT.MQ_TMT_MAP)
        self.data_dir = self.params.get('Utils').get('Data Folder')
        self.output_dir = self.params.get('Utils').get('Output Folder')
        # self.min_quant = self.params.get('Quantification).get('Minimum Quantity')

        self.queue.put(("stdout", f"{self.__class__.__name__} initialized."))

    def run(self):
        for sample in np.unique(self.metadata['sample_ID']):
            self.process_sample(sample)

    def process_sample(self, sample):
        sample_val_dir = self._locate_sample_dir(sample)
        if not sample_val_dir:
            self.queue.put(("stderr", f"{sample} validation directory not found"))
            return

        val_evidence_path = self.output_dir / 'MTP' / f'{sample}_Val_Evidence_Filtered_Stage_2.p'
        mtp_path = self.output_dir / 'MTP' / f'{sample}_MTP_Filtered_Stage_2.p'

        if not (os.path.exists(val_evidence_path) and os.path.exists(mtp_path)):
            self.queue.put(("stderr", f"Missing files for sample {sample}"))
            return

        with open(val_evidence_path, 'rb') as f:
            val_evidence = pickle.load(f)
        with open(mtp_path, 'rb') as f:
            mtp = pickle.load(f)

        ev_filter_seqs = np.unique(mtp[['DP Base Sequence', 'mistranslated sequence']].to_numpy())
        val_evidence = val_evidence[val_evidence['Sequence'].isin(ev_filter_seqs)].drop(['Frag.Type', 'Frag.Number'], axis=1).drop_duplicates()

        val_evidence = val_evidence.drop(['Raw file', 'Charge', 'MS/MS scan number'], axis=1).groupby('Sequence').sum()

        mtp = mtp[(mtp['DP Base Sequence'].isin(val_evidence.index)) & (mtp['mistranslated sequence'].isin(val_evidence.index))]
        mtp_df = val_evidence.loc[mtp['mistranslated sequence']]
        bp_df = val_evidence.loc[mtp['DP Base Sequence']]

        with np.errstate(divide='ignore'):
            if 'Reporter intensity corrected 1' in val_evidence.columns:
                quant = self._raas(mtp, mtp_df, bp_df, self.metadata, label_designation='TMT')
            else:
                quant = self._raas(mtp, mtp_df, bp_df, self.metadata, label_designation='Label-Free')

        quant = quant[np.isfinite(quant['Ratio']) & (~np.isnan(quant['Ratio']))]

        output_path = self.output_dir / 'MTP' / f'{sample}_MTP_Quant.p'
        with open(output_path, 'wb') as f:
            pickle.dump(quant, f)

        self.queue.put(("stdout", f"Quantification complete for sample: {sample}"))

    def _locate_sample_dir(self, sample):
        data_dir = Path(self.data_dir)
        for subdir in data_dir.rglob('*'):
            if subdir.is_dir() and subdir.name == f"{sample}_val":
                return subdir
        return None

    def _raas(self, mtp, mtp_df, bp_df, sample_df, label_designation):
        output_dict = {'DP Base Sequence' : mtp['DP Base Sequence'], 'mistranslated sequence' : mtp['mistranslated sequence'], 'aa subs' : mtp['aa subs']}

        if label_designation == 'Label-Free':
            ratios = np.log2(mtp_df['Intensity'].values / bp_df['Intensity'].values)
        elif label_designation == 'TMT':
            ratios = np.log10(mtp_df['Intensity'].values / bp_df['Intensity'].values)

            mtp_reporters = mtp_df.filter(regex=r'^(?!.*Normalised).*Reporter intensity corrected.*$', axis=1)
            bp_reporters = bp_df.filter(regex=r'^(?!.*Normalised).*Reporter intensity corrected.*$', axis=1)

            mtp_reporters_norm = mtp_df.filter(regex='Normalised Reporter intensity corrected', axis=1)
            bp_reporters_norm = bp_df.filter(regex='Normalised Reporter intensity corrected', axis=1)

            mtp_ratios = mtp_reporters.div(mtp_reporters.sum(axis=1).values, axis=0)
            bp_ratios = bp_reporters.div(bp_reporters.sum(axis=1).values, axis=0)

            mtp_distributed = mtp_ratios.mul(mtp_df['Intensity'], axis=0)
            bp_distributed = bp_ratios.mul(bp_df['Intensity'], axis=0)

            for tmt_plex in sample_df['MQ'].dropna().unique():
                channel = str(int(tmt_plex))
                output_dict[f'MTP.Plex.{channel}'] = mtp_distributed[f'Reporter intensity corrected {channel}'].values
                output_dict[f'BP.Plex.{channel}'] = bp_distributed[f'Reporter intensity corrected {channel}'].values
                output_dict[f'Ratio.Plex.{channel}'] = np.log10(
                    mtp_distributed[f'Reporter intensity corrected {channel}'].values /
                    bp_distributed[f'Reporter intensity corrected {channel}'].values
                )
                output_dict[f'MTP.Plex.{channel}.Norm.Sum'] = mtp_reporters_norm[f'Normalised Reporter intensity corrected {channel}'].values
                output_dict[f'BP.Plex.{channel}.Norm.Sum'] = bp_reporters_norm[f'Normalised Reporter intensity corrected {channel}'].values

        output_dict['MTP.Sum'] = mtp_df['Intensity'].values
        output_dict['BP.Sum'] = bp_df['Intensity'].values
        output_dict['Ratio'] = ratios

        return pd.DataFrame.from_dict(output_dict)
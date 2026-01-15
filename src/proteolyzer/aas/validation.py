import pandas as pd
import numpy as np
import os
import pickle
from multiprocessing import Queue
from typing import Dict, Union
from pathlib import Path
from . import utils
from .config import Config

CONFIG = Config()

class Validation:
    def __init__(self, params: Union[str, Dict], queue: Queue = None):
        if isinstance(params, str):
            self.params = utils._load_params(params)
        elif isinstance(params, dict):
            self.params = params
        else:
            raise ValueError("params must be a file path or dict")
        self.queue = queue if queue is not None else utils.NullQueue()

        # Input Parameters
        self.workflow = self.params.get('Utils').get('Workflow')
        self.label_setup = self.params.get('Utils').get('Labelling Setup')
        self.label_plex = self.params.get('Utils').get('Label Plex')
        self.metadata = pd.read_excel(self.params.get('Utils').get('Metadata File'), engine = 'openpyxl', usecols=['TMT plex', 'TMT channel', 'ParticipantID', 'Group', 'MQ', 'sample_name', 'sample_ID'])
        self.metadata['MQ'] = self.metadata['TMT channel'].map(CONFIG.TMT.MQ_TMT_MAP)
        self.data_dir = self.params.get('Utils').get('Data Folder')
        self.output_dir = self.params.get('Utils').get('Output Folder')

        # Validation thresholds
        self.validation_pep = self.params.get('Validation').get('Validation PEP')
        self.pif = self.params.get('Validation').get('PIF')
        self.frag_evidence = self.params.get('Validation').get('Fragment Evidence')
        
        self.queue.put(('stdout', f"{self.__class__.__name__} initialized."))


    def run(self):
        for sample in np.unique(self.metadata['sample_ID']):
            self.process_sample(sample)

    def process_sample(self, sample):
        sample_dir = self._locate_sample_dir(sample)
        if not sample_dir:
            self.queue.put(('stderr', f"{sample} not found in data directory."))
            return

        self.queue.put(('stdout', f"Processing sample: {sample}"))

        evidence_path = sample_dir / 'evidence.parquet'
        msms_path = sample_dir / 'msms.parquet'
        mtp_path = self.output_dir / 'MTP' / f'{sample}_MTP_Filtered_Stage_1.p'

        if not (os.path.exists(evidence_path) and os.path.exists(msms_path) and os.path.exists(mtp_path)):
            self.queue.put(('stderr', f"Missing data files for sample {sample}"))
            return

        evidence = pd.read_parquet(evidence_path, engine='fastparquet')
        msms = pd.read_parquet(msms_path, engine='fastparquet')
        with open(mtp_path, 'rb') as f:
            mtp = pickle.load(f)

        validated_mtp, val_evidence = self.mtp_validate(evidence, msms, mtp)

        output_path = self.output_dir / 'MTP' / f'{sample}_MTP_Filtered_Stage_2.p'
        with open(output_path, 'wb') as f:
            pickle.dump(validated_mtp, f)
        val_evidence_path = self.output_dir / 'MTP' / f'{sample}_Val_Evidence_Filtered_Stage_2.p'
        with open(val_evidence_path, 'wb') as f:
            pickle.dump(val_evidence, f)

        self.queue.put(('stdout', f"Validation complete for sample: {sample}"))

    def _locate_sample_dir(self, sample):
        data_dir = Path(self.data_dir)
        for subdir in data_dir.rglob('*'):
            if subdir.is_dir() and subdir.name == sample + '_val':
                return subdir
        return None

    def frags_containing_aas(self, row: pd.Series):
        seq = row['mistranslated sequence']
        pos = row['mistranslated aas positions']

        b_ions_pos = pos + 1
        y_ions_pos = len(seq) - pos

        b_ions_aas = np.arange(b_ions_pos, len(seq), dtype=int)
        y_ions_aas = np.arange(y_ions_pos, len(seq), dtype=int)

        return b_ions_aas, y_ions_aas

    def frag_count(self, row: pd.Series, frag_ev_merge: pd.DataFrame):
        b_ions = row['b_ions_aas']
        y_ions = row['y_ions_aas']
        seq = row['mistranslated sequence']
        scan = row['MS/MS scan number']
        charge = row['Charge']

        matches = frag_ev_merge[
            (frag_ev_merge['Raw file'] == row['Raw file']) &
            (frag_ev_merge['Sequence'] == seq) &
            (frag_ev_merge['MS/MS scan number'].isin(scan)) &
            (frag_ev_merge['Charge'] == charge)
        ]

        b_matches = matches[(matches['Frag.Type'] == 'b') & (matches['Frag.Number'].isin(b_ions))].shape[0]
        y_matches = matches[(matches['Frag.Type'] == 'y') & (matches['Frag.Number'].isin(y_ions))].shape[0]

        return b_matches + y_matches

    def mtp_validate(self, val_evidence: pd.DataFrame, val_msms: pd.DataFrame, mtp: pd.DataFrame):
        val_evidence = val_evidence[(val_evidence['PEP'] <= self.validation_pep) & (val_evidence['PIF'] >= self.pif)].copy()

        if 'Reporter intensity corrected 1' in val_evidence.columns:
            norm_reporters = val_evidence.filter(regex='Reporter intensity corrected').apply(lambda x: x / x.median(), axis=0)
            norm_reporters.columns = ['Normalised ' + col for col in norm_reporters.columns]
            val_evidence = pd.concat([val_evidence, norm_reporters], axis=1)
        else:
            val_evidence.loc[:, 'Intensity.Normalised'] = val_evidence['Intensity'] / val_evidence['Intensity'].median()

        filter_list = mtp[['Raw file', 'Charge', 'DP Base Sequence', 'mistranslated sequence']]
        filter_list = filter_list.melt(id_vars=['Raw file', 'Charge'], var_name='Type', value_name='Sequence').drop('Type', axis=1).drop_duplicates()

        val_evidence = val_evidence.merge(filter_list, how='inner', on=['Raw file', 'Charge', 'Sequence'])
        val_evidence = val_evidence.drop(['m/z', 'Mass', 'Mass error [ppm]', 'Retention time', 'PIF', 'PEP'], axis=1, errors='ignore')

        msms_ev_merge = val_msms.merge(val_evidence, how='inner', on=['Raw file', 'MS/MS scan number'])

        mtp.loc[:, 'b_ions_aas'], mtp.loc[:, 'y_ions_aas'] = zip(*mtp.apply(lambda x: self.frags_containing_aas(x), axis=1))
        mtp.loc[:, 'fragment_evidence'] = mtp.apply(lambda x: self.frag_count(x, msms_ev_merge), axis=1)

        # Group and sum fragment evidence across same MTP (ignore scan/raw/charge-level variation)
        mtp_grouped = mtp.groupby(['DP Base Sequence', 'mistranslated sequence', 'aa subs'])['fragment_evidence'].sum().reset_index()
        mtp_grouped = mtp_grouped[mtp_grouped['fragment_evidence'] > self.frag_evidence]

        return mtp_grouped, msms_ev_merge
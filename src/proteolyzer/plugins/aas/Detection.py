import pandas as pd
import numpy as np
import re
import os
import shutil
import ahocorasick
import pickle
from multiprocessing import Queue
from . import Utils
from Utils import NullQueue

class Detection:
    def __init__(self, params: dict, queue : Queue = None):
        self.params = params
        self.queue = queue if queue is not None else NullQueue()

        # General Params
        self.configuration = params.get('Configuration')
        self.label_setup = params.get('Labelling Setup')
        self.label_plex = int(params.get('Label Plex'))
        self.subs_ref = Utils.subs_ref_gen(self.label_setup)
        self.protease = params.get('Protease')
        self.translated_frames = params.get('Translated Frames Folder')
        self.mods_dict = Utils.gen_mod_dict(params.get('Modifications File'))
        self.metadata = pd.read_excel(params.get('Metadata File'), engine = 'openpyxl', usecols = ['TMT plex', 'TMT channel', 'ParticipantID' , 'Group', 'MQ', 'sample_name', 'sample_ID'])
        self.metadata['MQ'] = self.metadata['TMT channel'].map(Utils.MQ_TMT_dict)
        self.prot_fasta = params.get('Protein FASTA')
        self.data_dir = params.get('Data Folder')
        self.output_dir = params.get('Output Folder')
        # Detection Specific Params
        self.detection_pep = float(params.get('Detection PEP'))
        self.aa_sub_ppm = float(params.get('AA Substitution ppm'))
        self.ptm_ppm = float(params.get('PTM ppm'))
        self.pos_prob = float(params.get('Positional Probability Threshold'))
        self.cn_term_prob = float(params.get('C/n-term Modification Threshold'))

        Utils.prepare_output_dir(self.output_dir)
        Utils.txt_to_parquet(self.data_dir, self.label_setup, self.detection_pep, self.label_plex)
        
        self.queue.put(("stdout", f"{self.__class__.__name__} initialized."))

    def run(self):
        for sample in np.unique(self.metadata['sample_ID']):
            self.process_sample(sample)

    def process_sample(self, sample):
        sample_dir = self._locate_sample_dir(sample)
        if not sample_dir:
            self.queue.put(('stderr', f"Sample {sample} not found in data directory."))
            return

        self.queue.put(('stdout', f"Processing sample: {sample}"))
        all_pep_path = os.path.join(sample_dir, 'allPeptides.parquet')
        evidence_path = os.path.join(sample_dir, 'evidence.parquet')

        all_peps = pd.read_parquet(all_pep_path, engine='fastparquet')
        n_term_peps = np.load(os.path.dirname(all_pep_path) + '/n_term_peps.npy', allow_pickle=True)
        c_term_peps = np.load(os.path.dirname(all_pep_path) + '/c_term_peps.npy', allow_pickle=True)
        sample_df = self.find_potential_aas(all_peps, n_term_peps, c_term_peps)
        with open(self.output_dir / f'{sample}_Detection_DF.p', 'wb') as f:
            pickle.dump(sample_df, f)
        Utils.ptm_mtp_output(sample_df, sample, self.output_dir)

        with open(f'{self.output_dir}/MTP/{sample}_MTP.p', 'rb') as f:
            mtp = pickle.load(f)

        evidence = pd.read_parquet(evidence_path, engine='fastparquet')
        mtp = self.validate_mtp(evidence, mtp)

        with open(f'{self.output_dir}/MTP/{sample}_MTP_Filtered_Stage_1.p', 'wb') as f:
            pickle.dump(mtp, f)

        self.write_fasta(mtp, sample)


    def _locate_sample_dir(self, sample):
        for root, dirs, _ in os.walk(self.data_dir):
            if sample in dirs:
                return os.path.join(root, sample)
        return None

    def refine_localization_probabilities(self, modified_seq):
        modified_sites = [modified_seq[m.start()-1] for m in re.finditer(r'\(', modified_seq)]
        weights = [float(i) for i in re.findall(r'\(([^)]+)\)', modified_seq)]
        return len(modified_sites), modified_sites, weights

    def pep_cterm(self, modified_sequence):
        if modified_sequence[-1] == ')':
            prob = float(modified_sequence[:-1].split('(')[-1])
            return prob >= self.cn_term_prob
        return False

    def pep_nterm(self, modified_sequence):
        if modified_sequence[1] == '(':
            prob = float(modified_sequence[2:].split(')')[0])
            return prob >= self.cn_term_prob
        return False

    def get_aa_subs(self, row):
        candidate = row['DP candidate residues']
        pos_probs = row['DP positional probabilities']
        DP_deltam = row['DP Mass Difference']
        mtol = row['m/z'] * (self.aa_sub_ppm / 1e6)

        res_dict = self.subs_ref.get(candidate, {})
        for s in res_dict:
            delta_m = res_dict[s]
            if ((DP_deltam > 0 and delta_m > 0) or (DP_deltam < 0 and delta_m < 0)) and \
               (DP_deltam > delta_m - mtol) and (DP_deltam < delta_m + mtol) and (self.pos_prob < pos_probs):
                row['aa subs'] = s
                row['aa subs positional probability'] = pos_probs
                row['aa subs mass error (ppm)'] = abs(DP_deltam - delta_m)
                row['destination aa'] = s[-1]
        return row

    def get_mistranslated_seq(self, row):
        bp = row['DP Base Sequence']
        seq = row['DP Probabilities']
        sub = row['destination aa']
        pos = row['DP Positions']
        parsed_seq = np.array(list(seq))
        mtp = parsed_seq.copy()
        mtp[pos] = sub
        mtp = ''.join(mtp)
        mtp = re.sub(r'\([^)]*\)', '', mtp)

        position = int([i for i, x in enumerate(bp) if mtp[i] != x][0])
        row['mistranslated sequence'] = mtp
        row['mistranslated aas positions'] = position
        return row

    def find_PTMs(self, row):
        res = row['DP candidate residues']
        prob = row['DP positional probabilities']
        DP_deltam = row['DP Mass Difference']
        mtol = row['m/z'] * (self.ptm_ppm / 1e6)

        if row['Peptide N-term']:
            res = 'N-term'
        elif row['Peptide C-term']:
            res = 'C-term'

        filtered_mods = self.mods_dict.get(res, [])
        for mod in filtered_mods:
            modification, pos, delta_m = mod
            if (DP_deltam > delta_m - mtol) and (DP_deltam < delta_m + mtol):
                term_filter = (
                    pos == 'Anywhere' or
                    (pos == 'Protein N-term' and row['Protein N-term']) or
                    (pos == 'Any N-term' and row['Peptide N-term']) or
                    (pos == 'Protein C-term' and row['Protein C-term']) or
                    (pos == 'Any C-term' and row['Peptide C-term'])
                )
                if term_filter:
                    row['PTM'] = modification
                    row['PTM site'] = res
                    row['PTM positional probability'] = prob
                    row['PTM mass error [observed-expected] (ppm)'] = DP_deltam - delta_m
        return row

    def aho_corasick_search(self, translated_genome, patterns):
        A = ahocorasick.Automaton()
        results = {pattern: [] for pattern in patterns}
        for pattern in patterns:
            A.add_word(pattern, pattern)
        A.make_automaton()
        for end_index, pattern in A.iter(translated_genome):
            results[pattern].append(end_index - len(pattern) + 1)
        return {p: (bool(v), v) for p, v in results.items()}

    def aho_corasick_output_organize(self, result):
        return [k[1:] for k, v in result.items() if v[0]]

    def find_homologous_peptide(self, sample_df):
        cleavage_sites = Utils.protease_cleavage_sites[self.protease]
        all_mtps = sample_df['mistranslated sequence'][sample_df['mistranslated sequence'].notna()].unique()
        all_mtps = [prefix + seq for seq in all_mtps for prefix in cleavage_sites]

        for frame in range(1, 7):
            with open(self.translated_frames / f'W{frame}_aa_ambig.p', 'rb') as wf:
                w_aa = pickle.load(wf)
            with open(self.translated_frames / f's{frame}a_ambig.p', 'rb') as sf:
                s_aa = pickle.load(sf)
            w_aa_out = self.aho_corasick_search(w_aa, all_mtps)
            s_aa_out = self.aho_corasick_search(s_aa, all_mtps)
            matched = self.aho_corasick_output_organize(w_aa_out) + self.aho_corasick_output_organize(s_aa_out)
            col = f"{frame}-frame genome substring"
            sample_df[col] = False
            sample_df.loc[sample_df['mistranslated sequence'].isin(matched), col] = True
        return sample_df

    def evidence_ppm(self, evidence):
        me = [x for x in evidence['Mass error [ppm]'] if x > -1000]
        return np.mean(me), np.std(me)

    def posterior_aasub_prob(self, row, mean, std):
        pp = row['aa subs positional probability']
        merr = row['aa subs mass error (ppm)']
        coeff = 1 / (std * np.sqrt(2 * np.pi))
        exponent = np.exp(-0.5 * ((merr - mean) / std) ** 2)
        pdf = coeff * exponent

        return pdf * pp

    def q_val_calc(self, mtp):
        post_prob = mtp['Posterior subs probability'].to_numpy()
        pval = 1 - post_prob
        ranked_pval = np.sort(pval)
        cumsum = np.cumsum(ranked_pval)
        ranked_qval = np.array([x / (i+1) for i, x in enumerate(cumsum)])
        qval = ranked_qval[np.argsort(np.argsort(pval))]
        return qval

    def gen_metrics(self, mtp):
        qval = mtp['q-value'].to_numpy()
        thresh = np.max(qval)
        n_thresh = len([x for x in qval if x <= thresh])
        TP_thresh = np.floor((1 - thresh) * n_thresh)
        ref_model = np.logspace(-20, -1, num=100, base=10)
        metric_arr = np.zeros((len(ref_model), 8))
        for i, qt in enumerate(ref_model):
            TP = len([x for x in qval if x <= qt]) * (1 - qt)
            FP = len([x for x in qval if x <= qt]) * qt
            FN = TP_thresh - TP
            TN = len(qval) - TP - FN - FP
            P = TP / (TP + FP) if (TP + FP) > 0 else 0
            R = TP / (TP + FN) if (TP + FN) > 0 else 0
            F_score = (2 * P * R) / (P + R) if (P + R) > 0 else 0
            metric_arr[i, :] = qt, TP, FP, FN, TN, P, R, F_score
        return pd.DataFrame(metric_arr, columns=['q_threshold', 'TP', 'FP', 'FN', 'TN', 'Precision', 'Recall', 'F_score'])

    def mtp_filter(self, mtp, metric_df):
        max_F_idx = metric_df['F_score'].idxmax()
        q_thresh = metric_df.loc[max_F_idx, 'q_threshold']
        filtered_mtp = mtp[mtp['q-value'] <= q_thresh].reset_index(drop=True)
        return filtered_mtp

    def validate_mtp(self, evidence, mtp):
        mean, std = self.evidence_ppm(evidence)
        mtp['Posterior subs probability'] = mtp.apply(lambda x: self.posterior_aasub_prob(x, mean, std), axis=1)
        pp_sum = mtp.groupby(['Raw file', 'DP Base Sequence', 'Charge', 'DP PEP'])['Posterior subs probability'].sum().reset_index()
        pp_sum = pp_sum.rename({'Posterior subs probability': 'Posterior subs probability sum'}, axis=1)
        mtp = mtp.merge(pp_sum, on=['Raw file', 'DP Base Sequence', 'Charge', 'DP PEP'], how='left')
        mtp['Posterior subs probability'] = mtp['Posterior subs probability'] / mtp['Posterior subs probability sum'] * (1 - mtp['DP PEP'])
        mtp['q-value'] = self.q_val_calc(mtp)
        metrics = self.gen_metrics(mtp)
        return self.mtp_filter(mtp, metrics)

    def write_fasta(self, filtered_mtp, sample_name):
        # base_scans = filtered_mtp['DP Base Scan Number']
        # dp_scans = filtered_mtp['DP Mod Scan Number']
        base_seqs = filtered_mtp['DP Base Sequence']
        prots = filtered_mtp['Leading.Razor.DP.Protein']
        seqs = filtered_mtp['mistranslated sequence']
        aa_subs = filtered_mtp['aa subs'].str.replace(' to ', ':')
        sub_positions = filtered_mtp['mistranslated aas positions'].astype(int)
        output_fasta_path = self.output_dir / f'{sample_name}_validation.fasta'
        shutil.copy(self.prot_fasta, output_fasta_path)
        with open(output_fasta_path, 'a') as f:
            for i, seq in enumerate(seqs):
                # header = f">MTP|{i}__base{int(base_scans[i])}_DP{int(dp_scans[i])}"
                header = f">MTP|({base_seqs[i]})({sub_positions[i]})({aa_subs[i]})({prots[i]})"
                f.write(f"{header}\n{seq}\n")
        self.queue.put(('stdout', ('Validation Fasta Written.')))

    def find_potential_aas(self, dp_df, n_term_peps, c_term_peps):
        
        dp_df['Leading.Razor.DP.Protein'] = dp_df['DP Proteins'].str.split(';').str[0]
        dp_df['Con.DP.Protein'] = dp_df['DP Proteins'].str.contains('CON__', na = False)

        dp_df['DP Positions'] = dp_df.apply(lambda x: [i - 1 for i, a in enumerate(x['DP Probabilities']) if a == '('], axis=1)
        dp_df['count candidate residues per peptide'], dp_df['DP candidate residues'], dp_df['DP positional probabilities'] = zip(*dp_df.apply(lambda x: self.refine_localization_probabilities(x['DP Probabilities']), axis=1))
        dp_df['Protein N-term'] = dp_df.apply(lambda x: x['DP Base Sequence'] in n_term_peps, axis=1)
        dp_df['Protein C-term'] = dp_df.apply(lambda x: x['DP Base Sequence'] in c_term_peps, axis=1)
        dp_df['Peptide N-term'] = dp_df['DP Probabilities'].apply(self.pep_nterm)
        dp_df['Peptide C-term'] = dp_df['DP Probabilities'].apply(self.pep_cterm)
        dp_df = dp_df.explode(['DP candidate residues', 'DP positional probabilities', 'DP Positions'])
        self.queue.put(('progress', (1, 4)))

        dp_df[['aa subs', 'aa subs positional probability', 'aa subs mass error (ppm)', 'destination aa']] = None
        dp_df = dp_df.apply(lambda x: self.get_aa_subs(x), axis=1)
        dp_df[['mistranslated sequence', 'mistranslated aas positions']] = None
        dp_df = dp_df.apply(lambda x: self.get_mistranslated_seq(x) if x['aa subs'] else x, axis=1)
        self.queue.put(('progress', (2, 4)))

        dp_df[['PTM', 'PTM site', 'PTM positional probability', 'PTM mass error [observed-expected] (ppm)']] = None
        dp_df = dp_df.apply(lambda x: self.find_PTMs(x), axis=1)
        self.queue.put(('progress', (3, 4)))

        dp_df = self.find_homologous_peptide(dp_df)
        self.queue.put(('progress', (4, 4)))
        
        dp_df = dp_df[dp_df['Con.DP.Protein'] == False]
        dp_df = dp_df[dp_df['Reverse'].isna()]

        return dp_df
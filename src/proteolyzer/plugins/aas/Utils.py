import os
import pandas as pd
import numpy as np
import pickle
from pathlib import Path
from Config import FILES_NEEDED, COLSTOKEEP, PEP_THRESHOLD, AA_MW, TMT_MASS_SHIFT


class NullQueue:
    """A 'do-nothing' queue to replace the real Queue when multiprocessing isn't used."""
    def put(self, item):
        pass


def prepare_output_dir(output_dir : Path):
    parent_dir = output_dir.parent
    if not os.path.exists(parent_dir / 'Translation_Frames'):
        os.makedirs(parent_dir / 'Translation_Frames')
    if not os.path.exists(output_dir / 'PTM'):
        os.makedirs(output_dir / 'PTM')
    if not os.path.exists(output_dir / 'MTP'):
        os.makedirs(output_dir / 'MTP')


def column_mapping(df : pd.DataFrame, cols2keep : list):

    col_lookup = {col.lower(): col for col in df.columns}
    rename_map = {col_lookup[col.lower()]: col for col in cols2keep if col.lower() in col_lookup}
    df = df[list(rename_map.keys())].rename(columns=rename_map)
    df = df[[col for col in cols2keep if col in df.columns]]

    return df


## TEMPORARY ##

def txt_to_parquet(
    data_folder: str,
    label_designation: str,
    tmt_plex: int
):
    for root, dirs, _ in os.walk(data_folder):
        for dir in dirs:
            subdir = os.path.join(root, dir)
            files = os.listdir(subdir)

            if 'dependentPeptides.txt' in files:
                search_type = 'Detection'
            else:
                search_type = 'Validation'
            files_of_interest = FILES_NEEDED['MaxQuant'][search_type]

            for name in files:
                if name not in files_of_interest:
                    continue

                txt_path = os.path.join(subdir, name)
                parquet_path = os.path.join(subdir, os.path.splitext(name)[0] + '.parquet')

                if name == 'allPeptides.txt':
                    cols2keep = COLSTOKEEP['MaxQuant'][name].copy()
                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)
                    df = allpeptides_process(df)

                elif name == 'evidence.txt':
                    cols2keep = COLSTOKEEP['MaxQuant'][name].copy()
                    if label_designation == 'TMT' and '_val' in dir:
                        cols2keep += [f'Reporter intensity corrected {i}' for i in range(1, tmt_plex + 1)]

                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)
                    df = evidence_process(df)

                elif name == 'msms.txt':
                    cols2keep = COLSTOKEEP['MaxQuant'][name].copy()
                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)
                    df = msms_process(df)

                elif name == 'peptides.txt':
                    cols2keep = COLSTOKEEP['MaxQuant'][name].copy()
                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)
                    df = peptides_process(df)

                    n_term_peps = df[df['Amino acid after'] != '-']['Sequence'].to_numpy()
                    c_term_peps = df[df['Amino acid after'] == '-']['Sequence'].to_numpy()

                    np.save(os.path.join(subdir, 'n_term_peps.npy'), n_term_peps, allow_pickle=True)
                    np.save(os.path.join(subdir, 'c_term_peps.npy'), c_term_peps, allow_pickle=True)
                    print(f"Saved N-Terminus and C-Terminus Peptides to {subdir}")
                    continue

                df.reset_index(drop=True, inplace=True)
                df.to_parquet(parquet_path, engine='fastparquet')
                print(f"Converted {txt_path} to {parquet_path}")



def evidence_process(
                    evidence : pd.DataFrame
                    ):
    evidence = evidence[(evidence['Reverse'].isna()) & (evidence['Potential contaminant'].isna())]
    return evidence

def allpeptides_process(
                        allpeptides : pd.DataFrame
                        ):
    
    allpeptides['Leading.Razor.DP.Protein'] = allpeptides['DP Proteins'].str.split(';').str[0]
    allpeptides['Con.DP.Protein'] = allpeptides['DP Proteins'].str.contains('CON__', na = False)
    allpeptides = allpeptides[(allpeptides['Reverse'].isna()) & (allpeptides['Con.DP.Protein'] == False)]
    allpeptides = allpeptides.loc[~np.isnan(allpeptides['DP Mass Difference']) & (allpeptides['DP PEP'] <= PEP_THRESHOLD)]
    allpeptides = allpeptides[~allpeptides['MSMS Scan Numbers'].isna()]
    allpeptides['MSMS Scan Numbers'] = allpeptides['MSMS Scan Numbers'].str.split(';').apply(lambda x: list(map(int, x)))
    allpeptides.rename(columns={'MSMS Scan Numbers': 'MS/MS scan number'}, inplace=True)

    return allpeptides


def msms_process(
                msms : pd.DataFrame
                ):
    
    msms = msms[msms['Reverse'].isna()]
    msms = msms.rename({'Scan number' : 'MS/MS scan number'}, axis = 1)
    msms['Matches'] = msms['Matches'].str.split(';')
    msms = msms[msms['Matches'].notna()]
    msms['Matches'] = msms.apply(lambda x: [frag for frag in x['Matches'] if 'a' not in frag and '(' not in frag and 'NH3' not in frag and 'H2O' not in frag], axis = 1)
    msms['Frag.Type'] = msms.apply(lambda x: [frag[0] for frag in x['Matches']], axis = 1)
    msms['Frag.Number'] = msms.apply(lambda x: [int(frag[1:]) for frag in x['Matches']], axis = 1)
    msms = msms.drop('Matches', axis = 1).explode(['Frag.Type', 'Frag.Number'])

    return msms


def peptides_process(peptides : pd.DataFrame):
    peptides = peptides[(peptides['Reverse'].isna()) & (peptides['Potential contaminant'].isna())]

    peptides = peptides[
        (peptides['Start position'] == 1) |
        (peptides['Amino acid after'] == '-')
    ]
    return peptides


def subs_ref_gen(
                 label_designation : str
                 ):

    aa_masses = AA_MW.copy()
    if label_designation == 'TMT':
        aa_masses['K'] += TMT_MASS_SHIFT[11]  # assuming TMT11 for mass shift

    # list of amino acids
    AAs = 'ACDEFGHIKLMNPQRSTVWY'

    ## dictionary of AAS types and the theoretical mass shift of the substitution
    subs_dict = { i+' to '+j : aa_masses[j] - aa_masses[i] for i in aa_masses for j in aa_masses if i!=j}

    ## update subs_dict to account for the fact that I and L have the same mass
    del subs_dict['L to I']
    del subs_dict['I to L']
    subst_dict={}
    for k,v in subs_dict.items(): # unifies I and L
        if k[-1]!='I' and k[-1]!='L':
            subst_dict[k] = v
        elif k[-1]=='I':
            subst_dict[k+'/L']=v
            
    subs_dict = {}
    for a in AAs:
        subs_dict[a] = {k:v for k,v in subst_dict.items() if k[0]==a}

    return subs_dict


def ptm_mtp_output(
                  dp_df : pd.DataFrame,
                  sample_name : str,
                  output_dir : Path,
                  ):
    
    ptm_df = dp_df.copy()
    ptm_df = ptm_df[ptm_df['PTM'].notna()] 
    with open(output_dir / 'PTM' / f'{sample_name}_PTM.p', 'wb') as f:
        pickle.dump(ptm_df, f)

    # MTP are peptides with potential AAS that cannot be explained by PTM and do not have homologous sequence in genome
    mtp_df = dp_df.copy()
    mtp_df = mtp_df[mtp_df['PTM'].isna()] 
    mtp_df = mtp_df[~mtp_df['mistranslated sequence'].isnull()]
    mtp_df = mtp_df[~mtp_df.apply(lambda x: any([x['1-frame genome substring'], x['2-frame genome substring'], x['3-frame genome substring'],
                                                x['4-frame genome substring'], x['5-frame genome substring'], x['6-frame genome substring']]), axis = 1)]
    with open(output_dir / 'MTP' / f'{sample_name}_MTP.p', 'wb') as f:
        pickle.dump(mtp_df, f)


def gen_mod_dict(
                mod_path : str
                ):

    mod_df = pd.read_csv(mod_path, sep = ',', usecols = ['modification', 'position', 'site', 'mass shift'])
    all_mods_dict = {}
    for k, g in mod_df.groupby('site'):
        all_mods_dict[k] = g.drop('site', axis = 1).values.tolist()

    return all_mods_dict
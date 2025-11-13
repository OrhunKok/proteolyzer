import os
import pandas as pd
import numpy as np
import pickle
from pathlib import Path

codon_map = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

tmt_plex_map = {'126' : '1', '127N' : '2', '127C' : '3', '128N' : '4', '128C' : '5', '129N' : '6',
                '129C' : '7', '130N' : '8', '130C' : '9', '131N' : '10', '131C' : '11', '132N' : '12',
                '132C' : '13', '133N' : '14', '133C' : '15', '134N' : '16', '134C' : '17', '135N' : '18'}

MQ_TMT_dict = {'126' : 1,'127N' : 2, '127C' : 3, '128N' : 4, 
               '128C' : 5, '129N' : 6, '129C' : 7, '130N' : 8, 
               '130C' : 9, '131' : 10}

protease_cleavage_sites = {'Trypsin' : ['K', 'R', '*'], 'Lys-C' : ['K', '*'], 'Arg-C' : ['R', '*']}


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
import os
import numpy as np
import pandas as pd

def txt_to_parquet(
    data_folder: str,
    label_designation: str,
    pep: float,
    tmt_plex: int
):
    for root, dirs, _ in os.walk(data_folder):
        for dir in dirs:
            subdir = os.path.join(root, dir)
            files = os.listdir(subdir)

            if 'dependentPeptides.txt' in files:
                files_of_interest = ['evidence.txt', 'allPeptides.txt', 'peptides.txt']
            else:
                files_of_interest = ['evidence.txt', 'msms.txt']

            for name in files:
                if name not in files_of_interest:
                    continue

                txt_path = os.path.join(subdir, name)
                parquet_path = os.path.join(subdir, os.path.splitext(name)[0] + '.parquet')

                if name == 'allPeptides.txt':
                    cols2keep = [
                        'Raw file', 'Charge', 'm/z', 'Mass', 'Retention time',
                        'Intensity', 'DP Mass Difference', 'DP PEP', 'DP Decoy', 'DP Proteins',
                        'Reverse', 'DP Base Sequence', 'DP Probabilities', 'DP Positional Probability',
                        'DP Base Scan Number', 'DP Mod Scan Number', 'MSMS Scan Numbers'
                    ]
                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)

                    df['Leading.Razor.DP.Protein'] = df['DP Proteins'].str.split(';').str[0]
                    df['Con.DP.Protein'] = df['DP Proteins'].str.contains('CON__', na = False)
                    df = df[(df['Reverse'].isna()) & (df['Con.DP.Protein'] == False)]
                    df = df.loc[~np.isnan(df['DP Mass Difference']) & (df['DP PEP'] <= pep)]
                    df = df[~df['MSMS Scan Numbers'].isna()]
                    df['MSMS Scan Numbers'] = df['MSMS Scan Numbers'].str.split(';').apply(lambda x: list(map(int, x)))
                    df.rename(columns={'MSMS Scan Numbers': 'MS/MS scan number'}, inplace=True)

                elif name == 'evidence.txt':
                    cols2keep = [
                        'Raw file', 'Charge', 'm/z', 'Mass', 'Retention time', 'Reverse', 'Potential contaminant',
                        'Sequence', 'PIF', 'PEP', 'Mass error [ppm]', 'MS/MS scan number', 'Intensity'
                    ]
                    if label_designation == 'TMT' and '_val' in dir:
                        cols2keep += [f'Reporter intensity corrected {i}' for i in range(1, tmt_plex + 1)]

                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)

                    df = df[(df['Reverse'].isna()) & (df['Potential contaminant'].isna())]

                elif name == 'msms.txt':
                    cols2keep = ['Raw file', 'Scan number', 'Matches', 'Reverse']
                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)
                    df = msms_process(df)

                elif name == 'peptides.txt':
                    cols2keep = ['Sequence', 'Start position', 'Amino acid after', 'Amino acid before',
                                 'Reverse', 'Potential contaminant']
                    df = pd.read_csv(txt_path, delimiter='\t')
                    df = column_mapping(df, cols2keep)

                    df = df[(df['Reverse'].isna()) & (df['Potential contaminant'].isna())]

                    df = df[
                        (df['Start position'] == 1) |
                        (df['Amino acid after'] == '-')
                    ]

                    n_term_peps = df[df['Amino acid after'] != '-']['Sequence'].to_numpy()
                    c_term_peps = df[df['Amino acid after'] == '-']['Sequence'].to_numpy()

                    np.save(os.path.join(subdir, 'n_term_peps.npy'), n_term_peps, allow_pickle=True)
                    np.save(os.path.join(subdir, 'c_term_peps.npy'), c_term_peps, allow_pickle=True)
                    print(f"Saved N-Terminus and C-Terminus Peptides to {subdir}")
                    continue

                df.reset_index(drop=True, inplace=True)
                df.to_parquet(parquet_path, engine='fastparquet')
                print(f"Converted {txt_path} to {parquet_path}")



# def txt_to_parquet(
#                   data_folder : str,
#                   label_designation : str,
#                   pep : float,
#                   tmt_plex : int
#                   ):

#     for root, dirs, files in os.walk(data_folder):

#         if 'dependentPeptides.txt' in files:
#             files_of_interest = ['evidence.txt', 'allPeptides.txt', 'peptides.txt']

#         else:
#             files_of_interest = ['evidence.txt', 'msms.txt']

#         for name in files:
#             if name in files_of_interest:

#                 txt_path = os.path.join(root, name)
#                 parquet_path = os.path.join(root, os.path.splitext(name)[0] + '.parquet')

#                 if name == 'allPeptides.txt':

#                     cols2keep = ['Raw file', 'Charge', 'm/z', 'Mass', 'Retention time',
#                                 'Intensity', 'DP Mass Difference', 'DP PEP',
#                                 'DP Base Sequence', 'DP Probabilities', 'DP Positional Probability', 
#                                 'DP Base Scan Number', 'DP Mod Scan Number', 'MSMS Scan Numbers']
                    
#                     df = pd.read_csv(txt_path, delimiter = '\t', usecols = cols2keep)

#                     df = df.loc[~np.isnan(df['DP Mass Difference']) & (df['DP PEP'] <= pep), :]
#                     df = df[~df['MSMS Scan Numbers'].isna()]
#                     scans = df['MSMS Scan Numbers'].tolist()
#                     scans = [list(map(int, x.split(';'))) for x in scans]
#                     df['MSMS Scan Numbers'] = scans
#                     df = df.rename({'MSMS Scan Numbers' : 'MS/MS scan number'}, axis = 1)

#                 elif name == 'evidence.txt':

#                     cols2keep = ['Raw file', 'Charge', 'm/z', 'Mass', 'Retention time', 'Sequence', 'PIF', 'PEP', 'Mass error [ppm]', 'MS/MS scan number', 'Intensity']

#                     # Retrieving the necessary quantification columns
#                     if label_designation == 'TMT':
#                         cols2keep = cols2keep + ['Reporter intensity corrected ' + str(plex) for plex in range(1, tmt_plex + 1)]

#                     try:
#                         df = pd.read_csv(txt_path, delimiter = '\t', usecols = cols2keep)
#                     except ValueError as e:
#                         e.add_note('TMT reporter columns not found! Ensure correct data is being used, correct search settings are used, and correct tmt_plex is input.')
#                         raise

                    
#                 elif name == 'msms.txt':

#                     cols2keep = ['Raw file', 'Scan number', 'Matches']
#                     df = pd.read_csv(txt_path, delimiter = '\t', usecols = cols2keep)
#                     df = msms_process(df)

#                 elif name == 'peptides.txt':

#                     cols2keep = ['Sequence', 'Start position', 'Amino acid after', 'Amino acid before']
#                     df = pd.read_csv(txt_path, delimiter = '\t', usecols = cols2keep)

#                     # Also counting methionine as starting position

                      ### REVISIT ###
#                     df = df[(df['Start position'] == 1) | (df['Amino acid after'] == '-') | ((df['Start position'] == 2) & (df['Amino acid before'] == 'M'))]
                      ### REVISIT ###

#                     n_term_peps = df[df['Amino acid after'] != '-']['Sequence'].to_numpy()
#                     c_term_peps = df[df['Amino acid after'] == '-']['Sequence'].to_numpy()

#                     np.save(os.path.join(root, 'n_term_peps.npy'), n_term_peps, allow_pickle = True)
#                     np.save(os.path.join(root, 'c_term_peps.npy'), c_term_peps, allow_pickle = True)
#                     print(f"Saved N-Terminus and C-Terminus Peptides to {root}")
#                     continue

#                 df.index = range(0, len(df))
#                 df.to_parquet(parquet_path, engine = 'fastparquet')
                
#                 print(f"Converted {txt_path} to {parquet_path}")
                

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


def subs_ref_gen(
                 label_designation : str
                 ):

    if label_designation == 'Label Free':
        MW_dict = {
                    "G": 57.02147, "A" : 71.03712, "S" : 87.03203, "P" : 97.05277, "V" : 99.06842, "T" : 101.04768,
                    "I" : 113.08407, "L" : 113.08407, "N" : 114.04293, "D" : 115.02695, "Q" : 128.05858, "K" : 128.09497,
                    "E" : 129.0426,"M" : 131.04049, "H" : 137.05891, "F" : 147.06842, "R" : 156.10112, "C" : 160.030654, 
                    "Y" : 163.0633, "W" : 186.07932,
                  }

    elif label_designation == 'TMT':
        # Mass of lysine is adjusted to reflect TMT label
        MW_dict = {
                    "G": 57.02147, "A" : 71.03712, "S" : 87.03203, "P" : 97.05277, "V" : 99.06842, "T" : 101.04768,
                    "I" : 113.08407, "L" : 113.08407, "N" : 114.04293, "D" : 115.02695, "Q" : 128.05858, "K" : 357.257902, 
                    "E" : 129.0426,"M" : 131.04049, "H" : 137.05891, "F" : 147.06842, "R" : 156.10112, "C" : 160.030654, 
                    "Y" : 163.0633, "W" : 186.07932,
                  }

    # list of amino acids
    AAs = 'ACDEFGHIKLMNPQRSTVWY'

    ## dictionary of AAS types and the theoretical mass shift of the substitution
    subs_dict = { i+' to '+j : MW_dict[j] - MW_dict[i] for i in MW_dict for j in MW_dict if i!=j}

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
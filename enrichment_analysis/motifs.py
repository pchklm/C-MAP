import pandas as pd
import numpy as np
from Bio import motifs
from collections import defaultdict
from .constants import amino_acids, alphabet, alphabet_index, site_columns, site_columns_index, base_enzymes, base_enzyme_codes, base_enzyme_codes_without_P

def calculate_pssms(counts_by_code, background):
    '''
    Calculate the position specific scoring matrices for all enzyme candidates.

    args:
        counts_by_code: Dictionary containing all observed cleavages for all enzyme candidates.
        background: Dictionary with the total count of each amino acid.

    returns:
        pssms: List of all position specific scoring matrices for all enzyme candidates.
    '''

    pssms = defaultdict(list)

    for code in counts_by_code:

        site_counts = counts_by_code[code]
        counts_dict = {aa: list(site_counts[aa]) for aa in site_counts.columns}
        m = motifs.Motif(counts=counts_dict, alphabet=alphabet)
        m.background = background
        m.pseudocounts = 1
        pssm = m.pssm

        pssm_array = np.array([[pssm[aa][i] for aa in alphabet] + [0] for i in range(len(site_columns))])

        pssms[code] = pssm_array

    return pssms


def pssm_to_regex(pssms):
    '''
    Create a regex patterns from position specific scoring matrices.

    args:
        pssms: List of all position specific scoring matrices for all enzyme candidates.

    returns:
        regexes: Dictionary of regex patterns for all enzyme candidates.   
    '''

    regexes = defaultdict(list)

    for code, pssm in pssms.items():
        if code in base_enzyme_codes:
            regex = base_enzymes[code]["regex"]
            regexes[code] = regex
            continue
        regex=[]
        for i in site_columns_index:
            enriched_aa_list = [alphabet[j] for j in alphabet_index if pssm[i][j] > 1.68]

            if len(enriched_aa_list) == 0:
                depleted_aa_list = [("!"+alphabet[j]) for j in alphabet_index if pssm[i][j] < -1.68]
                if len(enriched_aa_list):
                    regex.append(depleted_aa_list)
                else:
                    regex.append(["X"])
            else:
                regex.append(enriched_aa_list)

        regexes[code] = regex
    
    return regexes


def analyze_enzymes(enzyme_df, background):
    '''
    Analyze all candidate enzymes, calculate position specific scoring matrices and create regex patterns

    args:
        enzyme_df: Pandas dataframe containing all enzyme candidates along with their observed cleavages.
        background: Dictionary with the total count of each amino acid.

    returns:
        pssms: List of all position specific scoring matrices for all enzyme candidates.
        regexes: Dictionary of regex patterns for all enzyme candidates.  
        code_to_name: Dicionary to map enzyme code to their real name.
    '''

    counts_by_code = defaultdict(lambda: pd.DataFrame(0, index=site_columns, columns=amino_acids))
    code_to_name = defaultdict(str)

    for _, row in enzyme_df.iterrows():

        code = row["code"]
        code_to_name[row["code"]] = row["enzyme_name"]

        for pos in site_columns:
            for aa in amino_acids:
                col_name = f"{pos}_{aa}"
                if col_name in row and pd.notna(row[col_name]):
                    counts_by_code[code].at[pos, aa] = row[col_name]

        if code in base_enzyme_codes_without_P:
            for aa in amino_acids:
                counts_by_code[code].at["Site_P1prime",aa] = background[aa]

    pssms = calculate_pssms(counts_by_code,background)
    regexes = pssm_to_regex(pssms)
    
    return pssms, regexes, code_to_name
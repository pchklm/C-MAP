import pandas as pd
import numpy as np
import math
from .constants import alphabet, aa_to_idx, site_columns_index, alphabet_index, alphabet_with_X
from .helper import normalize_background
from scipy.stats import norm
from collections import defaultdict
import time

def match_enzymes(df, trie, pssms, code_to_name, background):
    '''
    Match enzymes with observed cleavage while also calculating a p_value for each match.

    args:
        df: Pandas dataframe containing all observed cleavages along with their matched protein and metadata.
        trie: Search tree containing the regex patterns for all candidate enzymes.
        pssms: Dictionary containing the position specific scoring matrices for all candidate enzymes.
        code_to_name: Dicionary to map enzyme code to their real name.
        background: Dictionary with the total count of each amino acid.

    returns:
        Pandas dataframe containing all information for each cleavage.
    '''

    cleavage_sites, proteinIDs, enzymes, positions, p_values, samples = [], [], [], [], [], []
    
    df = df.reset_index(drop=True)

    mus, sigmas = precalculate_expected_p_values(pssms, background)

    n_windows = encode_window(df["n_term_cleavage_window"])
    c_windows = encode_window(df["c_term_cleavage_window"])

    n_codes, n_p_values = find_best_matches(n_windows, trie, pssms, mus, sigmas)
    c_codes, c_p_values = find_best_matches(c_windows, trie, pssms, mus, sigmas)

    for i, row in df.iterrows():

        if n_codes[i] is not None:
            cleavage_sites.append(row.n_term_cleavage_window)
            proteinIDs.append(row.proteinID)
            enzymes.append(code_to_name.get(n_codes[i], "unspecified cleavage"))
            positions.append(row.n_term_position)
            p_values.append(n_p_values[i])
            samples.append(row.Sample)

        if c_codes[i] is not None:
            cleavage_sites.append(row.c_term_cleavage_window)
            proteinIDs.append(row.proteinID)
            enzymes.append(code_to_name.get(c_codes[i], "unspecified cleavage"))
            positions.append(row.c_term_position)
            p_values.append(c_p_values[i])
            samples.append(row.Sample)
        
    result = pd.DataFrame({
        "cleavage_site": cleavage_sites,
        "protein_id": proteinIDs,
        "enzyme": enzymes,
        "position": positions,
        "p_value": p_values,
        "sample": samples
    })

    return result

def encode_window(series):
    return series.apply(lambda s: np.array([aa_to_idx.get(a, 20) for a in s], dtype=np.int32)).to_numpy()


def find_best_matches(windows, trie, pssms, mus, sigmas):
    all_codes = []
    all_pvals = []
    #start = time.perf_counter()
    for i, window in enumerate(windows):
        # if (i % 1000 == 0):
        #     current = time.perf_counter()
        #     print(f"time for last 1k cleavages: {current - start:.4f} seconds", "total:", i)
        #     start = current
        best_score = float("-inf")
        best_match = None

        candidates = trie.match("".join(alphabet_with_X[j] for j in window))

        for c in candidates:
            score = calculate_pssm_score(pssms[c], window)
            if score > best_score:
                best_score = score
                best_match = c
        if best_match is not None:
            all_codes.append(best_match)
            all_pvals.append(calculate_p_value(best_score, mus[best_match], sigmas[best_match]))
        else:
            all_codes.append("unspecified cleavage")
            all_pvals.append(None)

    return all_codes, all_pvals

def calculate_pssm_score(pssm, window):
    '''
    Calculate the pssm-score for a match.

    args:
        pssm: Position specific scoring matrix for the match.
        cleavage_site: String of the cleaved amino acid sequence.

    returns:
        score: Number indicating how good a match is.
    '''

    return pssm[site_columns_index, window].sum()

def precalculate_expected_p_values(pssms, background):

    mus = defaultdict(int)
    sigmas = defaultdict(int)

    relative_bg = normalize_background(background)
    bg_array = np.array([relative_bg[aa] for aa in alphabet], dtype=float)

    for code in pssms:
        pssm = pssms[code]
        score_matrix = np.array([[pssm[s][i] for i in alphabet_index] for s in site_columns_index])

        # Per-site expected score under background
        exp_per_site = score_matrix.dot(bg_array)

        # Per-site second moment -> E[s^2] = sum p(a) * s(a)^2
        e2_per_site = (score_matrix**2).dot(bg_array)
        var_per_site = e2_per_site - exp_per_site**2

        mu = exp_per_site.sum()
        sigma2 = var_per_site.sum()
        sigma = math.sqrt(sigma2) if sigma2 > 0 else 0.0

        mus[code] = mu
        sigmas[code] = sigma

    return mus, sigmas


def calculate_p_value(score, mu, sigma):
    '''
    Calculate the p_value for a match.

    args:
        pssm: Position specific scoring matrix for the match
        cleavage_site: String of the cleaved amino acid sequence.
        background_frequency: Dictionary containing background frequencies representing the probability of each amino acid at a random position.
    
    returns:
        p_value: Number for a match between 0 and 1 indicating how statistically significant the match is.
    '''

    # handle degenerate sigma
    if sigma == 0:
        # If sigma 0, the score doesn't vary under null; then p is 0 or 1
        p = 0.0 if score > mu else 1.0
        return p

    z = (score - mu) / sigma
    p_value = 1.0 - norm.cdf(z)
    return p_value
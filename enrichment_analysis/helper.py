import pandas as pd
from .constants import three_to_one

def search_function(input: str, list: list):
    if input == None:
        return list
    return [item for item in list if input in str(item)]

def convert_3to1(aa3: str):
    if aa3 is None:
        return "X"
    aa3 = aa3.capitalize()
    return three_to_one.get(aa3, "X")

def normalize_background(bg_counts: dict):
    total = sum(bg_counts.values())
    if total == 0:
        raise ValueError("Background counts sum to 0")
    bg_probs = {aa: count / total for aa, count in bg_counts.items()}
    return bg_probs

def counts_to_relative_motif(counts):
    '''
    Transform absolute counts of each amino acid for each position into a relative motif

    args:
        counts: Array of a dict containing absolute counts for each amino acid for a position

    returns:
        pd.Dataframe: Pandas dataframe with the relative frequency of each amino acid per site
    '''

    all_aas = set().union(*[d.keys() for d in counts])
    rows = []

    for position in counts:
        total = sum(position.values()) 
        row = {aa: (position[aa] / total if total > 0 else 0.0) for aa in all_aas}
        rows.append(row)
    
    return pd.DataFrame(rows, index=[-4,-3,-2,-1,1,2,3,4])
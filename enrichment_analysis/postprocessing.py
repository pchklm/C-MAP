import pandas as pd
from collections import defaultdict
from .helper import counts_to_relative_motif

def filter_results(results, proteinIDs = [], metadata_filter = []):
    '''
    Accumulate results for filter settings.

    args:
        results: Pandas dataframe containing all information for each cleavage.
        proteinIDs: A list of Protein IDs.
        metadata_filter: A List of Metadata filter strings.

    returns:
        Dictionary containing the wanted output data for the top k enzymes.
    '''
     
    mask = pd.Series(True, index=results.index)

    if len(proteinIDs) > 0:
        submask = results["protein_id"].apply(lambda x: any(proteinID in x for proteinID in proteinIDs))
        mask &= submask

    if len(metadata_filter) > 0:
        submask = results["sample"].apply(
            lambda samples: any(
                mf in sample for sample in samples for mf in metadata_filter)
            )
        mask &= submask              

    filtered_results = results[mask]

    return filtered_results


def group_theoretical_cleavages(df, proteinIDs):

    filtered_df = filter_results(df, proteinIDs)

    enzyme_summary = {}

    for enzyme, group in filtered_df.groupby("enzyme"):

        position_dicts = [defaultdict(int) for _ in range(8)]


        for cleavage_site in group["n_term_cleavage_window"]:
            if cleavage_site is not None:
                for i, aa in enumerate(cleavage_site):
                    position_dicts[i][aa] += 1

        for cleavage_site in group["c_term_cleavage_window"]:
            if cleavage_site is not None:
                for i, aa in enumerate(cleavage_site):
                    position_dicts[i][aa] += 1
        
        motif = counts_to_relative_motif(position_dicts)

        enzyme_summary[enzyme] = {
            "motif": motif,
        }

    return enzyme_summary



def group_results(df, proteinIDs, metadata_filter, k):
    '''
    Group enzymes and calculate their wanted output data.

    args:
        df: Filtered Pandas dataframe containing all information for each cleavage

    returns:
        Dictionary containing the wanted output data for the top k enzymes.
    '''

    filtered_df = filter_results(df, proteinIDs, metadata_filter)

    enzyme_counts = filtered_df["enzyme"].value_counts()
    enzyme_summary = {}

    if k is not None:
        top_enzymes = set(enzyme_counts.nlargest(k).index)
        filtered_df = filtered_df[filtered_df["enzyme"].isin(top_enzymes)]

    for enzyme, group in filtered_df.groupby("enzyme"):
        position_dicts = [defaultdict(int) for _ in range(8)]

        for cleavage_site in group["cleavage_site"]:
            for i, aa in enumerate(cleavage_site):
                position_dicts[i][aa] += 1

        mean_p = group["p_value"].mean()
        unique_positions = sorted(set(group["position"]))
        total_count = len(group)
        motif = counts_to_relative_motif(position_dicts)

        enzyme_summary[enzyme] = {
            "p_value": mean_p,
            "positions": unique_positions,
            "total_count": total_count,
            "motif": motif,
        }

    enzyme_summary = dict(
        sorted(enzyme_summary.items(), key=lambda x: x[1]["total_count"], reverse=True)
    )

    return enzyme_summary
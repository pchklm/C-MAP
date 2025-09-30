import json
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

# Example: load your data (replace with your filename if needed)

DEFAULT_COLORS = {
    'A': 'limegreen',     # alanine
    'R': 'darkorchid',    # arginine
    'N': 'mediumslateblue',  # asparagine
    'D': 'crimson',       # aspartic acid
    'C': 'gold',          # cysteine
    'Q': 'teal',          # glutamine
    'E': 'orangered',     # glutamic acid
    'G': 'deepskyblue',   # glycine
    'H': 'slategray',     # histidine
    'I': 'peru',          # isoleucine
    'L': 'darkorange',    # leucine
    'K': 'blueviolet',    # lysine
    'M': 'olive',         # methionine
    'F': 'firebrick',     # phenylalanine
    'P': 'sienna',        # proline
    'S': 'turquoise',     # serine
    'T': 'steelblue',     # threonine
    'W': 'indigo',        # tryptophan
    'Y': 'darkgoldenrod', # tyrosine
    'V': 'tomato',        # valine
    'B': 'lightgray',     # aspartic acid or asparagine
    'Z': 'gray',          # glutamic acid or glutamine
    'X': 'black',
}


with open("grouped_results.json") as f:
    grouped_results = json.load(f)

# Iterate over enzymes
motifnumber = 1
for enzyme, data in grouped_results.items():
    if "motif" not in data:
        continue  # skip if no motif
    
    motif = data["motif"]

    print(enzyme)
    
    # Convert motif dict → DataFram
    # df.index = df.index.astype(int)
    df = pd.DataFrame.from_dict(motif, orient="index").fillna(0)

    # Ensure index is numeric and sorted (positions)
    df.index = df.index.astype(int)
    df = df.sort_index()

    # Plot motif logo
    fig, ax = plt.subplots(figsize=(4, 3))  # wider, flatter
    logo = logomaker.Logo(df,
                        ax=ax,
                        color_scheme=DEFAULT_COLORS)

    # Titles and labels
    ax.set_xlabel("Position relative to cleavage site", fontsize=10)
    ax.set_ylabel("Frequency", fontsize=10)

    # Force ticks only at your positions (no 0)
    positions = df.index.tolist()
    ax.set_xticks(positions)
    ax.set_xticklabels(positions)

    # Style spines
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)

    plt.tight_layout()
    
    # Save as PDF
    out_file = f"motif_{motifnumber}.pdf"
    motifnumber += 1
    plt.savefig(out_file, format="pdf", bbox_inches="tight")
    plt.close()
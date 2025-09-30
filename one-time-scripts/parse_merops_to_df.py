import mysql.connector
import pandas as pd
from collections import defaultdict
from ..helper import convert_3to1

conn = mysql.connector.connect(
        host="localhost",
        user="merops",
        password="",
        database="merops"
    )

cursor = conn.cursor()

substrates = pd.read_sql_query("""
    SELECT code, Site_P4, Site_P3, Site_P2, Site_P1,
           Site_P4prime, Site_P3prime, Site_P2prime, Site_P1prime
    FROM Substrate_search
""", conn)

# Positions we care about
positions = ["Site_P4", "Site_P3", "Site_P2", "Site_P1",
             "Site_P4prime", "Site_P3prime", "Site_P2prime", "Site_P1prime"]

# Count amino acids per position per code
counts_by_code = defaultdict(lambda: {pos: defaultdict(int) for pos in positions})

for _, row in substrates.iterrows():
    code = row["code"]
    for pos in positions:
        aa = convert_3to1(row[pos])
        if pd.notna(aa):  # skip missing
            counts_by_code[code][pos][aa] += 1


names = pd.read_sql_query("""
    SELECT code, name
    FROM protein_name
    WHERE type = 'real'
""", conn)

code_to_name = dict(zip(names.code, names.name))

# Map code -> sequence_id
domain = pd.read_sql_query("""
    SELECT code, sequence_id
    FROM domain
    WHERE type = 'peptidase'
""", conn)

# Map sequence_id -> tax_id
sequence = pd.read_sql_query("""
    SELECT sequence_id, merops_taxonomy_id
    FROM sequence
""", conn)

# Map tax_id -> organism name
organism = pd.read_sql_query("""
    SELECT merops_taxonomy_id, name
    FROM organism_name 
    WHERE recommended_name = 'yes'
""", conn)

species_map = (
    domain.merge(sequence, on="sequence_id")
          .merge(organism, left_on="merops_taxonomy_id", right_on="merops_taxonomy_id")
)


species_map_grouped = (
    species_map.groupby("code")["name"]
               .apply(lambda x: sorted(set(x)))  # unique, sorted
               .to_dict()
)

rows = []
for code, positions_dict in counts_by_code.items():

    species_list = species_map_grouped.get(code, [])
    species_str = ", ".join(species_list) if species_list else None

    row = {"code": code,
           "enzyme_name": code_to_name.get(code, None),
           "species": species_str}
    
    observed_cleavages = 0

    # Flatten amino acid counts
    for pos, aa_counts in positions_dict.items():

        total_count = sum(aa_counts.values())  # total counts for this position
    
        if total_count > observed_cleavages:
            observed_cleavages = total_count

        for aa, count in aa_counts.items():
            row[f"{pos}_{aa}"] = count
    
    print(code, observed_cleavages)
    if observed_cleavages > 0:
        rows.append(row)

df_final = pd.DataFrame(rows)
df_final.to_parquet("enzyme_motifs.parquet", engine="pyarrow", index=False)
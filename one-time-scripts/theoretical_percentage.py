import pandas as pd

# --- load files ---
# your CSV with peptides (assuming the column with peptide sequences is called "sequence")
csv_df = pd.read_csv("theoretical_results.csv")

# MaxQuant peptides.txt (tab-separated, column "Sequence")
peptides_df = pd.read_csv("MaxQuant_peptides.txt", sep="\t")

# --- get sets of sequences ---
csv_peptides = set(csv_df["sequence"].dropna().str.strip())
mq_peptides = set(peptides_df["Sequence"].dropna().str.strip())

# --- calculate overlap ---
overlap = csv_peptides.intersection(mq_peptides)
percentage = (len(overlap) / len(mq_peptides)) * 100 if mq_peptides else 0

# --- print results ---
print(f"Total peptides in peptides.txt: {len(mq_peptides)}")
print(f"Total peptides in CSV: {len(csv_peptides)}")
print(f"Peptides also in CSV: {len(overlap)}")
print(f"Percentage overlap: {percentage:.2f}%")
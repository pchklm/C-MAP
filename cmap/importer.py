from pathlib import Path
from Bio import SeqIO
import pandas as pd

def import_peptides(filepath):

    all_cols = pd.read_csv(filepath, sep="\t", nrows=0).columns

    intensity_cols = [c for c in all_cols if c.startswith("Intensity")]
    # Keep Sequence + any column starting with "Intensity"
    usecols = ["Sequence"] + intensity_cols

    # Load MaxQuant peptides.txt
    df = pd.read_csv(filepath, sep="\t", usecols=usecols)
    

    if len(intensity_cols) == 1:
        # Only one intensity column → no samples
        result = df[["Sequence", intensity_cols[0]]].rename(
            columns={intensity_cols[0]: "Intensity"}
        )
        result["Sample"] = ""  # empty string for no sample
    else:
        # Multiple intensity columns → melt into long format
        df = df.drop(columns=["Intensity"], errors="ignore")
        intensity_cols.remove("Intensity")
        melted = df.melt(
            id_vars=["Sequence"],
            value_vars=intensity_cols,
            var_name="Sample",
            value_name="Intensity"
        )
        # Clean sample names ("Intensity " prefix removal)
        melted["Sample"] = melted["Sample"].str.replace(r"^Intensity\s*", "", regex=True)
        result = melted.reset_index(drop=True)

    return result[["Sequence", "Intensity", "Sample"]]

def import_fasta(file_path: Path):
    
    records = SeqIO.parse(file_path, "fasta")
    data = [{"id": rec.id, "sequence": str(rec.seq)} for rec in records]
    return pd.DataFrame(data, columns=["id", "sequence"])
from typing import Literal, get_args

# types

AminoAcid = Literal[
    "A",  # Alanine
    "R",  # Arginine
    "N",  # Asparagine
    "D",  # Aspartic acid
    "C",  # Cysteine
    "E",  # Glutamic acid
    "Q",  # Glutamine
    "G",  # Glycine
    "H",  # Histidine
    "I",  # Isoleucine
    "L",  # Leucine
    "K",  # Lysine
    "M",  # Methionine
    "F",  # Phenylalanine
    "P",  # Proline
    "S",  # Serine  
    "T",  # Threonine
    "W",  # Tryptophan
    "Y",  # Tyrosine
    "V"   # Valine
]

# constants

amino_acids = list(get_args(AminoAcid))
        
alphabet= "".join(x for x in amino_acids)

alphabet_with_X = alphabet + "X"

alphabet_index = range(len(alphabet))

aa_to_idx = {aa: i for i, aa in enumerate(alphabet)}
aa_to_idx["X"] = 20

three_to_one = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
    }

site_columns = [
        "Site_P4", "Site_P3", "Site_P2", "Site_P1",
        "Site_P1prime", "Site_P2prime", "Site_P3prime", "Site_P4prime"
    ]

site_columns_index = range(len(site_columns))

#Standard Enzymes with their regex patterns
    
base_enzymes = {
    "S01.151": {"name": "Trypsin", "regex": [['X'], ['X'], ['X'], ['K', 'R'], ['!P'], ['X'], ['X'], ['X']]},
    "S01.151/P": {"name": "Trypsin/P", "regex": [['X'], ['X'], ['X'], ['K', 'R'], ['X'], ['X'], ['X'], ['X']]},
    "S01.281": {"name": "Arg-C", "regex": [['X'], ['X'], ['X'], ['R'], ['!P'], ['X'], ['X'], ['X']]},
    "S01.281/P": {"name": "Arg-C/P", "regex": [['X'], ['X'], ['X'], ['R'], ['X'], ['X'], ['X'], ['X']]},
    "S01.269": {"name": "Glu-C", "regex": [['X'], ['X'], ['X'], ['D', 'E'], ['!P'], ['X'], ['X'], ['X']]},
    "S01.269/P": {"name": "Glu-C/P", "regex": [['X'], ['X'], ['X'], ['D', 'E'], ['X'], ['X'], ['X'], ['X']]},
    "A01.001": {"name": "PepsinA", "regex": [['X'], ['X'], ['X'], ['F', 'L', 'I'], ['!P'], ['X'], ['X'], ['X']]},
    "A01.001/P": {"name": "PepsinA/P", "regex": [['X'], ['X'], ['X'], ['F', 'L', 'I'], ['X'], ['X'], ['X'], ['X']]},
    "S01.280": {"name": "Lys-C", "regex": [['X'], ['X'], ['X'], ['K'], ['!P'], ['X'], ['X'], ['X']]},
    "S01.280/P": {"name": "Lys-C/P", "regex": [['X'], ['X'], ['X'], ['K'], ['X'], ['X'], ['X'], ['X']]},
    "M35.004": {"name": "Lys-N", "regex": [['X'], ['X'], ['X'], ['X'], ['K'], ['X'], ['X'], ['X']]},
    "S01.131": {"name": "Elastase", "regex": [['X'], ['X'], ['X'], ['A', 'L', 'I', 'V'], ['!P'], ['X'], ['X'], ['X']]},
    "S01.268": {"name": "Alpha-lytic protease", "regex": [['X'], ['X'], ['X'], ['T', 'A', 'S', 'V'], ['X'], ['X'], ['X'], ['X']]},
    "S09.001": {"name": "proline-endopeptidase", "regex": [['X'], ['X'], ['X'], ['P'], ['X'], ['X'], ['X'], ['X']]},
    "M72.001": {"name": "Asp-N", "regex": [['X'], ['X'], ['X'], ['X'], ['D'], ['X'], ['X'], ['X']]},
    "S01.001": {"name": "Chymotrypsin", "regex": [['X'], ['X'], ['X'], ['F', 'Y', 'W', 'L'], ['!P'], ['X'], ['X'], ['X']]},
    "S01.001/P": {"name": "Chymotrypsin/P", "regex": [['X'], ['X'], ['X'], ['F', 'Y', 'W', 'L'], ['X'], ['X'], ['X'], ['X']]}
    # "cyanogen-bromide": [['X'], ['X'], ['X'], ['M'], ['X'], ['X'], ['X'], ['X']], Chem
    # "Formic_acid": [['X'], ['X'], ['X'], ['D', 'N'], ['D', 'N'], ['X'], ['X'], ['X']], Chem
    # "2-iodobenzoate": [['X'], ['X'], ['X'], ['W'], ['X'], ['X'], ['X'], ['X']], Chem
}

base_enzyme_codes = base_enzymes.keys()

base_enzyme_codes_without_P = [c for c in base_enzyme_codes if c[-2:] == "/P"]
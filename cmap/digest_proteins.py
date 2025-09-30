import re
import pandas as pd

def motif_to_regex(motif):
    """
    Convert an 8-position motif definition into a regex pattern.
    motif = list of lists like [['X'], ['X'], ['X'], ['K'], ['!P'], ['X'], ['X'], ['X']]
    """
    # Map 'X' = any amino acid, '!X' = negated character
    def pos_to_regex(chars):
        
        if chars == ["X"]:
            return "."
        else:
            regex_parts = ""
            not_allowed = False
            for c in chars:
                if c.startswith("!"):
                    regex_parts+=c[1:]
                    not_allowed = True
                else:
                    regex_parts += c
            if not_allowed:
                return f"[^{regex_parts}]"
            else:
                 return f"[{regex_parts}]"

    # 8 characters → P4 P3 P2 P1 | P1′ P2′ P3′ P4′
    left = "".join(pos_to_regex(x) for x in motif[:4])   # before cleavage
    right = "".join(pos_to_regex(x) for x in motif[4:])  # after cleavage

    # Cleavage site is between P1 and P1′
    return f"(?<={left})(?={right})"



def digest_with_motif(sequence, motif, min_len=7, max_len=35):
    """
    Digest sequence with motif, return peptides + 8-aa cleavage window.
    """
    pattern = motif_to_regex(motif)
    sites = [0] + [m.start() for m in re.finditer(pattern, sequence)] + [len(sequence)]
    #print(pattern, sequence)
    #print(sites)
    
    peptides = []
    

    for site_index in range(len(sites)-1):
        start = sites[site_index]
        end_index = site_index+1
        end = sites[end_index]
        pep = sequence[start:end]
        while len(pep) <= max_len:
            if len(pep) >= min_len:
                if start == 0:
                    n_term_cleavage_window = None
                else:
                    n_term_cleavage_window = sequence[start-4:start+4]
                if end == len(sequence):
                    c_term_cleavage_window = None
                else:
                    c_term_cleavage_window = sequence[end-4:end+4]
                peptides.append({
                    "sequence": pep,
                    "n_term_cleavage_window": n_term_cleavage_window,
                    "c_term_cleavage_window": c_term_cleavage_window
                })
            end_index += 1
            if end_index > len(sites)-1:
                break
            end = sites[end_index]
            pep = sequence[start:end]
            
    return peptides

def digest_proteins(fasta_df, enzymes, min_len=7, max_len=35):
    all_peptides = []
    for _, row in fasta_df.iterrows():
        protein_id = row["id"]
        seq = row["sequence"]

        for enzyme in enzymes:
            #print(enzyme, enzymes[enzyme])

            peptides = digest_with_motif(seq, enzymes[enzyme], min_len, max_len)
            for pep in peptides:
                pep["protein_id"] = protein_id
                pep["enzyme"] = enzyme
                all_peptides.append(pep)

    return pd.DataFrame(all_peptides, columns=["protein_id", "sequence", "enzyme", "n_term_cleavage_window", "c_term_cleavage_window"])
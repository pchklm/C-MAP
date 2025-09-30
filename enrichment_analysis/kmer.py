from collections import defaultdict
from .constants import amino_acids

def build_kmer_index_and_background(fasta, k=6):
    '''
    Build a kmer index while also counting aminoacids to provide a background count of each amino acid.

    args:
        fasta: Fasta file containing protein id's and sequences.
        k: Number determining the length of the k-mers.

    returns:
        kmer_index: Dictionary mapping kmers to protein id's.
        protein_sequences: Dictionary mapping protein id's to protein sequences.
        background: Dictionary with the total count of each amino acid.
    '''

    kmer_index = defaultdict(list)
    protein_sequences = {}
    background = defaultdict(int, {aa: 1 for aa in amino_acids})

    for protein in fasta.itertuples():
        sequence = protein.sequence
        protein_sequences[protein.id] = sequence
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_index[kmer].append((protein.id, i))
            background[sequence[i]] += 1
        for j in sequence[-k:]:
            background[j] +=1

    return kmer_index, protein_sequences, background
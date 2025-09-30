from dataclasses import dataclass

from .helper import search_function
from .constants import alphabet
from .preprocessing import get_enzyme_df, get_filtered_enzyme_df, get_cleavage_sites
from .kmer import build_kmer_index_and_background
from .regex_trie import RegexTrie
from .motifs import analyze_enzymes
from .matching import match_enzymes
from .postprocessing import filter_results, group_results, group_theoretical_cleavages
from .digest_proteins import digest_proteins

import time

@dataclass
class CleavageEnrichmentAnalysis:
    _fasta = None
    _peptide_df = None
    _metadata = None

    use_standard_enzymes = True
    species = None
    enzymes = None
    theoretical_enzymes = None
    possible_species = None
    possible_enzymes = None

    _enzyme_df = None
    _kmer_index = None
    _protein_sequences = None
    _background = None
    _result = None
    _theoretical_result = None
    _calculated = False
    _theoretical_calculated = False

    def __post_init__(self):
        (self._enzyme_df,
         self.possible_species,
         self.possible_enzymes) = get_enzyme_df()
        
    def __setattr__(self, key, value):
        if ((key == "species" and self.species != value) or
            (key == "enzymes" and self.enzymes != value) or
            (key == "use_standard_enzymes" and self.use_standard_enzymes != value)):
            object.__setattr__(self, "_calculated", False)

        object.__setattr__(self, key, value)

    def set_fasta(self, fasta):
        starttime = time.perf_counter()
        self._fasta = fasta
        (self._kmer_index,
         self._protein_sequences,
         self._background) = build_kmer_index_and_background(self._fasta)
        endtime = time.perf_counter()
        print(f"Time to build k-mer index: {endtime - starttime:.4f} seconds")

        if self._peptide_df is not None:
            self._peptide_df = get_cleavage_sites(self._peptide_df, self._kmer_index, self._protein_sequences)
            
        
    def set_peptides(self, peptides):
        if self._fasta is not None:
            starttime = time.perf_counter()
            self._peptide_df = get_cleavage_sites(peptides, self._kmer_index, self._protein_sequences)
            endtime = time.perf_counter()
            print(f"Time to get cleavage sites: {endtime - starttime:.4f} seconds")
        else:
            self._peptide_df = peptides

    def get_results(self, proteinIDs, metadata_filter):
        if not self._calculated:
            self.calculate()
        return filter_results(self._result, proteinIDs, metadata_filter)
    
    def get_grouped_results(self, proteinIDs, metadata_filter, k = 3):
        if not self._calculated:
            self.calculate()
        return group_results(self._result, proteinIDs, metadata_filter, k)
    
    def get_theoretical_results(self, proteinIDs):
        if not self._theoretical_calculated:
            self.calculate_theoretical_peptides()

        return filter_results(self._theoretical_result, proteinIDs)

    def get_grouped_theoretical(self, proteinIDs):
        if not self._theoretical_calculated:
            self.calculate_theoretical_peptides()

        return group_theoretical_cleavages(self._theoretical_result, proteinIDs)

    def calculate(self):
        filtered_enzyme_df = get_filtered_enzyme_df(self._enzyme_df, self.use_standard_enzymes, self.species, self.enzymes)

        #calculate position sepecific scoring matrices and regexes for each enzyme
        starttime = time.perf_counter()
        pssms, regexs, code_to_name = analyze_enzymes(filtered_enzyme_df, self._background)
        endtime = time.perf_counter()
        print(f"Time to extract cleavage motifs: {endtime - starttime:.4f} seconds")

        # build Trie based on regexes
        trie = RegexTrie(alphabet)
        for code in regexs:
            regex = regexs[code]
            trie.insert(regex, code)
        
        #match enzymes for each cleavage
        starttime = time.perf_counter()
        self._result = match_enzymes(self._peptide_df, trie, pssms, code_to_name, self._background)
        endtime = time.perf_counter()
        print(f"Time to match enzymes: {endtime - starttime:.4f} seconds")
        self._calculated = True

    def calculate_theoretical_peptides(self):
        starttime = time.perf_counter()
        filtered_enzyme_df = get_filtered_enzyme_df(self._enzyme_df, False, None, self.theoretical_enzymes)
        _, regexes, _ = analyze_enzymes(filtered_enzyme_df, self._background)
        self._theoretical_result =  digest_proteins(self._fasta, regexes)
        endtime = time.perf_counter()
        print(f"Time to calculate theoretical cleavages: {endtime - starttime:.4f} seconds")
        self._theoretical_calculated = True

    def search_species(self, input):
        return search_function(input, self.possible_species)
    
    def search_enzymes(self, input):
        return search_function(input, self.possible_enzymes)
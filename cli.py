#!/usr/bin/env python
import click
import json
import pandas as pd
from enrichment_analysis.importer import import_peptides, import_fasta
from enrichment_analysis.cleavage_enrichment_analysis import CleavageEnrichmentAnalysis

@click.command()
@click.argument("peptide_file")
@click.argument("fasta_file")
@click.option(
    "--use-standard-enzymes", "-u", 
    is_flag = True,
    default=False, 
    help="Include a selection of common experimentally enzymes."
)
@click.option(
    "--species", "-s",
    default = None,
    help="Target species for the analysis. Example: -s 'homo sapiens'.")

@click.option(
    "--enzymes", "-e",
    default = [],
    multiple=True,
    help="One or more specific enzymes to use. Example: -e trypsin -e chymotrypsin."
)
@click.option(
    "--theoretical-enzymes", "-t",
    default = [],
    multiple=True,
    help="One or more specific enzymes to use. Example: -t trypsin -t chymotrypsin"
)
@click.option(
    "--proteins", "-p",
    default = [],
    multiple=True,
    help="One or more proteinIDs to filter the results."
)
@click.option(
    "--metadatafilter", "-m",
    default = [],
    multiple=True,
    help="List of included groups or sample names. Example: -m Sample1 -m Group42"
)
@click.option(
    "--top-k", "-k",
    "k",
    default = 3,
    help="Top-k enzymes to be included in grouped results."
)

def main(peptide_file, fasta_file, use_standard_enzymes, species, enzymes, theoretical_enzymes, proteins, metadatafilter, k):
    peptide_df = import_peptides(peptide_file)
    fasta = import_fasta(fasta_file)
    cmap = CleavageEnrichmentAnalysis()
    cmap.set_fasta(fasta)
    cmap.set_peptides(peptide_df)
    cmap.use_standard_enzymes = use_standard_enzymes
    cmap.species = species
    cmap.enzymes = list(enzymes)
    cmap.theoretical_enzymes = list(theoretical_enzymes)

    results: pd.DataFrame = cmap.get_results(proteins, metadatafilter)
    grouped_results = cmap.get_grouped_results(proteins, metadatafilter, k)
    theoretical_results: pd.DataFrame = cmap.get_theoretical_results(proteins)
    grouped_theoretical = cmap.get_grouped_theoretical(proteins)

    results.to_csv("results/results.csv")
    theoretical_results.to_csv("results/theoretical_results.csv")
    for enzyme, data in grouped_results.items():
        if "motif" in data and isinstance(data["motif"], pd.DataFrame):
            data["motif"] = data["motif"].to_dict(orient="index")
    with open("results/grouped_results.json","w") as fp:
        json.dump(grouped_results, fp)
    for enzyme, data in grouped_theoretical.items():
        if "motif" in data and isinstance(data["motif"], pd.DataFrame):
            data["motif"] = data["motif"].to_dict(orient="index")
    with open("results/grouped_theoretical_results.json","w") as fp:
        json.dump(grouped_theoretical, fp)

if __name__ == "__main__":
    main()
import mysql.connector
import pandas as pd
from collections import defaultdict

def get_enzyme_df(species, proteases):
    print(species, proteases)
    conn = mysql.connector.connect(
        host="localhost",
        user="cleavage_user",
        password="ILOVEJOSIE",
        database="merops"
    )
    cursor = conn.cursor()

    if (species != None):

        # 1. Get merops_taxonomy_ID for species
        taxonomy_query = f"""
            SELECT merops_taxonomy_id
            FROM organism_name
            WHERE LOWER(name) = '{species}';
        """
        taxonomy_df = pd.read_sql(taxonomy_query, conn)
        if taxonomy_df.empty:
            return pd.DataFrame()  # No match found
        taxonomy_id = taxonomy_df['merops_taxonomy_id'].iloc[0]

        # 2. Get sequence_ids for that taxonomy ID
        sequence_query = f"""
            SELECT sequence_id
            FROM sequence
            WHERE merops_taxonomy_id = '{taxonomy_id}';
        """
        sequence_df = pd.read_sql(sequence_query, conn)
        sequence_ids = tuple(sequence_df['sequence_id'].tolist())
        if not sequence_ids:
            return pd.DataFrame()

        # 3. Get codes from domain table for these sequence_ids
        if (proteases != None):
            specific_protease_query = [f"protein = '{p}'" for p in proteases]
            specific_protease_query =  " OR ".join(specific_protease_query)
        else:
            specific_protease_query = "FALSE"

        print(specific_protease_query)

        domain_query = f"""
            SELECT code
            FROM domain
            WHERE (sequence_id IN {sequence_ids} OR {specific_protease_query}) AND type = 'peptidase';
        """
        domain_df = pd.read_sql(domain_query, conn)
        codes = tuple(set(domain_df['code'].tolist()))

        # 4. Get final data from substrate_search
        final_query = f"""
            SELECT code, Uniprot, Protease, Site_P4, Site_P3, Site_P2, Site_P1,
                Site_P4prime, Site_P3prime, Site_P2prime, Site_P1prime,
                organism, Substrate_name
            FROM Substrate_search
            WHERE code IN {codes};
        """

        name_query = f"""
            SELECT code, name
            FROM protein_name
            WHERE code IN {codes} AND type = 'real'
        """
    else:

        final_query = f"""
            SELECT code, Uniprot, Site_P4, Site_P3, Site_P2, Site_P1,
                Site_P4prime, Site_P3prime, Site_P2prime, Site_P1prime
            FROM Substrate_search;
        """

        name_query = f"""
            SELECT code, name
            FROM protein_name
            WHERE type = 'real'
        """

    merops_df = pd.read_sql(final_query, conn)

    code_mapping_df = pd.read_sql(name_query, conn)

    code_to_proteins = defaultdict(list)

    for _, row in code_mapping_df.iterrows():
        code = row['code']
        protein = row['name']
        code_to_proteins[code] = protein
    

    cursor.close()
    conn.close()

    return merops_df, code_to_proteins
#!/usr/bin/env python3

'''
This mouse-specific script will determine the canonical transcript for each gene based on Ensembl 'is_canonical' annotation.
If no unambiguous canonical transcript can be derived from the annotation, the longest protein will be used.
Additionally, some dataframe reformatting is performed to resemble the human dataframe.
'''

import pandas as pd
import argparse

def load_MGI_data(ensembl_data, genemodel_data):
    """Loads MGI mouse data frames and combines relevant columns. 
    Quite a bit of column renaming is needed to get the same format as the human table.
    """

    ensembl_df = pd.read_csv(ensembl_data, header=None, sep="\t", index_col=False, dtype=str)
    gene_df = pd.read_csv(genemodel_data, header="infer", sep="\t", index_col=False, dtype=str)

    # manually select and rename important columns
    gene_df = gene_df.iloc[:,[2,0,1,3,5,10,11]]
    gene_df.columns = ["hgnc_symbol","hgnc_id","locus_group","approved_name","entrez_gene_id","ensembl_gene_id","chromosome"]

    ensembl_df = ensembl_df.iloc[:,[0,1,5,8]]
    ensembl_df.columns = ["hgnc_id", "hgnc_symbol", "ensembl_gene_id", "locus_type"]

    # merge relevant columns
    merged_df = gene_df.merge(ensembl_df[["hgnc_id", "locus_type"]], on="hgnc_id", how="right")

    # symbols should be uppercase
    merged_df["hgnc_symbol"] = merged_df["hgnc_symbol"].str.upper()
    
    return(merged_df)
    

def get_canonical_transcript_by_ensembl(transcript_info):
    """Sorts the dataframe on canonical (1/0) and then on protein length. 
    Then for each symbol retrieves the first gene.
    """

    transcripts = pd.read_csv(transcript_info, sep="\t", dtype=str)

    canonical_transcripts = transcripts.sort_values(["is_canonical", "protein_length"], ascending=False)\
        .groupby("hgnc_symbol").head(1)

    return(canonical_transcripts)


def add_canonical(merged_df, canonical_transcripts):
    """ Add canonical transcripts for each gene, and does some reformatting of the df
    """

    # add transcript ID per gene
    merged_df = merged_df.merge(canonical_transcripts[["gene_stable_id","transcript_stable_id"]], left_on="ensembl_gene_id", right_on="gene_stable_id", how="left")
    merged_df.drop(columns=["gene_stable_id"], inplace=True)
    merged_df.rename({
        "ensembl_gene_id": "ensembl_canonical_gene",
        "transcript_stable_id": "ensembl_canonical_transcript"
    }, axis=1, inplace=True)

    # add canonical transcript columns for overrides -- not sure if needed, but at least ot mimic the human file
    merged_df["genome_nexus_canonical_transcript"] = merged_df["ensembl_canonical_transcript"].copy()
    merged_df["mskcc_canonical_transcript"] = merged_df["ensembl_canonical_transcript"].copy()
    merged_df["uniprot_canonical_transcript"] = merged_df["ensembl_canonical_transcript"].copy()

    # ordering
    merged_df = merged_df[["hgnc_symbol","ensembl_canonical_gene",
    "ensembl_canonical_transcript","genome_nexus_canonical_transcript",
    "uniprot_canonical_transcript","mskcc_canonical_transcript",
    "hgnc_id","approved_name","locus_group","locus_type","entrez_gene_id",
    "chromosome"]]

    # add empty columns based on the human file, because the controller
    # will raise nullPointerExceptions if certain columns don't exist.
    merged_df = merged_df.reindex(columns=list(merged_df.columns)+["status","location_sortable","synonyms",
    "alias_name","previous_symbols","prev_name","gene_family","gene_family_id",
    "date_approved_reserved","date_symbol_changed","date_name_changed","date_modified",
    "vega_id","ucsc_id","accession_numbers","refseq_ids","ccds_id","uniprot_id",
    "pubmed_id","mgd_id","rgd_id","lsdb","cosmic","omim_id","mirbase","homeodb",
    "snornabase","bioparadigms_slc","orphanet","pseudogene.org","horde_id","merops",
    "imgt","iuphar","kznf_gene_catalog","mamit-trnadb","cd","lncrnadb","enzyme_id",
    "intermediate_filament_db","rna_central_ids"])
    
    return(merged_df)


def main(transcript_info, ensembl_data, genemodel_data, canonical_transcripts_per_symbol):
    mgi_data = load_MGI_data(ensembl_data, genemodel_data)
    canonical = get_canonical_transcript_by_ensembl(transcript_info)
    formatted_df = add_canonical(mgi_data, canonical)
    
    formatted_df.to_csv(canonical_transcripts_per_symbol, sep="\t", header=True, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("transcript_info",
                        help="tmp/ensembl_canonical_data.txt")
    parser.add_argument("ensembl_data",
                        help="common_input/mouse/MRK_ENSEMBL.rpt")
    parser.add_argument("genemodel_data",
                        help="common_input/mouse/MGI_Gene_Model_Coord.rpt")
    parser.add_argument("ensembl_biomart_canonical_transcripts_per_mgi",
                        help="tmp/ensembl_biomart_canonical_transcripts_per_mgi.txt")
    args = parser.parse_args()

    main(args.transcript_info, args.ensembl_data, args.genemodel_data, args.ensembl_biomart_canonical_transcripts_per_mgi)

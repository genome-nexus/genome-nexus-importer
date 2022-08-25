#!/usr/bin/env python3
"""Takes in a combined 2d and 3d hotspot file (grch37) and updates the transcript ids to their respective grch38 ids.

It does so by following a 2 step process:
1. Parse the file ../../data/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt
   and produce a map of hugo symbol X canonical transcript id.

2. Use this hugo symbol map and follow these rules to decide whether to REPLACE the transcript id or to DROP the hotspot row:
   2.1. get the grch37 transcript's translated sequence using (here using BRAF's transcripts as an example):
     curl 'https://grch37.rest.ensembl.org/sequence/id/ENST00000288602?type=protein' -H 'Content-type:text/plain'
    and do the same for the grch38 one using
     curl 'https://rest.ensembl.org/sequence/id/ENST00000646891?type=protein' -H 'Content-type:text/plain'
   2.2. if the sequences match, accept the grch38 transcript id, replacing the grch37 one with it in the hotspots file
   2.3. if there is a mismatch, remove the entry from the hotspots file. This hotspots entry will need a new dedicated analysis and validation,
        and cannot be lifted over as is.

Logs and outputs:
 - outputs the mapping from step 1
 - reports the total number of rows in original hotspots file, the number of rows replaced and the number of rows dropped
 - outputs a new updated combined 2d and 3d hotspot file (grch38)

Links to used API docs: 
 - https://rest.ensembl.org/documentation/info/sequence_id
 - https://grch38.genomenexus.org/swagger-ui.html 
 - https://www.cbioportal.org/api/swagger-ui/index.html
"""


import argparse
from xmlrpc.client import boolean
import requests, sys
import pandas as pd
 
ENSEMBL_GRCH38_SERVER = "https://rest.ensembl.org"
ENSEMBL_GRCH37_SERVER = "https://grch37.rest.ensembl.org"
CBIOPORTAL_SERVER = "https://www.cbioportal.org"
GENOME_NEXUS_GRCH38_SERVER = "https://grch38.genomenexus.org"


def get_translated_protein_sequence(ensembl_server: str, transcript_id: str) -> str:
    """ Returns the translated protein sequence for the given transcript id
    """
    api_url = "{0}/sequence/id/{1}?type=protein".format(ensembl_server, transcript_id) 
    response = requests.get(api_url, headers={ "Content-Type" : "text/plain"})
    if not response.ok:
       response.raise_for_status()
    return response.text


def get_gene_and_transcript_map(hugo_and_transcript_id_file_name: str) -> dict[str, str]:
    """ Returns a map of hugo symbol X canonical transcript id.
    """
    # parse file w/ pandas and build the result map:
    hugo_and_transcript_df = pd.read_csv(hugo_and_transcript_id_file_name, sep='\t')
    hugo_and_transcript_df.set_index("hgnc_symbol", inplace=True)
    hugo_and_transcript_map = hugo_and_transcript_df.to_dict(orient="index")
    return hugo_and_transcript_map


def validate_new_transcript_id(grch37_transcript_id: str,
                               proposed_grch38_transcript_id: str) -> boolean:
    """Compares the translated sequence of the grch37 transcript vs the grch38 one
       and returns True if the sequences match
    """
    grch37_sequence = get_translated_protein_sequence(ENSEMBL_GRCH37_SERVER, grch37_transcript_id)
    grch38_sequence = get_translated_protein_sequence(ENSEMBL_GRCH38_SERVER, proposed_grch38_transcript_id)
    return grch37_sequence == grch38_sequence

def generate_updated_grch38_hotspots_file(hugo_and_transcript_id_file_name: str,
                                          hotspots_2d_3d_grch37_file_name: str,
                                          isoform_override_source: str) -> None:
    hugo_and_transcript_map = get_gene_and_transcript_map(hugo_and_transcript_id_file_name)
    hotspots_grch37_df = pd.read_csv(hotspots_2d_3d_grch37_file_name, sep='\t')
    
    # iterate over the rows and replace with new transcript_id if valid:
    for index, row in hotspots_grch37_df.iterrows():
        hugo_symbol = row['hugo_symbol']
        grch37_transcript_id = row['transcript_id']
        proposed_grch38_transcript_id = hugo_and_transcript_map[hugo_symbol][isoform_override_source + "_canonical_transcript"]
        is_valid = validate_new_transcript_id(grch37_transcript_id,
                                              proposed_grch38_transcript_id)
        # if valid, replace transcript_id:
        if is_valid:
            row['transcript_id'] = proposed_grch38_transcript_id
        else:
            # drop any row where new transcript_id is not valid:
            hotspots_grch37_df.drop(index, inplace=True)
    
    # write the dataframe out as a new grch38 file:
    # TODO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hugo_and_transcript_id_file_name", default="../../data/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt",
                        type=str, help="hugo symbols with respective transcript ids mapping file")
    parser.add_argument("--isoform_override_source", default="mskcc",
                        type=str, help="which transcript override source is preferred. Possible values: [mskcc, uniprot, genome_nexus, ensembl]")
    parser.add_argument("--hotspots_2d_3d_grch37_file_name", default="../../data/grch37_ensembl92/export/hotspots_v2_and_3d.txt",
                        type=str, help="combined grch37 2D and 3D cancerhotspots data file")
    args = parser.parse_args()


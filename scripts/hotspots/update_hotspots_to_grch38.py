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
import logging
import math

ENSEMBL_GRCH38_SERVER = "https://rest.ensembl.org"
ENSEMBL_GRCH37_SERVER = "https://grch37.rest.ensembl.org"
CBIOPORTAL_SERVER = "https://www.cbioportal.org"
GENOME_NEXUS_GRCH38_SERVER = "https://grch38.genomenexus.org"

protein_sequence_cache = {}

def get_translated_protein_sequence(ensembl_server: str, transcript_id: str) -> str:
    """ Returns the translated protein sequence for the given transcript id.
        Returns None if the webservice failed with a 400 type code.
    """
    cache_key = "{0}_{1}".format(ensembl_server, transcript_id)
    if cache_key in protein_sequence_cache:
        return protein_sequence_cache[cache_key]
    api_url = "{0}/sequence/id/{1}?type=protein".format(ensembl_server, transcript_id) 
    response = requests.get(api_url, headers={ "Content-Type" : "text/plain"})
    if not response.ok:
        if response.status_code >= 400 and response.status_code < 500:
            return None
        else:
            response.raise_for_status()
    protein_sequence_cache[cache_key] = response.text
    return response.text


def get_gene_and_transcript_map(hugo_and_transcript_id_file_name: str) -> dict[str, str]:
    """ Returns a map of hugo symbol X canonical transcript id.
    """
    # parse file w/ pandas and build the result map:
    hugo_and_transcript_df = pd.read_csv(hugo_and_transcript_id_file_name, sep='\t', dtype=str)
    hugo_and_transcript_df.set_index("hgnc_symbol", inplace=True)
    hugo_and_transcript_map = hugo_and_transcript_df.to_dict(orient="index")
    return hugo_and_transcript_map


def new_transcript_id_is_valid(grch37_transcript_id: str,
                               proposed_grch38_transcript_id: str) -> boolean:
    """Compares the translated sequence of the grch37 transcript vs the grch38 one
       and returns True if the sequences match.
       Returns False if the sequences do not match or if the comparison was not
       possible (e.g. because the transcript id one either side is invalid)
    """
    grch37_sequence = get_translated_protein_sequence(ENSEMBL_GRCH37_SERVER, grch37_transcript_id)
    grch38_sequence = get_translated_protein_sequence(ENSEMBL_GRCH38_SERVER, proposed_grch38_transcript_id)
    return grch37_sequence == grch38_sequence


def get_new_hugo_symbol(old_hugo_symbol: str, hugo_and_transcript_map: str) -> str:
    """Looks for the give old_hugo_symbol in  "previous_symbols" or "synonyms" columns of the map values
       and when found, returns the corresponding key (the new hugo_symbol). If not found, returns None.
    """
    for new_hugo_symbol, other_columns in hugo_and_transcript_map.items():
        if type(other_columns['previous_symbols']) == str:
            previous_symbols = [x.strip() for x in other_columns['previous_symbols'].split(',')]
            if old_hugo_symbol in previous_symbols:
                logging.info("found {0} in 'previous_symbols", old_hugo_symbol)
                return new_hugo_symbol
        if type(other_columns['synonyms']) == str:
            synonyms = [x.strip() for x in other_columns['synonyms'].split(',')]
            if old_hugo_symbol in synonyms:
                logging.info("found {0} in 'synonyms", old_hugo_symbol)
                return new_hugo_symbol
    return None


def find_grch38_transcript_id(hugo_symbol: str, hugo_and_transcript_map: dict, grch38_isoform_override_source: str) -> str:
    """Tries to find the hugo symbol in the map's index or in the alias column. If found,
       returns the respective transcript id according to grch38_isoform_override_source.
       If not found, returns None.
    """
    if hugo_symbol not in hugo_and_transcript_map:
        # then try to find the row that has this symbol as one of the values in "previous_symbols" or "synonyms"
        hugo_symbol = get_new_hugo_symbol(hugo_symbol, hugo_and_transcript_map)
        if hugo_symbol is None:
            return None
          
    return hugo_and_transcript_map[hugo_symbol][grch38_isoform_override_source + "_canonical_transcript"]


def generate_updated_grch38_hotspots_info(grch38_hugo_and_transcript_id_file_name: str,
                                          grch38_isoform_override_source: str,
                                          grch37_hotspots_2d_3d_file_name: str) -> None:
    hugo_and_transcript_map = get_gene_and_transcript_map(grch38_hugo_and_transcript_id_file_name)
    grch37_hotspots_df = pd.read_csv(grch37_hotspots_2d_3d_file_name, sep='\t', dtype=str)
    
    # iterate over the rows and replace with new transcript_id if valid:
    for index, row in grch37_hotspots_df.iterrows():
        hugo_symbol = row['hugo_symbol']
        logging.info("Processing {0} ...", hugo_symbol)
        grch37_transcript_id = row['transcript_id']
        proposed_grch38_transcript_id = find_grch38_transcript_id(hugo_symbol, hugo_and_transcript_map, grch38_isoform_override_source)
        is_valid = proposed_grch38_transcript_id is not None and \
                   new_transcript_id_is_valid(grch37_transcript_id,
                                              proposed_grch38_transcript_id)
        # if valid, replace transcript_id:
        if is_valid:
            if row['transcript_id'] != proposed_grch38_transcript_id:
                logging.info('Updating transcript of {0} from {1} to {2}', hugo_symbol, row['transcript_id'], proposed_grch38_transcript_id)
                row['transcript_id'] = proposed_grch38_transcript_id
        else:
            # drop any row where new transcript_id is not valid:
            logging.info('Dropping entry for {0}', hugo_symbol)
            grch37_hotspots_df.drop(index, inplace=True)
    
    # return the updated dataframe:
    return grch37_hotspots_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--grch38_hugo_and_transcript_id_file_name", default="../../data/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt",
                        type=str, help="hugo symbols with respective transcript ids mapping file")
    parser.add_argument("--grch38_isoform_override_source", default="mskcc",
                        type=str, help="which transcript override source is preferred. Possible values: [mskcc, uniprot, genome_nexus, ensembl]")
    parser.add_argument("--grch37_hotspots_2d_3d_file_name", default="../../data/grch37_ensembl92/export/hotspots_v2_and_3d.txt",
                        type=str, help="combined grch37 2D and 3D cancerhotspots data file")
    args = parser.parse_args()


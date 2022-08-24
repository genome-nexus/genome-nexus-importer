#!/usr/bin/env python3
"""Takes in a combined 2d and 3d hotspot file (grch37) and updates the transcript ids to their respective grch38 ids.

It does so by following a 2 step process:
1. Take all genes from https://www.cbioportal.org/api/genes?direction=ASC&pageNumber=0&pageSize=10000000&projection=SUMMARY and
   query genomenexus on grch38 to retrieve the canonical transcript id for each gene (using isoformOverrideSource=mskcc):
    curl -X POST --header 'Content-Type: application/json' --header 'Accept: application/json' \
        -d '["TP53","PIK3CA","BRCA1"]' 'https://grch38.genomenexus.org/ensembl/canonical-transcript/hgnc?isoformOverrideSource=mskcc'
   and produce a map of hugo symbol X canonical transcript id.
2. Use this hugo symbol map and follow these rules to decide whether to REPLACE the transcript id or to DROP the hotspot row:
   2.1. get the grch37 transcript's translated sequence using (here using BRAF's transcripts as an example):
     curl 'https://grch37.rest.ensembl.org/sequence/id/ENST00000288602?type=protein' -H 'Content-type:text/x-fasta'
    and do the same for the grch38 one using
     curl 'https://rest.ensembl.org/sequence/id/ENST00000646891?type=protein' -H 'Content-type:text/x-fasta'
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
import requests, sys
 
ENSEMBL_GRCH38_SERVER = "https://rest.ensembl.org"
ENSEMBL_GRCH37_SERVER = "https://grch37.rest.ensembl.org"
CBIOPORTAL_SERVER = "https://www.cbioportal.org"
GENOME_NEXUS_GRCH38_SERVER = "https://grch38.genomenexus.org"


def get_translated_protein_sequence(ensembl_server: str, transcript_id: str) -> str:
    """ Returns the translated protein sequence for the given transcript id
    """
    api_url = "{0}/sequence/id/{1}?type=protein".format(ensembl_server, transcript_id) 
    response = requests.get(api_url, headers={ "Content-Type" : "text/x-fasta"})
    if not response.ok:
       response.raise_for_status()
    return response.text


def get_gene_and_transcript_map() -> dict[str, str]:
    """ Returns a map of hugo symbol X canonical transcript id.
    """
    hugo_symbols = get_cbio_genes()
    # call the canonical-transcript endpoint with all the hugo_symbols:
    api_url = "{0}/ensembl/canonical-transcript/hgnc?isoformOverrideSource=mskcc".format(GENOME_NEXUS_GRCH38_SERVER)
    symbols = '{ "symbols" : {0} }'.format(hugo_symbols)
    response = requests.post(api_url, headers={ "Content-Type" : "application/json", "Accept" : "application/json"}, data=symbols)
    if not response.ok:
       response.raise_for_status()
    # the response is a list with objects that contain many fields, including the following:
    #  {
    #     "transcriptId": "ENST00000269305",
    #     "hugoSymbols": [
    #       "TP53"
    #     ],
    #     "proteinLength": 393,
    #   },
    
    # TODO - build the result map

    return None


def get_cbio_genes() -> list[str]:
    api_url = "{0}/api/genes?direction=ASC&pageNumber=0&pageSize=10000000&projection=SUMMARY".format(CBIOPORTAL_SERVER) 
    response = requests.get(api_url)
    if not response.ok:
        response.raise_for_status()
    cbio_genes_json = response.json
    # the above is a list of objects like:
    #       {
    #        "entrezGeneId": 80119,
    #        "hugoGeneSymbol": "PIF1",
    #        "type": "protein-coding"
    #       },
    # return only the ones that have entrez_id > 0 (real ones):
    hugo_symbols =  [x.hugoGeneSymbol for x in cbio_genes_json if x.entrezGeneId > 0]
    return hugo_symbols


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hotspots_2d_3d_grch37", default="../../data/grch37_ensembl92/export/hotspots_v2_and_3d.txt",
                        type=str, help="combined grch37 2D and 3D cancerhotspots data file")
    args = parser.parse_args()


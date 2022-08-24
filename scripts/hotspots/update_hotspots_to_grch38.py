#!/usr/bin/env python3
"""Takes in a combined 2d and 3d hotspot file (grch37) and updates the transcript ids to their respective grch38 ids.

It does so by following a 2 step process:
1. Take all genes from https://www.cbioportal.org/api/genes?direction=ASC&pageNumber=0&pageSize=10000000&projection=SUMMARY and
   query genomenexus on grch38 to retrieve the canonical transcript id for each gene (using isoformOverrideSource=mskcc):
    curl -X POST --header 'Content-Type: application/json' --header 'Accept: application/json' \
        -d '["TP53","PIK3CA","BRCA1"]' 'https://grch38.genomenexus.org/ensembl/canonical-transcript/hgnc?isoformOverrideSource=mskcc'
   and produce a map of hugo symbol X canonical transcript id.
2. Use this hugo symbol map and follow these rules to decide whether to REPLACE the transcript id or to DROP the hotspot row:
   2.1. get the grch37 transcript's translated sequence using
    curl 'https://rest.ensembl.org/sequence/id/ENST00000288602?type=protein' -H 'Content-type:text/x-fasta'
    and do the same for the grch38 one
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


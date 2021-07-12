# Map transcript id to Uniprot id

### 1. Files and links
##### 1.1 Uniprot sequence FASTA (only have GRCh38)
- download "FASTA (canonical & isoform)" https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score
##### 1.2 Uniprot Cross-reference (only have GRCh38)
- add "Cross-reference (Ensembl)" in column and download "Tab-separated" https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score  
##### 1.3 GRCh37 Ensembl sequence FASTA
- http://ftp.ensembl.org/pub/grch37/release-104/fasta/homo_sapiens/pep/
##### 1.4 GRCh38 Ensembl sequence FASTA
- http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/pep/
##### 1.5 Cancer gene list
- https://www.oncokb.org/cancerGenes
##### 1.6 Biomart
- GRCh37: https://grch37.ensembl.org/biomart/martview/145410cb02e693b9da427f9c3dffe61f
- GRCh38: http://uswest.ensembl.org/biomart/martview/84262bd39b9c7d095757b6ce493a2a9f
##### 1.7. Previous mapping results
- GRCh37: https://docs.google.com/spreadsheets/d/14PN6RtFq_GTAu8OKNUyUNKJ_fj7OlcEhDjbMdapqqo8/edit#gid=0
- GRCh38: https://docs.google.com/spreadsheets/d/1slDx9zorUuA-xsmH1i9_6CIGjQB1GIgDlWDU5gw6f9I/edit#gid=0

### 2. Mapping
##### 2.1. Get all transcript ids - `df_transcript`
- columns: enst_id, ensp_id, ensembl_protein_length, ccds_id, uniprot_id
##### 2.2. Generate Uniprot sequence dictionary (1.1) - `sequence_to_uniprot_dict`
- key: sequence, value: [uniprot_ids]
##### 2.3. Generate Ensembl sequence dictionary (1.3 or 1.4) - `ensp_to_sequence_dict`
- key: ensp, value: sequence
##### 2.4. For every "ensp_id" in df_transcript, get "sequence_ensembl" from ensp_to_sequence_dict(2.3), then use sequence_ensembl as the key to get uniprot_ids list from sequence_to_uniprot_dict(2.2)
- add results to column: uniprot_id_with_isoform
##### 2.5. For every "ensp_id" in df_transcript, get uniprot id from biomart(1.6)
- add results to column: biomart_uniprot_id
##### 2.6. Compare "uniprot_id_with_isoform" with "uniprot_biomart". 
- First extract "uniprot_id" from "uniprot_id_with_isoform", which would be the substring before "-"(e.g. for "Q9Y3S1-3", we will compare "Q9Y3S1" with biomart uniprot, because biomart uniprot doesn't have isoform). If they match then return true in column "is_matched", otherwise return false.
- add results to column: is_matched
##### 2.7. Curation
- 2.7.1 If only one uniprot id in "uniprot_id_with_isoform"
    - use "uniprot_id_with_isoform" as final mapping id
- If multiple uniprot ids in "uniprot_id_with_isoform"
    - 2.7.2 if biomart uniprot id is one of the uniprot ids: use corresponding uniprot id (biomart uniprot id + isoform) 
    - 2.7.3 if no results got from above(), try to find uniprot sequence with the same length and levenshtein distance == 1
      - it's possible that multiple uniprot ids are found by sequence length, for now I don't see any of them, but it's better to handle it as well (TODO)
- 2.7.4 if no uniprot id in "uniprot_id_with_isoform", try to find uniprot sequence with the same length and levenshtein distance == 1
    - if multiple uniprot ids were mapping by sequence length, apply "2.7.2" to find a single mapping. If no single uniprot id could be mapped by "2.7.2", leave it empty
- 2.7.5 (Do we want this?)If this transcript was manually curated before (record as "step 2.1.2" and "step 2.2.1" in comment, in file 1.7)
    - what ever mapping result is, use previous mapping id
    - add results to column: final_mapping
##### 2.8. Mark columns that needs manual curation
- "is_matched" = True, "final_mapping" is empty: 
    - case 1: no ensp_id - no need for manual curation
    - case 2: both biomart and sequence matching have nothing returned - no need for manual curation
- "is_matched" = True, "final_mapping" is not empty: 
    - case 1: match successfully by sequence, and have only one uniprot sequence matching - no need for manual curation.
- "is_matched" = False, "final_mapping" is empty: 
    - case 1: couldn't get any mapping - need manual curation
- "is_matched" = False, "final_mapping" is not empty: 
    - case 1: have mapping results by one of the curation steps - no need for manual curation.
- exceptions:
    - have sequence matching(usually it's matching with one of the isoforms, not the canonical sequence), no biomart matching, on ensembl web page it has uniprot matching, but id matching by sequence and ensembl web page are different. e.g. ENST00000293826 ENSP00000293826: [ensembl](http://useast.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000248871;r=17:7549099-7561601;t=ENST00000293826) matches to A0A0A6YY99 and [biomart](http://uswest.ensembl.org/biomart/martview/6b3967a8b19ef08ea5e9574e4fd5972a?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version|hsapiens_gene_ensembl.default.feature_page.uniprotswissprot|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id|hsapiens_gene_ensembl.default.feature_page.ensembl_peptide_id_version&FILTERS=hsapiens_gene_ensembl.default.filters.ensembl_peptide_id."ENSP00000293826"&VISIBLEPANEL=resultspanel) matches with nothing, by sequence it's O43508-2
- add results to column: need_manual_curation
##### 2.9. Get all uniprot ids that couldn't map to a transcript
- Generate a dictionary with all uniprot ids that sussessfully mapping to a transcript - `uniprot_in_transcript_list`
- Generate a uniprot dataframe with all uniprot ids(1.1) - `df_all_uniprot_with_isoform`
    - columns: uniprot_id, gene 
- Go over all available uniprot ids and check if the uniprot id matches with a transcript
    - add results to column: matched 
- Go over all available uniprot ids and check if this uniprot id matches with it's canonical isoform (probably not useful)
    - add results to column: match_with_canonical
- Check if the uniprot id is a cancer gene, return true for cancer genes
    - add results to column: is_cancer_gene

##### 2.10. Find sequence difference for uniprot ids that couldn't match
- Generate a dictionary from transcript dataframe: biomart uniprot id to ensp - `uniprot_to_ensp_dict`.
- Find potential ensp ids for uniprot ids that couldn't match by lookup in uniprot_to_ensp_dict
    - add results to column: ensp_from_biomart
- For each ensp in column "ensp_from_biomart", find sequence difference between ensembl sequence and uniprot sequence. Because from biomart it's corresonding uniprot id doesn't have isoform, so we also check all isoforms associate with this uniprot id. If an isoform of uniprot sequence could match, we print out 
    - add retuslts to column: sequence_difference

### 3. Todo list
- cross reference on uniprot

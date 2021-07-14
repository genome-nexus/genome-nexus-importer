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
##### 2.7. Generate reviewed mapping dictionay(1.7) - `reviewed_mapping_dict`
- key: enst, value: uniprot_id
##### 2.8. Curation
| Mapping | | | | | | | | |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | :------: |
| Number of uniprot ids mapped by sequence | 0 | 0 | 1 | 1 | 1 | multiple | multiple | multiple |
| Number of uniprot ids from Biomart | 0 | 1 | 0 | 1 | 1 | 0 | 1 | 1 |
| Biomart uniprot is one of the uniprot sequence ids (or equals) | no | no | no | equals | no | no | yes (biomart is one of uniprot sequence mapping) | no |
| Find sequence mapping with the same sequence length and levenshtein distance == 1 | yes | yes | no | no | no | no | no | no |
| Return result | Mapping result from same sequence length and levenshtein distance == 1 | Mapping result from same sequence length and levenshtein distance == 1 | Uniprot sequence mapping result | Uniprot sequence mapping result | Uniprot sequence mapping result | No matching (return empty string) | Find biomart uniprot id in uniprot sequence list (have different isoforms) and return the one that has biomart id (with isoform if possible) | No matching (return empty string) |

- 2.7.1 If 0 uniprot ids in "uniprot_id_with_isoform"
    - Find sequence mapping with the same sequence length and levenshtein distance == 1
    - It's still possible that multiple ids will be found by the same sequence length and levenshtein distance == 1, for those cases:

| Mapping |  | | | | |
| ------ | ------ | ------ | ------ | ------ | ------ |
| Number of uniprot ids mapped by sequence | 0 | 1 | multiple | multiple | multiple |
| Number of uniprot ids from Biomart | na | na | 0 | 1 | 1 |
| Biomart uniprot is one of the uniprot sequence ids (or equals) | na | na | no | yes (biomart is one of uniprot seuqnce mapping) | no |
| Return result | No matching (return empty string) | Mapping result from same sequence length and levenshtein distance == 1 | No matching (return empty string) | Find biomart uniprot id in uniprot sequence list (have different isoforms) and return the one that has biomart id (with isoform if possible) | No matching (return empty string) |

- 2.7.2 If 1 uniprot id in "uniprot_id_with_isoform"
    - Return "uniprot_id_with_isoform"
- 2.7.3 If multiple uniprot ids in "uniprot_id_with_isoform"
    - 0 biomart uniprot id:
      - Return empty string
    - biomart uniprot id is one of the ids in "uniprot_id_with_isoform"
      - Return biomart + isoform (if possible)
    - biomart uniprot id is NOT one of the ids in "uniprot_id_with_isoform"
      - Return empty string
- 2.7.4 If no ids found from above, return results from `reviewed_mapping_dict`
- 2.7.5 If still have multiple ids for some reasons, return empty string
- add results to column: final_uniprot_id

### 3. Generate mapping
run Makefile:
```
make uniprot_mapping
```

### 4. Output
- `grchxx_ensemblxx_enst_to_uniprot_mapping_full.txt` is the mapping output with all columns. 
    - Columns: enst_id, ensp_id, ensembl_protein_length, ccds_id, biomart_uniprot_id, uniprot_id_with_isoform, is_matched, final_uniprot_id
- `grchxx_ensemblxx_enst_to_uniprot_mapping_id.txt` only contains enst_id and final_uniprot_id

### 5. Update files
- Homo_sapiens.grch37.pep.all.fa.gz: http://ftp.ensembl.org/pub/grch37/release-104/fasta/homo_sapiens/pep/
- Homo_sapiens.grch38.pep.all.fa.gz: http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/pep/
- uniprot_reviewed.fasta.gz: https://www.uniprot.org/uniprot/?query=+reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score, download FASTA (canonical + isoform)
- reviewed_map_grch37.tsv(no need to update): https://docs.google.com/spreadsheets/d/14PN6RtFq_GTAu8OKNUyUNKJ_fj7OlcEhDjbMdapqqo8/edit#gid=133083142
- reviewed_map_grch38.tsv(no need to update): https://docs.google.com/spreadsheets/d/1slDx9zorUuA-xsmH1i9_6CIGjQB1GIgDlWDU5gw6f9I/edit#gid=1399879714

### 6. Todo list
- cross reference on uniprot

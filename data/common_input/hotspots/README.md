# Transcripts changelog for grch37_ensembl111
The grch37_ensembl111 cancer hotspots transcripts are based on `isoform_overrides_at_mskcc_grch37.txt`
#### Need to change protein position in cancer hotspots:
- NF1: ENST00000358273 → ENST00000356175
  - New transcript doesn’t have 1370-1391 AA, all others are the same. We are good for location less than 1370, for location larger than 1370, the number needs to minus 21
  - R1870 changes to R1849
  - R2450 changes to R2429

#### Although protein sequences are different, no need to change any cancer hotspots
- CRLF2: ENST00000381567 → ENST00000400841
  - New transcript has one additional amino acid at the end. So the S128 remains unchanged.
- KRAS: ENST00000256078 → ENST00000311936
  - New transcript starts to have different protein sequence at 151, we don’t have cancer hotspots after 146, no changes needed
- SMARCA4: ENST00000344626 → ENST00000358026
  - Original transcript doesn’t have 1390 - 1421 AA, we are good for location less than 1390, all hotspots are less than 1390
- TGFBR2: ENST00000359013 → ENST00000295754
  - New transcript doesn’t have 2543-2535, all hotspots are less than 2543
#### Although protein sequences are different, the hotspots matches with new transcript, not matching with original transcript (original transcript was incorrect, these are changed on upsream Cancer Hotspots repo)
- ASXL2: ENST00000336112 → ENST00000435504
- PBRM1: ENST00000296302 → ENST00000394830
- RAC1: ENST00000348035 → ENST00000356142
#### Protein sequences are the same, no need to change hotspots
- DICER1: ENST00000343455 → ENST00000393063 (same protein change, ENST00000343455 has 27 exons, while ENST00000393063 has 28)
- HRAS: ENST00000311189 → ENST00000451590- (same protein change, ENST00000311189 has 6 exons while ENST00000451590 has 5)
- PIK3R1: ENST00000521381 → ENST00000274335 (same protein change, ENST00000521381 has 16 exons, ENST00000274335 has 15)
#### Genes don’t have transcripts in current hotspots file
- DOT1L: nan → ENST00000398665
- H3C13 (prev HIST2H3D): nan → ENST00000331491
- H3C6 (prev HIST1H3E): nan → ENST00000360408
- H3C8 (prev HIST1H3G): nan → ENST00000305910
- PAK5 (prev PAK7): nan → ENST00000353224
- PTPN11: nan → ENST00000351677

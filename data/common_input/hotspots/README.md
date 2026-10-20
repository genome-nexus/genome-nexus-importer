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

# grch38_ensembl111 changes
`grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt` is generated from `grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt` via `scripts/hotspots/liftover_hotspots_to_grch38.py` - a per-residue, alignment-based liftover (see "How to generate from scratch" below for the full pipeline).

Named `hotspots_v2_v3_and_3d.txt` (not `hotspots_v2_and_3d.txt`, the older name still used by some stale/frozen copies elsewhere in repo history) because the file now genuinely combines three cancerhotspots.org sources: the 2D `v2` dataset, the 2D `v3` dataset, and the 3D structural dataset - not just `v2` and 3D as the old name implied.

### Missing cancer hotspots in GRCh38
- WHR1 (STK19) D89 falls in the first 110 residues, which aren't present in the GRCh38 transcript at all (MANE is missing exon 1). 


# Gene name normalization
All gene names in `v2_multi_type_residue.txt`, `3d_hotspots.txt`, `isoform_overrides_at_mskcc_grch37.txt`, and `isoform_overrides_at_mskcc_grch38.txt` were checked against current HGNC symbols and updated where stale (e.g. `HIST1H3E` → `H3C6`, `WHSC1` → `NSD2`, `STK19` → `WHR1`, `TCEB1` → `ELOC`, and ~35 others). This matters because the liftover script normalizes a hotspot's gene symbol to its *current* HGNC name before looking up its GRCh38 transcript.

# How to generate `hotspots_v2_v3_and_3d.txt` from scratch

All commands below are run from `data/common_input/hotspots/`.

**1. Fetch the raw hotspot data** (all three are pinned to a specific commit in the `Makefile`'s `HOTSPOTS_RAW_URL`/`HOTSPOTS_V3_RAW_URL` - re-pin deliberately, not to a floating branch, since upstream branches get renamed/deleted):
```bash
make v2_multi_type_residue.txt v3_multi_type_residue.txt 3d_hotspots.txt
```
Note: `v3_multi_type_residue.txt` is reduced-precision - only genes with hotspots *not already in v2* get added, so delete and re-fetch it (and re-run step 2) whenever the upstream v3 branch moves.

**2. Combine v2 + v3 + 3D into one GRCh37 file:**
```bash
make hotspots_v2_v3_and_3d_grch37.txt
```
This also applies `isoform_overrides_v2.txt`/`isoform_overrides_3d.txt` (fills in transcript IDs that cancerhotspots.org itself left blank for certain genes - check its stderr output for "still missing transcript_id" warnings if you add new source data) and filters out low-confidence hotspots (mostly-truncating or 0-1 supporting mutations) into `removed_hotspots.txt`.

**3. Swap in MSKCC's preferred GRCh37 transcripts:**
```bash
make VERSION=grch37_ensembl111 ../../grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt
```
Runs `scripts/hotspots/overwrite_transcript_and_residue_to_mskcc_isoform.py`. Uses `isoform_overrides_at_mskcc_grch37.txt` (every candidate transcript per gene is tried and verified by alignment, not just renamed blindly - see "GRCh37 override step fixed..." above) plus `mskcc_isoform_residue_changes.txt` as a manual escape hatch for anything the automatic alignment can't resolve on its own. Check `grch37_override_review.txt` after running - any row there means no override candidate's sequence could be verified against that hotspot, so it was conservatively left on its original transcript; it should normally come out empty.

**4. Lift over to GRCh38:**
```bash
make VERSION=grch38_ensembl111 ../../grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt
```
Runs `scripts/hotspots/liftover_hotspots_to_grch38.py` - see "Transcripts verified correct for GRCh38" above for how it picks transcripts. Check `grch38_liftover_dropped.txt` after running; every dropped row has a reason (`deleted`, `grch37_ref_mismatch`, `aa_changed`, `no_grch38_transcript`, `no_sequence`).

**5. (Optional) Sync older builds:** `grch37_ensembl92`, `grch38_ensembl92`, and `grch38_ensembl95` don't get their own separate hotspot generation anymore - hotspot calls don't depend on which Ensembl release annotated the genome, so they just copy whichever of `grch37_ensembl111`/`grch38_ensembl111`'s file matches their genome build:
```bash
make ../../grch37_ensembl92/export/hotspots_v2_v3_and_3d.txt
make ../../grch38_ensembl92/export/hotspots_v2_v3_and_3d.txt
make ../../grch38_ensembl95/export/hotspots_v2_v3_and_3d.txt
```
Each of these re-invokes `make` with `VERSION` pinned to the right upstream build before copying, so they work regardless of whether `VERSION` was set for the outer invocation. (The old REST-API-based `update_hotspots_to_grch38.py` pathway that used to generate the `grch38_ensembl92`/`95` files independently has been removed - it's no longer used by the Makefile, and the script itself has been deleted.)

**Things to note / watch for when regenerating:**
- Both the GRCh37 override step and the GRCh38 liftover step share the same alignment/translation logic (`scripts/hotspots/hotspot_transcript_utils.py`) - if you need to change how a position gets translated between transcripts, fix it there once rather than in either script separately.
- Gene symbols must be current HGNC symbols in `v2_multi_type_residue.txt`, `3d_hotspots.txt`, and both `isoform_overrides_at_mskcc_grch3*.txt` files - an override filed under a stale symbol is silently invisible to the scripts, which normalize to the *current* symbol before looking anything up.
- A gene can legitimately have more than one row in `isoform_overrides_at_mskcc_grch3{7,8}.txt` (e.g. CDKN2A's two reading frames) - both scripts try all of them per hotspot row; don't collapse duplicates expecting only one to be "right."
- `grch37_override_review.txt` and `grch38_liftover_dropped.txt` are both empty-by-default expectations, not informational logs - anything in either one means a hotspot couldn't be automatically verified and needs a human to look at it (same as investigating TGFBR2/RAC1/CDKN2A/AKR7L/SYT7 above).
- `uniprot/input/Homo_sapiens.grch3{7,8}.pep.all.fa.gz` must be release-pinned (same Ensembl release as everything else in this build) - a stale cached FASTA silently produces wrong or missing liftover results with no error.
- The public MySQL mirror (`ensembldb.ensembl.org`, anonymous) is useful for investigating any gene flagged in `grch38_liftover_dropped.txt`: check `transcript.biotype` there before assuming a transcript doesn't exist for a gene - Ensembl sometimes reclassifies a gene's biotype (e.g. AKR7L's `polymorphic_pseudogene` → `protein_coding_LoF`) rather than removing it outright, and the "missing" transcript may just need a by-alignment (not by-raw-position) check against those less obvious candidates.

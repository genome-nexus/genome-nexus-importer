all: ../../$(VERSION)/export/hotspots_v2_and_3d.txt
.PHONY: all

HOTSPOTS_RAW_URL=https://raw.githubusercontent.com/cBioPortal/cancerhotspots/03b7523b2bc26178f8466f41ec943ad97b23c0cc
v2_multi_type_residue.txt:
	curl '$(HOTSPOTS_RAW_URL)/webapp/src/main/resources/data/v2_multi_type_residue.txt' | sed 's/^#//' > $@
3d_hotspots.txt:
	curl '$(HOTSPOTS_RAW_URL)/webapp/src/main/resources/data/3d_hotspots.txt' | sed 's/^#//' > $@

hotspots_v2_and_3d_grch37.txt: v2_multi_type_residue.txt 3d_hotspots.txt isoform_overrides_v2.txt isoform_overrides_3d.txt
	python ../../../scripts/hotspots/combine_2d_3d_add_mutation_type_counts_and_filter.py \
		--removed_hotspots removed_hotspots.txt \
		--override_unassigned_transcript_id_2d_hotspots $(word 3, $^) \
		--override_unassigned_transcript_id_3d_hotspots $(word 4, $^) \
		$(word 1, $^) $(word 2, $^)  > $@

# for grch38: take the grch37 version and port it to grch38:
hotspots_v2_and_3d_grch37.txt_grch38_ported.txt: hotspots_v2_and_3d_grch37.txt
	python ../../../scripts/hotspots/update_hotspots_to_grch38.py \
	--grch38_hugo_and_transcript_id_file_name ../../$(VERSION)/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt \
    --grch38_isoform_override_source mskcc \
	--grch37_hotspots_2d_3d_file_name $<

../../grch37_ensembl92/export/hotspots_v2_and_3d.txt: hotspots_v2_and_3d_grch37.txt
	cp $< $@

ENSEMBL_GENE_IDS=../../$(VERSION)/input/ensembl_biomart_geneids.txt
HGNC_FILE=../hgnc_v2024.10.1.txt
OVERRIDES_FILE=../isoform_overrides_at_mskcc_grch37.txt
RESIDUE_CHANGES_FILE=mskcc_isoform_residue_changes.txt
../../grch37_ensembl111/export/hotspots_v2_and_3d.txt: $(HGNC_FILE) $(OVERRIDES_FILE) $(RESIDUE_CHANGES_FILE) $(ENSEMBL_GENE_IDS) hotspots_v2_and_3d_grch37.txt
	python ../../../scripts/hotspots/overwrite_transcript_and_residue_to_mskcc_isoform.py \
	  --hgnc $(HGNC_FILE) \
	  --overrides $(OVERRIDES_FILE) \
	  --residue_changes $(RESIDUE_CHANGES_FILE) \
	  --biomart $(ENSEMBL_GENE_IDS) \
	  --hotspots hotspots_v2_and_3d_grch37.txt \
	  --out ../../grch37_ensembl111/export/hotspots_v2_and_3d.txt

../../grch38_ensembl92/export/hotspots_v2_and_3d.txt: hotspots_v2_and_3d_grch37.txt_grch38_ported.txt
	cp $< $@

../../grch38_ensembl95/export/hotspots_v2_and_3d.txt: hotspots_v2_and_3d_grch37.txt_grch38_ported.txt
	cp $< $@

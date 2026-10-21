all: ../../$(VERSION)/export/hotspots_v2_v3_and_3d.txt
.PHONY: all

HOTSPOTS_RAW_URL=https://raw.githubusercontent.com/cBioPortal/cancerhotspots/03b7523b2bc26178f8466f41ec943ad97b23c0cc
HOTSPOTS_V3_RAW_URL=https://raw.githubusercontent.com/cBioPortal/cancerhotspots/691b01e5e69bc35203d05b6b97dac81d4c2c907b
v2_multi_type_residue.txt:
	curl '$(HOTSPOTS_RAW_URL)/webapp/src/main/resources/data/v2_multi_type_residue.txt' | sed 's/^#//' > $@
v3_multi_type_residue.txt:
	curl '$(HOTSPOTS_V3_RAW_URL)/webapp/src/main/resources/data/v3_multi_type_residue.txt' | sed 's/^#//' > $@
3d_hotspots.txt:
	curl '$(HOTSPOTS_RAW_URL)/webapp/src/main/resources/data/3d_hotspots.txt' | sed 's/^#//' > $@

hotspots_v2_v3_and_3d_grch37.txt: v2_multi_type_residue.txt 3d_hotspots.txt isoform_overrides_v2.txt isoform_overrides_3d.txt v3_multi_type_residue.txt
	python ../../../scripts/hotspots/combine_2d_3d_add_mutation_type_counts_and_filter.py \
		--removed_hotspots removed_hotspots.txt \
		--override_unassigned_transcript_id_2d_hotspots $(word 3, $^) \
		--override_unassigned_transcript_id_3d_hotspots $(word 4, $^) \
		--hotspots_2d_v3 $(word 5, $^) \
		$(word 1, $^) $(word 2, $^)  > $@

ENSEMBL_GENE_IDS=../../$(VERSION)/input/ensembl_biomart_geneids.txt
HGNC_FILE=../hgnc_complete_set_oct_07_2025.txt
OVERRIDES_FILE=../isoform_overrides_at_mskcc_grch37.txt
RESIDUE_CHANGES_FILE=mskcc_isoform_residue_changes.txt
GRCH37_PEP_FASTA=../../uniprot/input/Homo_sapiens.grch37.pep.all.fa.gz
GRCH38_PEP_FASTA=../../uniprot/input/Homo_sapiens.grch38.pep.all.fa.gz
GRCH38_OVERRIDES_FILE=../isoform_overrides_at_mskcc_grch38.txt
MANE_FILE=../MANE.GRCh38.v1.2.summary.txt

../../grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt: $(HGNC_FILE) $(OVERRIDES_FILE) $(RESIDUE_CHANGES_FILE) $(ENSEMBL_GENE_IDS) $(GRCH37_PEP_FASTA) hotspots_v2_v3_and_3d_grch37.txt
	python ../../../scripts/hotspots/overwrite_transcript_and_residue_to_mskcc_isoform.py \
	  --hgnc $(HGNC_FILE) \
	  --overrides $(OVERRIDES_FILE) \
	  --residue_changes $(RESIDUE_CHANGES_FILE) \
	  --biomart $(ENSEMBL_GENE_IDS) \
	  --pep-fasta $(GRCH37_PEP_FASTA) \
	  --hotspots hotspots_v2_v3_and_3d_grch37.txt \
	  --out ../../grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt \
	  --review-out grch37_override_review.txt

# grch38_ensembl111 and later: per-residue liftover via sequence alignment, instead of the
# whole-gene drop-on-any-difference approach the old update_hotspots_to_grch38.py (REST-API
# based) used. See scripts/hotspots/liftover_hotspots_to_grch38.py for why.
../../grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt: $(HGNC_FILE) $(GRCH38_OVERRIDES_FILE) $(MANE_FILE) $(GRCH37_PEP_FASTA) $(GRCH38_PEP_FASTA) ../../grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt hotspots_v2_v3_and_3d_grch37.txt
	python ../../../scripts/hotspots/liftover_hotspots_to_grch38.py \
	  --grch37-hotspots ../../grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt \
	  --grch37-hotspots-original hotspots_v2_v3_and_3d_grch37.txt \
	  --grch37-pep-fasta $(GRCH37_PEP_FASTA) \
	  --grch38-pep-fasta $(GRCH38_PEP_FASTA) \
	  --grch38-overrides $(GRCH38_OVERRIDES_FILE) \
	  --mane-select $(MANE_FILE) \
	  --hgnc $(HGNC_FILE) \
	  --out ../../grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt \
	  --dropped-out grch38_liftover_dropped.txt

# Older builds (grch37_ensembl92, grch38_ensembl92, grch38_ensembl95) share the same hotspot
# data as the current grch37_ensembl111/grch38_ensembl111 builds rather than being
# regenerated via their own (now-removed) separate pathway - hotspot calls don't depend on
# which Ensembl release annotated the genome, so there's no reason to keep multiple
# diverging copies. Each target below re-invokes make with VERSION explicitly pinned to the
# right upstream build (grch37/grch38_ensembl111) before copying, rather than depending on
# that build's file directly - a direct prerequisite would make this target's own freshness
# check resolve the upstream rule's $(VERSION)-dependent prerequisites (ENSEMBL_GENE_IDS
# etc.) using *this* invocation's VERSION (unset, or whatever the caller passed), which
# breaks regardless of whether a rebuild is actually needed.
../../grch37_ensembl92/export/hotspots_v2_v3_and_3d.txt:
	$(MAKE) VERSION=grch37_ensembl111 ../../grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt
	cp ../../grch37_ensembl111/export/hotspots_v2_v3_and_3d.txt $@

../../grch38_ensembl92/export/hotspots_v2_v3_and_3d.txt:
	$(MAKE) VERSION=grch38_ensembl111 ../../grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt
	cp ../../grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt $@

../../grch38_ensembl95/export/hotspots_v2_v3_and_3d.txt:
	$(MAKE) VERSION=grch38_ensembl111 ../../grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt
	cp ../../grch38_ensembl111/export/hotspots_v2_v3_and_3d.txt $@

all: ../../export/hotspots_v2_and_3d.txt
.PHONY: all

HOTSPOTS_RAW_URL=https://raw.githubusercontent.com/cBioPortal/cancerhotspots/03b7523b2bc26178f8466f41ec943ad97b23c0cc
v2_multi_type_residue.txt:
	curl '$(HOTSPOTS_RAW_URL)/webapp/src/main/resources/data/v2_multi_type_residue.txt' | sed 's/^#//' > $@
3d_hotspots.txt:
	curl '$(HOTSPOTS_RAW_URL)/webapp/src/main/resources/data/3d_hotspots.txt' | sed 's/^#//' > $@

hotspots_v2_and_3d.txt: v2_multi_type_residue.txt 3d_hotspots.txt isoform_overrides_v2.txt isoform_overrides_3d.txt
	python ../../scripts/hotspots/combine_2d_3d_add_mutation_type_counts_and_filter.py \
		--removed_hotspots removed_hotspots.txt \
		--override_unassigned_transcript_id_2d_hotspots $(word 3, $^) \
		--override_unassigned_transcript_id_3d_hotspots $(word 4, $^) \
		$(word 1, $^) $(word 2, $^)  > $@

../../export/hotspots_v2_and_3d.txt: hotspots_v2_and_3d.txt
	cp $< $@

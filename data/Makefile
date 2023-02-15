# The function of this Makefile is to link all the scripts of the Genome Nexus Importer.
# Based on the Ensembl database version, it will download the nessecary datafiles and transforms them into the files imported by Genome Nexus MongoDB.
# Required parameters for this Makefile are VERSION (e.g. grch37_ensembl95), and GFF3_URL (ftp Ensembl GFF3).
#
# Example command:
# make all VERSION=grch38_ensembl95\
# GFF3_URL=ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.gff3.gz\
# QSIZE=1000\
# SPECIES=homo_sapiens
#
# The parameters QSIZE defines the size of the REST queries to Ensembl, and is not mandatory.
# The parameter SPECIES defines the pipeline that is run, which is by default homo_sapiens. 
# Setting SPECIES=mus_musculus will run the mouse data import pipeline which has a slightly different procedure.
# For more information check out the files README.md and docs/setup-genome-nexus-mouse.md
#
# some previously used URLs
# GFF3_URL=ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.gff3.gz
# GFF3_URL=ftp://ftp.ensembl.org/pub/grch37/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
# GFF3_URL=ftp://ftp.ensembl.org/pub/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh38.92.gff3.gz
# GFF3_URL=ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.gff3.gz
# GFF3_URL=ftp://ftp.ensembl.org/pub/release-95/gff3/mus_musculus/Mus_musculus.GRCm38.95.gff3.gz
# ------------------------------------------------------------------------------------------------------------------
SHELL:=/bin/bash

# default species is homo_sapiens unless overridden [mus_musculus]
SPECIES=homo_sapiens

# Specify the URL prefix for vcf2maf to retrieve `isoform_overrides_uniprot.txt`, `isoform_overrides_at_mskcc_grch37.txt` and `isoform_overrides_at_mskcc_grch38.txt`
VCF2MAF_RAW_URL=https://raw.githubusercontent.com/mskcc/vcf2maf/73369616c232b969a2b7b85e9ef9b84949bcc6a4

# dbPTM Base URL:
DB_PTM_URL=http://dbptm.mbc.nctu.edu.tw/download/experiment

# This is the folder to store intermediate files
TMP_DIR=$(VERSION)/tmp

# Ensembl REST query size. Lower this if Ensembl returns Timeout errors.
QSIZE=1000

# Genome build(grch37 or grch38). Use in Uniprot mapping
GENOME_BUILD=$(firstword $(subst _, ,$(VERSION)))

ifeq ($(GENOME_BUILD), grch38)
  MSKCC_ISOFORM_OVERRIDES_FILE_NAME=isoform_overrides_at_mskcc_grch38.txt
  ONCOKB_ISOFORM_OVERRIDES_FILE_NAME=isoform_overrides_oncokb_grch38.txt
else
  MSKCC_ISOFORM_OVERRIDES_FILE_NAME=isoform_overrides_at_mskcc_grch37.txt
  ONCOKB_ISOFORM_OVERRIDES_FILE_NAME=isoform_overrides_oncokb_grch37.txt
endif

# Generic rule to unzip prerequisite files
$(filter-out %.gz, %): %.gz
	gunzip $<

# Generate JSON content for signal-db mutations
signal/export/mutations.json.gz: signal/input/somatic_mutations_by_tumortype_merge.txt signal/input/mutations_cnv_by_tumortype_merge.txt signal/input/biallelic_by_tumortype_merge.txt signal/input/mutations_QCpass_by_tumortype_merge.txt signal/input/signaldb_all_variants_frequencies.txt signal/input/signaldb_msk_expert_review_variants.txt signal/input/signaldb_variants_by_cancertype_summary_statistics.txt
	python ../scripts/transform_signal_db_mutations.py $^ | gzip > $@

# Fetch CCDS files
common_input/CCDS2Sequence.current.txt:
	curl ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/22/CCDS2Sequence.current.txt > $@
common_input/CCDS2UniProtKB.current.txt:
	curl ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/22/CCDS2UniProtKB.current.txt > $@

# Fetch and extract PTM files
ptm/input/Acetylation.txt:
	curl $(DB_PTM_URL)/Acetylation.txt.gz | gunzip -c > $@
ptm/input/Amidation.txt:
	curl $(DB_PTM_URL)/Amidation.txt.gz | gunzip -c > $@
ptm/input/Carbamidation.txt:
	curl $(DB_PTM_URL)/Carbamidation.txt.gz | gunzip -c > $@
ptm/input/Carboxylation.txt:
	curl $(DB_PTM_URL)/Carboxylation.txt.gz | gunzip -c > $@
ptm/input/Citrullination.txt:
	curl $(DB_PTM_URL)/Citrullination.txt.gz | gunzip -c > $@
ptm/input/C-linkedGlycosylation.txt:
	curl $(DB_PTM_URL)/C-linkedGlycosylation.txt.gz | gunzip -c > $@
ptm/input/Crotonylation.txt:
	curl $(DB_PTM_URL)/Crotonylation.txt.gz | gunzip -c > $@
ptm/input/Formylation.txt:
	curl $(DB_PTM_URL)/Formylation.txt.gz | gunzip -c > $@
ptm/input/Gamma-carboxyglutamicAcid.txt:
	curl $(DB_PTM_URL)/Gamma-carboxyglutamicAcid.txt.gz | gunzip -c > $@
ptm/input/Glutarylation.txt:
	curl $(DB_PTM_URL)/Glutarylation.txt.gz | gunzip -c > $@
ptm/input/Glutathionylation.txt:
	curl $(DB_PTM_URL)/Glutathionylation.txt.gz | gunzip -c > $@
ptm/input/GPI-anchor.txt:
	curl $(DB_PTM_URL)/GPI-anchor.txt.gz | gunzip -c > $@
ptm/input/Hydroxylation.txt:
	curl $(DB_PTM_URL)/Hydroxylation.txt.gz | gunzip -c > $@
ptm/input/Lipoylation.txt:
	curl $(DB_PTM_URL)/Lipoylation.txt.gz | gunzip -c > $@
ptm/input/Malonylation.txt:
	curl $(DB_PTM_URL)/Malonylation.txt.gz | gunzip -c > $@
ptm/input/Methylation.txt:
	curl $(DB_PTM_URL)/Methylation.txt.gz | gunzip -c > $@
ptm/input/Myristoylation.txt:
	curl $(DB_PTM_URL)/Myristoylation.txt.gz | gunzip -c > $@
ptm/input/Neddylation.txt:
	curl $(DB_PTM_URL)/Neddylation.txt.gz | gunzip -c > $@
ptm/input/Nitration.txt:
	curl $(DB_PTM_URL)/Nitration.txt.gz | gunzip -c > $@
ptm/input/N-linkedGlycosylation.txt:
	curl $(DB_PTM_URL)/N-linkedGlycosylation.txt.gz | gunzip -c > $@
ptm/input/O-linkedGlycosylation.txt:
	curl $(DB_PTM_URL)/O-linkedGlycosylation.txt.gz | gunzip -c > $@
ptm/input/Oxidation.txt:
	curl $(DB_PTM_URL)/Oxidation.txt.gz | gunzip -c > $@
ptm/input/Palmitoylation.txt:
	curl $(DB_PTM_URL)/Palmitoylation.txt.gz | gunzip -c > $@
ptm/input/Phosphorylation.txt:
	curl $(DB_PTM_URL)/Phosphorylation.txt.gz | gunzip -c > $@
ptm/input/PyrrolidoneCarboxylicAcid.txt:
	curl $(DB_PTM_URL)/PyrrolidoneCarboxylicAcid.txt.gz | gunzip -c > $@
ptm/input/Pyruvate.txt:
	curl $(DB_PTM_URL)/Pyruvate.txt.gz | gunzip -c > $@
ptm/input/S-diacylglycerol.txt:
	curl $(DB_PTM_URL)/S-diacylglycerol.txt.gz | gunzip -c > $@
ptm/input/S-linkedGlycosylation.txt:
	curl $(DB_PTM_URL)/S-linkedGlycosylation.txt.gz | gunzip -c > $@
ptm/input/S-nitrosylation.txt:
	curl $(DB_PTM_URL)/S-nitrosylation.txt.gz | gunzip -c > $@
ptm/input/Succinylation.txt:
	curl $(DB_PTM_URL)/Succinylation.txt.gz | gunzip -c > $@
ptm/input/Sulfation.txt:
	curl $(DB_PTM_URL)/Sulfation.txt.gz | gunzip -c > $@
ptm/input/Sumoylation.txt:
	curl $(DB_PTM_URL)/Sumoylation.txt.gz | gunzip -c > $@
ptm/input/Ubiquitination.txt:
	curl $(DB_PTM_URL)/Ubiquitination.txt.gz | gunzip -c > $@

# Generate single PTM output from multiple input
ptm/export/ptm.json.gz: common_input/CCDS2UniProtKB.current.txt common_input/CCDS2Sequence.current.txt common_input/CCDS2Sequence.override.txt ptm/input
	python ../scripts/add_enst_id_to_ptm.py $^ | gzip > $@

# This will take a while. Only max 1000 transcripts can be retrieved per POST request. Temporary results are saved in
# $VERSION/tmp/transcript_info. This will make it possible to continue the process after the processes crashes, for
# example when the Ensembl API becomes unavailable due to too many requests. If Ensembl REST API returns timeout error,
# lower the query size (QSIZE).
$(TMP_DIR)/ensembl_canonical_data.txt: $(VERSION)/input/ensembl_biomart_geneids.txt
	python ../scripts/download_transcript_info_from_ensembl.py -q $(QSIZE) $< $@

$(TMP_DIR)/ensembl_biomart_transcripts.txt: $(TMP_DIR)/ensembl_canonical_data.txt
	csvcut -tc transcript_stable_id,gene_stable_id,hgnc_symbol,protein_stable_id,protein_length $< | csvsort -c transcript_stable_id,gene_stable_id,protein_stable_id,protein_length | uniq | csvformat -T > $@

# Specify the Ref genome and Ensembl Release in gff3 at top of this file
$(TMP_DIR)/$(SPECIES).gff3.gz:
	curl $(GFF3_URL) > $@

# Generate ensembl transcript info containing exons and UTRs (input & output file in script)
$(TMP_DIR)/ensembl_transcript_info.txt: $(TMP_DIR)/$(SPECIES).gff3.gz
	python ../scripts/transform_gff_to_tsv_for_exon_info_from_ensembl.py $^ $@

# Add HGNC symbols, exons, UTRs, PFAM domains and Uniprot id to Ensembl Transcript
$(TMP_DIR)/ensembl_biomart_transcripts.json.gz: $(TMP_DIR)/ensembl_biomart_transcripts.txt $(TMP_DIR)/ensembl_transcript_info.txt $(VERSION)/input/ensembl_biomart_pfam.txt $(VERSION)/input/ensembl_biomart_refseq.txt $(VERSION)/input/ensembl_biomart_ccds.txt uniprot/export/$(VERSION)_enst_to_uniprot_mapping_id.txt common_input/isoform_overrides_uniprot.txt common_input/$(MSKCC_ISOFORM_OVERRIDES_FILE_NAME) common_input/hgnc_complete_set_20221001.txt
	python ../scripts/add_domains_hugo_ccds_refseq_exon_info_uniprot_to_ensembl_transcript.py $^ $@

# for mouse a specific recipe without overrides
$(TMP_DIR)/ensembl_biomart_transcripts_mouse.json.gz: $(TMP_DIR)/ensembl_biomart_transcripts.txt $(TMP_DIR)/ensembl_transcript_info.txt $(VERSION)/input/ensembl_biomart_pfam.txt $(VERSION)/input/ensembl_biomart_refseq.txt $(VERSION)/input/ensembl_biomart_ccds.txt
	python ../scripts/build_transcript_json_mouse.py $^ $@

# give default/canonical geneid/transcript based on given hugo symbol takes
# about 50m to run (TODO: this can be easily optimized)
# isoform_overrides_genome_nexus.txt is made for genome nexus, others files are generated for vcf2maf
# Please note: we should keep hgnc_complete_set_20221001 in sync with https://github.com/cBioPortal/datahub-study-curation-tools/blob/master/gene-table-update/build-input-for-importer/hgnc_complete_set.txt
# isoform_overrides_oncokb_grch3*.txt is a list of OncoKB transcripts and genes that differ from msk_override, original GRCh38 file could be downloaded here: https://docs.google.com/spreadsheets/d/1ZZt8x0vvhrwL6VLQRzx3YE7XUOvEM-noQl7v6Tg7lCQ/edit#gid=0
$(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_hgnc.txt: $(TMP_DIR)/ensembl_canonical_data.txt common_input/hgnc_complete_set_20221001.txt common_input/isoform_overrides_uniprot.txt common_input/$(MSKCC_ISOFORM_OVERRIDES_FILE_NAME) common_input/isoform_overrides_genome_nexus.txt common_input/$(ONCOKB_ISOFORM_OVERRIDES_FILE_NAME) common_input/ignored_genes.txt
	python ../scripts/make_one_canonical_transcript_per_gene.py $^ $@

# mouse version. A different script is called that set the canonicals based on Ensembl lookup.
$(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_mgi.txt: $(TMP_DIR)/ensembl_canonical_data.txt common_input/mouse/MRK_ENSEMBL.rpt common_input/mouse/MGI_Gene_Model_Coord.rpt
	python ../scripts/make_canonical_transcript_mouse.py $^ $@

uniprot_mapping: $(VERSION)/export/ensembl_biomart_transcripts.json.gz uniprot/input/Homo_sapiens.$(GENOME_BUILD).pep.all.fa.gz uniprot/input/uniprot_reviewed.fasta.gz
	python ../scripts/enst_to_uniprot_mapping.py <(gunzip -c $(word 1, $^)) <(gunzip -c $(word 2, $^)) <(gunzip -c $(word 3, $^)) $(VERSION)
# vcf2maf canonical transcripts
common_input/isoform_overrides_uniprot.txt:
	curl '$(VCF2MAF_RAW_URL)/data/isoform_overrides_uniprot' | sed 's/^#//' > $@

common_input/isoform_overrides_at_mskcc_grch37.txt:
	curl '$(VCF2MAF_RAW_URL)/data/isoform_overrides_at_mskcc' | sed 's/^#//' | sed 's/dmp_refseq_id/refseq_id/' | sed 's/isoform_override/enst_id/g' > $@

common_input/isoform_overrides_at_mskcc_grch38.txt:
	curl '$(VCF2MAF_RAW_URL)/data/isoform_overrides_at_mskcc_grch38' | sed 's/^#//' | sed 's/dmp_refseq_id/refseq_id/' | sed 's/isoform_override/enst_id/g' > $@

# download PfamA.txt if not already available
common_input/pfamA.txt:
	curl "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/database_files/pfamA.txt.gz" | gunzip -c > $@

# download mouse data from MGI
common_input/mouse/MRK_ENSEMBL.rpt:
	wget http://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt -O $@

common_input/mouse/MGI_Gene_Model_Coord.rpt:
	wget http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt -O $@

$(VERSION)/input/ensembl_biomart_ccds.txt $(VERSION)/input/ensembl_biomart_geneids.txt $(VERSION)/input/ensembl_biomart_refseq.txt $(VERSION)/input/ensembl_biomart_pfam.txt:
	Rscript --vanilla ../scripts/retrieve_biomart_tables.R $(SPECIES) $(VERSION)/input/

# download OncoKB cancer genes list
# need to set ONCOKB_TOKEN first by "export ONCOKB_TOKEN="
common_input/oncokb_cancer_genes_list_from_API.json: 
	curl -X 'GET' "https://www.oncokb.org/api/v1/utils/cancerGeneList" -H "accept: application/json" -H "Authorization: Bearer $(ONCOKB_TOKEN)" | python -m json.tool > $@

common_input/oncokb_cancer_genes_list.txt: 
	curl -X 'GET' "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt" -H "accept: text/plain" -H "Authorization: Bearer $(ONCOKB_TOKEN)" > $@


# ClinVar version
# The latest version date number can be found on https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/ and https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
CLINVAR_VERSION=20210308
# download GRCh37 ClinVar VCF file from NCBI
clinvar/input/clinvar_grch37_input.vcf.gz:
	curl "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_${CLINVAR_VERSION}.vcf.gz" > $@
# generate GRCh37 ClinVar tsv file
clinvar/export/clinvar_grch37.txt.gz:
	gunzip -k ../data/clinvar/input/clinvar_grch37_input.vcf.gz && python ../scripts/transform_vcf_to_tsv.py ../data/clinvar/input/clinvar_grch37_input.vcf ../data/clinvar/export/clinvar_grch37.txt --fields chromosome,start_position,end_position,reference_allele,alternate_allele,clinvar_id,clnsig,clnsigconf && gzip ../data/clinvar/export/clinvar_grch37.txt && rm ../data/clinvar/input/clinvar_grch37_input.vcf
# download GRCh38 ClinVar VCF file from NCBI
clinvar/input/clinvar_grch38_input.vcf.gz:
	curl "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_${CLINVAR_VERSION}.vcf.gz" > $@
# generate GRCh38 ClinVar tsv file
clinvar/export/clinvar_grch38.txt.gz:
	gunzip -k ../data/clinvar/input/clinvar_grch38_input.vcf.gz && python ../scripts/transform_vcf_to_tsv.py ../data/clinvar/input/clinvar_grch38_input.vcf ../data/clinvar/export/clinvar_grch38.txt --fields chromosome,start_position,end_position,reference_allele,alternate_allele,clinvar_id,clnsig,clnsigconf && gzip ../data/clinvar/export/clinvar_grch38.txt && rm ../data/clinvar/input/clinvar_grch38_input.vcf

# Annotation sources version. When update any annotation sources version, we should update 'common_input/version_info.txt' and re-generate annotation_version.txt
$(VERSION)/export/annotation_version.txt:
	python ../scripts/annotation_version_file.py common_input/version_info.txt $(VERSION)/export/annotation_version.txt > $@

# create directories if not extistent
dirs: $(VERSION)/input $(VERSION)/export
$(VERSION)/input $(VERSION)/export $(VERSION)/tmp:
	mkdir -p $(VERSION)/input && mkdir -p $(VERSION)/export

# check if variables are defined
input:
	@echo SPECIES
	test $(SPECIES)
	@echo GFF3_URL
	test $(GFF3_URL)
	@echo VERSION
	test $(VERSION)

all: input dirs $(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_hgnc.txt $(TMP_DIR)/ensembl_biomart_transcripts.json.gz common_input/pfamA.txt common_input/hotspots_v2_and_3d.txt
	cp $(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_hgnc.txt $(TMP_DIR)/ensembl_biomart_transcripts.json.gz common_input/pfamA.txt common_input/hotspots_v2_and_3d.txt common_input/oncokb_cancer_genes_list_from_API.json $(VERSION)/export/

# For mouse various steps of this Makefile are not applicable. Therefore this is a mouse specific recipe.
mouse: input dirs $(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_mgi.txt $(TMP_DIR)/ensembl_biomart_transcripts_mouse.json.gz common_input/pfamA.txt
	cp $(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_mgi.txt $(TMP_DIR)/ensembl_biomart_transcripts_mouse.json.gz common_input/pfamA.txt $(VERSION)/export/ && mv $(VERSION)/export/ensembl_biomart_canonical_transcripts_per_mgi.txt $(VERSION)/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt && mv $(VERSION)/export/ensembl_biomart_transcripts_mouse.json.gz $(VERSION)/export/ensembl_biomart_transcripts.json.gz

.PHONY: all mouse
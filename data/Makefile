# Specify version, consistent with the folder in `data/` created for this version
VERSION=grch38_ensembl95
# GFF3_URL=ftp://ftp.ensembl.org/pub/grch37/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
# GFF3_URL=ftp://ftp.ensembl.org/pub/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh38.92.gff3.gz
GFF3_URL=ftp://ftp.ensembl.org/pub/release-95/gff3/homo_sapiens/Homo_sapiens.GRCh38.95.gff3.gz

# Specify the URL prefix for vcf2maf to retrieve `isoform_overrides_uniprot.txt` and `isoform_overrides_at_mskcc.txt`
VCF2MAF_RAW_URL=https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.16

# dbPTM Base URL:
DB_PTM_URL=http://dbptm.mbc.nctu.edu.tw/download/experiment

# This is the folder to store intermediate files
TMP_DIR=$(VERSION)/tmp

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
ptm/export/ptm.json.gz: common_input/CCDS2UniProtKB.current.txt common_input/CCDS2Sequence.current.txt ptm/input
	python3 ../scripts/add_enst_id_to_ptm.py $^ | gzip > $@

# This will take a while. Only 1000 transcripts can be retrieved per POST request. Temporary results are saved in
# $VERSION/tmp/transcript_info . This will make it possible to continue the process after the processes crashes, for
# example when the Ensembl API becomes unavailable due to too many requests.
$(TMP_DIR)/ensembl_biomart_geneids.transcript_info.txt: $(VERSION)/input/ensembl_biomart_geneids.txt
	python3 ../scripts/download_transcript_info_from_ensembl.py $< $@

$(TMP_DIR)/ensembl_biomart_transcripts.txt: $(TMP_DIR)/ensembl_biomart_geneids.transcript_info.txt
	csvcut -tc transcript_stable_id,gene_stable_id,hgnc_symbol,protein_stable_id,protein_length $< | csvsort -c transcript_stable_id,gene_stable_id,protein_stable_id,protein_length | uniq | csvformat -T > $@

# Specify the Ref genome and Ensembl Release in gff3 at top of this file
$(TMP_DIR)/Homo_sapiens.gff3.gz:
	curl $(GFF3_URL) > $@

# Generate ensembl transcript info containing exons and UTRs (input & output file in script)
$(TMP_DIR)/ensembl_transcript_info.txt: $(TMP_DIR)/Homo_sapiens.gff3.gz
	python3 ../scripts/transform_gff_to_tsv_for_exon_info_from_ensembl.py $^ $@

# Add HGNC symbols, exons, UTRs and PFAM domains to Ensembl Transcript
$(TMP_DIR)/ensembl_biomart_transcripts.json.gz: $(TMP_DIR)/ensembl_biomart_transcripts.txt $(TMP_DIR)/ensembl_transcript_info.txt $(VERSION)/input/ensembl_biomart_pfam.txt $(VERSION)/input/ensembl_biomart_refseq.txt $(VERSION)/input/ensembl_biomart_ccds.txt common_input/isoform_overrides_uniprot.txt common_input/isoform_overrides_at_mskcc.txt
	python3 ../scripts/add_domains_hugo_ccds_refseq_exon_info_to_ensembl_transcript.py $^ $@

# give default/canonical geneid/transcript based on given hugo symbol takes
# about 50m to run (TODO: this can be easily optimized)
# isoform_overrides_genome_nexus.txt is made for genome nexus, others files are generated for vcf2maf
$(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_hgnc.txt: $(TMP_DIR)/ensembl_biomart_geneids.transcript_info.txt common_input/hgnc_complete_set_20190122.txt common_input/isoform_overrides_uniprot.txt common_input/isoform_overrides_at_mskcc.txt common_input/isoform_overrides_genome_nexus.txt
	python3 ../scripts/make_one_canonical_transcript_per_gene.py $^ $@

# vcf2maf canonical transcripts
common_input/isoform_overrides_uniprot.txt:
	curl '$(VCF2MAF_RAW_URL)/data/isoform_overrides_uniprot_from_biomart_91' | sed 's/^#//' > $@

common_input/isoform_overrides_at_mskcc.txt:
	curl '$(VCF2MAF_RAW_URL)/data/isoform_overrides_at_mskcc' | sed 's/^#//' | sed 's/dmp_refseq_id/refseq_id/' | sed 's/isoform_override/enst_id/g' > $@

all: $(TMP_DIR)/ensembl_biomart_canonical_transcripts_per_hgnc.txt $(TMP_DIR)/ensembl_biomart_transcripts.json.gz common_input/pfamA.txt common_input/hotspots_v2_and_3d.txt
	cp $^ $(VERSION)/export/

.PHONY: all

#!/bin/bash
set -e

# GRCh37 ClinVar
curl "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" | gunzip -c > ../data/clinvar/input/clinvar_grch37.vcf
python transform_clinvar_vcf_to_txt.py ../data/clinvar/input/clinvar_grch37.vcf ../data/clinvar/input/clinvar_grch37.txt

gzip ../data/clinvar/input/clinvar_grch37.txt
rm ../data/clinvar/input/clinvar_grch37.vcf


# GRCh38 ClinVar
curl "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" | gunzip -c > ../data/clinvar/input/clinvar_grch38.vcf
python transform_clinvar_vcf_to_txt.py ../data/clinvar/input/clinvar_grch38.vcf ../data/clinvar/input/clinvar_grch38.txt

gzip ../data/clinvar/input/clinvar_grch38.txt
rm ../data/clinvar/input/clinvar_grch38.vcf
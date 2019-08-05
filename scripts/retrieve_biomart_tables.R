#!/usr/bin/env Rscript --vanilla

# This script will download various data tables from Ensembl BioMart,
# which serves as input for the Genome Nexus Import pipeline.

library(biomaRt)

args = commandArgs(trailingOnly=TRUE)

# set working dir to the correct genome/version input dir
species <- args[1]
path <- args[2]

stopifnot(species %in% c('homo_sapiens', 'mus_musculus'))

# set species to short name for biomart
species <- ifelse(species=='homo_sapiens', 'hsapiens', 'mmusculus')

# select mart
# listEnsembl()

ensembl <- useMart(biomart='ensembl', dataset=paste0(species, '_gene_ensembl'))

# list datasets and attributes
# listDatasets(ensembl)
# listAttributes(ensembl)

# if species is mouse, get IDs from MGI
attributes <- c('ensembl_gene_id', 'ensembl_transcript_id', 'hgnc_symbol', 'hgnc_id')
if (species=='mmusculus') attributes[3:4] <- c('mgi_symbol', 'mgi_id')

genes <- getBM(attributes=attributes, mart=ensembl)
pfam <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'pfam', 'pfam_start', 'pfam_end'), mart=ensembl)
refseq <- getBM(attributes=c('ensembl_transcript_id', 'refseq_mrna'), mart=ensembl)
ccds <- getBM(attributes=c('ensembl_transcript_id', 'ccds'), mart=ensembl)

write.table(pfam, paste0(path, '/ensembl_biomart_pfam.txt'), na='', sep='\t', quote=FALSE, row.names=FALSE, 
  col.names=c('Gene stable ID', 'Transcript stable ID', 'Gene name', 'Pfam domain ID', 'Pfam domain start', 'Pfam domain end'))
write.table(genes, paste0(path, '/ensembl_biomart_geneids.txt'), na='', sep='\t', quote=FALSE,  row.names=FALSE, 
  col.names=c('Gene stable ID', 'Transcript stable ID', 'HGNC symbol', 'HGNC ID'))
write.table(refseq, paste0(path, '/ensembl_biomart_refseq.txt'), na='', sep='\t', quote=FALSE,  row.names=FALSE, 
  col.names=c('Transcript stable ID', 'RefSeq mRNA ID'))
write.table(ccds, paste0(path, '/ensembl_biomart_ccds.txt'), na='', sep='\t', quote=FALSE,  row.names=FALSE, 
  col.names=c('Transcript stable ID', 'CCDS ID'))
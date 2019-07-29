# Setting up a Genome Nexus for mouse

Genome Nexus is adapted to use mouse genomes. Loading mouse data requires using an alternative pipeline. The script will take as a starting point several Ensembl BioMart tables, and will download additional information from Ensembl and the [mouse genome informatics website (MGI)](http://www.informatics.jax.org/downloads/reports/index.html). The mouse pipeline will format the files to resemble the human data files, but then of course with mouse information.

There are several human-specific annotation tables, such as the cancer hotspots and canonical transcript overrides by MSK or Uniprot. Those will not be used. For each mouse gene it will query all transcripts through the REST webservice, and determine the gene's canonical transcript based on Ensembl annotation, and alternatively on protein length.

## Downloading Ensembl BioMart data
Ensembl data must be downloaded from BioMart as a starting point for data annotation. The data may be downloaded using the described [manual steps](https://github.com/genome-nexus/genome-nexus-importer#2-manually-download-ensembl-biomart-files), or using the R script `scripts/retrieve_biomart_tables.R`.

First make sure the directories are created and environment variables are set. For this example we will use Ensembl 95.
```
cd ~/git/genome-nexus-importer
export REF_ENSEMBL_GENOME=grcm38_ensembl95
export SPECIES=mus_musculus

mkdir -p data/${REF_ENSEMBL_GENOME}/input
mkdir -p data/${REF_ENSEMBL_GENOME}/export
mkdir -p data/${REF_ENSEMBL_GENOME}/tmp
```

##### If downloading the annotation tables through the Ensembl BioMart website:
- Download the appropriate files [as described](https://github.com/genome-nexus/genome-nexus-importer#2-manually-download-ensembl-biomart-files), <b>however</b>:
-- Select the columns `mgi_symbol` and `mgi_id` instead of `hgnc_symbol` and `hgnc_id`
-- Rename the columns `hgnc_symbol` and `hgnc_id` to `mgi_symbol` and `mgi_id`
##### If using the script `scripts/retrieve_biomart_tables.R` to download annotation tables:
- Open the R script in your favorite editor
- Make sure the `species` variable is assigned `"mmusculus"`
- Adjust the `setwd` command to appropriate `data/<refgenome_ensemblversion>/input` location
- Run the commands in an R shell

## Running the pipeline
A mouse-specific recipe is implemented in the pipeline (Makefile) that is responsible for the data transformation.

Run the pipeline:
```
cd data/
make mouse \
  VERSION=grcm38_ensembl95 \
  GFF3_URL=ftp://ftp.ensembl.org/pub/release-95/gff3/mus_musculus/Mus_musculus.GRCm38.95.gff3.gz \
  SPECIES=mus_musculus \
  QSIZE=100
```

The QSIZE parameter was introduced because it takes significantly longer to query the Ensembl database for non-human species. When using a query size of 1000 genes at a time (default), timeout errors may occur. In the case of a timeout error, the QSIZE can be lowered.
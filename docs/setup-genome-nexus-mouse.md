# Setting up a Genome Nexus for mouse

Genome Nexus is adapted to use mouse genomes. Loading mouse data requires using an alternative pipeline. The script will take as a starting point several Ensembl BioMart tables, and will download additional information from Ensembl and the [mouse genome informatics website (MGI)](http://www.informatics.jax.org/downloads/reports/index.html). The mouse pipeline will format the files to resemble the human data files, but then of course with mouse information.

There are several human-specific annotation tables, such as the cancer hotspots and canonical transcript overrides by MSK or Uniprot. Those will not be used. For each mouse gene it will query all transcripts through the REST webservice, and determine the gene's canonical transcript based on Ensembl annotation, and alternatively on protein length.

## Running the pipeline
A mouse-specific recipe is implemented in the pipeline (Makefile) that is responsible for the data transformation.

Now call the `mouse` recipe and make sure to override the `SPECIES` parameter for mus_musculus.
```bash
cd data/
make mouse \
  VERSION=grcm38_ensembl95 \
  GFF3_URL=ftp://ftp.ensembl.org/pub/release-95/gff3/mus_musculus/Mus_musculus.GRCm38.95.gff3.gz \
  SPECIES=mus_musculus \
  QSIZE=100
```

The QSIZE parameter was introduced because it takes significantly longer to query the Ensembl database for non-human species. When using a query size of 1000 genes at a time (default), timeout errors may occur. In the case of a timeout error, the QSIZE can be lowered.

## Directly importing data into MongoDB
If you want to import data directly, as described [https://github.com/genome-nexus/genome-nexus-importer#using-the-mongo-database](here), please make sure to also set the `SPECIES` env variable to mus_musculus.
```bash
MONGO_URI="mongodb://127.0.0.1:27017/annotator"
REF_ENSEMBL_VERSION="grcm38_ensembl95"
SPECIES="mus_musculus"
./scripts/import_mongo.sh
```
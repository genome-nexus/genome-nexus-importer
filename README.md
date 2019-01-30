# Genome Nexus Data Importer
This repo contains data to set up the mongo
database for genome nexus, as well as scripts to generate new data. When setting
up a database container for Genome Nexus, it is not required to generate new
data. Data for several reference genomes and Ensembl releases is available in
the `data` folder.

## Using the mongo database

### Using docker container
There's a mongo docker container that has all the data imported. You can use
the docker compose file in the genome nexus repo itself to start both the web
app and the database: [genome
nexus](https://github.com/genome-nexus/genome-nexus).

### Directly import to mongo database
Run the script [scripts/import_mongo.sh](scripts/import_mongo.sh) to import processed data files into a running 
database. When running this script, please specify:
- `MONGO_URI`: Mongo database address, for example `mongodb://127.0.0.1:27017/annotator`.
- `REF_ENSEMBL_VERSION`: Reference Genome and Ensembl Release, for example `grch37_ensembl92` or `grch38_ensembl92`. Files are imported from 
`data/<refgenome_ensemblversion>/export/`.

Example:
```bash
MONGO_URI="mongodb://127.0.0.1:27017/annotator"
REF_ENSEMBL_VERSION="grch37_ensembl92"
./scripts/import_mongo.sh
```

## Generating data
This repository contains a pipeline to retrieve data for a specified reference genome and Ensembl build. Generated data is saved in:
```
data/
```

To generated data for a different reference genome and Ensembl Release, follow the instructions below.
1. Create folders in Git repository
2. Manually download Ensembl BioMart files
3. Run data transformation pipeline
4. Verify data
5. Commit, Push and create Pull Request

### 1. Create folders in Git repository
- `data/<refgenome_ensemblversion>/input`
- `data/<refgenome_ensemblversion>/export`
- `data/<refgenome_ensemblversion>/tmp`

The `input` and `export` folders are tracked by Git, while the `tmp` folder contains the intermediate files and is not 
tracked by Git.

### 2. Manually download Ensembl BioMart files

- GRCh37 / hg19: https://grch37.ensembl.org/biomart
- GRCh38 / hg38: https://www.ensembl.org/biomart

http://www.ensembl.org/Help/ArchiveList contains archived BioMart releases.

##### PFAM endpoint
Ensembl Biomart file is required by the PFAM endpoint. In order to download this file
follow these steps:

1. Go to the Ensemble BioMart website of the reference genome and Ensembl Release of interest.
2. Select `Ensemble Genes` from the `Database` dropdown menu.
3. Select `Human Genes` from the `Dataset` dropdown menu.
4. Click on `Attributes`, and select these ones:
Gene stable ID, Transcript stable Id, Protein stable Id, Gene name, Pfam domain ID, Pfam domain start, Pfam domain end.
5. Click on `Results`, and export all results to a `TSV` file.
6. Save the downloaded file as `ensembl_biomart_pfam.txt` in `data/<refgenome_ensemblversion>/input`.

##### Ensembl endpoint 
1. Go to the Ensemble BioMart website of the reference genome and Ensembl Release of interest.
2. Select `Ensemble Genes` from the `Database` dropdown menu.
3. Select `Human Genes` from the `Dataset` dropdown menu.
4. Click on `Attributes`, and select these ones:
Gene stable ID, Transcript stable Id, HGNC Symbol, HGNC ID
5. Click on `Results`, and export all results to a `TSV` file.
6. Save the downloaded file as `ensembl_biomart_geneids.txt` in `data/<refgenome_ensemblversion>/input`.

##### Ensembl endpoint RefSeq identifiers
1. Go to the Ensemble BioMart website of the reference genome and Ensembl Release of interest.
2. Select `Ensemble Genes` from the `Database` dropdown menu.
3. Select `Human Genes` from the `Dataset` dropdown menu.
4. Click on `Attributes`, and select these ones:
Transcript stable Id, RefSeq mRNA ID
5. Click on `Results`, and export all results to a `TSV` file.
6. Save the downloaded file as `ensembl_biomart_refseq.txt` in `data/<refgenome_ensemblversion>/input`.

##### Ensembl endpoint CCDS identifiers
1. Go to the Ensemble BioMart website of the reference genome and Ensembl Release of interest.
2. Select `Ensemble Genes` from the `Database` dropdown menu.
3. Select `Human Genes` from the `Dataset` dropdown menu.
4. Click on `Attributes`, and select these ones:
Transcript stable Id, CCDS ID
5. Click on `Results`, and export all results to a `TSV` file.
6. Save the downloaded file as `ensembl_biomart_ccds.txt` in `data/<refgenome_ensemblversion>/input`.

### Run data transformation pipeline
To run the pipeline that transforms all intermediate data, run the following. This takes _roughly two hours_ to complete.

```bash
cd data
make all
```

Please make sure you have the requirements in `requirements.txt` installed. If the pipeline crashes somewhere, for example when the Ensembl REST API is down, sometimes an empty file is created. To continue the pipeline, remove the empty file and run `make all` again.

##### Canonical transcripts
During this process, every transcript in `ensembl_biomart_geneids.txt` is assessed to be either canonical or not, by 
querying the Ensembl REST API. This takes a while, because only 1000 transcripts can be queried at a time. Progress can 
be viewed by inspecting the temporary files created in  `data/<refgenome_ensemblversion>/tmp/transcript_info`. Gene 
source file `ensembl_biomart_geneids.txt` contains about _224596_ transcripts, so the pipeline will save about _225_ 
of these files.

### Verify data
To verify the pipeline ran data for the correct reference genome, you can verify the exon coordinates in
`export/ensembl_biomart_transcripts.json.gz`. Select an Ensembl Exon ID, query it on Ensembl GRCh38 or GRCh37, select
the gene, select transcript, and select 'Exons'. This will display all the exons of the transcript and their genomic
coordinates.

### Commit, Push and create Pull Request
When new data has been created, create a PR to Genome-Nexus to add this data to the master branch.

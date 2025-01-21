# Genome Nexus Data Importer
This repo contains data to set up the mongo
database for genome nexus, as well as scripts to generate new data. When setting
up a database container for Genome Nexus, it is not required to generate new
data. Data for several reference genomes and Ensembl releases is available in
the `data` folder.

The parent image is Bitnami's mongodb image, found [here](https://hub.docker.com/r/bitnami/mongodb/). 

Current installed MongoDB version: **4.0.12**

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
export MONGO_URI="mongodb://127.0.0.1:27017/annotator"
export REF_ENSEMBL_VERSION="grch37_ensembl92"
./scripts/import_mongo.sh
```

## Generating data
This repository contains a pipeline to retrieve data for a specified reference genome and Ensembl build. Generated data is saved in:
```
data/
```

To generated data for a different reference genome and Ensembl Release, follow the instructions below.
The main driver of the data loading pipeline is the Makefile found in `data/`. It will download relevant tables from Ensembl BioMart, Pfam, HGNC, and transforms to the proper format for MongoDB.

The Makefile will create and fill the directories
- `data/<refgenome_ensemblversion>/input`: Input tables retrieved from Ensembl Biomart
- `data/<refgenome_ensemblversion>/export`: Pipeline output, used by MongoDB.
- `data/<refgenome_ensemblversion>/tmp`: Temp files.

The `input` and `export` folders are tracked by Git, while the `tmp` folder contains the intermediate files and is not 
tracked by Git.

#### Dependencies
There are few dependencies of this pipeline for either python or R.
The python dependencies can be installed from the file `requirements.txt`:
```bash
cd scripts
pip install -r requirements.txt
```
For R there is only the dependency on the biomaRt library.
```bash
R -e "source('https://bioconductor.org/biocLite.R'); biocLite('biomaRt')"
```

#### Updating versions
##### data/<ref_genome_ensembl_version>/input
Delete old input files, new files will be generated while running `Makefile`. 
If you want to update data to a specific verison, check [Ensembl Archieves](https://useast.ensembl.org/info/website/archives/index.html) to find the corresponding url and replace `host_url` in the `retrieve_biomart_tables.R`. Please note only GRCh38 supports versioning url, GRCh37 can only point to `https://grch37.ensembl.org/` which uses the latest version.

##### Onocokb isoform overrides
Delete old files, new files will be downloaded by `download_oncokb_isoform_overrides.py` while running `Makefile`. Please make sure to set `ONCOKB_VERSION=` in `Makefile`. Also this version need to be updated in `/data/common_input/version_info.txt`.

##### Clinvar
Delete old input and output files, new files can be genrated by manually run:
```
make clinvar/input/clinvar_grch37_input.vcf.gz
make clinvar/export/clinvar_grch37.txt.gz
make clinvar/input/clinvar_grch38_input.vcf.gz
make clinvar/export/clinvar_grch38.txt.gz

``` 
The latest version date number can be found on [https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/) and [https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/), please make sure to provide correct `CLINVAR_VERSION=` in `Makefile` (for exampel: `CLINVAR_VERSION=20250106`). This new version number also needs to be added in `/data/common_input/version_info.txt`.

##### HGNC
`/data/common_inpu/hgnc_xxxx.txt` should keep in sync with [https://github.com/cBioPortal/datahub-study-curation-tools/blob/master/gene-table-update/build-input-for-importer/hgnc_complete_set.txt](https://github.com/cBioPortal/datahub-study-curation-tools/blob/master/gene-table-update/build-input-for-importer/hgnc_complete_set.txt). When updating HGNC file, download new file from the github, rename to `hgnc_xxxx.txt` (`xxxx` indicates the version, for example v2024.10.1). Please also update this new file name in `Makefile` where the files is used, and update the version in `/data/common_input/version_info.txt`.

#### Running
Run the import pipeline using the command below. This will take a few hours to complete.
```bash
cd data
make all \
VERSION=grch37_ensembl92 \
GFF3_URL=ftp://ftp.ensembl.org/pub/grch37/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
```

To change the reference genome to build the data files for, change the
`VERSION` and `GFF3_URL` variables accordingly (examples are in the Makefile).

If the pipeline crashes, for example when the Ensembl REST API is down, sometimes an empty file is created. To continue the pipeline, remove the empty file and run `make all` again.

Additionally, mouse data can be processed to build a database for mouse. This is described [here](docs/setup-genome-nexus-mouse.md).

##### Canonical transcripts
During this process, every transcript in `data/<refgenome_ensemblversion>/input/ensembl_biomart_geneids.txt` is assessed to be either canonical or not, by querying the Ensembl REST API. This takes a while, because a maximum of 1000 transcripts can be queried at a time. Progress can be viewed by inspecting the temporary files created in  `data/<refgenome_ensemblversion>/tmp/transcript_info`. Gene source file `ensembl_biomart_geneids.txt` contains about _224596_ transcripts, so the pipeline will save about _225_ of these files.

When the REST API is slow for whatever reason, the server can return a timeout error. When that happens, the `QSIZE` parameter can be used to reduce query size (e.g. to 100 transcripts at a time).
```
make all \
VERSION=grch37_ensembl92 \
GFF3_URL=ftp://ftp.ensembl.org/pub/grch37/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz \
QSIZE=100
```

### Verify data
To verify the pipeline ran data for the correct reference genome, you can verify the exon coordinates in
`export/ensembl_biomart_transcripts.json.gz`. Select an Ensembl Exon ID, query it on Ensembl GRCh38 or GRCh37, select
the gene, select transcript, and select 'Exons'. This will display all the exons of the transcript and their genomic
coordinates.

### Update hotspots

The hotspots files are updated by a slightly different `make` process. To update the files,
run the following command (replacing `VERSION=grch38_ensembl95` with your build of interest):
```
cd data/common_input/hotspots
make all VERSION=grch38_ensembl95
```

### Commit, Push and create Pull Request
When new data has been created, create a PR to Genome-Nexus to add this data to the master branch.

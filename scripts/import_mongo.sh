#!/bin/bash
set -e

# Set default if ENV variables are not set
MONGO_URI=${MONGO_URI:-"mongodb://127.0.0.1:27017/annotator"}
REF_ENSEMBL_VERSION=${REF_ENSEMBL_VERSION:-"grch37_ensembl92"}
SPECIES=${SPECIES:-"homo_sapiens"}
AWS_S3_BUCKET_NAME=${AWS_S3_BUCKET_NAME:-"genome-nexus-mongo-dump"}
AWS_S3_BUCKET_DUMP_DIR=${AWS_S3_BUCKET_DUMP_DIR:-"dump/genome-nexus"}
SKIP_MONGODUMP_RESTORE=${SKIP_MONGODUMP_RESTORE:-"0"}

echo "MONGO_URI:" ${MONGO_URI}
echo "REF_ENSEMBL_VERSION:" ${REF_ENSEMBL_VERSION}
echo "SPECIES:" ${SPECIES}
echo "AWS_BUCKET_NAME:" ${AWS_S3_BUCKET_NAME}
echo "AWS_S3_BUCKET_DUMP_DIR:" ${AWS_S3_BUCKET_DUMP_DIR}

# Get name of latest mongodump, extract contents, and run mongorestore
#  mongo dump directory path is set in Dockerfile
if [ -z "${AWS_BUCKET_NAME}" ] || [ -z "${AWS_S3_BUCKET_DUMP_DIR}" ] ; then
    SKIP_MONGODUMP_RESTORE="1"
fi

if [ "${SKIP_MONGODUMP_RESTORE}" == "0" ] ; then
    cd /mongodump
    wget https://${AWS_S3_BUCKET_NAME}.s3.amazonaws.com/${AWS_S3_BUCKET_DUMP_DIR}/latest_mongo_dump.txt
    LATEST_MONGODUMP_FILENAME=$(cat latest_mongo_dump.txt)
    wget https://${AWS_S3_BUCKET_NAME}.s3.amazonaws.com/${AWS_S3_BUCKET_DUMP_DIR}/$LATEST_MONGODUMP_FILENAME
    tar -xzvf $LATEST_MONGODUMP_FILENAME
    mongorestore --uri ${MONGO_URI} --verbose dump
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

import() {
    collection=$1
    file=$2
    extraoptions=$3

    mongoimport --uri ${MONGO_URI} --drop --collection $collection $extraoptions --file $file
}

if [[ ! -d "${DIR}/../data/${REF_ENSEMBL_VERSION}" ]]; then
	echo "Can't find directory for given reference genome and ensembl release: ${DIR}/../data/"${REF_ENSEMBL_VERSION}
	exit
fi

##TODO: get this config from some JSON file, so both bash and Java can read it
import ensembl.biomart_transcripts <(gunzip -c ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/ensembl_biomart_transcripts.json.gz) '--type json'
import ensembl.canonical_transcript_per_hgnc ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt '--type tsv --headerline'
import pfam.domain ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/pfamA.txt '--type tsv --headerline'
import ptm.experimental <(gunzip -c ${DIR}/../data/ptm/export/ptm.json.gz) '--type json'

# Exit if species is not homo_sapiens. Next import steps are human specific
[[ "$SPECIES" == "homo_sapiens" ]] || exit 0

echo "Executing human-specific import steps"

import hotspot.mutation ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/hotspots_v2_and_3d.txt '--type tsv --headerline --mode upsert --upsertFields hugo_symbol,residue,type,tumor_count'
import signal.mutation <(gunzip -c ${DIR}/../data/signal/export/mutations.json.gz) '--type json'

# import oncokb cancer genes list
import oncokb.gene ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/oncokb_cancer_genes_list_from_API.json '--type json --jsonArray'


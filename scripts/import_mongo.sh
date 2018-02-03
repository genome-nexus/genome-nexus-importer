#!/bin/bash
set -e

# use given first argument as MONGO_URI if starting with mongodb://
if [[ $1 == mongodb://* ]]
then
    MONGO_URI=$1
else
    MONGO_URI=mongodb://127.0.0.1:27017/annotator
fi
echo ${MONGO_URI}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

import() {
    collection=$1
    file=$2
    extraoptions=$3

    mongoimport --uri ${MONGO_URI} --drop --collection $collection $extraoptions --file $file
}

#TODO: get this config from some JSON file, so both bash and Java can read it
import ensembl.biomart_transcripts ${DIR}/../export/ensembl_biomart_transcripts.json '--type json'
import ensembl.canonical_transcript_per_hgnc ${DIR}/../export/ensembl_biomart_canonical_transcripts_per_hgnc.txt '--type tsv --headerline'
import pfam.domain ${DIR}/../export/pfamA.txt '--type tsv --headerline'

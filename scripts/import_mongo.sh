#!/bin/bash
set -e

# Set default if ENV variables are not set
MONGO_URI=${MONGO_URI:-"mongodb://127.0.0.1:27017/annotator"}
REF_ENSEMBL_VERSION=${REF_ENSEMBL_VERSION:-"grch37_ensembl92"}
SPECIES=${SPECIES:-"homo_sapiens"}
MUTATIONASSESSOR=${MUTATIONASSESSOR:-"false"}

echo "MONGO_URI:" ${MONGO_URI}
echo "REF_ENSEMBL_VERSION:" ${REF_ENSEMBL_VERSION}
echo "SPECIES:" ${SPECIES}
echo "MUTATIONASSESSOR:" ${MUTATIONASSESSOR}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

import() {
    collection=$1
    file=$2
    extraoptions=$3

    mongoimport --uri ${MONGO_URI} --collection $collection $extraoptions --file $file
}

if [[ ! -d "${DIR}/../data/${REF_ENSEMBL_VERSION}" ]]; then
	echo "Can't find directory for given reference genome and ensembl release: ${DIR}/../data/"${REF_ENSEMBL_VERSION}
	exit
fi

##TODO: get this config from some JSON file, so both bash and Java can read it
import ensembl.biomart_transcripts <(gunzip -c ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/ensembl_biomart_transcripts.json.gz) '--drop --type json'
import ensembl.canonical_transcript_per_hgnc ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt '--drop --type tsv --headerline'
import pfam.domain ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/pfamA.txt '--drop --type tsv --headerline'
import ptm.experimental <(gunzip -c ${DIR}/../data/ptm/export/ptm.json.gz) '--drop --type json'

# Exit if species is not homo_sapiens. Next import steps are human specific
[[ "$SPECIES" == "homo_sapiens" ]] || exit 0

echo "Executing human-specific import steps"

import hotspot.mutation ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/hotspots_v2_and_3d.txt '--drop --type tsv --headerline --mode upsert --upsertFields hugo_symbol,residue,type,tumor_count'
import signal.mutation <(gunzip -c ${DIR}/../data/signal/export/mutations.json.gz) '--drop --type json'

# import oncokb cancer genes list
import oncokb.gene ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/oncokb_cancer_genes_list_from_API.json '--drop --type json --jsonArray'

# import ClinVar
if [[ ${REF_ENSEMBL_VERSION} == *"grch37"* ]]; then
    import clinvar.mutation <(gunzip -c ${DIR}/../data/clinvar/export/clinvar_grch37.txt.gz) '--drop --type tsv --headerline --columnsHaveTypes --parseGrace autoCast'
elif [[ ${REF_ENSEMBL_VERSION} == *"grch38"* ]]; then
    import clinvar.mutation <(gunzip -c ${DIR}/../data/clinvar/export/clinvar_grch38.txt.gz) '--drop --type tsv --headerline --columnsHaveTypes --parseGrace autoCast'
fi

# import mutation assessor
declare -a mutation_assessor_files=(
    "https://genome-nexus-static-data.s3.amazonaws.com/mutationassessor_v4_1.tsv.gz"
    "https://genome-nexus-static-data.s3.amazonaws.com/mutationassessor_v4_2.tsv.gz"
    "https://genome-nexus-static-data.s3.amazonaws.com/mutationassessor_v4_3.tsv.gz"
    "https://genome-nexus-static-data.s3.amazonaws.com/mutationassessor_v4_4.tsv.gz"
)
for url in "${mutation_assessor_files[@]}"
do
    filename=$(basename $url)
    # Download file and extract it
    echo "Downloading $filename"
    curl $url -o ${DIR}/../data/common_input/${filename}
    echo "Download completed."

    echo "Extracting $filename"
    gunzip ${DIR}/../data/common_input/${filename}
    mutation_assessor_tsv_file="${filename%.gz}"

    echo "Transforming $mutation_assessor_tsv_file"
    # Rename the columns
    if [ -f "${DIR}/../data/common_input/${mutation_assessor_tsv_file}" ]; then
        echo "File exists: ${DIR}/../data/common_input/${mutation_assessor_tsv_file}"
    else
        echo "File not found: ${DIR}/../data/common_input/${mutation_assessor_tsv_file}"
        exit 1
    fi
    sed 's/uniprotId\tSV\thgvspShort\tF_score\tF_impact\tMSA\tMAV/uniprotId\tsv\thgvspShort\tf_score\tf_impact\tmsa\tmav/' ${DIR}/../data/common_input/$mutation_assessor_tsv_file
    # Add a new column "_id" (uniprotId,hgvspShort)
    awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "_id",$0; next} {print $1","$3,$0}' ${DIR}/$mutation_assessor_tsv_file > ${DIR}/../data/common_input/processed_$mutation_assessor_tsv_file

    # Import the data into MongoDB
    echo "Importing $mutation_assessor_tsv_file"
    import mutation_assessor.annotation ${DIR}/../data/common_input/processed_$mutation_assessor_tsv_file "--type tsv --headerline"

    rm ${DIR}/../data/common_input/processed_$mutation_assessor_tsv_file
    rm ${DIR}/../data/common_input/$mutation_assessor_tsv_file
done

# import annotation sources version
import version ${DIR}/../data/${REF_ENSEMBL_VERSION}/export/annotation_version.txt '--drop --type tsv --headerline'
#!/bin/bash
set -e

# Set default if ENV variables are not set
MONGO_URI=${MONGO_URI:-"mongodb://127.0.0.1:27017/annotator"}
REF_ENSEMBL_VERSION=${REF_ENSEMBL_VERSION:-"grch37_ensembl92"}
SPECIES=${SPECIES:-"homo_sapiens"}

echo "MONGO_URI:" ${MONGO_URI}
echo "REF_ENSEMBL_VERSION:" ${REF_ENSEMBL_VERSION}
echo "SPECIES:" ${SPECIES}

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
echo "Downloading Mutation assessor data"

curl http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz -o ${DIR}/../data/common_input/MA_scores_rel3_hg19_full.tar.gz

echo "Extracting Mutation assessor data"
tar -xvf ${DIR}/../data/common_input/MA_scores_rel3_hg19_full.tar.gz
rm ${DIR}/../data/common_input/MA_scores_rel3_hg19_full.tar.gz

echo "Transforming Mutation assessor data"
sed -i -e 's/"Mutation","RefGenome variant","Gene","Uniprot","Info","Uniprot variant","Func. Impact","FI score"/_id,rgaa,gene,uprot,info,var,F_impact,F_score/g' MA_scores_rel3_hg19_full/MA_scores_rel3_hg19_chr*
sed -i -e 's/hg19,//g' MA_scores_rel3_hg19_full/MA_scores_rel3_hg19_chr*

echo "Importing Mutation assessor data"
for filename in MA_scores_rel3_hg19_full/*.csv; do import mutation_assessor.annotation $filename '--type csv --headerline' && rm $filename; done
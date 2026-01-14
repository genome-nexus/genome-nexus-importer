#!/bin/bash
set -e

# ---------- Config & defaults ----------
MONGO_URI=${MONGO_URI:-"mongodb://127.0.0.1:27017/annotator"}
# Explicit DB name if you set it; else infer from URI path (fallback to 'annotator')
MONGO_DB=${MONGO_DB:-$(echo "$MONGO_URI" | awk -F/ '{print $4}' | awk -F\? '{print $1}')}
MONGO_DB=${MONGO_DB:-annotator}

REF_ENSEMBL_VERSION=${REF_ENSEMBL_VERSION:-"grch37_ensembl111"}
SPECIES=${SPECIES:-"homo_sapiens"}
MUTATIONASSESSOR=${MUTATIONASSESSOR:-"false"}
SLIM_MODE=${SLIM_MODE:-"false"}

# Comma-separated list of collections to force-update even if they already exist
COLLECTION_UPDATE_LIST_RAW=${COLLECTION_UPDATE_LIST:-""}

echo "MONGO_URI: ${MONGO_URI}"
echo "MONGO_DB: ${MONGO_DB}"
echo "REF_ENSEMBL_VERSION: ${REF_ENSEMBL_VERSION}"
echo "SPECIES: ${SPECIES}"
echo "MUTATIONASSESSOR: ${MUTATIONASSESSOR}"
echo "SLIM_MODE: ${SLIM_MODE}"
echo "COLLECTION_UPDATE_LIST: ${COLLECTION_UPDATE_LIST_RAW}"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# ---------- Helpers ----------
# Normalize update list into an array (trim spaces)
IFS=',' read -r -a COLLECTION_UPDATE_LIST_ARR <<< "$(echo "$COLLECTION_UPDATE_LIST_RAW" | tr -d '[:space:]')"

in_update_list() {
  local target="$1"
  for c in "${COLLECTION_UPDATE_LIST_ARR[@]}"; do
    [[ -n "$c" && "$c" == "$target" ]] && return 0
  done
  return 1
}

collection_exists() {
  local col="$1"
  # Use the mongosh shell to check if the collection exists in MONGO_DB
  # Returns "1" if present, "0" otherwise.
  local js="db.getSiblingDB('${MONGO_DB}').getCollectionNames().indexOf('${col}') !== -1 ? '1' : '0'"
  mongosh "$MONGO_URI" --quiet --eval "$js"
}

should_import() {
  local col="$1"
  local exists
  exists="$(collection_exists "$col")"
  if [[ "$exists" == "1" ]]; then
    if in_update_list "$col"; then
      echo "[INFO] $col exists but is in COLLECTION_UPDATE_LIST → will import/update."
      return 0
    else
      echo "[INFO] $col exists and is NOT in COLLECTION_UPDATE_LIST → skipping import."
      return 1
    fi
  else
    echo "[INFO] $col does not exist → will import."
    return 0
  fi
}

import() {
  local collection="$1"
  local file="$2"
  local extraoptions="${3:-}"

  mongoimport --uri "$MONGO_URI" --collection "$collection" $extraoptions --file "$file"
}

import_if_needed() {
  local collection="$1"
  local file="$2"
  local extraoptions="${3:-}"

  if should_import "$collection"; then
    import "$collection" "$file" "$extraoptions"
  fi
}

# ---------- Pre-flight ----------
if [[ ! -d "${DIR}/../data/${REF_ENSEMBL_VERSION}" ]]; then
  echo "Can't find directory for given reference genome and ensembl release: ${DIR}/../data/${REF_ENSEMBL_VERSION}"
  exit 1
fi

# ---------- Core imports ----------
import_if_needed "ensembl.biomart_transcripts" <(gunzip -c "${DIR}/../data/${REF_ENSEMBL_VERSION}/export/ensembl_biomart_transcripts.json.gz") '--drop --type json'
import_if_needed "ensembl.canonical_transcript_per_hgnc" "${DIR}/../data/${REF_ENSEMBL_VERSION}/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt" '--drop --type tsv --headerline'

# Exit early if slim mode (only import ensembl.biomart_transcripts, ensembl.canonical_transcript_per_hgnc, and version)
if [[ "$SLIM_MODE" == "true" ]]; then
  echo "[INFO] SLIM_MODE enabled. Only keep ensembl.biomart_transcripts, ensembl.canonical_transcript_per_hgnc, and version."
  # Import annotation sources version with only VEP and HGNC rows
  echo "[INFO] Filtering annotation_version.txt to only include VEP and HGNC for slim mode."
  import version <(awk 'NR==1 || /^VEP\t/ || /^HGNC\t/' "${DIR}/../data/${REF_ENSEMBL_VERSION}/export/annotation_version.txt") '--drop --type tsv --headerline'
  # Stop running the script afterwards
  exit 0
fi

# Exit early if non-human
if [[ "$SPECIES" != "homo_sapiens" ]]; then
  echo "[INFO] Non-human species detected. Skipping human-specific imports."
  exit 0
fi

echo "[INFO] Executing human-specific import steps"

import_if_needed "pfam.domain" "${DIR}/../data/${REF_ENSEMBL_VERSION}/export/pfamA.txt" '--drop --type tsv --headerline'
import_if_needed "ptm.experimental" <(gunzip -c "${DIR}/../data/ptm/export/ptm.json.gz") '--drop --type json'

import_if_needed "hotspot.mutation" "${DIR}/../data/${REF_ENSEMBL_VERSION}/export/hotspots_v2_and_3d.txt" '--drop --type tsv --headerline --mode upsert --upsertFields hugo_symbol,residue,type,tumor_count'
import_if_needed "signal.mutation" <(gunzip -c "${DIR}/../data/signal/export/mutations.json.gz") '--drop --type json'

# import oncokb cancer genes list
import_if_needed "oncokb.gene" "${DIR}/../data/${REF_ENSEMBL_VERSION}/export/oncokb_cancer_genes_list_from_API.json" '--drop --type json --jsonArray'

# import ClinVar
if should_import "clinvar.mutation"; then
  if [[ ${REF_ENSEMBL_VERSION} == *"grch37"* ]]; then
    curl https://genome-nexus-static-data.s3.us-east-1.amazonaws.com/clinvar_grch37.txt.gz -o "${DIR}/../data/common_input/clinvar_grch37.txt.gz"
    import "clinvar.mutation" <(gunzip -c "${DIR}/../data/common_input/clinvar_grch37.txt.gz") '--drop --type tsv --headerline --columnsHaveTypes --parseGrace autoCast'
  elif [[ ${REF_ENSEMBL_VERSION} == *"grch38"* ]]; then
    curl https://genome-nexus-static-data.s3.us-east-1.amazonaws.com/clinvar_grch38.txt.gz -o "${DIR}/../data/common_input/clinvar_grch38.txt.gz"
    import "clinvar.mutation" <(gunzip -c "${DIR}/../data/common_input/clinvar_grch38.txt.gz") '--drop --type tsv --headerline --columnsHaveTypes --parseGrace autoCast'
  fi
else
  echo "[INFO] Skipping ClinVar download/import."
fi

# import mutation assessor
if [[ "${MUTATIONASSESSOR}" == "true" ]]; then
  if should_import "mutation_assessor.annotation"; then
    declare -a mutation_assessor_files=(
      "https://genome-nexus-static-data.s3.us-east-1.amazonaws.com/mutationassessor_v4_1.tsv.gz"
      "https://genome-nexus-static-data.s3.us-east-1.amazonaws.com/mutationassessor_v4_2.tsv.gz"
      "https://genome-nexus-static-data.s3.us-east-1.amazonaws.com/mutationassessor_v4_3.tsv.gz"
      "https://genome-nexus-static-data.s3.us-east-1.amazonaws.com/mutationassessor_v4_4.tsv.gz"
    )
    first_file=true
    for url in "${mutation_assessor_files[@]}"; do
      filename=$(basename "$url")
      echo "Downloading $filename"
      curl "$url" -o "${DIR}/../data/common_input/${filename}"

      echo "Extracting $filename"
      gunzip "${DIR}/../data/common_input/${filename}"
      mutation_assessor_tsv_file="${filename%.gz}"

      echo "Transforming $mutation_assessor_tsv_file"
      sed -i 's/uniprotId\tSV\thgvspShort\tF_score\tF_impact\tMSA\tMAV/uniprotId\tsv\thgvspShort\tf_score\tf_impact\tmsa\tmav/' "${DIR}/../data/common_input/$mutation_assessor_tsv_file"
      awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "_id",$0; next} {print $1","$3,$0}' "${DIR}/../data/common_input/$mutation_assessor_tsv_file" > "${DIR}/../data/common_input/processed_$mutation_assessor_tsv_file"

      echo "Importing $mutation_assessor_tsv_file"
      
      import_opts="--type tsv --headerline"
      if [ "$first_file" = true ]; then
        import_opts="--drop $import_opts"
        first_file=false
      fi

      import "mutation_assessor.annotation" "${DIR}/../data/common_input/processed_$mutation_assessor_tsv_file" "$import_opts"

      rm "${DIR}/../data/common_input/processed_$mutation_assessor_tsv_file" "${DIR}/../data/common_input/$mutation_assessor_tsv_file"
    done
  else
    echo "[INFO] Skipping MutationAssessor (collection up-to-date and not forced)."
  fi
fi

# Always import annotation sources version
import version "${DIR}/../data/${REF_ENSEMBL_VERSION}/export/annotation_version.txt" '--drop --type tsv --headerline'

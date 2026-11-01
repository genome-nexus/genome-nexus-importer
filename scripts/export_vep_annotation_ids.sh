#!/bin/bash
# Export every _id from a vep.annotation collection into filtered, sharded key
# files, in one streaming pass. Resumable: re-running continues after the last
# key already written.
#
# Usage: export_vep_annotation_ids.sh <mongo-uri> <outdir> [max-keys]
#
# Filters mirror CachedVariantAnnotationFetcher.isValidId in genome-nexus, plus
# two measured corruption classes. Dropped counts are reported so that
# kept + dropped can be reconciled against estimatedDocumentCount().
set -euo pipefail

URI="${1:?usage: $0 <mongo-uri> <outdir> [max-keys]}"
OUTDIR="${2:?usage: $0 <mongo-uri> <outdir> [max-keys]}"
MAX_KEYS="${3:-0}"
PER_SHARD="${PER_SHARD:-250000}"
HERE="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUTDIR"

# Resume: pick up after the last key in the highest-numbered shard, and keep
# numbering where it stopped so shard files stay uniform.
RESUME_AFTER=""
START_INDEX=0
LAST_SHARD=$(ls "$OUTDIR"/ids.* 2>/dev/null | sort | tail -1 || true)
if [ -n "$LAST_SHARD" ]; then
  RESUME_AFTER=$(tail -1 "$LAST_SHARD")
  COMPLETE=$(ls "$OUTDIR"/ids.* | wc -l | tr -d ' ')
  START_INDEX=$(( (COMPLETE - 1) * PER_SHARD + $(wc -l < "$LAST_SHARD" | tr -d ' ') ))
  echo "resuming after '$RESUME_AFTER' at index $START_INDEX" >&2
fi

export RESUME_AFTER MAX_KEYS

# SCAN_CMD exists so the filter/shard half can be tested without a mongod;
# it defaults to the real covered scan. mongosh --quiet so only print() output
# reaches stdout.
SCAN_CMD="${SCAN_CMD:-mongosh \"$URI\" --quiet --file \"$HERE/scan_vep_annotation_ids.js\"}"
eval "$SCAN_CMD" \
  | awk -v dir="$OUTDIR" -v per="$PER_SHARD" -v n="$START_INDEX" -v n0="$START_INDEX" '
      # Corruption seen in prod: keys with leading/trailing whitespace (e.g.
      # " ENST00000447712.2:p.G703S") and with a leading quote. Trailing-quote
      # variants survive these and are deferred later by the converter.
      /^[[:space:]]|[[:space:]]$/ { bad["whitespace"]++; next }
      /^"/               { bad["quote"]++;   next }
      !/:/               { bad["nocolon"]++; next }   # isValidId
      /-|undefined|g\.0/ { bad["invalid"]++; next }   # isValidId
      # JsonMongoRepositoryImpl refuses to cache keys >= 900 bytes, so storing
      # an annotation under one is pointless.
      length($0) >= 900  { bad["toolong"]++; next }
      {
        # f == "" covers a resume landing mid-shard, where n % per != 0 and the
        # file handle would otherwise never be opened. Append rather than
        # truncate so a partial shard is continued instead of overwritten.
        if (f == "" || n % per == 0) {
          if (f != "") close(f)
          f = sprintf("%s/ids.%05d", dir, int(n/per))
        }
        print >> f
        n++
      }
      END {
        if (f != "") close(f)
        for (k in bad) { printf "dropped %s %d\n", k, bad[k] > "/dev/stderr"; d += bad[k] }
        printf "kept %d  dropped %d  total %d\n", n - n0, d, (n - n0) + d > "/dev/stderr"
      }'

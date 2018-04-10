#!/bin/bash
# Downloads necessary import and script files to import into a mongo database
# given a github url to the export/ folder. Outputs temp directory where the
# data has been downloaded
#
# Usage e.g.:
#
# ./download_files_from_github_url.sh https://github.com/genome-nexus/genome-nexus-importer/tree/master/export
#
# or if using file remotely:
#
# bash <(curl https://raw.githubusercontent.com/genome-nexus/genome-nexus-importer/master/scripts/download_files_from_github_url.sh) https://github.com/genome-nexus/genome-nexus-importer/tree/master/export

# halt on error
set -e

GITHUB_EXPORT_URL=$1 # https://github.com/genome-nexus/genome-nexus-importer/tree/master/export

if [[ "$(uname)" = "Darwin" ]]; then
    CHMOD="gchmod"
else
    CHMOD="chmod"
fi

# create temp dir
TMPDIR=$(mktemp -d)
$CHMOD -t ${TMPDIR}
$CHMOD g+rwX -R ${TMPDIR}
cd $TMPDIR

EXPORT_FILES=$(curl "$GITHUB_EXPORT_URL" | grep -oE '/genome-nexus[^ "]*(.txt|.json|.json.gz)')
mkdir -p export
cd export

for f in $EXPORT_FILES; do
    curl -GLO "https://github.com/$f" -d 'raw=true'
done

cd ..
mkdir scripts
cd scripts
curl -GLO "${GITHUB_EXPORT_URL/tree/blob}/../scripts/import_mongo.sh" -d 'raw=true'
echo $TMPDIR

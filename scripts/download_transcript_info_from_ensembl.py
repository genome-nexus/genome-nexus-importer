#!/usr/bin/env python3

'''
This script uses the table of gene IDs downloaded from Ensembl BioMart, and for each
transcript determines whether it is the canonical transcript for its gene, its protein
(translation) stable ID, and its protein length.

For homo_sapiens/GRCh37 or GRCh38 (--genome-build + --ensembl-release given), this is
done by querying Ensembl's public, anonymous, read-only MySQL mirror
(ensembldb.ensembl.org) for is_canonical/transcript version/protein stable ID, and by
downloading the release-pinned protein FASTA (Homo_sapiens.GRCh3x.pep.all.fa.gz) for
protein length. This keeps every file in the pipeline pinned to the same Ensembl
release and avoids the Ensembl REST API entirely, which (a) isn't release-pinned - it
always answers from whatever is currently live, which can silently drift from an older
pinned release - and (b) has historically been unreliable/down.

For anything else (e.g. mouse), this falls back to the original Ensembl REST /lookup/id
API, batched in blocks of <qsize>, since there's no verified MySQL equivalent for it yet.
'''

import pandas as pd
import numpy as np
import requests
import subprocess
import urllib.request
import gzip
import sys
import os
import argparse
import re

MYSQL_HOST = "ensembldb.ensembl.org"
MYSQL_PORT = "3306"
MYSQL_USER = "anonymous"

ASSEMBLY_SUFFIX = {"grch37": "37", "grch38": "38"}
# GRCh37's core database is permanently frozen at Ensembl release 75 - later release
# numbers (as used for the repackaged FTP protein FASTA, see pep_fasta_url below) don't
# correspond to a real core database for it.
FROZEN_DB_RELEASE = {"grch37": "75"}


# --- MySQL + FTP path (homo_sapiens / grch37 or grch38) ---------------------------

def run_mysql_query(database, sql):
    result = subprocess.run(
        ["mysql", "-h", MYSQL_HOST, "-P", MYSQL_PORT, "-u", MYSQL_USER, "-N", "-B", database, "-e", sql],
        capture_output=True, text=True, check=True,
    )
    return result.stdout


def database_exists(database):
    result = subprocess.run(
        ["mysql", "-h", MYSQL_HOST, "-P", MYSQL_PORT, "-u", MYSQL_USER, "-N", "-B",
         "-e", f"SHOW DATABASES LIKE '{database}';"],
        capture_output=True, text=True, check=True,
    )
    return result.stdout.strip() == database


def fetch_canonical_and_protein_ids(genome_build, ensembl_release):
    db_release = FROZEN_DB_RELEASE.get(genome_build, ensembl_release)
    database = f"homo_sapiens_core_{db_release}_{ASSEMBLY_SUFFIX[genome_build]}"

    if not database_exists(database):
        sys.stderr.write(f"No core database '{database}' on {MYSQL_HOST} for {genome_build} release {ensembl_release}.\n")
        sys.exit(1)

    print(f"Querying {database} on {MYSQL_HOST} for canonical/translation info...")
    sql = """
        SELECT t.stable_id, t.version, IF(t.transcript_id = g.canonical_transcript_id, 1, 0),
               IFNULL(tl.stable_id, '')
        FROM transcript t
        JOIN gene g ON g.gene_id = t.gene_id
        LEFT JOIN translation tl ON tl.transcript_id = t.transcript_id
        WHERE t.is_current = 1;
    """
    records = {}
    for row in run_mysql_query(database, sql).strip('\n').split('\n'):
        if not row:
            continue
        stable_id, version, is_canonical, protein_stable_id = row.split('\t')
        records[stable_id] = {
            'transcript_id_version': version,
            'is_canonical': is_canonical == '1',
            'protein_stable_id': protein_stable_id if protein_stable_id else np.nan,
        }
    return records


def pep_fasta_url(genome_build, ensembl_release):
    if genome_build == "grch37":
        return (f"https://ftp.ensembl.org/pub/grch37/release-{ensembl_release}/fasta/"
                f"homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz")
    return (f"https://ftp.ensembl.org/pub/release-{ensembl_release}/fasta/"
            f"homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz")


def fetch_protein_lengths(genome_build, ensembl_release, tmp_dir):
    url = pep_fasta_url(genome_build, ensembl_release)
    local_path = os.path.join(tmp_dir, os.path.basename(url))
    if not os.path.exists(local_path):
        print(f"Downloading {url} ...")
        urllib.request.urlretrieve(url, local_path)
    else:
        print(f"Using cached {local_path}")

    lengths_by_transcript = {}
    with gzip.open(local_path, 'rt') as f:
        transcript_id = None
        length = 0
        for line in f:
            if line.startswith('>'):
                if transcript_id is not None:
                    lengths_by_transcript[transcript_id] = length
                transcript_id = None
                for field in line.split():
                    if field.startswith('transcript:'):
                        transcript_id = field.split(':', 1)[1].split('.')[0]
                length = 0
            else:
                length += len(line.strip())
        if transcript_id is not None:
            lengths_by_transcript[transcript_id] = length
    return lengths_by_transcript


def lookup_transcripts_via_mysql(gene_info, tmp_dir, genome_build, ensembl_release):
    canonical_and_protein = fetch_canonical_and_protein_ids(genome_build, ensembl_release)
    protein_lengths = fetch_protein_lengths(genome_build, ensembl_release, tmp_dir)

    def lookup(transcript_stable_id):
        info = canonical_and_protein.get(transcript_stable_id, {})
        protein_stable_id = info.get('protein_stable_id', np.nan)
        protein_length = protein_lengths.get(transcript_stable_id, np.nan)
        return pd.Series(
            ['1' if info.get('is_canonical') else '0', info.get('transcript_id_version', np.nan),
             protein_stable_id, protein_length],
            index=['is_canonical', 'transcript_id_version', 'protein_stable_id', 'protein_length']
        )

    return gene_info['transcript_stable_id'].apply(lookup)


# --- Ensembl REST fallback (anything other than homo_sapiens grch37/grch38) --------

def request_transcript_ids(transcripts):
    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    transcripts_formatted = '", "'.join(transcripts)
    data = '{"expand": 1, "format":"full", "ids": ["%s"] }' % transcripts_formatted

    try:
        r = requests.post(server + ext, headers=headers, data=data)
    except requests.exceptions.ConnectionError:
        sys.stderr.write('Connection error when trying to query Ensembl API')
        sys.exit(1)

    if not r.ok:
        r.raise_for_status()
        sys.exit(1)

    return r.json()


def get_transcript_info(transcript, ensembl_transcript_response):
    # Attempt to parse API response. Sometimes the transcript does not have an API response, probably because the API is
    # on a newer Ensembl release than the input files. Restarting the pipeline will attempt to continue.
    try:
        data = ensembl_transcript_response.get(transcript, {})
        is_canonical = data.get('is_canonical', np.nan)

        transcript_id_version = str(data.get('version')) if 'version' in data else np.nan

        translation = data.get('Translation', {})
        protein_stable_id = translation.get('id', np.nan)

        # ensure protein_length is integer (or NaN)
        protein_length = translation.get('length', np.nan)
        protein_length = (
            int(protein_length) if not pd.isna(protein_length) else np.nan
        )

    except (TypeError, KeyError):
        is_canonical = np.nan
        transcript_id_version = np.nan
        protein_stable_id = np.nan
        protein_length = np.nan

    return pd.Series(
        [is_canonical, transcript_id_version, protein_stable_id, protein_length],
        index=['is_canonical', 'transcript_id_version', 'protein_stable_id', 'protein_length']
    )


def get_rest_jobs(tmp_dir, ngenes):
    '''Reads the tmp dir and determines which indexes still need processing.
    '''

    # determine which indices are already processed. Start and stop are derived from filenames.
    filename_re = re.compile('transcript_info_([0-9]*)-([0-9]*)\.txt$')
    tmp_files = [fn for fn in os.listdir(tmp_dir) if re.match(filename_re, fn) != None]

    ngenes_set = set(range(0, ngenes))
    passed = set()
    for fn in tmp_files:
        match = re.search(filename_re, fn)
        start, end = int(match.group(1)), int(match.group(2))
        passed.update(range(start, end))

    return ngenes_set.difference(passed)


def lookup_transcripts_via_rest(gene_info, tmp_dir, jobs, query_size):
    '''Loops through gene IDs and looks them up in Ensembl if needed.
    Results are merged into a data frame.
    '''

    last_job = len(gene_info)

    transcripts_all = gene_info['versioned_transcript_id']
    transcript_info = pd.DataFrame({'is_canonical': [], 'transcript_id_version': [], 'protein_stable_id': [], 'protein_length': []})

    # Iterate over transcripts
    for low_index in range(0, last_job, query_size):
        high_index = min(low_index + query_size, last_job)
        tmp_file = os.path.join(tmp_dir, 'transcript_info_%s-%s.txt' % (low_index, high_index))

        # Check whether we already processed all indices between start and stop, else redo this chunk.
        if len(set(range(low_index, high_index)).intersection(jobs)) > 0:
            print('Retrieving %s-%s of %s transcripts from Ensembl' % (low_index, high_index, last_job))
            transcripts_chunk = transcripts_all[low_index:high_index]

            # Request, decode and save info for transcripts
            decoded = request_transcript_ids(transcripts_chunk.tolist())
            transcript_info_chunk = transcripts_chunk.apply(lambda x: get_transcript_info(x, decoded))
            transcript_info_chunk.to_csv(tmp_file, sep='\t', index=True)

    # read all files and append to final table, but keep index
    filename_re = re.compile('transcript_info_([0-9]*)-([0-9]*)\.txt$')
    tmp_files = [fn for fn in os.listdir(tmp_dir) if re.match(filename_re, fn) != None]
    for fn in tmp_files:
        transcript_info_chunk = pd.read_csv(os.path.join(tmp_dir, fn), sep='\t', index_col=0)
        transcript_info = pd.concat([transcript_info, transcript_info_chunk], ignore_index=False)

    # Prepare transcript info table for merging. Sort and remove duplicate indices.
    transcript_info.sort_index(inplace=True)
    transcript_info = transcript_info.loc[~transcript_info.index.duplicated(keep='first')]
    transcript_info['is_canonical'].fillna(False, inplace=True)
    transcript_info['is_canonical'] = transcript_info['is_canonical'].astype(bool)
    transcript_info['is_canonical'].replace({True: '1', False: '0'}, inplace=True)

    return transcript_info


def main(ensembl_biomart_geneids, ensembl_canonical_data, query_size, genome_build, ensembl_release):
    gene_info = pd.read_csv(ensembl_biomart_geneids, sep='\t', dtype=str)
    gene_info.columns = [c.lower().replace(' ', '_') for c in gene_info.columns]

    # Create temporary directory to save files for each transcripts file
    tmp_dir = os.path.join(os.path.dirname(ensembl_canonical_data), 'transcript_info/')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    if genome_build in ASSEMBLY_SUFFIX and ensembl_release:
        transcript_info = lookup_transcripts_via_mysql(gene_info, tmp_dir, genome_build, ensembl_release)
        transcript_info.index = gene_info.index
    else:
        ngenes = len(gene_info)
        jobs = get_rest_jobs(tmp_dir, ngenes)
        transcript_info = lookup_transcripts_via_rest(gene_info, tmp_dir, jobs, query_size)

    # check whether the total number of jobs is correct
    assert(len(transcript_info.index) == len(gene_info.index))
    # check whether the order of transcripts is incremental and nonredundant
    assert(list(transcript_info.index) == list(range(0, len(transcript_info.index))))

    # merge with gene IDs, save
    gene_transcript_info = pd.concat([gene_info, transcript_info], axis=1, sort=False)
    gene_transcript_info.to_csv(ensembl_canonical_data, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ensembl_biomart_geneids",
                        help="input/ensembl_biomart_geneids.txt")
    parser.add_argument("ensembl_canonical_data",
                        help="tmp/ensembl_canonical_data.txt")
    parser.add_argument("-q", "--querysize",
                        help="REST fallback only: the number of Ensembl IDs submitted per POST request",
                        default=1000,
                        type=int)
    parser.add_argument("--genome-build",
                        help="grch37 or grch38 - queries Ensembl's MySQL mirror directly instead of "
                             "the REST API. Omit (or pass anything else, e.g. grcm38 for mouse) to "
                             "fall back to the original REST-based lookup.",
                        default=None)
    parser.add_argument("--ensembl-release",
                        help="Ensembl release number to pin to, required together with --genome-build",
                        default=None)
    args = parser.parse_args()

    main(args.ensembl_biomart_geneids, args.ensembl_canonical_data, args.querysize,
         args.genome_build, args.ensembl_release)

import pandas as pd
import numpy as np
import requests
import sys
import os
import argparse


def request_transcript_ids(transcripts, inputfile):
    # Prepare API call
    if 'grch37' in inputfile:
        server = "https://grch37.rest.ensembl.org"
    else:
        server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    transcripts_formatted = '", "'.join(transcripts)
    data = '{"expand": 1, "format":"full", "ids": ["%s"] }' % transcripts_formatted

    # Perform API call
    try:
        r = requests.post(server+ext, headers=headers, data=data)
    except requests.exceptions.ConnectionError:
        sys.stderr.write('Connection error when trying to query Ensembl API')
        sys.exit(1)

    # Check if return JSON is ok
    if not r.ok:
        r.raise_for_status()
        sys.exit(1)

    return r.json()


def get_transcript_info(transcript, ensembl_transcript_response):
    # Attempt to parse API response. Sometimes the transcript does not have an API response, probably because the API is
    # on a newer Ensembl release than the input files. Restarting the pipeline will attempt to continue.
    try:
        is_canonical = ensembl_transcript_response[transcript]['is_canonical']

        try:
            protein_stable_id = ensembl_transcript_response[transcript]['Translation']['id']
        except KeyError:
            protein_stable_id = np.nan

        try:
            # store as string to prevent integer -> float
            protein_length = str(ensembl_transcript_response[transcript]['Translation']['length'])
        except KeyError:
            protein_length = np.nan

    except TypeError:
        # print('Transcript %s has no response from Ensembl API, perhaps API is on newer Ensembl Release and ID is '
        #       'deprecated.' % transcript)
        is_canonical = np.nan
        protein_stable_id = np.nan
        protein_length = np.nan

    return pd.Series(
        [is_canonical, protein_stable_id, protein_length],
        index='is_canonical protein_stable_id protein_length'.split()
    )


def main(ensembl_biomart_geneids, ensembl_biomart_geneids_transcript_info):
    gene_info = pd.read_csv(ensembl_biomart_geneids, sep='\t')
    gene_info.columns = [c.lower().replace(' ', '_') for c in gene_info.columns]
    # print('Retrieving transcript information per gene to retrieve:\n'
    #       '- whether transcript is canonical\n'
    #       '- protein ID\n'
    #       '- protein length')
    # print('Can only retrieve 1000 transcripts per POST request, see '
    #       'https://github.com/Ensembl/ensembl-rest/wiki/POST-Requests')

    # Create temporary directory to save files for each 1000 transcripts
    tmp_dir = ensembl_biomart_geneids.replace('input/ensembl_biomart_geneids.txt', 'tmp/transcript_info')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    transcripts_all = gene_info['transcript_stable_id']
    transcript_info = pd.DataFrame()

    for low_index in range(0, len(gene_info), 1000):
        high_index = min(low_index + 1000, len(gene_info))
        tmp_file = os.path.join(tmp_dir, 'transcript_info_%s-%s.txt' % (low_index, high_index))

        # Check if this file has been made before
        if not os.path.isfile(tmp_file):
            # print('Retrieving %s-%s of %s transcripts from Ensembl' % (low_index, high_index, len(gene_info)))
            transcripts_chunk = transcripts_all[low_index:high_index]

            # Request, decode and save info for 1000 transcripts
            decoded = request_transcript_ids(transcripts_chunk.tolist(), sys.argv[1])
            transcript_info_chunk = transcripts_chunk.apply(lambda x: get_transcript_info(x, decoded))
            transcript_info_chunk.to_csv(tmp_file, sep='\t', index=False)
        # Open previously created transcript info
        else:
            # print('Found tmp file for %s-%s transcripts' % (low_index, high_index))
            transcript_info_chunk = pd.read_table(tmp_file, sep='\t', dtype=str)
        transcript_info = transcript_info.append(transcript_info_chunk, ignore_index=False)

    # Prepare transcript info table for merging
    transcript_info.reset_index(drop=True, inplace=True)
    transcript_info['is_canonical'].fillna('0', inplace=True)
    transcript_info['is_canonical'].replace({'1.0': '1', '0.0': '0'}, inplace=True)

    gene_transcript_info = pd.concat([gene_info, transcript_info], axis=1, sort=False)
    gene_transcript_info.to_csv(ensembl_biomart_geneids_transcript_info, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ensembl_biomart_geneids",
                        help="input/ensembl_biomart_geneids.txt")
    parser.add_argument("ensembl_biomart_geneids_transcript_info",
                        help="tmp/ensembl_biomart_geneids.transcript_info")
    args = parser.parse_args()

    main(args.ensembl_biomart_geneids, args.ensembl_biomart_geneids_transcript_info)

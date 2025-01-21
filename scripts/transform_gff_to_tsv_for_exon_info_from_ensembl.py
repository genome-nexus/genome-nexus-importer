#!/usr/bin/env python3
"""
Copyright (c) 2018 The Hyve B.V.
This code is licensed under the GNU Affero General Public License (AGPL),
version 3, or (at your option) any later version.
"""

import pandas as pd
import gzip
import argparse


def main(gff_file, ensembl_transcript_info):
    """Transform GFF3 file to TSV for exons and UTRs. Output is written to
    stdout, as the Makefile wrapper expects. Input file is GFF3 file, see Makefile"""

    # Dataframe to append the transcript information
    transcript_info = pd.DataFrame(columns=['transcript_id', 'type', 'id', 'start', 'end', 'rank', 'strand', 'version'])
    # list to append transcript information per row
    rows = []

    # Open gff file and read lines, when line contains transcript information, extract this information
    with gzip.open(gff_file, 'rt') as gff:
        for line in gff:
            if line[0] != '#':
                list_line = line.strip('\n').split('\t')
                entry_dict = {}

                # Extract UTRs and exons
                if len(list_line) > 1 and list_line[2] in ['exon', 'five_prime_UTR', 'three_prime_UTR']:

                    meta_info = list_line[8].split(';')
                    entry_dict['transcript_id'] = meta_info[0].split(':')[1]
                    entry_dict['type'] = list_line[2]
                    entry_dict['start'] = list_line[3]
                    entry_dict['end'] = list_line[4]
                    strand = list_line[6]
                    # Convert plus strand into 1 and minus strand into -1
                    if strand == '+':
                        entry_dict['strand'] = '1'
                    elif strand == '-':
                        entry_dict['strand'] = '-1'
                    else:
                        entry_dict['strand'] = ''

                    if list_line[2] == 'exon':
                        entry_dict['id'] = meta_info[5].split('=')[1]
                        entry_dict['rank'] = meta_info[6].split('=')[1]
                        entry_dict['version'] = meta_info[7].split('=')[1]
                    else:
                        entry_dict['id'] = ''
                        entry_dict['rank'] = ''
                        entry_dict['version'] = ''
                    rows.append(entry_dict)

    # By first appending it to a list and only adding it to a DF once, performance is greatly improved
    transcript_info = pd.concat([transcript_info, pd.DataFrame(rows)], ignore_index=True, sort=False)
    transcript_info.to_csv(ensembl_transcript_info, sep='\t', index=False)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_file",
                        help="tmp/annotations.gff3.gz")
    parser.add_argument("ensembl_transcript_info",
                        help="tmp/ensembl_transcript_info.txt")
    args = parser.parse_args()

    main(args.gff_file, args.ensembl_transcript_info)

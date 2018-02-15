#!/usr/bin/env python3
"""Add domain as nested property to transcript. Same for hugo symbol and exon
info. Output resulting JSON"""

import argparse
import pandas as pd
import numpy as np
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ensembl_biomart_transcripts", default="../data/ensembl_biomart_transcripts.txt", type=str, help="Ensembl Biomart Transcript info")
    parser.add_argument("ensembl_biomart_pfam", default="../data/ensembl_biomart_pfam_grch37.p13.txt", type=str, help="Ensembl Biomart PFAM domain info")
    parser.add_argument("ensembl_exon_info", default="../data/ensembl_exon_info.txt", type=str, help="Ensembl Exon info extracted from GFF file")
    args = parser.parse_args()

    transcripts = pd.read_csv(args.ensembl_biomart_transcripts, sep='\t')

    # make nested object hugo symbols
    transcripts = transcripts.set_index('transcript_stable_id')
    def get_hgnc_symbol(transcript_id):
        hgnc_symbols = transcripts.loc[transcript_id]
        if hgnc_symbols.ndim == 1:
            if pd.isnull(hgnc_symbols.hgnc_symbol):
                return hgnc_symbols.hgnc_symbol
            else:
                return [hgnc_symbols.hgnc_symbol]
        else:
            return list(hgnc_symbols.hgnc_symbol.values)
    hgnc_symbol_list = transcripts.index.drop_duplicates().map(get_hgnc_symbol)
    # make one row per transcript_stable_id by removing hgnc_symbol
    unique_transcripts = transcripts.copy().reset_index()
    del unique_transcripts['hgnc_symbol']
    unique_transcripts = unique_transcripts.drop_duplicates()
    # should only be one row per transcript now
    assert(len(unique_transcripts) == len(transcripts.index.drop_duplicates()))
    # should have the same order after dropping duplicates as
    # hgnc_symbol_list
    assert((unique_transcripts.transcript_stable_id != transcripts.index.drop_duplicates()).sum() == 0)
    unique_transcripts['hgnc_symbols'] = hgnc_symbol_list
    unique_transcripts = unique_transcripts.set_index("transcript_stable_id")

    # add nested exons
    exons = pd.read_csv(args.ensembl_exon_info, sep="\t")
    def get_exon_dict(x):
        d = x.to_dict(orient='record')
        for rec in d:
            # remove key that's already in the index of group by
            del rec['transcriptId']
        return d
    exon_dicts = exons.groupby('transcriptId').apply(get_exon_dict)
    def get_exons_for_transcript(x):
        try:
            return exon_dicts.loc[x]
        except KeyError:
            return np.nan
    unique_transcripts["exons"] = unique_transcripts.index.map(get_exons_for_transcript)

    # add nested domains
    trans_domain =  pd.read_csv(args.ensembl_biomart_pfam, sep='\t')
    trans_domain.columns = [c.lower().replace(' ','_') for c in trans_domain.columns]
    domain_grouped = trans_domain.groupby("transcript_stable_id")
    def get_domain_for_transcript(x):
        try:
            # Fix for: https://github.com/pandas-dev/pandas/issues/9466
            # use groupby to get indices, loc is slow for non-unique
            return trans_domain.iloc[domain_grouped.groups[x]]['pfam_domain_id pfam_domain_start pfam_domain_end'.split()]
        except KeyError:
            return np.nan
    unique_transcripts["domains"] = unique_transcripts.index.map(get_domain_for_transcript)


    # print records as json
    unique_transcripts.reset_index().to_json(sys.stdout, orient='records', lines=True)

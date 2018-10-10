#!/usr/bin/env python3
"""Add domain as nested property to transcript. Same for hugo symbol and exon
info. Output resulting JSON"""

import pandas as pd
import numpy as np
import argparse
import sys


def add_nested_hgnc(transcripts):
    """ Make nested object HGNC symbols per transcript"""

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
    assert (len(unique_transcripts) == len(transcripts.index.drop_duplicates()))
    # should have the same order after dropping duplicates as
    # hgnc_symbol_list
    assert (0 == (unique_transcripts.transcript_stable_id != transcripts.index.drop_duplicates()).sum())
    unique_transcripts['hgnc_symbols'] = hgnc_symbol_list
    return unique_transcripts.set_index("transcript_stable_id")


def add_nested_transcript_info(transcripts, transcript_info):
    """ Make nested object with exons and UTR per transcript. """

    def get_list_of_info_dicts(transcript_group):
        list_of_info_dicts = transcript_group.to_dict(orient='record')
        for info_dict in list_of_info_dicts:
            # remove key that's already in the index of group by
            del info_dict['transcript_id']

            # For UTRs, remove rank
            if info_dict['type'] in ['five_prime_UTR', 'three_prime_UTR']:
                del info_dict['id']
                del info_dict['rank']
                del info_dict['version']

            if info_dict['type'] == 'exon':
                del info_dict['type']

        return list_of_info_dicts

    # Split in exons and UTRs
    exon_info = transcript_info.loc[transcript_info.type == 'exon']
    utr_info = transcript_info.loc[transcript_info.type.isin(['five_prime_UTR', 'three_prime_UTR'])]

    # Per transcriptID, take all exons and create a list of exons dictionaries. Save all these lists in a series.
    series_of_lists_of_exon_dicts = exon_info.groupby('transcript_id').apply(get_list_of_info_dicts)
    series_of_lists_of_utr_dicts = utr_info.groupby('transcript_id').apply(get_list_of_info_dicts)

    # Add a list of exon dictionaries to every transcript
    transcripts['exons'] = transcripts.index.map(series_of_lists_of_exon_dicts)
    transcripts['utrs'] = transcripts.index.map(series_of_lists_of_utr_dicts)

    return transcripts


def add_nested_pfam_domains(transcripts, pfam_domains):
    """ Add nested PFAM domains"""
    pfam_domains.columns = [c.lower().replace(' ', '_') for c in pfam_domains.columns]
    domain_grouped = pfam_domains.groupby("transcript_stable_id")

    def get_domain_for_transcript(x):
        try:
            # Fix for: https://github.com/pandas-dev/pandas/issues/9466
            # use groupby to get indices, loc is slow for non-unique
            return pfam_domains.iloc[domain_grouped.groups[x]][
                'pfam_domain_id pfam_domain_start pfam_domain_end'.split()]
        except KeyError:
            return np.nan

    transcripts["domains"] = transcripts.index.map(get_domain_for_transcript)
    return transcripts

def add_refseq(transcripts, refseq):
    """Add one refseq id for each transcript. There can be multiple. Pick
    highest number transcript id in that case."""
    refseq.columns = [c.lower().replace(' ','_') for c in refseq.columns]
    refseq_grouped = refseq.groupby("transcript_stable_id")

    def get_refseq_for_transcript(x):
        try:
            # Fix for: https://github.com/pandas-dev/pandas/issues/9466
            # use groupby to get indices, loc is slow for non-unique
            return refseq.iloc[refseq_grouped.groups[x]]['refseq_mrna_id'].sort_values(ascending=False).values[0]
        except KeyError:
            return np.nan

    transcripts['refseq_mrna_id'] = transcripts.index.map(get_refseq_for_transcript)
    return transcripts

def add_ccds(transcripts, ccds):
    """Add one ccds id for each transcript. There is only one per transcript."""
    ccds.columns = [c.lower().replace(' ','_') for c in ccds.columns]
    ccds = ccds[~pd.isnull(ccds["ccds_id"])]
    # assume each transcript has only one CCDS
    assert(any(ccds["transcript_stable_id"].duplicated()) == False)
    ccds = ccds.set_index("transcript_stable_id")

    def get_ccds_for_transcript(x):
        try:
            return ccds.loc[x].values[0]
        except KeyError:
            return np.nan

    transcripts["ccds_id"] = transcripts.index.map(get_ccds_for_transcript)
    return transcripts


def main(ensembl_biomart_transcripts, ensembl_transcript_info, ensembl_biomart_pfam, ensembl_biomart_refseq, ensembl_biomart_ccds):

    # Read input and set index column
    transcripts = pd.read_csv(ensembl_biomart_transcripts, sep='\t')
    transcript_info = pd.read_csv(ensembl_transcript_info, sep="\t")
    pfam_domains = pd.read_csv(ensembl_biomart_pfam, sep='\t')
    transcripts.set_index('transcript_stable_id', inplace=True)

    # import refseq
    refseq = pd.read_csv(ensembl_biomart_refseq, sep="\t")
    transcripts = add_refseq(transcripts, refseq)

    # import ccds
    ccds = pd.read_csv(ensembl_biomart_ccds, sep="\t")
    transcripts = add_ccds(transcripts, ccds)

    # Add nested HGNC, exons and PFAM domains
    transcripts = add_nested_hgnc(transcripts)
    transcripts = add_nested_transcript_info(transcripts, transcript_info)
    transcripts = add_nested_pfam_domains(transcripts, pfam_domains)

    # print records as json
    transcripts.reset_index().to_json(sys.stdout, orient='records', lines=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ensembl_biomart_transcripts",
                        default="../data/ensembl_biomart_transcripts.txt",
                        help="Ensembl Biomart Transcript info")

    parser.add_argument("ensembl_transcript_info",
                        default="../data/ensembl_transcript_info.txt.gz",
                        help="Ensembl transcript info extracted from GFF file")

    parser.add_argument("ensembl_biomart_pfam",
                        default="../data/ensembl_biomart_pfam_grch37.p13.txt",
                        help="Ensembl Biomart PFAM domain info")
    parser.add_argument("ensembl_biomart_refseq",
                        default="../data/ensembl_biomart_refseq_grch37.p13.txt",
                        help="Ensembl Biomart RefSeq info")
    parser.add_argument("ensembl_biomart_ccds",
                        default="../data/ensembl_biomart_ccds_grch37.p13.txt",
                        help="Ensembl Biomart CCDS info")

    args = parser.parse_args()
    main(args.ensembl_biomart_transcripts, args.ensembl_transcript_info, args.ensembl_biomart_pfam, args.ensembl_biomart_refseq, args.ensembl_biomart_ccds)

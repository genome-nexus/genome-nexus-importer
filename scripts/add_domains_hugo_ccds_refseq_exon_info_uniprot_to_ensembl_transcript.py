#!/usr/bin/env python3
"""Add domain as nested property to transcript. Same for hugo symbol and exon
info. Output resulting JSON"""

import pandas as pd
import numpy as np
import argparse


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


def add_refseq(transcripts, refseq, isoform_overrides_uniprot, isoform_overrides_mskcc):
    """Add one refseq id for each transcript. There can be multiple. Pick
    highest number transcript id in that case."""
    refseq.columns = [c.lower().replace(' ', '_') for c in refseq.columns]
    refseq_grouped = refseq.groupby("transcript_stable_id")

    def get_refseq_for_transcript(x):
        # use previously assigned uniprot refseq id
        try:
            override_refseq_id = isoform_overrides_uniprot.loc[x]['refseq_id']
            if hasattr(override_refseq_id, 'split'):
                return override_refseq_id.split('.')[0]
        except KeyError:
            pass

        # use previously assigned mskcc refseq id
        try:
            override_refseq_id = isoform_overrides_mskcc.loc[x]['refseq_id']
            if hasattr(override_refseq_id, 'split'):
                return override_refseq_id.split('.')[0]
        except KeyError:
            pass

        try:
            # Fix for: https://github.com/pandas-dev/pandas/issues/9466
            # use groupby to get indices, loc is slow for non-unique
            return refseq.iloc[refseq_grouped.groups[x]]['refseq_mrna_id'].sort_values(ascending=False).values[0]
        except KeyError:
            return np.nan

    transcripts['refseq_mrna_id'] = transcripts.index.map(get_refseq_for_transcript)
    return transcripts


def add_ccds(transcripts, ccds, isoform_overrides_uniprot, isoform_overrides_mskcc):
    """Add one ccds id for each transcript. There is only one per transcript."""
    ccds.columns = [c.lower().replace(' ', '_') for c in ccds.columns]
    ccds = ccds[~pd.isnull(ccds["ccds_id"])]
    # assume each transcript has only one CCDS
    assert(any(ccds["transcript_stable_id"].duplicated()) == False)
    ccds = ccds.set_index("transcript_stable_id")

    def get_ccds_for_transcript(x):
        try:
            override_ccdsid = isoform_overrides_uniprot.loc[x]['ccds_id']
            if hasattr(override_ccdsid, 'split'):
                return override_ccdsid.split('.')[0]
        except KeyError:
            pass

        try:
            override_ccdsid = isoform_overrides_mskcc.loc[x]['ccds_id']
            if hasattr(override_ccdsid, 'split'):
                return override_ccdsid.split('.')[0]
        except KeyError:
            pass

        try:
            return ccds.loc[x].values[0]
        except KeyError:
            return np.nan

    transcripts["ccds_id"] = transcripts.index.map(get_ccds_for_transcript)
    return transcripts


def add_uniprot(transcripts, uniprot):
    """Add one Uniprot id for each transcript. There is only one per transcript."""
    uniprot.columns = [u.lower().replace(' ', '_') for u in uniprot.columns]
    uniprot = uniprot[~pd.isnull(uniprot["reviewed_uniprot_accession"])]
    # assume each transcript has only one Uniprot ID
    assert(any(uniprot["enst_id"].duplicated()) == False)
    uniprot = uniprot.set_index("enst_id")

    def get_uniprot_for_transcript(x):
        try:
            return uniprot.loc[x].values[0]
        except KeyError:
            return np.nan

    transcripts["uniprot_id"] = transcripts.index.map(get_uniprot_for_transcript)
    return transcripts

def main(ensembl_biomart_transcripts,
         ensembl_transcript_info,
         ensembl_biomart_pfam,
         ensembl_biomart_refseq,
         ensembl_biomart_ccds,
         enst_to_uniprot,
         isoform_overrides_uniprot,
         isoform_overrides_mskcc,
         ensembl_biomart_transcripts_json):

    # Read input and set index column
    transcripts = pd.read_csv(ensembl_biomart_transcripts, sep='\t')
    transcript_info = pd.read_csv(ensembl_transcript_info, sep="\t")
    pfam_domains = pd.read_csv(ensembl_biomart_pfam, sep='\t')
    transcripts.set_index('transcript_stable_id', inplace=True)

    # import refseq
    refseq = pd.read_csv(ensembl_biomart_refseq, sep="\t")
    isoform_overrides_uniprot = pd.read_csv(isoform_overrides_uniprot, sep="\t").set_index('enst_id')
    isoform_overrides_mskcc = pd.read_csv(isoform_overrides_mskcc, sep="\t").set_index('enst_id')
    transcripts = add_refseq(transcripts, refseq, isoform_overrides_uniprot, isoform_overrides_mskcc)

    # import ccds
    ccds = pd.read_csv(ensembl_biomart_ccds, sep="\t")
    transcripts = add_ccds(transcripts, ccds, isoform_overrides_uniprot, isoform_overrides_mskcc)

    # Add nested HGNC, exons and PFAM domains
    transcripts = add_nested_hgnc(transcripts)
    transcripts = add_nested_transcript_info(transcripts, transcript_info)
    transcripts = add_nested_pfam_domains(transcripts, pfam_domains)

    # Add Uniprot id
    enst_to_uniprot_map = pd.read_csv(enst_to_uniprot, sep='\t')
    transcripts = add_uniprot(transcripts, enst_to_uniprot_map)

    # print records as json
    transcripts.reset_index().to_json(ensembl_biomart_transcripts_json,
                                      orient='records', lines=True, compression='gzip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ensembl_biomart_transcripts",
                        help="tmp/ensembl_biomart_transcripts.txt")
    parser.add_argument("ensembl_transcript_info",
                        help="tmp/ensembl_transcript_info.txt")
    parser.add_argument("ensembl_biomart_pfam",
                        help="input/ensembl_biomart_pfam.txt")
    parser.add_argument("ensembl_biomart_refseq",
                        help="input/ensembl_biomart_refseq.txt")
    parser.add_argument("ensembl_biomart_ccds",
                        help="input/ensembl_biomart_ccds.txt")
    parser.add_argument("enst_to_uniprot",
                        help="../uniprot/export/grch37_ensembl92_enst_to_uniprot_mapping_id.txt or ../uniprot/export/grch38_ensembl92_enst_to_uniprot_mapping_id.txt")
    parser.add_argument("vcf2maf_isoform_overrides_uniprot",
                        help="common_input/isoform_overrides_uniprot.txt")
    parser.add_argument("vcf2maf_isoform_overrides_mskcc",
                        help="common_input/isoform_overrides_at_mskcc.txt")
    parser.add_argument("ensembl_biomart_transcripts_json",
                        help="tmp/ensembl_biomart_transcripts.json.gz")

    args = parser.parse_args()
    main(args.ensembl_biomart_transcripts,
         args.ensembl_transcript_info,
         args.ensembl_biomart_pfam,
         args.ensembl_biomart_refseq,
         args.ensembl_biomart_ccds,
         args.enst_to_uniprot,
         args.vcf2maf_isoform_overrides_uniprot,
         args.vcf2maf_isoform_overrides_mskcc,
         args.ensembl_biomart_transcripts_json)

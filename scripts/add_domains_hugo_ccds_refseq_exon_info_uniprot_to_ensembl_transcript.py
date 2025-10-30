#!/usr/bin/env python3
"""Add domain as nested property to transcript. Same for hugo symbol and exon
info. Output resulting JSON"""

import pandas as pd
import numpy as np
import argparse

def main_transcript_id(enst_or_versioned: str) -> str:
    return enst_or_versioned.split(".", 1)[0] if isinstance(enst_or_versioned, str) else enst_or_versioned

def normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [c.lower().replace(" ", "_") for c in df.columns]
    # unify transcript_id_version column name from multiple sources
    if "transcript_id" in df.columns:
        df.rename(columns={"transcript_id": "transcript_stable_id"}, inplace=True)
    return df

def add_nested_hgnc(transcripts, hgnc_df):
    """ Make nested object HGNC symbols per transcript"""

    def get_hgnc_symbol(transcript_id, hgnc_dict):
        hgnc_symbols = transcripts.loc[transcript_id]
        if hgnc_symbols.ndim == 1:
            if pd.isnull(hgnc_symbols.hgnc_symbol):
                return hgnc_symbols.hgnc_symbol
            else:
                if (hgnc_symbols.hgnc_symbol in hgnc_dict):
                    return [hgnc_dict[hgnc_symbols.hgnc_symbol]]
                else:
                    return [hgnc_symbols.hgnc_symbol]
        else:
            symbol_list = []
            for symbol in hgnc_symbols.hgnc_symbol.values:
                if (symbol in hgnc_dict):
                    symbol_list.append(hgnc_dict[symbol])
                else:
                    symbol_list.append(symbol)
            return symbol_list

    hgnc_dict = dict()
    for index, row in hgnc_df.iterrows():
        for symbol in row['prev_symbol'].split('|'):
            hgnc_dict[symbol] = index

    hgnc_symbol_list = transcripts.index.drop_duplicates().map(lambda transcript_id: get_hgnc_symbol(transcript_id, hgnc_dict))
    # make one row per transcript_stable_id by removing hgnc_symbol
    unique_transcripts = transcripts.copy().reset_index()
    del unique_transcripts['hgnc_symbol']
    # remove duplicated transcripts
    unique_transcripts = unique_transcripts.drop_duplicates(subset=['versioned_transcript_id'], keep='first')
    # should only be one row per transcript now
    assert (len(unique_transcripts) == len(transcripts.index.drop_duplicates()))
    # should have the same order after dropping duplicates as
    # hgnc_symbol_list
    assert (0 == (unique_transcripts.versioned_transcript_id != transcripts.index.drop_duplicates()).sum())
    unique_transcripts['hgnc_symbols'] = hgnc_symbol_list
    return unique_transcripts.set_index("transcript_stable_id")


def add_nested_transcript_info(transcripts, transcript_info):
    """ Make nested object with exons and UTR per transcript. """

    def get_list_of_info_dicts(transcript_group):
        list_of_info_dicts = transcript_group.to_dict(orient='records')
        for info_dict in list_of_info_dicts:
            info_dict.pop('transcript_stable_id', None)

            # For UTRs, remove rank fields present only in Ensembl payload
            if info_dict.get('type') in ['five_prime_utr', 'three_prime_utr']:
                info_dict.pop('id', None)
                info_dict.pop('rank', None)
                info_dict.pop('version', None)

            if info_dict.get('type') == 'exon':
                info_dict.pop('type', None)

        return list_of_info_dicts

    # Split in exons and UTRs
    exon_info = transcript_info.loc[transcript_info.type.str.lower() == 'exon']
    utr_info = transcript_info.loc[transcript_info.type.str.lower().isin(['five_prime_utr', 'three_prime_utr'])]

    # Group by main transcript id
    series_of_lists_of_exon_dicts = exon_info.groupby('transcript_stable_id', dropna=False).apply(get_list_of_info_dicts)
    series_of_lists_of_utr_dicts  = utr_info.groupby('transcript_stable_id', dropna=False).apply(get_list_of_info_dicts)

    transcripts['exons'] = transcripts.index.map(series_of_lists_of_exon_dicts)
    transcripts['utrs']  = transcripts.index.map(series_of_lists_of_utr_dicts)
    return transcripts


def add_nested_pfam_domains(transcripts, pfam_domains):
    """ Add nested PFAM domains"""
    domain_grouped = pfam_domains.groupby("versioned_transcript_id", dropna=False)

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
    refseq = normalize_cols(refseq)
    refseq_grouped = refseq.groupby("versioned_transcript_id", dropna=False)

    def get_refseq_for_transcript(enst_versioned):
        bare = main_transcript_id(enst_versioned)
        # Maybe retire uniprot?
        for ov in (isoform_overrides_mskcc, isoform_overrides_uniprot):
            try:
                override_refseq_id = ov.loc[bare, 'refseq_id']
                if isinstance(override_refseq_id, str) and override_refseq_id:
                    return override_refseq_id.split('.')[0]
            except KeyError:
                pass
        # from biomart table (by ENST.X)
        try:
            return refseq.iloc[refseq_grouped.groups[enst_versioned]]['refseq_mrna_id'].sort_values(ascending=False).values[0]
        except KeyError:
            return np.nan

    transcripts['refseq_mrna_id'] = transcripts.index.map(get_refseq_for_transcript)
    return transcripts


def add_ccds(transcripts, ccds, isoform_overrides_uniprot, isoform_overrides_mskcc):
    """Add one ccds id for each transcript. There is only one per transcript."""
    ccds = normalize_cols(ccds)
    ccds = ccds[~pd.isnull(ccds["ccds_id"])]
    # assume each transcript ENST.X has only one CCDS
    assert not any(ccds["versioned_transcript_id"].duplicated())
    ccds = ccds.set_index("versioned_transcript_id")

    def get_ccds_for_transcript(enst_versioned):
        bare = main_transcript_id(enst_versioned)
        # Maybe retire uniprot?
        for ov in (isoform_overrides_mskcc, isoform_overrides_uniprot):
            try:
                override_ccdsid = ov.loc[bare, 'ccds_id']
                if isinstance(override_ccdsid, str) and override_ccdsid:
                    return override_ccdsid.split('.')[0]
            except KeyError:
                pass
        try:
            return ccds.loc[enst_versioned].values[0]
        except KeyError:
            return np.nan

    transcripts["ccds_id"] = transcripts.index.map(get_ccds_for_transcript)
    return transcripts


def add_uniprot(transcripts, uniprot):
    """Add one Uniprot id for each transcript.
       Prefer versioned join if file provides versioned transcript ID,
       else fall back to bare ENST join."""
    uniprot = normalize_cols(uniprot)
    uniprot = uniprot[~pd.isnull(uniprot["final_uniprot_id"])]

    if "versioned_transcript_id" in uniprot.columns:
        # prefer versioned mapping
        assert not any(uniprot["versioned_transcript_id"].duplicated())
        uniprot = uniprot.set_index("versioned_transcript_id")
        def get_up(enst_versioned):
            try:
                return uniprot.loc[enst_versioned].values[0]
            except KeyError:
                return np.nan
        transcripts["uniprot_id"] = transcripts.index.map(get_up)
    else:
        # legacy: only enst_id (bare)
        assert not any(uniprot["enst_id"].duplicated())
        uniprot = uniprot.set_index("enst_id")
        def get_up(enst_versioned):
            try:
                return uniprot.loc[main_transcript_id(enst_versioned)].values[0]
            except KeyError:
                return np.nan
        transcripts["uniprot_id"] = transcripts.index.map(get_up)

    return transcripts

def main(ensembl_biomart_transcripts,
         ensembl_transcript_info,
         ensembl_biomart_pfam,
         ensembl_biomart_refseq,
         ensembl_biomart_ccds,
         enst_to_uniprot,
         isoform_overrides_uniprot,
         isoform_overrides_mskcc,
         hgnc_symbol_set,
         ensembl_biomart_transcripts_json
         ):

    # Read & normalize
    transcripts = normalize_cols(pd.read_csv(ensembl_biomart_transcripts, sep="\t", dtype=str))
    transcript_info = normalize_cols(pd.read_csv(ensembl_transcript_info, sep="\t", dtype=str))
    pfam_domains = normalize_cols(pd.read_csv(ensembl_biomart_pfam, sep="\t", dtype=str))

    # Make sure we have full version and bare id in transcripts
    # Expect transcripts to already include transcript_id_version (ENST.X) and transcript_stable_id (bare)
    if "versioned_transcript_id" not in transcripts.columns:
        raise RuntimeError("transcripts file missing versioned_transcript_id")

    # Use versioned id as the index for all downstream joins that support versions
    transcripts.set_index("versioned_transcript_id", inplace=True, drop=True)

    # import refseq
    refseq = pd.read_csv(ensembl_biomart_refseq, sep="\t")
    isoform_overrides_uniprot = pd.read_csv(isoform_overrides_uniprot, sep="\t").set_index('enst_id')
    isoform_overrides_mskcc = pd.read_csv(isoform_overrides_mskcc, sep="\t").set_index('enst_id')
    transcripts = add_refseq(transcripts, refseq, isoform_overrides_uniprot, isoform_overrides_mskcc)

    # import ccds
    ccds = pd.read_csv(ensembl_biomart_ccds, sep="\t")
    transcripts = add_ccds(transcripts, ccds, isoform_overrides_uniprot, isoform_overrides_mskcc)

    # Add nested PFAM domains
    transcripts = add_nested_pfam_domains(transcripts, pfam_domains)

    # Add nested HGNC, exons and 
    hgnc_df = pd.read_csv(hgnc_symbol_set, sep='\t', usecols = ['symbol', 'prev_symbol'], index_col=0).dropna()
    transcripts = add_nested_hgnc(transcripts, hgnc_df)

    # transcripts.index is main id from now
    transcripts = add_nested_transcript_info(transcripts, transcript_info)

    # Add Uniprot id
    enst_to_uniprot_map = pd.read_csv(enst_to_uniprot, sep='\t')
    transcripts = add_uniprot(transcripts, enst_to_uniprot_map)

    transcripts = transcripts.copy()
    # print records as json
    transcripts.drop(columns=['versioned_transcript_id'], inplace=True)

    # show 1 instead of 1.0
    transcripts["transcript_id_version"] = transcripts["transcript_id_version"].astype(float).astype(int).astype(str)
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
                        help="../uniprot/export/grch37_ensembl92_enst_to_uniprot_mapping_id.txt or ../uniprot/export/grch38_ensembl92_enst_to_uniprot_mapping_id.txt or ../uniprot/export/grch38_ensembl95_enst_to_uniprot_mapping_id.txt")
    parser.add_argument("vcf2maf_isoform_overrides_uniprot",
                        help="common_input/isoform_overrides_uniprot.txt")
    parser.add_argument("vcf2maf_isoform_overrides_mskcc",
                        help="common_input/mskcc_isoform_overrides_grch37.txt")
    parser.add_argument("hgnc_symbol_set", help="common_input/hgnc_v2024.10.1.txt")
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
         args.hgnc_symbol_set,
         args.ensembl_biomart_transcripts_json
         )

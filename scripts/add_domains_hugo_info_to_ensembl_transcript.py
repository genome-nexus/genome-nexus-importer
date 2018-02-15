"""Add domain as nested property to transcript. Same for hugo symbol. Output resulting JSON"""

import pandas as pd
import numpy as np
import sys


if __name__ == "__main__":
    transcripts = pd.read_csv("../data/ensembl_biomart_transcripts.txt", sep='\t')

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
    unique_transcripts['hgnc_symbol'] = hgnc_symbol_list

    # add nested domain
    trans_domain =  pd.read_csv("../data/ensembl_biomart_pfam_grch37.p13.txt", sep='\t')
    trans_domain.columns = [c.lower().replace(' ','_') for c in trans_domain.columns]
    domain_grouped = trans_domain.groupby("transcript_stable_id")
    unique_transcripts = unique_transcripts.set_index("transcript_stable_id")
    def get_domain_for_transcript(x):
        try:
            # Fix for: https://github.com/pandas-dev/pandas/issues/9466
            # use groupby to get indices, loc is slow for non-unique
            return trans_domain.iloc[domain_grouped.groups[x]]['pfam_domain_id pfam_domain_start pfam_domain_end'.split()]
        except KeyError:
            return np.nan
    unique_transcripts["domains"] = unique_transcripts.index.map(get_domain_for_transcript)
    unique_transcripts.reset_index().to_json(sys.stdout, orient='records', lines=True)

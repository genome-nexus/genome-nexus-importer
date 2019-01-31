#!/usr/bin/env python3
import pandas as pd
import numpy as np
import itertools
import argparse


def get_overrides_transcript(overrides_tables, ensembl_table, hgnc_symbol, hgnc_canonical_genes):
    """Find canonical transcript id for given hugo symbol. Overrides_tables is
    a list of different override tables"""
    for overrides in overrides_tables:
        try:
            # corner case when there are multiple overrides for a given gene symbol
            if overrides.loc[hgnc_symbol].ndim > 1:
                transcript = overrides.loc[hgnc_symbol].isoform_override.values[0]
            else:
                transcript = overrides.loc[hgnc_symbol].isoform_override
            return transcript
        except KeyError:
            pass
    else:
        # get ensembl canonical version otherwise
        return get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(ensembl_table, hgnc_symbol, hgnc_canonical_genes, 'transcript_stable_id')


def pick_canonical_longest_transcript_from_ensembl_table(ensembl_rows, field):
    """Get canonical transcript id with largest protein length or if there is
    no such thing, pick biggest gene id"""
    return ensembl_rows.sort_values('is_canonical protein_length gene_stable_id'.split(), ascending=False)[field].values[0]

def get_ensembl_canonical(ensembl_rows, field):
    if (ensembl_rows.ndim == 1 and len(ensembl_rows) == 0) or (ensembl_rows.ndim == 2 and len(ensembl_rows) == 0):
        return np.nan
    elif ensembl_rows.ndim == 1:
        return ensembl_rows[field]
    elif ensembl_rows.ndim == 2:
        if len(ensembl_rows) == 1:
            return ensembl_rows[field].values[0]
        else:
            return pick_canonical_longest_transcript_from_ensembl_table(ensembl_rows, field)


def get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(ensembl_table, hgnc_symbol, hgnc_canonical_genes, field):
    """Determine canonical transcript based on hgnc mappings to ensembl id.
    If not possible use ensembl's data."""
    try:
        hgnc_gene_rows = hgnc_canonical_genes.loc[hgnc_symbol]
    except KeyError:
        raise(Exception("Unknown hugo symbol"))

    if hgnc_gene_rows.ndim == 1:
        result = get_ensembl_canonical(ensembl_table[ensembl_table.gene_stable_id == hgnc_gene_rows.ensembl_gene_id], field)
        if pd.isnull(result):
            # use ensembl data to find canonical transcript if nan is found
            # there's actually 222 of these (see notebook)
            try:
                return get_ensembl_canonical(ensembl_table.loc[hgnc_symbol], field)
            except KeyError:
                return np.nan
        else:
            return result
    else:
        raise(Exception("One hugo symbol expected in hgnc_canonical_genes"))

def lowercase_set(x):
    return set({i.lower() for i in x})

def ignore_rna_gene(x):
    return set({i for i in x if not i.startswith('rn') and not i.startswith('mir') and not i.startswith('linc')})

def ignore_certain_genes(x):
    ignore_genes = {'fam25hp', 'rmrpp1', 'hla-as1', 'pramef16', 'arid4b-it1',
    'ctglf8p', 'htr4-it1', 'nicn1-as1', 'pramef3', 'c9orf38', 'tbc1d4-as1',
    'rmrpp3', 'magi2-it1', 'rmrpp2', 'rmrpp4', 'ercc6-pgbd3', 'tbce-as1', 'hpvc1', 'fam231a', 'fam231b', 'fam231c', 'fam231d'}

    return set({i for i in x if i not in ignore_genes})


def main(ensembl_biomart_geneids_transcript_info,
         hgnc_complete_set,
         isoform_overrides_uniprot,
         isoform_overrides_at_mskcc,
         isoform_overrides_genome_nexus,
         ensembl_biomart_canonical_transcripts_per_hgnc):

    # input files
    tsv = pd.read_csv(ensembl_biomart_geneids_transcript_info, sep='\t', dtype={'is_canonical':bool})
    tsv = tsv.drop_duplicates()
    uniprot = pd.read_csv(isoform_overrides_uniprot, sep='\t')\
        .rename(columns={'enst_id':'isoform_override'})\
        .set_index('gene_name'.split())
    mskcc = pd.read_csv(isoform_overrides_at_mskcc, sep='\t')\
        .rename(columns={'enst_id':'isoform_override'})\
        .set_index('gene_name'.split())
    custom = pd.read_csv(isoform_overrides_genome_nexus, sep='\t')\
        .rename(columns={'enst_id':'isoform_override'})\
        .set_index('gene_name'.split())
    hgnc = pd.read_csv(hgnc_complete_set, sep='\t', dtype=object)
    hgnc = hgnc[hgnc['name'] != 'entry withdrawn'].copy()
    hugos = hgnc['symbol'].unique()
    hgnc = hgnc.set_index('symbol')
    # assume each row has approved symbol
    assert(len(hugos) == len(hgnc))

    # only test the cancer genes for oddities (these are very important)
    cgs = set(pd.read_csv('common_input/oncokb_cancer_genes_list_20170926.txt',sep='\t')['Hugo Symbol'])
    # each cancer gene stable id should have only one associated cancer gene symbol
    assert(tsv[tsv.hgnc_symbol.isin(cgs)].groupby('gene_stable_id').hgnc_symbol.nunique().sort_values().nunique() == 1)
    # each transcript stable id always belongs to only one gene stable id
    assert(tsv.groupby('transcript_stable_id').gene_stable_id.nunique().sort_values().nunique() == 1)

    # create hgnc_symbol to gene id mapping
    # ignore hugo symbols from ensembl data dump (includes prev symbols and alias_symbol)
    syns = hgnc.alias_symbol.str.strip('"').str.split("|").dropna()
    syns = set(itertools.chain.from_iterable(syns))
    prev_symbol = hgnc.prev_symbol.str.strip('"').str.split("|").dropna()
    prev_symbol = set(itertools.chain.from_iterable(prev_symbol))

    # there is overlap between symbols, alias_symbol and previous symbols
    # therefore use logic in above order when querying
    # assert(len(syns.intersection(prev_symbol)) == 0) #329
    # assert(len(set(hugos).intersection(syns)) == 0) # 495
    # assert(len(set(hugos).intersection(prev_symbol)) == 0) #227

    # all cancer genes and hugo symbols in ensembl data dump should be
    # contained in hgnc approved symbols and alias_symbol
    # c12orf9 is only in sanger's cancer gene census and has been withdrawn
    assert(len(lowercase_set(set(cgs)) - set(['c12orf9']) - lowercase_set(set(hugos).union(syns).union(prev_symbol))) == 0)
    no_symbols_in_hgnc = lowercase_set(tsv.hgnc_symbol.dropna().unique()) - lowercase_set(set(hugos).union(syns).union(prev_symbol))
    assert(len(ignore_certain_genes(ignore_rna_gene(no_symbols_in_hgnc))) == 0)

    tsv = tsv.set_index('hgnc_symbol')
    # for testing use
    # NSD3 replaces WHSC1L1
    # AATF has uniprot canonical transcript, not hgnc ensembl gene id, but
    # there are multiple in ensembl data dump
    # hugos = ['KRT18P53', 'NSD3', 'AATF']
    one_transcript_per_hugo_symbol = pd.Series(hugos).apply(lambda x:
        pd.Series(
            [
                get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(tsv, x, hgnc, 'gene_stable_id'),
                get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(tsv, x, hgnc, 'transcript_stable_id'),
                get_overrides_transcript([custom], tsv, x, hgnc),
                get_overrides_transcript([uniprot, custom], tsv, x, hgnc),
                get_overrides_transcript([mskcc, uniprot, custom], tsv, x, hgnc),
            ],
            index="""
            ensembl_canonical_gene
            ensembl_canonical_transcript
            genome_nexus_canonical_transcript
            uniprot_canonical_transcript
            mskcc_canonical_transcript
            """.split()
        )
    )
    one_transcript_per_hugo_symbol.index = hugos
    one_transcript_per_hugo_symbol.index.name = 'hgnc_symbol'

    # merge in other hgnc fields
    merged = pd.merge(one_transcript_per_hugo_symbol.reset_index(), hgnc.reset_index(), left_on='hgnc_symbol', right_on='symbol')
    del merged['symbol']
    del merged['ensembl_gene_id']
    merged.to_csv(ensembl_biomart_canonical_transcripts_per_hgnc, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ensembl_biomart_geneids_transcript_info",
                        help="tmp/ensembl_biomart_geneids.transcript_info.txt")
    parser.add_argument("hgnc_complete_set",
                        help="common_input/hgnc_complete_set_20190122.txt")
    parser.add_argument("isoform_overrides_uniprot",
                        help="common_input/isoform_overrides_uniprot.txt")
    parser.add_argument("isoform_overrides_at_mskcc",
                        help="common_input/isoform_overrides_at_mskcc.txt")
    parser.add_argument("isoform_overrides_genome_nexus",
                        help="common_input/isoform_overrides_genome_nexus.txt")
    parser.add_argument("ensembl_biomart_canonical_transcripts_per_hgnc",
                        help="tmp/ensembl_biomart_canonical_transcripts_per_hgnc.txt")
    args = parser.parse_args()

    main(args.ensembl_biomart_geneids_transcript_info,
         args.hgnc_complete_set,
         args.isoform_overrides_uniprot,
         args.isoform_overrides_at_mskcc,
         args.isoform_overrides_genome_nexus,
         args.ensembl_biomart_canonical_transcripts_per_hgnc)

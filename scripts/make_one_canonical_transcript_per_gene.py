#!/usr/bin/env python3
import pandas as pd
import numpy as np
import itertools
import argparse

def normalize_version_string(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none"}:
        return None
    if s.endswith(".0"):
        s = s[:-2]
    return s if s.isdigit() else None

def split_enst_id(enst):
    # split enst id, e.g. ENST00000123456.7 -> 'ENST00000123456','7'
    if enst is None or pd.isna(enst):
        return None, None
    s = str(enst).strip()
    if "." in s:
        b, v = s.split(".", 1)
        return b, normalize_version_string(v)
    return s, None

def find_transcript_version(transcript_info_df, base_enst):
    if not base_enst:
        return np.nan
    rows = transcript_info_df[transcript_info_df['transcript_stable_id'] == base_enst]
    if rows.empty:
        return np.nan
    # Prefer transcript_id_version, else parse from versioned_transcript_id
    if 'transcript_id_version' in rows:
        version = rows['transcript_id_version'].apply(normalize_version_string).dropna()
    else:
        version = pd.Series([], dtype=object)
    if version.empty and 'versioned_transcript_id' in rows:
        version = rows['versioned_transcript_id'].astype(str).str.extract(r'\.(\d+)$')[0].dropna()
    uniq = sorted({v for v in version if isinstance(v, str) and v.isdigit()}, key=lambda x: int(x))
    return uniq[-1] if uniq else np.nan

def get_overrides_transcript(overrides_tables, ensembl_table, ensembl_table_indexed_by_gene_stable_id, hgnc_symbol, hgnc_canonical_genes, overrides_table_names):
    """Find canonical transcript id for given hugo symbol. Overrides_tables is
    a list of different override tables"""
    for index in range(len(overrides_tables)):
        overrides = overrides_tables[index]
        try:
            # corner case when there are multiple overrides for a given gene symbol
            if overrides.loc[hgnc_symbol].ndim > 1:
                transcript = overrides.loc[hgnc_symbol].isoform_override.values[0]
                isoform_override = overrides_table_names[index]
            else:
                transcript = overrides.loc[hgnc_symbol].isoform_override
                isoform_override = overrides_table_names[index]
            # return transcript and isoform_override_explanation
            return transcript, isoform_override
        except KeyError:
            pass
    # get ensembl canonical version otherwise
    return get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(ensembl_table, ensembl_table_indexed_by_gene_stable_id, hgnc_symbol, hgnc_canonical_genes, 'transcript_stable_id')

def pick_canonical_longest_transcript_from_ensembl_table(ensembl_rows, field):
    """Get canonical transcript id with largest protein length or if there is
    no such thing, pick biggest gene id"""
    return ensembl_rows.sort_values('is_canonical protein_length gene_stable_id'.split(), ascending=False)[field].values[0], "ensembl longest"

def get_ensembl_canonical(ensembl_rows, field):

    if (ensembl_rows.ndim == 1 and len(ensembl_rows) == 0) or (ensembl_rows.ndim == 2 and len(ensembl_rows) == 0):
        return np.nan, np.nan
    elif ensembl_rows.ndim == 1:
        # this gene onlyhas one transcript from ensembl_table
        return ensembl_rows[field], "ensembl only one transcript"
    elif ensembl_rows.ndim == 2:
        # there are multiple transcripts found by hugo_symbol or ensembl_gene_id
        # we need to pick the longest transcript in this case
        return pick_canonical_longest_transcript_from_ensembl_table(ensembl_rows, field)

def get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(ensembl_table, ensembl_table_indexed_by_gene_stable_id, hgnc_symbol, hgnc_canonical_genes, field):
    """Determine canonical transcript based on hgnc mappings to ensembl id.
    If not possible use ensembl's data.
    ensembl_table is the same as ensembl_table_indexed_by_gene_stable_id
    But ensembl_table has hugo_symbol as index
    ensembl_table_indexed_by_gene_stable_id has gene_stable_id as index
    Adding ensembl_table_indexed_by_gene_stable_id to improve performance since this is the most time consuming part
    """
    try:
        # try to find hgnc_symbol from hgnc_canonical_genes
        # it should only have one row returned, otherwise raise an exception
        hgnc_gene_rows = hgnc_canonical_genes.loc[hgnc_symbol]
    except KeyError:
        raise(Exception("Unknown hugo symbol"))
    if hgnc_gene_rows.ndim == 1:
        try:
            # from the only one record found in hgnc_canonical_genes, we could have the ensembl_gene_id from the record
            # find transcripts from ensembl_table (same as ensembl_table_indexed_by_gene_stable_id) by searching ensembl_gene_id
            # here we use ensembl_table_indexed_by_gene_stable_id instead to increase performance
            # get_ensembl_canonical returns one transcript from all transcripts found by ensembl_gene_id
            return get_ensembl_canonical(ensembl_table_indexed_by_gene_stable_id.loc[hgnc_gene_rows.ensembl_gene_id], field)
        except KeyError:
            # if couldn't find any transcripts by ensembl_gene_id, switch to searching by hgnc_symbol
            # there's actually 222 of these (see notebook)
            try:
                return get_ensembl_canonical(ensembl_table.loc[hgnc_symbol], field)
            except KeyError:
                return np.nan, np.nan
    else:
        raise(Exception("One hugo symbol expected in hgnc_canonical_genes"))

def get_transcript_id_and_explanation(transcript_info_df, transcript_info_indexed_by_gene_stable_id, hugo_symbol, hgnc_df, mskcc, uniprot, custom):
    ensembl_canonical_gene, ensembl_canonical_gene_explanation = get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(transcript_info_df, transcript_info_indexed_by_gene_stable_id, hugo_symbol, hgnc_df, 'gene_stable_id')
    ensembl_canonical_transcript, ensembl_canonical_transcript_explanation = get_ensembl_canonical_transcript_id_from_hgnc_then_ensembl(transcript_info_df, transcript_info_indexed_by_gene_stable_id, hugo_symbol, hgnc_df, 'versioned_transcript_id')
    genome_nexus_overrides_transcript, genome_nexus_overrides_transcript_explanation = get_overrides_transcript([custom], transcript_info_df, transcript_info_indexed_by_gene_stable_id, hugo_symbol, hgnc_df, ["genome nexus isoform override"])
    uniprot_overrides_transcript, uniprot_overrides_transcript_explanation = get_overrides_transcript([uniprot], transcript_info_df, transcript_info_indexed_by_gene_stable_id, hugo_symbol, hgnc_df, ["uniprot isoform override"])
    mskcc_overrides_transcript, mskcc_overrides_transcript_explanation = get_overrides_transcript([mskcc], transcript_info_df, transcript_info_indexed_by_gene_stable_id, hugo_symbol, hgnc_df, ["mskcc isoform override"])

    ensembl_transcript_id_base, ensembl_transcript_id_version = split_enst_id(ensembl_canonical_transcript)
    gn_transcript_id_base,  gn_transcript_id_version  = split_enst_id(genome_nexus_overrides_transcript)
    uniprot_transcript_id_base, uniprot_transcript_id_version = split_enst_id(uniprot_overrides_transcript)
    mskcc_transcript_id_base, mskcc_transcript_id_version = split_enst_id(mskcc_overrides_transcript)

    if not ensembl_transcript_id_version:
        ensembl_transcript_id_version = find_transcript_version(transcript_info_df, ensembl_transcript_id_base)
    if not gn_transcript_id_version:
        gn_transcript_id_version = find_transcript_version(transcript_info_df, gn_transcript_id_base)
    if not uniprot_transcript_id_version:
        uniprot_transcript_id_version = find_transcript_version(transcript_info_df, uniprot_transcript_id_base)
    if not mskcc_transcript_id_version:
        mskcc_transcript_id_version = find_transcript_version(transcript_info_df, mskcc_transcript_id_base)

    try:
        if 'note' in mskcc.columns:
            m = mskcc.loc[hugo_symbol]
            note = m['note'].iloc[0] if isinstance(m, pd.DataFrame) else m.get('note')
            if isinstance(note, str) and note.strip():
                mskcc_overrides_transcript_explanation = f"{mskcc_overrides_transcript_explanation}; {note.strip()}"
    except KeyError:
        pass

    return pd.Series([
                ensembl_canonical_gene,
                ensembl_transcript_id_base,
                ensembl_transcript_id_version,
                ensembl_canonical_transcript_explanation,
                gn_transcript_id_base,
                gn_transcript_id_version,
                genome_nexus_overrides_transcript_explanation,
                uniprot_transcript_id_base,
                uniprot_transcript_id_version,
                uniprot_overrides_transcript_explanation,
                mskcc_transcript_id_base,
                mskcc_transcript_id_version,
                mskcc_overrides_transcript_explanation
            ],
            index="""
            ensembl_canonical_gene
            ensembl_canonical_transcript
            ensembl_canonical_transcript_version
            ensembl_canonical_transcript_explanation
            genome_nexus_canonical_transcript
            genome_nexus_canonical_transcript_version
            genome_nexus_canonical_transcript_explanation
            uniprot_canonical_transcript
            uniprot_canonical_transcript_version
            uniprot_canonical_transcript_explanation
            mskcc_canonical_transcript
            mskcc_canonical_transcript_version
            mskcc_canonical_transcript_explanation
            """.split()
        )

def lowercase_set(x):
    return set({i.lower() for i in x})


def ignore_rna_gene(x):
    return set({i for i in x if not i.startswith('rn') and not i.startswith('mir') and not i.startswith('linc')})


def ignore_certain_genes(x, ignored_genes_file_name):
    with open(ignored_genes_file_name) as ignored_genes_file:
        ignore_genes = [gene.rstrip() for gene in ignored_genes_file]
        return set({i for i in x if i not in ignore_genes})


def main(ensembl_biomart_geneids_transcript_info,
         hgnc_complete_set,
         isoform_overrides_uniprot,
         isoform_overrides_at_mskcc,
         isoform_overrides_genome_nexus,
         ignored_genes_file_name,
         ensembl_biomart_canonical_transcripts_per_hgnc):
    # input files
    transcript_info_df = pd.read_csv(ensembl_biomart_geneids_transcript_info, sep='\t', dtype={'is_canonical':bool, 'transcript_id_version': object})
    if 'transcript_id_version' in transcript_info_df.columns:
        transcript_info_df['transcript_id_version'] = transcript_info_df['transcript_id_version'].apply(normalize_version_string)
    transcript_info_df = transcript_info_df.drop_duplicates()
    uniprot = pd.read_csv(isoform_overrides_uniprot, sep='\t')\
        .rename(columns={'enst_id':'isoform_override'})\
        .set_index('gene_name'.split())
    mskcc = pd.read_csv(isoform_overrides_at_mskcc, sep='\t')\
        .rename(columns={'enst_id':'isoform_override'})\
        .set_index('gene_name'.split())
    custom = pd.read_csv(isoform_overrides_genome_nexus, sep='\t')\
        .rename(columns={'enst_id':'isoform_override'})\
        .set_index('gene_name'.split())
    hgnc_df = pd.read_csv(hgnc_complete_set, sep='\t', dtype=object)

    # Convert new column names to old stable column names. If this is not done properly, Genome Nexus and any other
    # downstream applications break
    # TODO: Update Genome Nexus to accept the latest HGNC column names so that remapping is not necessary.
    column_name_mapping = {'name': 'approved_name',
                           'symbol': 'approved_symbol',
                           'prev_symbol': 'previous_symbols',
                           'alias_symbol': 'synonyms',
                           'location': 'chromosome',
                           'entrez_id': 'entrez_gene_id',
                           'ena': 'accession_numbers',
                           'refseq_accession': 'refseq_ids',
                           'uniprot_ids': 'uniprot_id',
                           'ensembl_id': 'ensembl_gene_id'}
    hgnc_df.rename(columns=column_name_mapping, inplace=True)
    hgnc_df = hgnc_df[hgnc_df['approved_name'] != 'entry withdrawn'].copy()
    hugos = hgnc_df['approved_symbol'].unique()
    hgnc_df = hgnc_df.set_index('approved_symbol')
    # assume each row has approved symbol
    assert(len(hugos) == len(hgnc_df))

    # only test the cancer genes for oddities (these are very important)
    cgs = set(pd.read_csv('common_input/oncokb_cancer_genes_list.txt',sep='\t')['Hugo Symbol'])
    # each cancer gene stable id should have only one associated cancer gene symbol
    assert(transcript_info_df[transcript_info_df.hgnc_symbol.isin(cgs)].groupby('gene_stable_id').hgnc_symbol.nunique().sort_values().nunique() == 1)
    # each transcript stable id always belongs to only one gene stable id
    assert(transcript_info_df.groupby('transcript_stable_id').gene_stable_id.nunique().sort_values().nunique() == 1)

    # create hgnc_symbol to gene id mapping
    # ignore hugo symbols from ensembl data dump (includes prev symbols and synonyms)
    syns = hgnc_df.synonyms.str.strip('"').str.split("|").dropna()
    syns = set(itertools.chain.from_iterable(syns))
    previous_symbols = hgnc_df.previous_symbols.str.strip('"').str.split("|").dropna()
    previous_symbols = set(itertools.chain.from_iterable(previous_symbols))

    # there is overlap between symbols, synonyms and previous symbols
    # therefore use logic in above order when querying
    # assert(len(syns.intersection(previous_symbols)) == 0) #329
    # assert(len(set(hugos).intersection(syns)) == 0) # 495
    # assert(len(set(hugos).intersection(previous_symbols)) == 0) #227

    # all cancer genes and hugo symbols in ensembl data dump should be
    # contained in hgnc approved symbols and synonyms
    # c12orf9 is only in sanger's cancer gene census and has been withdrawn
    assert(len(lowercase_set(set(cgs)) - set(['c12orf9']) - lowercase_set(set(hugos).union(syns).union(previous_symbols))) == 0)
    no_symbols_in_hgnc = lowercase_set(transcript_info_df.hgnc_symbol.dropna().unique()) - lowercase_set(set(hugos).union(syns).union(previous_symbols))
    new_genes = ignore_certain_genes(ignore_rna_gene(no_symbols_in_hgnc),ignored_genes_file_name)
    if len(new_genes) != 0:
        print('------ New genes need to be added into ignored_genes.txt ------\n' +
        '------ Start of new genes list ------\n' +
        '\n'.join(new_genes) + '\n'
        '------ End of new genes list ------\n')
    assert(len(new_genes) == 0)

    transcript_info_df = transcript_info_df.set_index('hgnc_symbol').sort_index()
    transcript_info_indexed_by_gene_stable_id = transcript_info_df
    # duplicate a temp column to use as index
    transcript_info_indexed_by_gene_stable_id["gene_stable_id_temp"] = transcript_info_df['gene_stable_id']
    transcript_info_indexed_by_gene_stable_id = transcript_info_indexed_by_gene_stable_id.set_index('gene_stable_id_temp').sort_index()

    # for testing use
    # NSD3 replaces WHSC1L1
    # AATF has uniprot canonical transcript, not hgnc ensembl gene id, but
    # there are multiple in ensembl data dump
    # hugos = ['KRT18P53', 'NSD3', 'AATF']

    # TODO Optimize this part, as this part takes most time
    one_transcript_per_hugo_symbol = pd.Series(hugos).apply(lambda x: get_transcript_id_and_explanation(transcript_info_df, transcript_info_indexed_by_gene_stable_id, x, hgnc_df, mskcc, uniprot, custom ))
    one_transcript_per_hugo_symbol.index = hugos
    one_transcript_per_hugo_symbol.index.name = 'hgnc_symbol'

    # merge in other hgnc fields
    merged = pd.merge(one_transcript_per_hugo_symbol.reset_index(), hgnc_df.reset_index(), left_on='hgnc_symbol', right_on='approved_symbol')
    del merged['approved_symbol']
    del merged['ensembl_gene_id']

    # Replace '|' to ', ' to be in the correct format for Genome Nexus
    # TODO: Update Genome Nexus to accept the latest HGNC format so that replacement is not necessary.
    merged = merged.astype(str).replace({'\\|': ', '}, regex=True)

    # Write file
    merged.to_csv(ensembl_biomart_canonical_transcripts_per_hgnc, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ensembl_biomart_geneids_transcript_info",
                        help="tmp/ensembl_biomart_geneids.transcript_info.txt")
    parser.add_argument("hgnc_complete_set",
                        help="common_input/hgnc_v2024.10.1.txt")
    parser.add_argument("isoform_overrides_uniprot",
                        help="common_input/isoform_overrides_uniprot.txt")
    parser.add_argument("isoform_overrides_at_mskcc",
                        help="common_input/isoform_overrides_mskcc_grch37.txt or common_input/isoform_overrides_at_mskcc_grch38.txt")
    parser.add_argument("isoform_overrides_genome_nexus",
                        help="common_input/isoform_overrides_genome_nexus_grch37.txt or common_input/isoform_overrides_genome_nexus_grch38.txt")
    parser.add_argument("ignored_genes_file_name",
                        help="common_input/ignored_genes.txt")
    parser.add_argument("ensembl_biomart_canonical_transcripts_per_hgnc",
                        help="tmp/ensembl_biomart_canonical_transcripts_per_hgnc.txt")
    args = parser.parse_args()

    main(args.ensembl_biomart_geneids_transcript_info,
         args.hgnc_complete_set,
         args.isoform_overrides_uniprot,
         args.isoform_overrides_at_mskcc,
         args.isoform_overrides_genome_nexus,
         args.ignored_genes_file_name,
         args.ensembl_biomart_canonical_transcripts_per_hgnc)

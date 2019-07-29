#!/usr/bin/env python3

import pandas as pd
import argparse

def exons_per_transcript(exons):
    '''Builds a nested data frame from exon file for JSON output
    Structure: { transcript_id: { exons: [ { id: 'abc', start: 123, end: 456, rank: 1, strand: 1, version: 'x'}, {...} ] } }
    '''

    all_exons = {}

    for i, row in enumerate(exons.itertuples(index=False)):
        if i % 100000 == 0:
            print("Processed %i exons..." %i)

        transcript_id, exontype, exonid, start, end, rank, strand, version = row
        if exontype == 'exon':
            exon = {'id': exonid, 'start': int(start), 'end': int(end), 'rank': rank, 'strand': int(strand), 'version': version}
            all_exons.setdefault(transcript_id, {}).setdefault('exons', []).append(exon)
        elif exontype in ['five_prime_UTR', 'three_prime_UTR']:
            utr = {'start': int(start), 'end': int(end), 'strand': int(strand)}
            all_exons.setdefault(transcript_id, {}).setdefault('utrs', []).append(utr)
    
    return(pd.DataFrame.from_dict(all_exons, orient='index'))


def pfam_domains_per_transcript(pfam):
  
    all_domains = {}
    for i, row in enumerate(pfam.itertuples(index=False)):
        gene_id, transcript_id, symbol, pfam_id, start, end = row
        
        # Skip if no domain
        if pd.isna(pfam_id):
            continue

        all_domains.setdefault(transcript_id, {}).setdefault('domains', [])\
                                .append({'pfam_domain_id': pfam_id, 'pfam_domain_start': int(start), 'pfam_domain_end': int(end)})

    return(pd.DataFrame.from_dict(all_domains, orient='index'))


def combine_tables(transcripts, refseq, exons, pfam, ccds):
    '''Merges the relevant tables on transcript ID
    '''
 
    table = transcripts.copy()

    table.index.set_names(['transcript_id'])

    # put symbols into a list, and drop the original column
    table['hgnc_symbols'] = pd.Series(table['hgnc_symbol'].map(lambda x: [x] if not pd.isna(x) else []))
    table = table.drop(['hgnc_symbol'], axis=1)

    print('Merging RefSeq IDs...')
    table = table.merge(refseq, left_index=True, right_index=True, how='left').rename(columns={'RefSeq mRNA ID': 'refseq_mrna_id'})
    print('Merging CCDS IDs...')
    table = table.merge(ccds, left_index=True, right_index=True, how='left').rename(columns={'CCDS ID': 'ccds_id'})
    print('Merging exons...')
    table = table.merge(exons, left_index=True, right_index=True, how='left')
    print('Merging Pfam domains...')
    table = table.merge(pfam, left_index=True, right_index=True, how='left')

    return(table.reset_index().rename(columns={'index': 'transcript_stable_id'}))

def main(ensembl_biomart_transcripts,
         ensembl_transcript_info,
         ensembl_biomart_pfam,
         ensembl_biomart_refseq,
         ensembl_biomart_ccds,
         ensembl_biomart_transcripts_json):

    # Read input and set index column
    transcripts_df = pd.read_csv(ensembl_biomart_transcripts, sep='\t', index_col=0).sort_index()
    ccds_df = pd.read_csv(ensembl_biomart_ccds, sep='\t', index_col=0).sort_index()
    refseq_df = pd.read_csv(ensembl_biomart_refseq, sep='\t', index_col=0).sort_index()
    exons_df = pd.read_csv(ensembl_transcript_info, sep='\t')
    pfam_df = pd.read_csv(ensembl_biomart_pfam, sep='\t')

    #combined = combine_tables(transcripts, refseq, exons, pfam, ccds)
    exons = exons_per_transcript(exons_df).sort_index()
    pfam = pfam_domains_per_transcript(pfam_df).sort_index()

    # merge all tables
    merged = combine_tables(transcripts_df, refseq_df, exons, pfam, ccds_df)

    # print records as json
    merged.to_json(ensembl_biomart_transcripts_json,
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
    parser.add_argument("ensembl_biomart_transcripts_json",
                        help="tmp/ensembl_biomart_transcripts.json.gz")

    args = parser.parse_args()
    main(args.ensembl_biomart_transcripts,
         args.ensembl_transcript_info,
         args.ensembl_biomart_pfam,
         args.ensembl_biomart_refseq,
         args.ensembl_biomart_ccds,
         args.ensembl_biomart_transcripts_json)

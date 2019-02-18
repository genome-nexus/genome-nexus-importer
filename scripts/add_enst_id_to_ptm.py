import pandas as pd
import os
import argparse
import sys


def index_ccds_ids_by_uniprot(ccds_to_uniprot_df):
    # trim the trailing version information (-1, -2, etc.) since we don't have versioning in the PTM files
    ccds_to_uniprot_df['UniProtKB_short'] = ccds_to_uniprot_df['UniProtKB'].apply(lambda u: u.split('-')[0])
    # use the new trimmed UniProtKB column for grouping
    return ccds_to_uniprot_df.groupby('UniProtKB_short')["#ccds"].apply(set).to_dict()


def index_enst_by_ccds_ids(ccds_to_sequence_df):
    filtered = ccds_to_sequence_df[ccds_to_sequence_df['nucleotide_ID'].str.contains("ENST")]
    grouped = filtered.groupby("#ccds")
    return grouped["nucleotide_ID"].apply(set).to_dict()


def find_enst_by_ccds (ccds_id, ccds_to_enst_dict):
    try:
        return ccds_to_enst_dict[ccds_id]
    except KeyError:
        return 'NA'


def find_enst_by_uniprot(uniprot_id, ccds_to_uniprot_dict, ccds_to_enst_dict):
    try:
        flat_enst_id_list = []

        # add all into a flat list
        for enst_ids in map(lambda ccds_id: find_enst_by_ccds(ccds_id, ccds_to_enst_dict),
                            ccds_to_uniprot_dict[uniprot_id]):
            for enst_id in enst_ids:
                flat_enst_id_list.append(enst_id)

        # remove duplicates (if any)
        flat_enst_id_set = set(flat_enst_id_list)

        # remove NA
        flat_enst_id_set = set(filter(lambda enst_id: enst_id != 'NA', flat_enst_id_set))

        return flat_enst_id_set
    except KeyError:
        return set([])


# Convert comma (or colon, or semi-colon) separated PubMedID column into a proper list of strings
def parse_pubmed_ids(pubmed_ids):
    # clean up & split by possible delimiters
    pubmed_id_list = str(pubmed_ids).replace('doi:', '').replace(':', ';').replace(",", ";").split(';')
    # filter out empty/invalid string values
    pubmed_id_list = list(filter(lambda pubmed_id: len(pubmed_id.strip()) > 0, pubmed_id_list))

    # return a set of strings (or an empty set)
    if len(pubmed_id_list) > 0:
        return set(pubmed_id_list)
    else:
        return set([])


# for each ptm file map UniprotKB to ENST, and add a new EnsemblTranscript column
def add_enst_column_to_ptm_files(ccds_to_uniprot_dict, ccds_to_enst_dict, ptm_input_dir):
    frames = []

    # read and process all files under the directory
    for ptm_file in os.listdir(ptm_input_dir):
        ptm_df = pd.read_csv(f'{ptm_input_dir}/{ptm_file}',
                             sep='\t',
                             names=["uniprot_entry", "uniprot_accession", "position", "type", "pubmed_ids", "sequence"])
        # parse PubMed ids
        ptm_df['pubmed_ids'] = ptm_df.apply(lambda row: parse_pubmed_ids(row["pubmed_ids"]), axis=1)
        # add EnsemblTrascript info
        ptm_df['ensembl_transcript_ids'] = ptm_df.apply(
            lambda row: find_enst_by_uniprot(row["uniprot_accession"], ccds_to_uniprot_dict, ccds_to_enst_dict),
            axis=1)
        frames.append(ptm_df)

    # combine all frames and output a single PTM file
    pd.concat(frames).to_json(sys.stdout, orient='records', lines=True)


def main(ccds_to_uniprot, ccds_to_sequence, ptm_input_dir):
    # parse ccds mapping files
    ccds_to_uniprot_df = pd.read_csv(ccds_to_uniprot, sep='\t')
    ccds_to_sequence_df = pd.read_csv(ccds_to_sequence, sep='\t')

    # create dictionaries
    ccds_to_uniprot_dict = index_ccds_ids_by_uniprot(ccds_to_uniprot_df)
    ccds_to_enst_dict = index_enst_by_ccds_ids(ccds_to_sequence_df)

    # add ENST to PTM files
    add_enst_column_to_ptm_files(ccds_to_uniprot_dict, ccds_to_enst_dict, ptm_input_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ccds_to_uniprot",
                        help="common_input/CCDS2UniProtKB.current.txt")
    parser.add_argument("ccds_to_sequence",
                        help="common_input/CCDS2Sequence.current.txt")
    parser.add_argument("ptm_input_dir",
                        help="ptm/input")
    args = parser.parse_args()

    main(args.ccds_to_uniprot, args.ccds_to_sequence, args.ptm_input_dir)

import pandas as pd
import argparse
import sys

VARIANT_COUNT_POSTFIX = "_variant_count"
TUMOR_TYPE_COUNT_POSTFIX = "_tumortype_count"


def generate_count_list(mutation_row, tumor_types):
    return [
        {
            "tumor_type": tumor_type,
            "tumor_type_count": mutation_row[tumor_type + TUMOR_TYPE_COUNT_POSTFIX],
            "variant_count": mutation_row[tumor_type + VARIANT_COUNT_POSTFIX]
        } for tumor_type in tumor_types
    ]


def extract_tumor_types_from_col_names(mutations_df):
    # get columns (keys) ending with count postfixes
    count_keys = list(filter(lambda colName: colName.endswith(TUMOR_TYPE_COUNT_POSTFIX) or
                                             colName.endswith(VARIANT_COUNT_POSTFIX),
                             mutations_df.keys()))
    return set([
        count_key.replace(TUMOR_TYPE_COUNT_POSTFIX, "").replace(VARIANT_COUNT_POSTFIX, "") for count_key in count_keys
    ])


def process_data_frame(mutations_df, mutation_status):
    # add the list for tumor type counts
    mutations_df["counts_by_tumor_type"] = mutations_df.apply(
        lambda row: generate_count_list(row, extract_tumor_types_from_col_names(mutations_df)), axis=1
    )
    mutations_df["Mutation_Status"] = mutation_status
    # pick only the columns we need
    df = mutations_df[["Hugo_Symbol",
                       "Chromosome",
                       "Start_Position",
                       "End_Position",
                       "Reference_Allele",
                       "Alternate_Allele",
                       "Mutation_Status",
                       "classifier_pathogenic_final",
                       "penetrance",
                       "counts_by_tumor_type"]]
    # rename columns for JSON format
    df.columns = ["hugo_gene_symbol",
                  "chromosome",
                  "start_position",
                  "end_position",
                  "reference_allele",
                  "variant_allele",
                  "mutation_status",
                  "pathogenic",
                  "penetrance",
                  "counts_by_tumor_type"]
    return df


def merge_mutations(germline_mutations_df, biallelic_mutations_df, qc_pass_mutations_df):
    # generate ID for each data frame
    generate_id(germline_mutations_df)
    generate_id(biallelic_mutations_df)
    generate_id(qc_pass_mutations_df)

    # create a dictionary for germline mutations
    biallelic_mutations_index = biallelic_mutations_df.set_index("id").to_dict('index')
    qc_pass_mutations_index = qc_pass_mutations_df.set_index("id").to_dict('index')

    # add counts_by_tumor_type from biallelic and qc_pass into germline
    add_counts_by_tumor_type(germline_mutations_df, biallelic_mutations_index, "biallelic_counts_by_tumor_type")
    add_counts_by_tumor_type(germline_mutations_df, qc_pass_mutations_index, "qc_pass_counts_by_tumor_type")

    # drop the id columns from the original data frames
    germline_mutations_df.drop(["id"], axis=1)
    biallelic_mutations_df.drop(["id"], axis=1)
    qc_pass_mutations_df.drop(["id"], axis=1)


def add_counts_by_tumor_type(target_df, source_index, col_name):
    target_df[col_name] = target_df.apply(lambda row: get_counts_from_index(row["id"], source_index), axis=1)


def get_counts_from_index(row_id, source_index):
    try:
        return source_index[row_id]["counts_by_tumor_type"]
    except KeyError:
        return None


def generate_id(df):
    df["id"] = df["chromosome"].map(str) + "_" + \
               df["start_position"].map(str) + "_" + \
               df["end_position"].map(str) + "_" + \
               df["reference_allele"].map(str) + "_" + \
               df["variant_allele"].map(str)


def main(input_somatic, input_germline, input_biallelic, input_qc_pass):
    # parse mutations file
    somatic_mutations_df = pd.read_csv(input_somatic, sep='\t')
    germline_mutations_df = pd.read_csv(input_germline, sep='\t')
    biallelic_mutations_df = pd.read_csv(input_biallelic, sep='\t')
    qc_pass_mutations_df = pd.read_csv(input_qc_pass, sep='\t')
    # process original input
    somatic_mutations_df = process_data_frame(somatic_mutations_df, "somatic")
    germline_mutations_df = process_data_frame(germline_mutations_df, "germline")
    biallelic_mutations_df = process_data_frame(biallelic_mutations_df, "germline")
    qc_pass_mutations_df = process_data_frame(qc_pass_mutations_df, "germline")
    # merge counts from biallelic and qc_pass files into main germline mutations df
    merge_mutations(germline_mutations_df, biallelic_mutations_df, qc_pass_mutations_df)
    # convert processed data frames to JSON format
    somatic_mutations_df.to_json(sys.stdout, orient='records', lines=True)
    germline_mutations_df.to_json(sys.stdout, orient='records', lines=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_somatic",
                        help="insight/somatic_mutations_by_tumortype_merge.txt")
    parser.add_argument("input_germline",
                        help="insight/mutations_cnv_by_tumortype_merge.txt")
    parser.add_argument("input_biallelic",
                        help="insight/biallelic_by_tumortype_merge.txt")
    parser.add_argument("input_qc_pass",
                        help="insight/mutations_QCpass_by_tumortype_merge.txt")
    args = parser.parse_args()
    main(args.input_somatic, args.input_germline, args.input_biallelic, args.input_qc_pass)

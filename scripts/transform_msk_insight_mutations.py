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


def main(input_somatic, input_germline):
    # parse mutations file
    somatic_mutations_df = pd.read_csv(input_somatic, sep='\t')
    germline_mutations_df = pd.read_csv(input_germline, sep='\t')
    # process original input
    somatic_mutations_df = process_data_frame(somatic_mutations_df, "somatic")
    germline_mutations_df = process_data_frame(germline_mutations_df, "germline")
    # convert processed data frames to JSON format
    somatic_mutations_df.to_json(sys.stdout, orient='records', lines=True)
    germline_mutations_df.to_json(sys.stdout, orient='records', lines=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("inputSomatic",
                        help="insight/somatic_mutations_by_tumortype_merge_fixed.txt")
    parser.add_argument("inputGermline",
                        help="insight/mutations_by_tumortype_merge.txt")
    args = parser.parse_args()
    main(args.inputSomatic, args.inputGermline)

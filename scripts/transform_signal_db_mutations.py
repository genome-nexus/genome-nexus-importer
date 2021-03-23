import pandas as pd
import argparse
import sys

VARIANT_COUNT_POSTFIX = "_variant_count"
TUMOR_TYPE_COUNT_POSTFIX = "_tumortype_count"

GENERAL_POPULATION_F_PREFIX = "f_"
GENERAL_POPULATION_N_PREFIX = "n_"
SIG_PREFIX = "Sig."


def generate_count_list(mutation_row, tumor_types):
    return [
        {
            "tumor_type": tumor_type,
            "tumor_type_count": mutation_row[tumor_type + TUMOR_TYPE_COUNT_POSTFIX],
            "variant_count": mutation_row[tumor_type + VARIANT_COUNT_POSTFIX]
        } for tumor_type in tumor_types
    ]


def generate_stats_by_prefix(mutation_row, col_names, prefix):
    return {
        col_name.replace(prefix, "").lower(): mutation_row[col_name] for col_name in col_names
    }


def generate_general_population_stats(mutation_row, count_col_names, freq_col_names):
    # add the lists for count and frequency stats
    return {
        "counts": generate_stats_by_prefix(
            mutation_row,
            count_col_names,
            GENERAL_POPULATION_N_PREFIX
        ),
        "frequencies": generate_stats_by_prefix(
            mutation_row,
            freq_col_names,
            GENERAL_POPULATION_F_PREFIX
        )
    }


def generate_tumor_type_stats(mutation_row, sig_col_names):
    return {
        "tumor_type": mutation_row["Proposed_level"],
        "n_cancer_type_count": mutation_row["n_cancer_type_count"],
        "f_cancer_type_count": mutation_row["f_cancer_type_count"],
        "f_biallelic": mutation_row["f_biallelic"],
        "age_at_dx": mutation_row["age_at_dx"],
        "tmb": mutation_row["tmb"],
        "msi_score": mutation_row["msi_score"],
        "n_with_sig": mutation_row["n_with_sig"],
        "signatures": generate_stats_by_prefix(mutation_row, sig_col_names, SIG_PREFIX),
        "hrd_score": {
            "lst": mutation_row["lst"],
            "ntelomeric_ai": mutation_row["ntelomeric_ai"],
            "fraction_loh": mutation_row["fraction_loh"],
        },
        "n_germline_homozygous": mutation_row["n_germline_homozygous"]
    }


def extract_tumor_types_from_col_names(mutations_df):
    # get columns (keys) ending with count postfixes
    count_keys = list(filter(lambda colname: colname.endswith(TUMOR_TYPE_COUNT_POSTFIX) or
                                             colname.endswith(VARIANT_COUNT_POSTFIX),
                             mutations_df.keys()))
    return set([
        count_key.replace(TUMOR_TYPE_COUNT_POSTFIX, "").replace(VARIANT_COUNT_POSTFIX, "") for count_key in count_keys
    ])


def extract_col_names_by_prefix(mutations_df, stats_prefix):
    # get columns (keys) starting with population stats prefixes
    return list(filter(lambda colname: colname.startswith(stats_prefix), mutations_df.keys()))


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


def process_all_variant_freq_df(mutations_df, mutation_status):
    mutations_df["Mutation_Status"] = mutation_status
    mutations_df["general_population_stats"] = mutations_df.apply(
        lambda row: generate_general_population_stats(
            row,
            extract_col_names_by_prefix(mutations_df, GENERAL_POPULATION_N_PREFIX),
            extract_col_names_by_prefix(mutations_df, GENERAL_POPULATION_F_PREFIX),
        ), axis=1
    )
    # pick only the columns we need
    df = mutations_df[["Hugo_Symbol",
                       "Chromosome",
                       "Start_Position",
                       "End_Position",
                       "Reference_Allele",
                       "Alternate_Allele",
                       "Mutation_Status",
                       "general_population_stats",
                       "n_germline_homozygous"]]
    # rename columns for JSON format
    df.columns = ["hugo_gene_symbol",
                  "chromosome",
                  "start_position",
                  "end_position",
                  "reference_allele",
                  "variant_allele",
                  "mutation_status",
                  "general_population_stats",
                  "n_germline_homozygous"]
    return df


def process_variants_by_cancertype_summary_df(mutations_df, mutation_status):
    mutations_df["Mutation_Status"] = mutation_status
    mutations_df["stats_by_tumor_type"] = mutations_df.apply(
        lambda row: generate_tumor_type_stats(
            row,
            extract_col_names_by_prefix(mutations_df, SIG_PREFIX)
        ), axis=1
    )
    # pick only the columns we need
    df = mutations_df[["Hugo_Symbol",
                       "Chromosome",
                       "Start_Position",
                       "End_Position",
                       "Reference_Allele",
                       "Alternate_Allele",
                       "Mutation_Status",
                       "stats_by_tumor_type"]]
    # rename columns for JSON format
    df.columns = ["hugo_gene_symbol",
                  "chromosome",
                  "start_position",
                  "end_position",
                  "reference_allele",
                  "variant_allele",
                  "mutation_status",
                  "stats_by_tumor_type"]
    return df


def create_variants_by_cancertype_summary_index(mutations_df):
    # all other column names except stats_by_tumor_type
    col_names = list(filter(lambda colname: colname != "stats_by_tumor_type", mutations_df.keys()))
    # keep the first for every column except stats_by_tumor_type, they are the same
    aggregator = {col_name: 'first' for col_name in col_names}
    # aggregate stats_by_tumor_type as a list
    aggregator['stats_by_tumor_type'] = lambda s: list(s)
    return mutations_df.groupby('id').agg(aggregator).set_index("id").to_dict('index')


def process_msk_expert_review_df(mutations_df, mutation_status):
    mutations_df["Mutation_Status"] = mutation_status
    mutations_df["MSK_Expert_Review"] = True
    # pick only the columns we need
    df = mutations_df[["Hugo_Symbol",
                       "Chromosome",
                       "Start_Position",
                       "End_Position",
                       "Reference_Allele",
                       "Alternate_Allele",
                       "Mutation_Status",
                       "MSK_Expert_Review"]]
    # rename columns for JSON format
    df.columns = ["hugo_gene_symbol",
                  "chromosome",
                  "start_position",
                  "end_position",
                  "reference_allele",
                  "variant_allele",
                  "mutation_status",
                  "msk_expert_review"]
    return df


def merge_mutations(germline_mutations_df,
                    biallelic_mutations_df,
                    qc_pass_mutations_df,
                    all_variants_freq_df,
                    msk_expert_review_df,
                    variants_by_cancertype_summary_df):
    # generate ID for each data frame
    generate_id(germline_mutations_df)
    generate_id(biallelic_mutations_df)
    generate_id(qc_pass_mutations_df)
    generate_id(all_variants_freq_df)
    generate_id(msk_expert_review_df)
    generate_id(variants_by_cancertype_summary_df)

    # create a dictionary for germline mutations
    biallelic_mutations_index = biallelic_mutations_df.set_index("id").to_dict('index')
    qc_pass_mutations_index = qc_pass_mutations_df.set_index("id").to_dict('index')
    all_variants_freq_index = all_variants_freq_df.set_index("id").to_dict('index')
    msk_expert_review_index = msk_expert_review_df.set_index("id").drop_duplicates().to_dict('index')
    variants_by_cancertype_summary_index = create_variants_by_cancertype_summary_index(
        variants_by_cancertype_summary_df)

    # add data from other sources into the main germline data frame
    add_column_data(germline_mutations_df,
                    biallelic_mutations_index,
                    "biallelic_counts_by_tumor_type",
                    "counts_by_tumor_type")
    add_column_data(germline_mutations_df,
                    qc_pass_mutations_index,
                    "qc_pass_counts_by_tumor_type",
                    "counts_by_tumor_type")
    add_column_data(germline_mutations_df,
                    all_variants_freq_index,
                    "general_population_stats",
                    "general_population_stats")
    add_column_data(germline_mutations_df,
                    all_variants_freq_index,
                    "n_germline_homozygous",
                    "n_germline_homozygous")
    add_column_data(germline_mutations_df,
                    msk_expert_review_index,
                    "msk_expert_review",
                    "msk_expert_review",
                    default_value=False)
    add_column_data(germline_mutations_df,
                    variants_by_cancertype_summary_index,
                    "stats_by_tumor_type",
                    "stats_by_tumor_type")

    # drop the id columns from the original data frames
    germline_mutations_df.drop(columns=["id"], inplace=True)
    biallelic_mutations_df.drop(columns=["id"], inplace=True)
    qc_pass_mutations_df.drop(columns=["id"], inplace=True)
    all_variants_freq_df.drop(columns=["id"], inplace=True)
    msk_expert_review_df.drop(columns=["id"], inplace=True)
    variants_by_cancertype_summary_df.drop(columns=["id"], inplace=True)


def add_column_data(target_df, source_index, target_col_name, source_col_name, default_value=None):
    target_df[target_col_name] = target_df.apply(
        lambda row: pick_column_from_index(row["id"], source_index, source_col_name, default_value), axis=1
    )


def pick_column_from_index(row_id, source_index, source_col_name, default_value=None):
    try:
        return source_index[row_id][source_col_name]
    except KeyError:
        return default_value


def generate_id(df):
    df["id"] = df["chromosome"].map(str) + "_" + \
               df["start_position"].map(str) + "_" + \
               df["end_position"].map(str) + "_" + \
               df["reference_allele"].map(str) + "_" + \
               df["variant_allele"].map(str)


# workaround for mixed NAs and integers
# see https://stackoverflow.com/questions/39666308/pd-read-csv-by-default-treats-integers-like-floats
def fix_na_values(df):
    df[["Start_Position", "End_Position"]] = df[["Start_Position", "End_Position"]].fillna(-1).astype(int)


def parse_file(input, sep):
    # make sure that chromosome is always parsed as string
    df = pd.read_csv(input, sep=sep, dtype={"Chromosome": object})
    fix_na_values(df)
    return df


def main(input_somatic,
         input_germline,
         input_biallelic,
         input_qc_pass,
         input_all_variants_freq,
         input_msk_expert_review,
         input_variants_by_cancertype_summary):
    # parse mutation files
    somatic_mutations_df = parse_file(input_somatic, sep='\t')
    germline_mutations_df = parse_file(input_germline, sep='\t')
    biallelic_mutations_df = parse_file(input_biallelic, sep='\t')
    qc_pass_mutations_df = parse_file(input_qc_pass, sep='\t')
    all_variants_freq_df = parse_file(input_all_variants_freq, sep='\t')
    msk_expert_review_df = parse_file(input_msk_expert_review, sep='\t')
    variants_by_cancertype_summary_df = parse_file(input_variants_by_cancertype_summary, sep='\t')
    # process original input
    somatic_mutations_df = process_data_frame(somatic_mutations_df, "somatic")
    germline_mutations_df = process_data_frame(germline_mutations_df, "germline")
    biallelic_mutations_df = process_data_frame(biallelic_mutations_df, "germline")
    qc_pass_mutations_df = process_data_frame(qc_pass_mutations_df, "germline")
    all_variants_freq_df = process_all_variant_freq_df(all_variants_freq_df, "germline")
    msk_expert_review_df = process_msk_expert_review_df(msk_expert_review_df, "germline")
    variants_by_cancertype_summary_df = process_variants_by_cancertype_summary_df(variants_by_cancertype_summary_df,
                                                                                  "germline")
    # merge everything into the main germline mutations data frame
    merge_mutations(germline_mutations_df,
                    biallelic_mutations_df,
                    qc_pass_mutations_df,
                    all_variants_freq_df,
                    msk_expert_review_df,
                    variants_by_cancertype_summary_df)
    # convert processed data frames to JSON format
    somatic_mutations_df.to_json(sys.stdout, orient='records', lines=True)
    germline_mutations_df.to_json(sys.stdout, orient='records', lines=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_somatic",
                        help="signal/somatic_mutations_by_tumortype_merge.txt")
    parser.add_argument("input_germline",
                        help="signal/mutations_cnv_by_tumortype_merge.txt")
    parser.add_argument("input_biallelic",
                        help="signal/biallelic_by_tumortype_merge.txt")
    parser.add_argument("input_qc_pass",
                        help="signal/mutations_QCpass_by_tumortype_merge.txt")
    parser.add_argument("input_all_variants_freq",
                        help="signal/signaldb_all_variants_frequencies.txt")
    parser.add_argument("input_msk_expert_review",
                        help="signal/signaldb_msk_expert_review_variants.txt")
    parser.add_argument("input_variants_by_cancertype_summary",
                        help="signal/signaldb_variants_by_cancertype_summary_statistics.txt")
    args = parser.parse_args()
    main(args.input_somatic,
         args.input_germline,
         args.input_biallelic,
         args.input_qc_pass,
         args.input_all_variants_freq,
         args.input_msk_expert_review,
         args.input_variants_by_cancertype_summary)

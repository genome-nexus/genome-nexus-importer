import argparse
import requests
import pandas as pd
import os
import re

# check genome build (only support grch37 and grch38)
def validate_genome_build(genome_build_string):
    genome_build_values = genome_build_string.split('_')
    return (len(genome_build_values) == 2 and (genome_build_values[0].lower() == 'grch37' or genome_build_values[0].lower() == 'grch38'))

def get_annotation_sources_version_df(annotation_sources_version_file, genome_build):
    annotation_sources_version_df = pd.read_csv(annotation_sources_version_file, sep='\t')
    return annotation_sources_version_df[annotation_sources_version_df['genome_build'] == genome_build].drop('genome_build', axis=1)

def main(annotation_sources_version_file, annotation_version_output_file):
    # read genome build from environment variables
    try:
        genome_build_string = os.environ["VERSION"]
        if (validate_genome_build(genome_build_string) == True):
            genome_build = genome_build_string.split('_')[0]
            annotation_sources_version_df = get_annotation_sources_version_df(annotation_sources_version_file, genome_build)
            annotation_sources_version_df.to_csv(annotation_version_output_file, sep="\t", header=True, index=False)
        else:
            raise Exception("Invalid 'VERSION', only grch37 and grch38 are supported")
    except KeyError:
        print("No 'VERSION' variable found, please export a valid 'VERSION'. (e.g. export VERSION=grch38_ensembl95)")
    except Exception as e: 
        print(str(e))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("annotation_sources_version_file", help="../data/common_input/version_info.txt")
    parser.add_argument("annotation_version_output_file", help="../data/grch37_ensembl92/export/annotation_version.txt or ../data/grch38_ensembl95/export/annotation_version.txt or ../data/grch38_ensembl92/export/annotation_version.txt")
    args = parser.parse_args()
    main(args.annotation_sources_version_file, args.annotation_version_output_file)
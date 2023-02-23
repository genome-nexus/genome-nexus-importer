import argparse
import requests
import pandas as pd

def main(oncokb_isoform_overrides_file):
    # check if genome is grch37 (hg19), drop grch38 columns if reference_genome is grch37
    if 'grch37' in oncokb_isoform_overrides_file:
        drop_columns_name = 'grch38'
    else:
        drop_columns_name = 'grch37'

    url ='https://www.oncokb.org/api/v1/utils/allCuratedGenes?includeEvidence=false'
    oncokb_df = pd.json_normalize(requests.get(url).json())

    oncokb_df.drop(list(oncokb_df.filter(regex =drop_columns_name)), axis = 1, inplace = True)
    column_name_mapping = {"background": "background",
                            "entrezGeneId": "entrez_gene_id",
                            "grch37Isoform": "enst_id",
                            "grch37RefSeq": "ref_seq",
                            "grch38Isoform": "enst_id",
                            "grch38RefSeq": "ref_seq",
                            "highestResistancLevel": "highest_resistanc_level",
                            "highestResistanceLevel": "highest_resistance_level",
                            "highestSensitiveLevel": "highest_sensitive_level",
                            "hugoSymbol": "hugo_symbol",
                            "oncogene": "oncogene",
                            "summary": "summary",
                            "tsg": "tsg"}
    oncokb_df.rename(columns=column_name_mapping, inplace=True)
    oncokb_df.to_csv(oncokb_isoform_overrides_file, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("oncokb_isoform_overrides_file",
                        help="isoform_overrides_oncokb_grch37.txt or isoform_overrides_oncokb_grch38.txt")
    args = parser.parse_args()

    main(args.oncokb_isoform_overrides_file)

import argparse
import requests
import pandas as pd

def main(reference_genome):
    # check if genome is grch37 (hg19), drop grch38 columns if reference_genome is grch37
    drop_columns_name = 'grch38' if reference_genome == 'grch37' else 'grch37'

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
    oncokb_df.to_csv('common_input/isoform_overrides_oncokb_' + reference_genome + '.txt', sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_genome",
                        help="grch37 or grch38")
    args = parser.parse_args()

    main(args.reference_genome)

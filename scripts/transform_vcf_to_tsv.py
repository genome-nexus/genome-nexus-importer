import argparse
import pandas as pd 
from cyvcf2 import VCF

def __main__():
   parser = argparse.ArgumentParser()
   parser.add_argument('input_vcf', help='Input VCF file, for example:../data/clinvar/input/clinvar_grch37.vcf')
   parser.add_argument('out_tsv', help='Output file, for example:../data/clinvar/tsv/clinvar_grch37.tsv')
   parser.add_argument('--fields', help='Select fields and separate by comma. Avalivable fields: chromosome, start_position, end_position, reference_allele, alternate_allele, clinvar_id, quality, filter, af_esp, af_exac, af_tgp, alleleid, clndisdb, clndisdbincl, clndn, clndnincl, clnhgvs, clnrevstat, clnsig, clnsigconf, clnsigincl, clnvc, clnvcso, clnvi, dbvarid, geneinfo, mc, origin, rs, ssr', required=False, default=None)
   args = parser.parse_args()
   vcf2tsv(args.input_vcf, args.out_tsv, args.fields)

def get_genomic_location(variant):
   genomic_location = []
   genomic_location.append(str(variant.CHROM))
   if len(variant.ALT) == 1 and len(variant.ALT[0]) > 0 and len(variant.REF) > 0:
      if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:
         # SNP
         # 1	877523	C	G -> 1,877523,877523,C,G   
         # start position
         genomic_location.append(str(variant.POS))
         # end position
         genomic_location.append(str(variant.POS))
         # ref
         genomic_location.append(variant.REF)
         # var
         genomic_location.append(variant.ALT[0])
      elif len(variant.REF) > 1 and len(variant.ALT[0]) == 1 and variant.REF[0:1] == variant.ALT[0][0:1]:
         # DEL
         # 1 3342785 AAACGGT	A -> 1,3342786,3342791,AACGGT,-
         # start position
         genomic_location.append(str(variant.POS + 1))
         # end position
         genomic_location.append(str(variant.end))
         # ref
         genomic_location.append(variant.REF[1:])
         # var
         genomic_location.append('-')
      elif len(variant.REF) == 1 and len(variant.ALT[0]) > 1 and variant.REF[0:1] == variant.ALT[0][0:1]:
         # INS 
         # 1	6529194	C	CTCT -> 1,6529194,6529195,-,TCT
         # start position
         genomic_location.append(str(variant.POS))
         # end position 
         genomic_location.append(str(variant.POS + 1))
         # ref
         genomic_location.append('-')
         # var
         genomic_location.append(variant.ALT[0][1:])
      else:
         # DELINS
         # 1	12062157	AG	CT -> 1,12062157,12062158,AG,CT
         # 1	173797451	T	CC -> 1,173797451,173797451,T,CC
         # start position
         genomic_location.append(str(variant.POS))
         # end position
         genomic_location.append(str(variant.end))
         # ref
         genomic_location.append(variant.REF)
         # var
         genomic_location.append(variant.ALT[0])
      # skip if it's multiple variant allele(ALT:G,T), or no variant allele(ALT:.), see examples in VCF format specification: https://samtools.github.io/hts-specs/VCFv4.1.pdf
   return genomic_location

def parse_INFO_column(vcf, column_types):
   info_columns_header = []
   for h in vcf.header_iter():
      header_element = h.info()
      if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
         if header_element['HeaderType'] == 'INFO':
            column_types[header_element['ID']] = header_element['Type']
            info_columns_header.append(header_element['ID'])
   return sorted(info_columns_header)

def vcf2tsv(input_vcf, out_tsv, fields):
   vcf = VCF(input_vcf, gts012 = True)
   fixed_columns_header = ['chromosome','start_position','end_position','reference_allele','alternate_allele','clinvar_id','quality','filter']
   column_types = {}
   
   info_columns_header = parse_INFO_column(vcf, column_types)
   header = fixed_columns_header + list(map(str.lower, info_columns_header))
   # select what fields will be included
   # will keep all fields if not specified
   fields_list = fields.lower().split(",") if fields != None else map(str.lower, header)
   
   tsv_data = []
   for variant in vcf:
      # add genomic location fields to row element
      tsv_fields = get_genomic_location(variant)
      # only keep data with full genomic location (e.g. skip when ref is missing)
      if (len(tsv_fields)) > 1:
         # Variant ID
         tsv_fields.append(str(variant.ID))
         
         # Quality
         if variant.QUAL is None:
            tsv_fields.append(None)
         else:
            tsv_fields.append(str("{0:.2f}".format(variant.QUAL)))

         # Filter
         if variant.FILTER is None:
            tsv_fields.append(None)
         else:
            tsv_fields.append(str(variant.FILTER))
         
         # INFO
         variant_info = variant.INFO
         for info_field in sorted(info_columns_header):
            if type(variant_info.get(info_field)) is list or type(variant_info.get(info_field)) is tuple:
               tsv_fields.append(",".join(str(n) for n in variant_info.get(info_field)))
            else:
               if variant_info.get(info_field) is None:
                  tsv_fields.append(None)
               else:
                  if column_types[info_field] == 'Float':
                     tsv_fields.append(str("{0:.5f}".format(variant_info.get(info_field))))
                  else:
                     tsv_fields.append(str(variant_info.get(info_field)))
         tsv_data.append(tsv_fields)
      else:
         continue
   # write tsv data into DataFrame, then drop unused columns
   df = pd.DataFrame(tsv_data, columns=header)
   for column in df.columns:
      if column.lower() in fields_list:
         continue
      else:
         df = df.drop(column, axis=1)
   # set column type
   df.rename(columns=lambda x: x + ".string()" if x == "chromosome" else x + ".auto()", inplace=True)
   df.to_csv(out_tsv, sep="\t", index=False)

if __name__=="__main__": __main__()
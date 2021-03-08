import argparse
from cyvcf2 import VCF
import subprocess
import gzip

def __main__():
   parser = argparse.ArgumentParser()
   parser.add_argument('clinvar_vcf', help='../data/clinvar/input/clinvar_grch37.vcf')
   parser.add_argument('out_tsv', help='../data/clinvar/tsv/clinvar_grch37.tsv')
   # parser.add_argument('compress_output', help='set compress_output will compress output, compress by default', default=True)
   args = parser.parse_args()
   vcf2tsv(args.clinvar_vcf, args.out_tsv)
   # vcf2tsv('/Users/lix2/Documents/GitHub/genome-nexus-importer/data/clinvar/input/clinvar_grch37.vcf.gz', '/Users/lix2/Documents/GitHub/genome-nexus-importer/data/clinvar/tsv/clinvar_grch37.tsv')

def get_genomic_location(variant):
   genomic_location = []
   genomic_location.append(variant.CHROM)
   if len(variant.ALT) == 1:
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
      elif len(variant.REF) > 1 and len(variant.ALT[0]) == 1:
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
      elif len(variant.REF) == 1 and len(variant.ALT[0]) > 1:
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
         # start position
         genomic_location.append(str(variant.POS))
         # end position
         genomic_location.append(str(variant.end))
         # ref
         genomic_location.append(variant.REF)
         # var
         genomic_location.append(variant.ALT[0])
   else:
      # multiple variant allele, or no variant allele
      # start position
      genomic_location.append(str('-1'))
      # end position
      genomic_location.append(str('-1'))
      # ref
      genomic_location.append(variant.REF)
      # var
      genomic_location.append(",".join(str(n) for n in variant.ALT))
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

def vcf2tsv(clinvar_vcf, out_tsv):
   vcf = VCF(clinvar_vcf, gts012 = True)
   out = open(out_tsv,'w')
   fixed_columns_header = ['Chromosome','Start_Position','End_Position','Reference_Allele','Alternate_Allele','ClinVar_ID','Quality','Filter']
   column_types = {}
   
   info_columns_header = parse_INFO_column(vcf, column_types)
   header = fixed_columns_header + info_columns_header
   out.write(str('\t'.join(header)) + '\n')
   
   for variant in vcf:
      # add genomic location fields to row element
      tsv_fields = get_genomic_location(variant)

      # ClinVar ID
      tsv_fields.append(str(variant.ID))
      
      # Quality
      if variant.QUAL is None:
         tsv_fields.append('')
      else:
         tsv_fields.append(str("{0:.2f}".format(variant.QUAL)))

      # Filter
      if variant.FILTER is None:
         tsv_fields.append('')
      else:
         tsv_fields.append(str(variant.FILTER))
      
      # INFO
      variant_info = variant.INFO
      for info_field in sorted(info_columns_header):
         if type(variant_info.get(info_field)) is list or type(variant_info.get(info_field)) is tuple:
            tsv_fields.append(",".join(str(n) for n in variant_info.get(info_field)))
         else:
            if variant_info.get(info_field) is None:
               tsv_fields.append('')
            else:
               if column_types[info_field] == 'Float':
                  tsv_fields.append(str("{0:.5f}".format(variant_info.get(info_field))))
               else:
                  tsv_fields.append(str(variant_info.get(info_field)))

      out.write('\t'.join(tsv_fields) + '\n')
   out.close()
   # print(compress_output)
   # if compress_output == True:
   #    subprocess.run('gzip -f ' + str(out_tsv), shell=True)
if __name__=="__main__": __main__()


   
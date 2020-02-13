
import pandas as pd
import sys

gwas_summary_file_hg37 = sys.argv[1]
gwas_summary_file_liftedOver_to_hg38 = sys.argv[2]
dbSNP150_file = sys.argv[3]
gwas_summary_folder = sys.argv[4]
gwas_summary_filename_mod=sys.argv[5]

# example to pass an argument variable to a string 
# gwas_file = pd.read_csv(filepath_or_buffer='/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/%s/tmp3.txt' %infile,sep='\t')

gwas_file = pd.read_csv(filepath_or_buffer=gwas_summary_file_hg37, sep='\t')
crossmap_file = pd.read_csv(filepath_or_buffer=gwas_summary_file_liftedOver_to_hg38, sep='\t')
gwas_file_liftedOver_hg38 = pd.merge(left=gwas_file, right=crossmap_file, left_on="variant_", right_on="variant_hg37", how='left')	
dbSNP150_hg38=pd.read_csv(filepath_or_buffer=dbSNP150_file, header=None,sep='\t')
gwas_file_liftedOver_hg38_dbSNP150 = pd.merge(left=gwas_file_liftedOver_hg38, right=dbSNP150_hg38, left_on="variant_hg38", right_on=0, how='left')	
gwas_file_liftedOver_hg38_dbSNP150.to_csv(path_or_buf=f"{gwas_summary_folder}/{gwas_summary_filename_mod}_hg38.tsv", sep='\t', index=False)


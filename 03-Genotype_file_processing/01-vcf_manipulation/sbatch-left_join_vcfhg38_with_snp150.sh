#!/bin/bash

set -e

#cat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.vcf | grep -v ^#  | awk '{print $1"_"$2"_"$4"_"$5 "\t" $0}'> GTEx.hg38.noHash.vcf

# hg38_avsnp150_SNPonly.txt - downloaded dbSNP150 SNP only table from annovar

cat hg38_avsnp150_SNPonly.txt | awk '{print "chr"$1"_"$2"_"$4"_"$5 "\t" $6}' > hg38_avsnp150_SNPonlyV2.txt

cat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.vcf | grep '#' >  GTEx.hg38.rsIDs.Hash.vcf

awk 'NR==FNR{a[$1]=$2;next}NR!=FNR{c=$1; if(c in a){print $0 "\t" a[c]}}' hg38_avsnp150_SNPonlyV2.txt GTEx.hg38.noHash.vcf | awk '{$4=$NF}1' | awk -v OFS="\t" 'NF{NF-=1}1' | cut -f2- > GTEx.hg38.rsIDs.noHash.vcf

cat GTEx.hg38.rsIDs.noHash.vcf >> GTEx.hg38.rsIDs.Hash.vcf

rm -f GTEx.hg38.rsIDs.noHash.vcf

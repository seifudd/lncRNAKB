#!/bin/bash

set -e

awk '/#/{print;next}NR==FNR{a[$2]=$1;next}NR!=FNR{c=$3; if (c in a){print $0,a[c]}}' gtexID.snpID.txt GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.vcf | awk '{$3=$NF}1' | awk 'NF{NF-=1}1'  > GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.rsIDs.vcf
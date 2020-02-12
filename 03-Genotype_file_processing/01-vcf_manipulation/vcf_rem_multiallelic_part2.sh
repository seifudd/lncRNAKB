#!/bin/bash

set -e

awk '/#/{print; next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38n.vcf > GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.vcf
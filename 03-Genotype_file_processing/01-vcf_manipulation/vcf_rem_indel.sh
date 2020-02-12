#!/bin/bash

set -e

module load vcftools
vcftools --gzvcf GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz  --remove-indels --recode --recode-INFO-all --out GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz
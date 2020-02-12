#!/bin/bash

set -e

module load vcftools
vcf-stats GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
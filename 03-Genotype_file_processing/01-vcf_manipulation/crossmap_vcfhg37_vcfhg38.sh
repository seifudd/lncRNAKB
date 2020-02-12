#!/bin/bash

set -e

# hg38.nochr.fa - from UCSC Genome Browser

module load crossmap
CrossMap.py vcf GRCh37_to_GRCh38.chain.gz GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.vcf hg38.nochr.fa GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38n.vcf



#####################


#cat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.vcf | head -1000 | grep ^# > GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.mod.vcf

#cat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.vcf | grep -v ^# | head -2000 | awk '{print "chr"$0}' >> GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.mod.vcf

#crossmap vcf GRCh37_to_GRCh38.chain.gz GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.mod.vcf hg38.fa > GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.withgenome38.vcf





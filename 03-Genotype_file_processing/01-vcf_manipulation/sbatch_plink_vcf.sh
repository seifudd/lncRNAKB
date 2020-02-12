#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=plinkvcf \
		--partition=norm \
		--time=24:00:00 \
		--mem=128 \
		--gres=lscratch:150 \
		--cpus-per-task=4 \
		--error=plink_vcf.stderr \
		--output=plink_vcf.stdout \
		/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/plink_vcf.sh
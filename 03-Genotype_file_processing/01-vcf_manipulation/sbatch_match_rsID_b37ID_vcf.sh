#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=match_rsID_b37ID_vcf \
		--partition=norm \
		--time=72:00:00 \
		--mem=128 \
		--gres=lscratch:150 \
		--cpus-per-task=4 \
		--error=match_rsID_b37ID_vcf.stderr \
		--output=match_rsID_b37ID_vcf.stdout \
		/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/match_rsID_b37ID_vcf.sh
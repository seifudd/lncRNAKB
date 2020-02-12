#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=crossmap_vcf \
			--partition=norm \
			--time=96:00:00 \
			--mem=128g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=crossmap_vcf.stderr \
			--output=crossmap_vcf.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/crossmap_vcfhg37_vcfhg38.sh
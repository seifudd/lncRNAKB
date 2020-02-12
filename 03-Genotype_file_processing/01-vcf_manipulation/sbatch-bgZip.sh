#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=bgZip \
			--partition=norm \
			--time=48:00:00 \
			--mem=128g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=bgZip.stderr \
			--output=bgZip.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_manipulation_ks/bgZip.sh
#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=remMissingMulti \
			--partition=norm \
			--time=96:00:00 \
			--mem=128g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=remMissingMulti.stderr \
			--output=remMissingMulti.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_rem_missing_multiallelic.sh
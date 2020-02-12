#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=remMultiPart2 \
			--partition=norm \
			--time=96:00:00 \
			--mem=128g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=remMultiPart2.stderr \
			--output=remMultiPart2.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_rem_multiallelic_part2.sh
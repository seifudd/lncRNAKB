#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=vcfstats \
			--partition=norm \
			--time=24:00:00 \
			--mem=64g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=vcfstats.stderr \
			--output=vcfstats.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_stats.sh
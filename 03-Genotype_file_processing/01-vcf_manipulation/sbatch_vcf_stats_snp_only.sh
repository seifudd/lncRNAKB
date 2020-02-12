#!/bin/bash

cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=vcfstats_SNPOnly \
			--partition=norm \
			--time=3-00:00:00 \
			--mem=156g \
			--gres=lscratch:200 \
			--cpus-per-task=4 \
			--error=vcfstats_snpOnly.stderr \
			--output=vcfstats_snpOnly.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_stats_snp_only.sh
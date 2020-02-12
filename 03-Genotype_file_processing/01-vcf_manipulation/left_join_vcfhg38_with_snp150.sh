#!/bin/bash



cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1

sbatch --job-name=leftjoin \
		--cpus-per-task=16 \
		--time=4:00:00 \
		--mem=128g \
		--error=leftjoinvcfhg38.stderr \
		--output=leftjoinvcfhg38.stdout \
		/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_manipulation_ks/sbatch-left_join_vcfhg38_with_snp150.sh
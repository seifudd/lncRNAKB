#!/bin/bash

function do_vcf_subset (){
	
	tissue=$1
	tissue_subset_file=$2
	vcfDir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1"
	outputDir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks"
	mkdir $outputDir/$tissue
	echo $tissue
	sbatch --job-name=vcf_${tissue}_subset \
			--partition=norm \
			--time=96:00:00 \
			--mem=128g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=$outputDir/vcf_subset.$tissue.stderr \
			--output=$outputDir/vcf_subset.$tissue.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks/vcf_subsetks.sh \
			$tissue \
			$tissue_subset_file \
			$vcfDir \
			$outputDir    
}

# do_vcf_subset Adipose_Tissue Adipose_Tissue.363.subset.txt
do_vcf_subset Salivary_Gland Salivary_Gland.63.subset.txt

# while read tissue tissuebysubsetfile; do
#	do_vcf_subset ${tissue} ${tissuebysubsetfile}
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/subset_files/00-by_tissue_subset_files_list_redo.txt"


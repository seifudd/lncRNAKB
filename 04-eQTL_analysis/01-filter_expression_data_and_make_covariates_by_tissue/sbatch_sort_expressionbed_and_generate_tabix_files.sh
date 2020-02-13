#!/bin/bash

function do_sort_expressionbed_and_generate_tabix_for_fasteqtl (){
	tissue=$1
	expressiondatafile=$2
	inputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS"
	outputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS"
	echo $tissue
	sbatch --job-name=${tissue}.sort.bed.tabix \
			--partition=norm \
			--time=5:00:00 \
			--mem=64g \
			--gres=lscratch:100 \
			--cpus-per-task=4 \
			--error=$outputdir/$tissue.sort.bed.tabix.stderr \
			--output=$outputdir/$tissue.sort.bed.tabix.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/sort_expressionbed_and_generate_tabix_files.sh \
			$expressiondatafile \
			$tissue \
			$inputdir \
			$outputdir  
}

# do_sort_expressionbed_and_generate_tabix_for_fasteqtl Bladder Bladder.9samples.15597pcs.13098lincs.gene.CPM.GTExId.txt

while read tissue expressiondatafile; do
 	do_sort_expressionbed_and_generate_tabix_for_fasteqtl ${tissue} ${expressiondatafile}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/list_of_filtered_cpm_files_by_tissue.txt"


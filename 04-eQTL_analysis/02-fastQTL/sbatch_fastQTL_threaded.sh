#!/bin/bash

function do_fasteqtl_threaded (){
	tissue=$1
	expressiondatafile=$2
	chunks=100
	threads=25
	inputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS"
	vcfinputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_FS"
	outputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS"
	window=500000
	permutations=0
	echo "TISSUE=$tissue"
	mkdir -p $outputdir/$tissue
	sbatch --job-name=${tissue}.fastqtl.threaded \
			--partition=norm \
			--time=5:00:00 \
			--mem=100g \
			--gres=lscratch:100 \
			--cpus-per-task=$threads \
			--error=$outputdir/$tissue.chr$chr.$window.$permutations.fastqtl.threaded.stderr \
			--output=$outputdir/$tissue.chr$chr.$window.$permutations.fastqtl.threaded.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded/fastQTL_threaded.sh \
			$inputdir \
			$vcfinputdir \
			$expressiondatafile \
			$tissue \
			$window \
			$permutations \
			$chunks \
			$threads \
			$outputdir 
}

while read tissue expressiondatafile; do
		do_fasteqtl_threaded ${tissue} ${expressiondatafile}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/list_of_filtered_cpm_files_by_tissue.txt"



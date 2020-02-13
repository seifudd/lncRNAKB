#!/bin/bash

function do_fasteqtl_threaded_permutations (){
	tissue=$1
	expressiondatafile=$2
	chunks=30
	threads=100
	inputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6"
	vcfinputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks"
	outputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded_permutations"
	window=500000
	permutations1=1000
	permutations2=10000
	echo "TISSUE=$tissue"
	mkdir $outputdir/$tissue
	sbatch --job-name=${tissue}.perm.fastqtl.threaded \
			--partition=largemem \
			--time=12:00:00 \
			--mem=1000g \
			--gres=lscratch:100 \
			--cpus-per-task=$threads \
			--error=$outputdir/$tissue.$window.$permutations1.$permutations2.fastqtl.threaded.stderr \
			--output=$outputdir/$tissue.$window.$permutations1.$permutations2.fastqtl.threaded.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded_permutations/fastQTL_threaded_permutations.sh \
			$inputdir \
			$vcfinputdir \
			$expressiondatafile \
			$tissue \
			$window \
			$permutations1 \
			$permutations2 \
			$chunks \
			$threads 
}

while read tissue expressiondatafile; do
		do_fasteqtl_threaded_permutations ${tissue} ${expressiondatafile}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6/cpm_filteredbycount_filteredbytpm_file_list_redo_perm.txt"

# Prostate	Prostate.101samples.15931pcs.14442lincs.gene.CPM.GTExId.txt
# Pancreas	Pancreas.167samples.14210pcs.9359lincs.gene.CPM.GTExId.txt

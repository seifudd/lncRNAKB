#!/bin/bash
#$ -cwd

function do_calculate_FDR () {
	scriptdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS"
	tissue=$1
	covariatefilename=$2
	numcpus=8
	tissue_with_samplesize=`echo $covariatefilename | sed 's/.covariate.sorted.txt//g'`
	echo "TISSUE=$covariatefilename"
	echo "TISSUE=$tissue_with_samplesize"
	sbatch  --job-name="$tissue" \
		--partition=norm \
		--time=3:00:00 \
		--mem=64g \
		--cpus-per-task=$numcpus \
		--error="$scriptdir/$tissue.calculate.fdr.slurm.err.txt" \
		--output="$scriptdir/$tissue.calculate.fdr.slurm.out.txt" \
		$scriptdir/sbatch-calculate_FDR.sh $scriptdir $tissue $tissue_with_samplesize
}

# do_calculate_FDR Adipose_Tissue Adipose_Tissue.363samples.covariate.sorted.txt

while read tissue covariatedatafilename; do
 	do_calculate_FDR ${tissue} ${covariatedatafilename}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS/list_of_covariate_files_by_tissue.txt"


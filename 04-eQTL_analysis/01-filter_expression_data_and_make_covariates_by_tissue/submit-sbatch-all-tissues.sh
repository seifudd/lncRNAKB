#!/bin/bash
#$ -cwd

function do_make_covariate_filtered_tpm_count_cpm_files () {
	scriptdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/"
	tissue=$1
	tissue_ID_R1_R2=$2
	numcpus=8
	echo "TISSUE=$tissue"
	mkdir $scriptdir/$tissue
	sbatch  --job-name="$tissue" \
		--partition=norm \
		--time=8:00:00 \
		--mem=64g \
		--gres=lscratch:150 \
		--cpus-per-task=$numcpus \
		--error="$scriptdir/$tissue.slurm.err.txt" \
		--output="$scriptdir/$tissue.slurm.out.txt" \
		$scriptdir/sbatch-make-covariate-filtered-tpm-count-cpm-files.sh $scriptdir $tissue
}

# do_make_covariate_filtered_tpm_count_cpm_files Adipose_Tissue	Adipose_Tissue_ID_R1_R2.txt

function do_PEER () {
	scriptdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS"
	tissue=$1
#	expressiondatafilename=$2
	covariatedatafilename=$2
#	numcpus=4
	echo "TISSUE=$tissue"
#	mkdir $scriptdir/$tissue
#	sbatch  --job-name="$tissue" \
#		--partition=norm \
#		--time=5:00:00 \
#		--mem=32g \
#		--gres=lscratch:150 \
#		--cpus-per-task=$numcpus \
#		--error="$scriptdir/$tissue.peer.slurm.err.txt" \
#		--output="$scriptdir/$tissue.peer.slurm.out.txt" \
sh		$scriptdir/run_PEER.sh \
		$covariatedatafilename	\
		$tissue \
		$scriptdir 
}

# do_make_covariate_filtered_tpm_count_cpm_files Adipose_Tissue	Adipose_Tissue_ID_R1_R2.txt
# do_PEER Adipose_Tissue Adipose_Tissue.363samples.15024pcs.10575lincs.gene.CPM.GTExId.txt

while read tissue covariatedatafilename; do
	do_PEER ${tissue} ${covariatedatafilename}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/list_of_covariate_files_by_tissue.txt"

# do_make_covariate_filtered_tpm_count_cpm_files ${tissue} ${tissue_ID_R1_R2}
# /data/NHLBI_BCB/Fayaz/21-GTEx-data/01-by-tissue-files-and-ids-rnaseq-only/histological_type_s.txt

# do_PEER ${tissue} ${expressiondatafilename}
# /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/list_of_filtered_cpm_files_by_tissue.txt


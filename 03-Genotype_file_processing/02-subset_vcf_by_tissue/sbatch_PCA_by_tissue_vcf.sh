#!/bin/bash
#$ -cwd

function do_PCA_by_tissue_vcf () {
	scriptdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files/"
	vcfbytissuedir="/data//NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks"
	covbytissueoutdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6"
	tissue=$1
	tissuebysubsetfile=$2
	tissuesubsetfilename=`echo $tissuebysubsetfile | sed 's/.txt//'`
	numcpus=8
	echo "TISSUE=$tissue"
#	mkdir $scriptdir/$tissue
	sbatch  --job-name="${tissue}.PCA" \
		--partition=norm \
		--time=24:00:00 \
		--mem=200g \
		--gres=lscratch:150 \
		--cpus-per-task=$numcpus \
		--error="$scriptdir/$tissue.PCA.slurm.err.txt" \
		--output="$scriptdir/$tissue.PCA.slurm.out.txt" \
		$scriptdir/PCA_by_tissue_vcf.sh $tissue $tissuesubsetfilename $vcfbytissuedir $covbytissueoutdir
}

# do_PCA_by_tissue_vcf Adipose_Tissue Adipose_Tissue.363.subset.txt

while read tissue tissuebysubsetfile; do
	do_PCA_by_tissue_vcf ${tissue} ${tissuebysubsetfile}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/subset_files/00-by_tissue_subset_files_list_redo.txt"


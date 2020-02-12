#!/bin/bash

function do_MAF_calc_filter_subset () {
	tissue=$1
	tissuebysubsetfile=$2
#	vcfsubsetdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks/blood"
	vcfsubsetdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks"
	covbytissueoutdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6"
	echo $tissue
	sbatch --job-name=$tissue.vcf_calcMAF \
			--partition=norm \
			--time=12:00:00 \
			--mem=64g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=$vcfsubsetdir/$tissue.slurm.calcMAF_filter_subset.stderr \
			--output=$vcfsubsetdir/$tissue.slurm.calcMAF_filter_subset.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks/MAF_calc_filter_subset.sh \
			$tissue \
			$tissuebysubsetfile \
			$vcfsubsetdir \
			$covbytissueoutdir

}

# do_MAF_calc_filter_subset Adipose_Tissue Adipose_Tissue.363.subset.txt
# do_MAF_calc_filter_subset Bladder Bladder.9.subset.txt
do_MAF_calc_filter_subset Salivary_Gland Salivary_Gland.63.subset.txt

# while read tissue tissuebysubsetfile; do
#	do_MAF_calc_filter_subset ${tissue} ${tissuebysubsetfile}
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/subset_files/00-by_tissue_subset_files_list_redo.txt"


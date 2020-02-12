#!/bin/bash
#$ -cwd

function do_update_rsIDs_of_tissue_VCF_files_GTEx () {
#	cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1
	
	inputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_ks"
	outputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_FS"
	tissue=$1
	bytissue_vcf_subset_file=$2

	echo TISSUE=$tissue
	mkdir -p $outputdir/$tissue

	sbatch --job-name=${tissue}.leftjoin \
			--cpus-per-task=32 \
			--time=4:00:00 \
			--mem=128g \
			--error=${tissue}.leftjoinvcfhg38.stderr \
			--output=${tissue}.leftjoinvcfhg38.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_FS/subsettedvcf_import_rsIDs.sh $tissue \
			$inputdir \
			$bytissue_vcf_subset_file \
			$outputdir		
}

# do_update_rsIDs_of_tissue_VCF_files_GTEx Salivary_Gland Salivary_Gland.63.subset.GTEx.subset.MAF0.05.FINAL.sorted.vcf.gz

while read tissue tissuebysubsetfile; do
	do_update_rsIDs_of_tissue_VCF_files_GTEx ${tissue} ${tissuebysubsetfile}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_FS/list_of_vcf_files_by_tissue.txt"



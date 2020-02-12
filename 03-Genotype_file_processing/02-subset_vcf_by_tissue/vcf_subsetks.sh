#!/bin/bash

set -e #stop script if a command fails
function fail {
echo “FAIL: $@“ >&2
exit 1 # signal failure
}

tissue=$1
tissue_subset_file=$2
vcfDir=$3
outputDir=$4

module load bcftools

function do_vcf_subset {
	bcftools view -S \
	$vcfDir/subset_files/$tissue_subset_file \
	$vcfDir/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.rsIDs.vcf > $outputDir/$tissue/$tissue.GTEx.subset.vcf

#	echo "TISSUE=$tissue"
#	echo "TISSUE_SUBSET_DIR=$tissue_subset_file"

}

do_vcf_subset


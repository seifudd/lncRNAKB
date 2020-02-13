#!/bin/sh -l
#$ -cwd

scriptdir=$1
tissue=$2
tissue_with_samplesize=$3

echo "TISSUE=$tissue"
echo "TISSUE.SAMPLE.SIZE=$tissue_with_samplesize"

cd $scriptdir
module load R
Rscript ./calculate_FDR.R $tissue $tissue_with_samplesize


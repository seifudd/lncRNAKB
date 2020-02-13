#!/bin/sh -l
#$ -cwd

scriptdir=$1
tissue=$2
cpmexpressionfile=$3

echo "TISSUE=$tissue"
echo "CPMFILE=$cpmexpressionfile"

cd $scriptdir
module load R
Rscript ./runcemitools.R $tissue $cpmexpressionfile

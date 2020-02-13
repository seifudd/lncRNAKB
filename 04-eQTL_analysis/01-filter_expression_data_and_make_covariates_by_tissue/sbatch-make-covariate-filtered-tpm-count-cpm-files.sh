#!/bin/sh -l
#$ -cwd

scriptdir=$1
tissue=$2
tissueSRA=`echo $tissue | sed 's/_/ /g'`

echo "TISSUE=$tissue"
echo "TISSUESRA=$tissueSRA"

cd $scriptdir
module load R
Rscript ./make_covariate_filtered_tpm_count_cpm_files.R $tissue $tissueSRA


#!/bin/bash

data_dir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb-v6"
#output_dir="/data/NHLBI_BCB/Abhilash/05_Tissue_Specificity/data"

for i in `find $data_dir -type f -name "*.gene.TPM.with.annotation.filtered.txt"`; do
	echo $i
	echo `wc -l $i`
	tissue=`echo $i | cut -d"/" -f 7`
	sed 's/\"//g' $i | cut -d" " -f 1,6- > "/data/NHLBI_BCB/Abhilash/09_PCA/data_v6/${tissue}_TPM.txt"
done

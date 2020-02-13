#!/bin/bash

inputdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK"

while read -r ukbbpheno ukphenodownloadlink
do
	mkdir -p "$inputdir/$ukbbpheno"
	cd "$inputdir/$ukbbpheno"
	$ukphenodownloadlink
	filename=`echo $ukphenodownloadlink | awk '{print $4}'`
	newfilename=`echo $filename | sed 's/.bgz/.gz/'`
	mv -f $filename $newfilename
	gunzip $newfilename

done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/ukbb_phenotypes_only_names_and_downloadlinks_Ncases50000_for_SMR_lncRNAKb.txt"


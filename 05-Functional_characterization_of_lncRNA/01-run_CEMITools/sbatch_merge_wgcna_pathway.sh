#!/bin/bash



function do_merge () {
	tissue=$1
	dir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo"
	echo $tissue
	cd /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo
	sbatch --job-name=merge \
			--partition=norm \
			--time=1:00:00 \
			--mem=128g \
			--gres=lscratch:150 \
			--cpus-per-task=4 \
			--error=/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/merge.stderr \
			--output=/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/merge.stdout \
			/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/merge_wgcna_pathway.sh \
			$tissue \
			$dir 
}

while read tissue; do
		do_merge ${tissue}
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/tissue2.txt"
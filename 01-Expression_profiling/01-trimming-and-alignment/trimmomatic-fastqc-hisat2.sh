#!/bin/bash 

function do_trimmomatic_fastqc_hisat () {

tissue=$1
tissue_ID_R1_R2=$2
basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"
bytissuefilesdir="$basedir/00-by-tissue-files-and-ids"
datadir="$basedir/dbGaP-13516/fastq.gzs"
fastqcdir="$basedir/01-fastqc-gtex-trimmomatic"
hisatdir="$basedir/03-trimmomatic-alignment-hisat2"
mkdir "$hisatdir/$tissue"
reference="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/grch38_tran"
numcpus=4

while read SAMPLE READ1 READ2; do
    	echo $tissue,$SAMPLE
	mkdir "$hisatdir/$tissue/$SAMPLE"
	sbatch  --job-name="$tissue.$SAMPLE" \
		--partition=norm \
		--time=24:00:00 \
		--mem=64g \
		--gres=lscratch:150 \
		--cpus-per-task=$numcpus \
		--error="$hisatdir/$tissue/$SAMPLE.slurm.err.txt" \
		--output="$hisatdir/$tissue/$SAMPLE.slurm.out.txt" \
		"$hisatdir/sbatch-trimmomatic-fastqc-hisat2.sh $SAMPLE $datadir $READ1 $READ2 $reference $hisatdir $tissue $numcpus $fastqcdir"
done < "$bytissuefilesdir/$tissue_ID_R1_R2"

}

#do_trimmomatic_fastqc_hisat Heart	Heart_ID_R1_R2.txt
do_trimmomatic_fastqc_hisat Colon Colon_ID_R1_R2.txt

# basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"

# while read tissue tissue_ID_R1_R2; do
# 	do_trimmomatic_fastqc_hisat ${tissue} ${tissue_ID_R1_R2}
# done < "$basedir/01-fastqc-gtex/histological_type_s.txt"


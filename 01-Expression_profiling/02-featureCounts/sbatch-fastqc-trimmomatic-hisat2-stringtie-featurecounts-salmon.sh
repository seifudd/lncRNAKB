#!/bin/bash 

function do_fastqc_trimmomatic_hisat_stringtie_featurecounts_salmon () {
	tissue=$1
	tissue_ID_R1_R2=$2
	basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"
	bytissuefilesdir="$basedir/01-by-tissue-files-and-ids-rnaseq-only"
	datadir="$basedir/dbGaP-13516/fastq.gzs"
#	outdir="$basedir/05-featurecounts-lncRNAKb"
	outdir="$basedir/05-featurecounts-lncRNAKb-v6"
#	mkdir "$outdir/$tissue"
	reference="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/grch38_tran"
	numcpus=4

#	while read SAMPLE READ1 READ2; do
#	    	echo $tissue,$SAMPLE
	    	echo $tissue
#		mkdir "$outdir/$tissue/$SAMPLE"
		SAMPLE="null"
		READ1="null"
		READ2="null"
#		sbatch  --job-name="$tissue.$SAMPLE.fc" \
#			--partition=norm \
#			--time=24:00:00 \
#			--mem=64g \
#			--gres=lscratch:150 \
#			--cpus-per-task=$numcpus \
#			--error="$outdir/$tissue/$SAMPLE.$tissue.slurm.err.txt" \
#			--output="$outdir/$tissue/$SAMPLE.$tissue.slurm.out.txt" \
#		 	$outdir/fastqc-trimmomatic-hisat2-stringtie-featurecounts-salmon.sh $SAMPLE $datadir $READ1 $READ2 $reference $outdir $numcpus $tissue $bytissuefilesdir/$tissue_ID_R1_R2
#		sbatch  --job-name="$tissue.merge.fc" \
#			--partition=norm \
#			--time=4:00:00 \
#			--mem=64g \
#			--cpus-per-task=$numcpus \
#			--error="$outdir/$tissue/$tissue.fc.merge.slurm.err.txt" \
#			--output="$outdir/$tissue/$tissue.fc.merge.slurm.out.txt" \
sh		 	$outdir/fastqc-trimmomatic-hisat2-stringtie-featurecounts-salmon.sh $SAMPLE $datadir $READ1 $READ2 $reference $outdir $numcpus $tissue $bytissuefilesdir/$tissue_ID_R1_R2
#	done < "$bytissuefilesdir/$tissue_ID_R1_R2"
}

# To Run individual tissues example
# do_fastqc_trimmomatic_hisat_stringtie_featurecounts_salmon Heart	Heart_ID_R1_R2.txt
# do_fastqc_trimmomatic_hisat_stringtie_featurecounts_salmon Blood	Blood_ID_R1_R2.txt

basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"

while read tissue tissue_ID_R1_R2; do
 	do_fastqc_trimmomatic_hisat_stringtie_featurecounts_salmon ${tissue} ${tissue_ID_R1_R2}
done < "$basedir/01-by-tissue-files-and-ids-rnaseq-only/histological_type_s.txt"


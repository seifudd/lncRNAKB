#!/bin/bash 

function do_fastqc_trimmomatic_hisat_stringtie_featurecounts_salmon () {
	########################################################################################################################
	tissue=$1
	tissue_R1_R2_filelist=$2
	basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"
	datadir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/fastq.gzs"
	scriptdir="$basedir/05-salmon"
	outdir="$basedir/05-salmon"
	reference="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/grch38_tran"
	numcpus=4
	transcript_reference="/data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.primary.assembly.fa"

#	while read SAMPLE READ1 READ2; do
#	while read SAMPLE; do
#	    	echo ${tissue}.${SAMPLE}
#	    	echo ${tissue},${tissue_R1_R2_filelist}
#		mkdir -p "$outdir/$tissue"
#		mkdir -p "$outdir/$tissue/$SAMPLE"
		SAMPLE="null"
		READ1="null"
		READ2="null"
#		tissue="null"
#		tissue_R1_R2_filelist="null"
#		--error="$outdir/build_salmon_index.slurm.err.txt" \
#		--output="$outdir/build_salmon_index.slurm.out.txt" \
#		--error="$outdir/$tissue/${tissue}.${SAMPLE}.salmon.slurm.err.txt" \
#		--output="$outdir/$tissue/${tissue}.${SAMPLE}.salmon.slurm.out.txt" \
#		sbatch  --job-name="build_salmon_index" \
#		sbatch  --job-name="$tissue.${SAMPLE}.salmon" \
		sbatch  --job-name="$tissue.${SAMPLE}.merge.salmon" \
			--error="$outdir/$tissue/${tissue}.${SAMPLE}.merge.salmon.slurm.err.txt" \
			--output="$outdir/$tissue/${tissue}.${SAMPLE}.merge.salmon.slurm.out.txt" \
			--time=6:00:00 \
			--partition=norm \
			--mem=200g \
			--cpus-per-task=$numcpus \
			"$scriptdir/fastqc-trimmomatic-hisat2-stringtie-featurecounts-salmon.sh" ${SAMPLE} $datadir $READ1 $READ2 $reference $outdir $numcpus $transcript_reference $tissue ${tissue_R1_R2_filelist}
#	done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-by-tissue-files-and-ids-rnaseq-only/${tissue_R1_R2_filelist}"
	########################################################################################################################
}

do_fastqc_trimmomatic_hisat_stringtie_featurecounts_salmon Fallopian_Tube Fallopian_Tube_ID_R1_R2.txt

# while read tissue tissue_R1_R2_filelist
# do
# 	do_fastqc_trimmomatic_hisat_stringtie_featurecounts_salmon $tissue $tissue_R1_R2_filelist
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-by-tissue-files-and-ids-rnaseq-only/histological_type_s.txt"


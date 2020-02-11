#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

SAMPLE=$1
DATAPATH=$2
READ1=$3
READ2=$4
REF=$5
hisatdir=$6
tissue=$7
numcpus=$8
fastqcdir=$9

module load hisat
module load samtools
module load trimmomatic
module load fastqc

date
java -jar $TRIMMOJAR PE \
            -threads $numcpus \
            "$DATAPATH/$READ1" \
            "$DATAPATH/$READ2" \
	    -baseout "/lscratch/${SLURM_JOBID}/${SAMPLE}.fastq.gz" \
            ILLUMINACLIP:"/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
            MINLEN:50

echo "trimmomatic done"
date

DATAPATH_TRIMMED_READS="/lscratch/${SLURM_JOBID}"

# fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

mkdir -p "$fastqcdir/$tissue/$SAMPLE"

fastqc -o "$fastqcdir/$tissue/$SAMPLE"  \
  --nogroup \
  "$DATAPATH_TRIMMED_READS/${SAMPLE}_1P.fastq.gz"  \
  "$DATAPATH_TRIMMED_READS/${SAMPLE}_2P.fastq.gz"  \
  || fail "fastqc failed"


echo "fastqc done"
date

hisat2 	-p $numcpus \
		-x $REF/genome_tran \
		--downstream-transcriptome-assembly \
		-1 "$DATAPATH_TRIMMED_READS/${SAMPLE}_1P.fastq.gz" \
		-2 "$DATAPATH_TRIMMED_READS/${SAMPLE}_2P.fastq.gz" \
		--rg-id $SAMPLE --rg SM:$SAMPLE \
| samtools view -h -f 3 -O SAM - \
| perl -nle  'print if m/^@(?:[A-Z]{2})\s|\bNH:i:1\b/' \
| samtools sort -@ $numcpus \
	-o "$hisatdir/$tissue/$SAMPLE/$SAMPLE.unique.bam" \
	-T /lscratch/${SLURM_JOB_ID}/${SAMPLE}_chunk -
   
samtools index "$hisatdir/$tissue/$SAMPLE/$SAMPLE.unique.bam"
#		-S "$hisatdir/$tissue/$SAMPLE/$SAMPLE.out.sam" \
# samtools sort -@ 8 -o "$hisatdir/$tissue/$SAMPLE/$SAMPLE.out.sorted.bam" "$hisatdir/$tissue/$SAMPLE/$SAMPLE.out.sam"

# exclude unmapped (-F 4) : mapped
# proper mapped -f 3
# exclude multiple-mapped (grep -v "XS:") 
# NH:i:n # to get the uniq : n=1

# samtools view -H  "$hisatdir/$tissue/$SAMPLE/$SAMPLE.out.sorted.bam" > "$hisatdir/$tissue/$SAMPLE/$SAMPLE.header.sam"
# samtools view -f 3 "$hisatdir/$tissue/$SAMPLE/$SAMPLE.out.sorted.bam" \
#     | grep  -w "NH:i:1" | cat "$hisatdir/$tissue/$SAMPLE/$SAMPLE.header.sam" - \
#     | samtools view -b - > "$hisatdir/$tissue/$SAMPLE/$SAMPLE.unique.bam"

# samtools index "$hisatdir/$tissue/$SAMPLE/$SAMPLE.unique.bam"

# rm "$hisatdir/$tissue/$SAMPLE/$SAMPLE.header.sam"
# rm "$hisatdir/$tissue/$SAMPLE/$SAMPLE.out.sam"
echo "hisat2 done"
date

# rm -f "$DATAPATH/$READ1"
# rm -f "$DATAPATH/$READ2"
# date

# python extract_data.py --tissue_name $tissue
# echo 'Data extraction done.'
# date

echo 'Changing Permissions of folders.'

chmod -R 777 "$hisatdir/$tissue/$SAMPLE/"
chmod -R 777 "$fastqcdir/$tissue/$SAMPLE/"

date



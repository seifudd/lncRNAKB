#!/bin/bash

set -e   

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

SAMPLE=$1
DATAPATH=$2
READ1=$3
READ2=$4
OUTPUTDIR=$5
tissue=$6

module load fastqc

# fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

mkdir -p "$OUTPUTDIR/$tissue/$SAMPLE" 

fastqc -o "$OUTPUTDIR/$tissue/$SAMPLE"  \
  $DATAPATH/$READ1  \
  $DATAPATH/$READ2  \
  || fail "fastqc failed"

date
echo "fastqc done"


#!/bin/bash 

module load fastqc

function do_fastqc () {
tissue=$1
tissue_ID_R1_R2=$2
basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"
datadir="$basedir/dbGaP-13516/fastq.gzs/"
fastqcdir="$basedir/01-fastqc-gtex"
mkdir "$fastqcdir/$tissue"
while read SAMPLE READ1 READ2; do
    	echo $tissue,$SAMPLE
	sbatch  --job-name="$tissue.$SAMPLE" \
		--partition=norm \
		--time=2:00:00 \
		--mem=10g \
		--cpus-per-task=2 \
		--error="$fastqcdir/$tissue/$SAMPLE.slurm.err.txt" \
		--output="$fastqcdir/$tissue/$SAMPLE.slurm.out.txt" \
		"$fastqcdir/sbatch-fastqc.sh $SAMPLE $datadir $READ1 $READ2 $fastqcdir $tissue"
done < "$fastqcdir/$tissue_ID_R1_R2"
}

# while read tissue tissue_ID_R1_R2; do
#	do_fastqc ${tissue} ${tissue_ID_R1_R2}
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-fastqc-gtex/histological_type_s.txt"

function do_unzip () {
tissue=$1
tissue_ID_R1_R2=$2
basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"
datadir="$basedir/dbGaP-13516/fastq.gzs/"
fastqcdir="$basedir/01-fastqc-gtex"
	# unzip them
	while read SAMPLE READ1 READ2; do
		echo $tissue,$SAMPLE
		unzip "$fastqcdir/$tissue/$SAMPLE/${SAMPLE}_1_fastqc.zip" -d $fastqcdir/$tissue/$SAMPLE/
		unzip "$fastqcdir/$tissue/$SAMPLE/${SAMPLE}_2_fastqc.zip" -d $fastqcdir/$tissue/$SAMPLE/
		wait
	done < "$fastqcdir/$tissue_ID_R1_R2"
}

# while read tissue tissue_ID_R1_R2; do
#	do_unzip ${tissue} ${tissue_ID_R1_R2}
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-fastqc-gtex/histological_type_s.txt"

function get_fastqc_stats
{

basedir="/data/NHLBI_BCB/Fayaz/21-GTEx-data"
datadir="$basedir/dbGaP-13516/fastq.gzs/"
fastqcdir="$basedir/01-fastqc-gtex"

rm -f tmp1

while read SAMPLE READ1 READ2; do
  	cd $fastqc_dir/$i/${i}_R1_001_fastqc/
  	csplit fastqc_data.txt '/>>/' {*}
  	seqn=`grep "Overrepresented sequences" xx* | head -1 | cut -f 1 -d':'`
	cat $seqn | sed '1,2d' | awk '{if($3>0) {print $0} }' | cut -f 3 | awk '{sum+=$NF} END {print sum}'
	



  cat $seqn | sed '1,2d' | awk '{print $0}' | cut -f 1,4 | cut -f 1 -d',' >> ../../tmp1

  cd $fastqc_dir/$i/${i}_R2_001_fastqc/
  csplit fastqc_data.txt '/>>/' {*}
  seqn=`grep "Overrepresented sequences" xx* | head -1 | cut -f 1 -d':'`
  cat $seqn | sed '1,2d' | awk '{print $0}' | cut -f 1,4 | cut -f 1 -d',' >> ../../tmp1
done < "$fastqcdir/$tissue_ID_R1_R2"

cd $fastqc_dir
count=1
for i in `cat tmp11 | cut -f 1 | sort | uniq ` ; do
  echo "> overrep $count" >> tmp2
  echo $i >> tmp2
   (( count += 1 ))
done

rm -f tmp1

}

# while read tissue tissue_ID_R1_R2; do
#	do_fastqc ${tissue} ${tissue_ID_R1_R2}
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-fastqc-gtex/histological_type_s.txt"



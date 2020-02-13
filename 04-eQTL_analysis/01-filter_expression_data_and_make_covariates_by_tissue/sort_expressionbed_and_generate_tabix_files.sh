#!/bin/bash

module load samtools

expressiondatafile=$1
tissue=$2
inputdir=$3
outputdir=$4

expressiondatafilename=`echo $expressiondatafile | sed 's/.txt//g'`

cat $inputdir/$tissue/$expressiondatafile | sed '1,1d' > $outputdir/$tissue/$tissue.temp
cat $outputdir/$tissue/$tissue.temp | cut -f6- > $outputdir/$tissue/$tissue.expression
cat $outputdir/$tissue/$tissue.temp | cut -f2-4 > $outputdir/$tissue/$tissue.chr_start_end
cat $outputdir/$tissue/$tissue.temp | cut -f1 > $outputdir/$tissue/$tissue.gene_id
paste $outputdir/$tissue/$tissue.chr_start_end $outputdir/$tissue/$tissue.gene_id $outputdir/$tissue/$tissue.expression > $outputdir/$tissue/$tissue.temp2

cat $inputdir/$tissue/$expressiondatafile | head -n1 |sed 's/.chr/chr/g' |sed 's/\./-/g' | awk '{print "gene_id""\t"$0}' | awk '{print "#"$2"\t"$3"\t"$4"\t"$1}' > $outputdir/$tissue/$tissue.temp3
cat $inputdir/$tissue/$expressiondatafile | head -n1 |sed 's/.chr/chr/g' |sed 's/\./-/g' | awk '{print "gene_id""\t"$0}' | cut -f6- > $outputdir/$tissue/$tissue.temp4
paste $outputdir/$tissue/$tissue.temp3 $outputdir/$tissue/$tissue.temp4 > $outputdir/$tissue/$tissue.temp5
cat $outputdir/$tissue/$tissue.temp5 $outputdir/$tissue/$tissue.temp2 > $outputdir/$tissue/$tissue.temp6

sort -k1,1 -k2,2n $outputdir/$tissue/$tissue.temp6 > $outputdir/$tissue/$tissue.temp7
cat $outputdir/$tissue/$tissue.temp7 | grep -v "chrY" > $outputdir/$tissue/$tissue.temp8
bgzip -c $outputdir/$tissue/$tissue.temp8 > $inputdir/$tissue/$expressiondatafilename.sorted.bed.gz
tabix -f -p bed $inputdir/$tissue/$expressiondatafilename.sorted.bed.gz

rm -f $outputdir/$tissue/$tissue.temp
rm -f $outputdir/$tissue/$tissue.expression
rm -f $outputdir/$tissue/$tissue.chr_start_end
rm -f $outputdir/$tissue/$tissue.gene_id
rm -f $outputdir/$tissue/$tissue.temp2
rm -f $outputdir/$tissue/$tissue.temp3
rm -f $outputdir/$tissue/$tissue.temp4
rm -f $outputdir/$tissue/$tissue.temp5
rm -f $outputdir/$tissue/$tissue.temp6
rm -f $outputdir/$tissue/$tissue.temp7
rm -f $outputdir/$tissue/$tissue.temp8


#!/bin/bash

set -e

module load samtools
module load picard

tissue=$1
inputdir=$2
bytissue_vcf_subset_file=$3
outputdir=$4

# zcat $inputdir/$tissue/$bytissue_vcf_subset_file | grep -v ^#  | awk '{print $1"_"$2"_"$4"_"$5 "\t" $0}'> $outputdir/$tissue/$tissue.noHash.vcf
#cat hg38_avsnp150_SNPonly.txt | awk '{print "chr"$1"_"$2"_"$4"_"$5 "\t" $6}' > hg38_avsnp150_SNPonlyV3.txt
# zcat $inputdir/$tissue/$bytissue_vcf_subset_file |  grep '#' >  $outputdir/$tissue/$tissue.Hash.vcf
# awk 'NR==FNR{a[$1]=$2;next}NR!=FNR{c=$1; if(c in a){print $0 "\t" a[c]}}' /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/hg38_avsnp150_SNPonlyV3.txt $outputdir/$tissue/$tissue.noHash.vcf | awk '{$4=$NF}1' | awk -v OFS="\t" 'NF{NF-=1}1' | cut -f2- > $outputdir/$tissue/$tissue.noHashV2.vcf
# cat $outputdir/$tissue/$tissue.noHashV2.vcf >> $outputdir/$tissue/$tissue.Hash.vcf
# awk 'NR==FNR{a[$1];next}!($1 in a)' /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/hg38_avsnp150_SNPonlyV3.txt $outputdir/$tissue/$tissue.noHash.vcf | cut -f2- >> $outputdir/$tissue/$tissue.Hash.vcf
# rm -f $outputdir/$tissue/$tissue.noHashV2.vcf $outputdir/$tissue/$tissue.noHash.vcf
bytissue_vcf_subset_file_without_gz_extension=`echo $bytissue_vcf_subset_file | sed 's/.gz//'`
zcat $outputdir/$tissue/$bytissue_vcf_subset_file > $outputdir/$tissue/$bytissue_vcf_subset_file_without_gz_extension
bytissue_vcf_subset_file_without_vcf_extension=`echo $bytissue_vcf_subset_file_without_gz_extension | sed 's/.vcf//'`
java -Xmx100g -XX:ParallelGCThreads=16 -jar $PICARDJARPATH/picard.jar SortVcf I=$outputdir/$tissue/$bytissue_vcf_subset_file_without_gz_extension O=$outputdir/$tissue/${bytissue_vcf_subset_file_without_vcf_extension}.by.picard.vcf
# mv -f $outputdir/$tissue/$tissue.Hash.vcf  $outputdir/$tissue/$bytissue_vcf_subset_file
rm -f $outputdir/$tissue/$bytissue_vcf_subset_file
rm -f $outputdir/$tissue/${bytissue_vcf_subset_file}.tbi
bgzip -c $outputdir/$tissue/${bytissue_vcf_subset_file_without_vcf_extension}.by.picard.vcf > $outputdir/$tissue/$bytissue_vcf_subset_file
tabix -f -p vcf $outputdir/$tissue/$bytissue_vcf_subset_file



chr1	mitranscriptome	transcript	77506	79261	1000.0	-	.	tcat "lncrna"; gene_id "G000006"; tss_id "TSS000008"; uce "FALSE"; transcript_id "T000008"; tstatus "unannotated"; tgenic "intergenic"; func_name_final "NA";
chr1	mitranscriptome	exon	77506	78400	1000.0	-	.	exon_number "0"; tcat "lncrna"; gene_id "G000006"; tss_id "TSS000008"; uce "FALSE"; transcript_id "T000008"; tstatus "unannotated"; tgenic "intergenic"; func_name_final "NA";
chr1	mitranscriptome	exon	78906	79261	1000.0	-	.	exon_number "1"; tcat "lncrna"; gene_id "G000006"; tss_id "TSS000008"; uce "FALSE"; transcript_id "T000008"; tstatus "unannotated"; tgenic "intergenic"; func_name_final "NA";
chr1	mitranscriptome	transcript	157150	158081	1000.0	.	.	tcat "lncrna"; gene_id "G000008"; tss_id "TSS000017"; uce "FALSE"; transcript_id "T000026"; tstatus "annotated"; tgenic "intergenic"; func_name_final "NA";
chr1	mitranscriptome	exon	157150	158081	1000.0	.	.	exon_number "0"; tcat "lncrna"; gene_id "G000008"; tss_id "TSS000017"; uce "FALSE"; transcript_id "T000026";

gtf="mitranscriptome.hg38.v2.lncrna.gtf"
db="mitranscriptome"

#### transcripts row number # 175259
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR }' > tmp1
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR-1 }' | sed '1,1d' > tmp2
cat $gtf | wc -l >> tmp2
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 | sed 's/; /\n/g' | grep gene_id | sed 's/gene_id "//g' | sed 's/"//g' > tmp3
paste tmp1 tmp2 tmp3 > chr_bp_genename_transcrips_all


# now for each gene: # 63505
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 | sed 's/; /\n/g' | grep gene_id | sed 's/gene_id "//g' | sed 's/"//g' |  uniq | sort | uniq > gene_names

run_it(){
m=$1
n=$2
cat gene_names | sed -n "$m,${n}p"   | while read i; do
  st=`cat $gtf  | grep -n "\"$i\";" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "\"$i\";" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep    "\"$i\";" | cut -f 1 | head -1`
  strand=`cat $gtf     | grep    "\"$i\";" | cut -f 7 | head -1`
  gene_id=`cat $gtf    | grep -n "\"$i\";" | cut -f 9 | sed 's/; /\n/g;' | sort -r | uniq | grep gene_id | cut -f 2 -d' ' | sed 's/"//g'`
#  gene_alias=`cat $gtf | grep -n "$i\";" | cut -f 9 | sed 's/;/\n/g;' | sort -r | uniq | grep gene_alias | cut -f 2 -d'=' | cut -f 1 -d'.' |   sort | uniq | tr '\r\n' ' ' | sed 's/ /;gene_alias=/g' | awk '{print "gene_alias="$0}' | rev | cut -c 13- | rev`
  mkdir -p ${db}_genes/$chr
  echo -e "$chr\t$db\tgene\t$st\t$end\t.\t$strand\t.\tgene_id \"$gene_id\";$gene_alias" > ${db}_genes/$chr/${chr}_${st}_$end
  
  cat chr_bp_genename_transcrips_all | grep -w $i |  while read a b c; do 
    cat $gtf |  sed -n "$a,${b}p" >> ${db}_genes/$chr/${chr}_${st}_$end
  done
done
}

  cat chr_bp_genename_transcrips_all | grep -w $i |  while read a b c; do 
     echo $a  $b $c
     cat $gtf |  sed -n "$a,${b}p" 
  done


# loop it
for i in {0..14}; do
   let m=($i*4000)+1 ; let n=($i+1)*4000;   run_it $m $n &
done
run_it 60001 63505 &
wait

for i in `ls ${db}_genes ` ; do ls ${db}_genes/$i/* | cut -f 3 -d'/' >> gene_bp_all; done
# 63489 #16 less than gene_names. were they duplicates? :

run_it(){
m=$1
n=$2
o=$3
cat gene_names | sed -n "$m,${n}p"   | while read i; do
  st=`cat $gtf | grep -n "\"$i\"" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "\"$i\"" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep  "\"$i\"" | cut -f 1 | head -1`
  strand=`cat $gtf | grep  "\"$i\"" | cut -f 7 | head -1`
  echo $i
  echo -e "$chr\t$st\t$end\t$strand\t$i" >> chr_bp_gene_$o
done
}
# loop
for i in {0..14}; do
   let m=($i*4000)+1 ; let n=($i+1)*4000;   run_it $m $n $i &
done
run_it 60001 63505 15 &
wait
cat chr_bp_gene_* | sort | uniq > tmp1 ; mv tmp1 chr_bp_gene
cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1>1) print $2"\t"$3"\t"$4}'  > duplicates_chr_bp 
wc -l duplicates_chr_bp  #16
cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1<2) print $2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n > tmp1
cat tmp1 | while read a b c; do 
   cat chr_bp_gene | grep -P "$a\t$b\t$c" >> chr_bp_gene_dedup
done
wc -l chr_bp_gene*
  63505 chr_bp_gene
  63473 chr_bp_gene_dedup

################################## no need here

mkdir -p ${db}_genes_attr
run_it(){
m=$1
n=$2
cat chr_bp_gene_dedup | sed -n "$m,${n}p"   | while read chr b c e f ; do
  mkdir -p ${db}_genes_attr/$chr
  cat ${db}_genes/$chr/${chr}_${b}_${c} | head -1 | cut -f 9 | sed 's/;/\n/g' | sed 's/=/\t/g' > ${db}_genes_attr/$chr/${chr}_${b}_${c}
done
}

# loop it


###################################

#chr1	mitranscriptome	transcript	832	72894	1000.0	+	.	tcat "lncrna"; gene_id "G004230"; tss_id "TSS008698"; transcript_id "T018741";  

mkdir -p ${db}_trans
run_it(){
m=$1
n=$2
cat chr_bp_gene_dedup | sed -n "$m,${n}p"   | while read chr b c1 e1 f1 ; do
  mkdir -p ${db}_trans/$chr
   c="${chr}_${b}_${c1}"
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR }' > tmp11_$c
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR-1 }' | sed '1,1d' > tmp22_$c
   cat ${db}_genes/$chr/$c | wc -l >> tmp22_$c
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print $1"_"$4"_"$5"\t"$9 }' \
                           | sed 's/transcript_id /\t/g'  | cut -f 1,3 | sed 's/;/\t/g' | cut -f 1,2 | sed 's/\t/_/g' | sed 's/:/_/g'  | sed 's/"//g' > tmp_tname_$c
   paste tmp11_$c tmp22_$c tmp_tname_$c > tmp_tall_$c
   mkdir -p  ${db}_trans/$chr/$c
   cat tmp_tall_$c | while read d e f; do
     cat ${db}_genes/$chr/$c | sed -n "$d,${e}p" > ${db}_trans/$chr/$c/$f
   done
   rm -f tmp11_$c tmp22_$c tmp_tname_$c tmp_tall_$c
done
}

# loop it
for i in {0..14}; do
   let m=($i*4000)+1 ; let n=($i+1)*4000;   run_it $m $n $i &
done
run_it 60001 63473 15 &
wait


############################ no need here

mkdir -p ${db}_trans_attr
run_it(){
m=$1
n=$2
cat chr_bp_gene_dedup | sed -n "$m,${n}p"  | while read chr b c1 e1 f1 ; do
   mkdir -p ${db}_trans_attr/$chr
   c="${chr}_${b}_${c1}"
   for i in `ls ${db}_trans/$chr/$c`; do
     mkdir -p ${db}_trans_attr/$chr/$c/
     cat ${db}_trans/$chr/$c/$i | head -1 | cut -f 9 | sed 's/;/\n/g' | sed 's/=/\t/g' > ${db}_trans_attr/$chr/$c/$i
   done
done
}

# loop it


############## 


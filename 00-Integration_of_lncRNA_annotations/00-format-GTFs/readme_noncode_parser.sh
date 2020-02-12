
chr1	Cufflinks	transcript	11872	14412	0	+	.	gene_id "NONHSAG000001.2"; transcript_id "NONHSAT000002.2"; FPKM "0"; exon_number 3;
chr1	Cufflinks	exon	11872	12227	0	+	.	gene_id "NONHSAG000001.2"; transcript_id "NONHSAT000002.2"; FPKM "0"; exon_number 3;
chr1	Cufflinks	exon	12613	12721	0	+	.	gene_id "NONHSAG000001.2"; transcript_id "NONHSAT000002.2"; FPKM "0"; exon_number 3;
chr1	Cufflinks	exon	13225	14412	0	+	.	gene_id "NONHSAG000001.2"; transcript_id "NONHSAT000002.2"; FPKM "0"; exon_number 3;
chr1	Cufflinks	transcript	11874	14409	0	+	.	gene_id "NONHSAG000001.2"; transcript_id "NONHSAT000003.2"; FPKM "0"; exon_number 4;
chr1	Cufflinks	exon	11874	12227	0	+	.	gene_id "NONHSAG000001.2"; transcript_id "NONHSAT000003.2"; FPKM "0"; exon_number 4;
chr1	Cufflinks	exon	12595	12721	0	+	.	gene_id "NONHSAG000001.2"; transcript_id "NONHSAT000003.2"; FPKM "0"; exon_number 4;

gtf="NONCODEv5_human_hg38_lncRNA.gtf"
db="NONCODE"

#### transcripts row number
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR }' > tmp1
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR-1 }' | sed '1,1d' > tmp2
cat $gtf | wc -l >> tmp2
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 | sed 's/;/\t/g' | cut -f 1 | sed 's/gene_id "//g' | sed 's/"//g' > tmp3
paste tmp1 tmp2 tmp3 > chr_bp_genename_transcrips_all


# now for each gene: #96308
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 | sed 's/;/\t/g' | cut -f 1 | sed 's/gene_id "//g' | sed 's/"//g' |  uniq | sort | uniq > gene_names

run_it(){
m=$1
n=$2
cat gene_names | sed -n "$m,${n}p"   | while read i; do
  st=`cat $gtf  | grep -n "$i\";" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "$i\";" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep    "$i\";" | cut -f 1 | head -1`
  strand=`cat $gtf     | grep    "$i\";" | cut -f 7 | head -1`
  gene_id=`cat $gtf    | grep -n "$i\";" | cut -f 9 | sed 's/;/\n/g;' | sort -r | uniq | grep gene_id | cut -f 2 -d' ' | sed 's/"//g'`
#  gene_alias=`cat $gtf | grep -n "$i\";" | cut -f 9 | sed 's/;/\n/g;' | sort -r | uniq | grep gene_alias | cut -f 2 -d'=' | cut -f 1 -d'.' |   sort | uniq | tr '\r\n' ' ' | sed 's/ /;gene_alias=/g' | awk '{print "gene_alias="$0}' | rev | cut -c 13- | rev`
  mkdir -p ${db}_genes/$chr
  echo -e "$chr\t$db\tgene\t$st\t$end\t.\t$strand\t.\tgene_id \"$gene_id\";$gene_alias" > ${db}_genes/$chr/${chr}_${st}_$end
  
  cat chr_bp_genename_transcrips_all | grep -w $i |  while read a b c; do 
    cat $gtf |  sed -n "$a,${b}p" >> ${db}_genes/$chr/${chr}_${st}_$end
  done
done
}

# loop it
for i in {0..14}; do
   let m=($i*6000)+1 ; let n=($i+1)*6000;   run_it $m $n &
done
run_it 84001 96308 &
wait

for i in `ls ${db}_genes | head -1` ; do ls ${db}_genes/$i/* | cut -f 3 -d'/' >> gene_bp_all; done
# 95944 #364 less than gene_names. were they duplicates? :

cat gene_names | sed '1,53657d' | while read i; do
  st=`cat $gtf | grep -n "\"$i\"" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "\"$i\"" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep  "\"$i\"" | cut -f 1 | head -1`
  strand=`cat $gtf | grep  "\"$i\"" | cut -f 7 | head -1`
  echo $i
  echo -e "$chr\t$st\t$end\t$strand\t$i" >> chr_bp_gene #96308
done
cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1>1) print $2"\t"$3"\t"$4}'  > duplicates_chr_bp #364

# remove #364x2=728 rows duplicates
cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1<2) print $2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n  > tmp1 #95580
cat tmp1 | while read a b c; do 
   cat chr_bp_gene | grep -P "$a\t$b\t$c" >> chr_bp_gene_dedup
done
wc -l chr_bp_gene*
  96308 chr_bp_gene
  95580 chr_bp_gene_dedup


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
for i in {0..14}; do
   let m=($i*6000)+1 ; let n=($i+1)*6000;   run_it $m $n &
done
run_it 84001 MPMPMP &
wait


###################################

chr1	11869	14412	+	NONHSAG000001.2
chr1	14407	29370	-	NONHSAG000002.2
chr=chr1
b=14407
c1=29370
e1="-"
f1=NONHSAG000002.2


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
                           | sed 's/; transcript_id /\t/g'  | cut -f 1,3 | sed 's/; /\t/g' | sed 's/"//g' | cut -f 1,2 | sed 's/\t/_/g' | sed 's/:/_/g' > tmp_tname_$c
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
   let m=($i*6000)+1 ; let n=($i+1)*6000;   run_it $m $n &
done
run_it 90001 95580 &
wait


############################ # no need

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



chr1	CAFE	transcript	29554	31097	.	+	.	gene_id "BIG-gene-chr1-4"; transcript_id "BIG-trx-chr1-4.1"; tag "known_lncRNA,ENST00000473358.1";
chr1	CAFE	exon	29554	30039	.	+	.	gene_id "BIG-gene-chr1-4"; transcript_id "BIG-trx-chr1-4.1"; tag "known_lncRNA,ENST00000473358.1";
chr1	CAFE	exon	30564	30667	.	+	.	gene_id "BIG-gene-chr1-4"; transcript_id "BIG-trx-chr1-4.1"; tag "known_lncRNA,ENST00000473358.1";
chr1	CAFE	exon	30976	31097	.	+	.	gene_id "BIG-gene-chr1-4"; transcript_id "BIG-trx-chr1-4.1"; tag "known_lncRNA,ENST00000473358.1";
chr1	CAFE	transcript	30267	31109	.	+	.	gene_id "BIG-gene-chr1-4"; transcript_id "BIG-trx-chr1-4.2"; tag "known_lncRNA,ENST00000469289.1";
chr1	CAFE	exon	30267	

gtf="BIGTranscriptome_lncRNA_catalog.hg38.gtf"
db="BIGTranscriptome"

#### transcripts row number #26591
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR }' > tmp1
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR-1 }' | sed '1,1d' > tmp2
cat $gtf | wc -l >> tmp2
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 | sed 's/;/\t/g' | cut -f 1 | sed 's/gene_id "//g' | sed 's/"//g' > tmp3
paste tmp1 tmp2 tmp3 > chr_bp_genename_transcrips_all


# now for each gene: #14090
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 | sed 's/;/\t/g' | cut -f 1 | sed 's/gene_id "//g' | sed 's/"//g' |  uniq | sort | uniq > gene_names

run_it(){
m=$1
n=$2
cat gene_names | sed -n "$m,${n}p"   | while read i; do
  st=`cat $gtf  | grep -n "$i\";" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "$i\";" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep    "$i\";" | cut -f 1 | head -1`
  strand=`cat $gtf     | grep    "$i\";" | cut -f 7 | head -1`
  gene_id=`cat $gtf    | grep -n "$i\";" | cut -f 9 | sed 's/; /\n/g;' | sort -r | uniq | grep gene_id | cut -f 2 -d' ' | sed 's/"//g'`
#  gene_alias=`cat $gtf | grep -n "$i\";" | cut -f 9 | sed 's/;/\n/g;' | sort -r | uniq | grep gene_alias | cut -f 2 -d'=' | cut -f 1 -d'.' |   sort | uniq | tr '\r\n' ' ' | sed 's/ /;gene_alias=/g' | awk '{print "gene_alias="$0}' | rev | cut -c 13- | rev`
  mkdir -p ${db}_genes/$chr
  echo -e "$chr\t$db\tgene\t$st\t$end\t.\t$strand\t.\tgene_id \"$gene_id\";$gene_alias" > ${db}_genes/$chr/${chr}_${st}_$end
  
  cat chr_bp_genename_transcrips_all | grep -w $i |  while read a b c; do 
    cat $gtf |  sed -n "$a,${b}p" >> ${db}_genes/$chr/${chr}_${st}_$end
  done
done
}

# loop it
for i in {0..1}; do
   let m=($i*5000)+1 ; let n=($i+1)*5000;   run_it $m $n &
done
run_it 10001 14090 &
wait

for i in `ls ${db}_genes ` ; do ls ${db}_genes/$i/* | cut -f 3 -d'/' >> gene_bp_all; done
# 14087 #3 less than gene_names. were they duplicates? :

cat gene_names  | while read i; do
  st=`cat $gtf | grep -n "\"$i\"" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "\"$i\"" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep  "\"$i\"" | cut -f 1 | head -1`
  strand=`cat $gtf | grep  "\"$i\"" | cut -f 7 | head -1`
  echo $i
  echo -e "$chr\t$st\t$end\t$strand\t$i" >> chr_bp_gene
done
cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1>1) print $2"\t"$3"\t"$4}'  > duplicates_chr_bp
cat duplicates_chr_bp | sort -k1,1 -k2,2n > tmp1
cat tmp1
chr1	103414879	103425465
chr1	143877731	143885076
chr9	41092123	41095600

for i in `cat tmp1 | cut -f 2`; do grep $i chr_bp_gene; done
chr1	103414879	103425465	+	BIG-gene-chr1-1911
chr1	103414879	103425465	+	BIG-gene-chr1-1913
chr1	143877731	143885076	-	BIG-gene-chr1-2200
chr1	143877731	143885076	-	BIG-gene-chr1-2427
chr9	41092123	41095600	-	BIG-lncRNA-1006
chr9	41092123	41095600	-	BIG-lncRNA-1008

#exclude BIG-gene-chr1-1913, BIG-gene-chr1-2427, BIG-lncRNA-1008
cat chr_bp_gene | grep -v BIG-gene-chr1-1913 | grep -v  BIG-gene-chr1-2427 | grep -v  BIG-lncRNA-1008 > chr_bp_gene_dedup
wc -l chr_bp_gene*
  14090 chr_bp_gene
  14087 chr_bp_gene_dedup

################################## 

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
for i in {0..1}; do
   let m=($i*5000)+1 ; let n=($i+1)*5000;   run_it $m $n &
done
run_it 10001 14087 &
wait


###################################

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
                           | sed 's/;transcript_id=/\t/g'  | cut -f 1,3 | sed 's/;/\t/g' | cut -f 1,2 | sed 's/\t/_/g' | sed 's/:/_/g' > tmp_tname_$c
   paste tmp11_$c tmp22_$c tmp_tname_$c > tmp_tall_$c
   mkdir -p  ${db}_trans/$chr/$c
   cat tmp_tall_$c | while read d e f; do
     cat ${db}_genes/$chr/$c | sed -n "$d,${e}p" > ${db}_trans/$chr/$c/$f
   done
   rm -f tmp11_$c tmp22_$c tmp_tname_$c tmp_tall_$c
done
}

# loop it
for i in {0..1}; do
   let m=($i*5000)+1 ; let n=($i+1)*5000;   run_it $m $n &
done
run_it 10001 14087 &
wait


############################ 

mkdir -p ${db}_trans_attr
run_it(){
m=$1
n=$2
cat chr_bp_gene_dedup | sed -n "$m,${n}p"  | while read chr b c1 e1 f1 ; do
   mkdir -p ${db}_trans_attr/$chr
   c="${chr}_${b}_${c1}"
   for i in `ls ${db}_trans/$chr/$c`; do
     mkdir -p ${db}_trans_attr/$chr/$c/
     cat ${db}_trans/$chr/$c/$i | head -1 | cut -f 9 | sed 's/; /\n/g' | sed 's/ /\t/g' > ${db}_trans_attr/$chr/$c/$i
   done
done
}

# loop it
for i in {0..1}; do
   let m=($i*5000)+1 ; let n=($i+1)*5000;   run_it $m $n &
done
run_it 10001 14087 &
wait

############## 


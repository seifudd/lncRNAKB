
# https://lncipedia.org/download

#wget  https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc_hg38.gtf
#wget  https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc_hg38.gff
wget  https://lncipedia.org/downloads/lncipedia_5_2/full-database/lncipedia_5_2_hg38.gtf
wget  https://lncipedia.org/downloads/lncipedia_5_2/full-database/lncipedia_5_2_hg38.gff

Go with all
extract geneID > get the start and end > thats the record

gtf="lncipedia_5_2_hg38.gff"
db="lncipedia"

#### transcripts row number
echo "2" > tmp1
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="lnc_RNA") print NR+1 }' >> tmp1
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="lnc_RNA") print NR }' > tmp2
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="lnc_RNA") print }' | cut -f 9 | sed 's/;gene_id=/\t/g' | cut -f 2 | sed 's/;/\t/g' | cut -f 1  > tmp3
paste tmp1 tmp2 tmp3 > chr_bp_genename_transcrips_all


# now for each gene:
cat $gtf | cut -f 9 | sed 's/;gene_id=/\t/g' | cut -f 2 | sed 's/;/\t/g' | cut -f 1 | uniq | sort | uniq > gene_names

run_it(){
m=$1
n=$2
cat gene_names | sed -n "$m,${n}p"   | while read i; do
  st=`cat $gtf | grep -n "$i;" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "$i;" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep  "$i;" | cut -f 1 | head -1`
  strand=`cat $gtf | grep  "$i;" | cut -f 7 | head -1`
  gene_id=`cat $gtf | grep -n "$i;" | cut -f 9 | sed 's/;/\n/g;' | sort -r | uniq | grep gene_id | cut -f 2 -d'='`
  gene_alias=`cat $gtf | grep -n "$i;" | cut -f 9 | sed 's/;/\n/g;' | sort -r | uniq | grep gene_alias | cut -f 2 -d'=' | cut -f 1 -d'.' |   sort | uniq | tr '\r\n' ' ' | sed 's/ /;gene_alias=/g' | awk '{print "gene_alias="$0}' | rev | cut -c 13- | rev`
  mkdir -p ${db}_genes/$chr
  echo -e "$chr\t$db\tgene\t$st\t$end\t.\t$strand\t.\tGENE_NAME=$gene_id;$gene_alias" > ${db}_genes/$chr/${chr}_${st}_$end
  
  cat chr_bp_genename_transcrips_all | grep -w $i |  while read a b c; do 
    cat $gtf | head -$b | tail -1 >> ${db}_genes/$chr/${chr}_${st}_$end
    let d=$b-1
    cat $gtf |  sed -n "$a,${d}p" >> ${db}_genes/$chr/${chr}_${st}_$end
  done
done
}

# loop it
for i in {0..14}; do
   let m=($i*4000)+1 ; let n=($i+1)*4000;   run_it $m $n &
done
run_it 56001 56947 &
wait

for i in {1..22} X Y ; do ls lncipedia_genes/chr$i/* | cut -f 3 -d'/' >> gene_bp_all; done
# 56585 #362 less than gene_names. were they duplicates? :
gtf="lncipedia_5_2_hg38.gff"
cat gene_names | while read i; do
  st=`cat $gtf | grep -n "$i;" | cut -f 4 | sort -n | head -1`
  end=`cat $gtf | grep -n "$i;" | cut -f 5 | sort -n | tail -1`
  chr=`cat $gtf | grep  "$i;" | cut -f 1 | head -1`
  strand=`cat $gtf | grep  "$i;" | cut -f 7 | head -1`
  echo $i
  echo -e "$chr\t$st\t$end\t$strand\t$i" >> chr_bp_gene
done
cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1>1) print}'
cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1>1) print $2"\t"$3"\t"$4}'  > duplicates_chr_bp

cat chr_bp_gene | cut -f -3 | sort | uniq -c | awk '{if($1<2) print $2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n > tmp1
cat tmp1 | while read a b c; do 
   cat chr_bp_gene | grep -P "$1\t$2\t$3" >> chr_bp_gene_dedup
done
wc -l chr_bp_gene*
  56947 chr_bp_gene
  56224 chr_bp_gene_dedup

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
for i in {0..14}; do
   let m=($i*4000)+1 ; let n=($i+1)*4000;   run_it $m $n &
done
run_it 56001 56224 &
wait


###################################
mkdir -p ${db}_trans
run_it(){
m=$1
n=$2
cat chr_bp_gene_dedup | sed -n "$m,${n}p"   | while read chr b c1 e1 f1 ; do
  mkdir -p ${db}_trans/$chr
   c="${chr}_${b}_${c1}"
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="lnc_RNA") print NR }' > tmp11_$c
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="lnc_RNA") print NR-1 }' | sed '1,1d' > tmp22_$c
   cat ${db}_genes/$chr/$c | wc -l >> tmp22_$c
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="lnc_RNA") print $1"_"$4"_"$5"\t"$9 }' \
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
for i in {0..14}; do
   let m=($i*4000)+1 ; let n=($i+1)*4000;   run_it $m $n &
done
run_it 56001 56224 &
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
     cat ${db}_trans/$chr/$c/$i | head -1 | cut -f 9 | sed 's/;/\n/g' | sed 's/=/\t/g' > ${db}_trans_attr/$chr/$c/$i
   done
done
}

# loop it
for i in {0..14}; do
   let m=($i*4000)+1 ; let n=($i+1)*4000;   run_it $m $n &
done
run_it 56001 56224 &
wait


############## 


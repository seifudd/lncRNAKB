
gtf="FANTOM_CAT.lv3_robust.only_lncRNA.hg38.gtf"
db="fantom"

cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="gene") print NR }' > tmp1
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="gene") print NR-1 }' | sed '1,1d' > tmp2
cat $gtf | wc -l >> tmp2
cat $gtf | awk -v FS='\t' -v OFS='\t'  '{if($3=="gene") print $1"_"$4"_"$5 }' > tmp_name
paste tmp* > chr_bp_gene_all

cat chr_bp_gene_all | cut -f 3 | wc
cat chr_bp_gene_all | cut -f 3 | sort | uniq | wc
cat chr_bp_gene_all | cut -f 3 | sort | uniq -c | sort | tail
"
chr10_75398581_75412482 :  removed rm 365394	365406	chr10_75398581_75412482
chr12_53991834_54000698: rm 288555	288573	chr12_53991834_54000698
chr19_54424186_54449430: rm 197121	197131	chr19_54424186_54449430
chr2_5601735_5692160: rm 178751	178770	chr2_5601735_5692160
chr5_50964893_50970255: rm 288171	288180	chr5_50964893_50970255
"


mkdir -p ${db}_genes
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p ${db}_genes/$chr
  cat $gtf | sed -n "$a,${b}p" > ${db}_genes/$chr/$c
done
}

# loop it
for i in {0..14}; do
   let m=($i*3000)+1 ; let n=($i+1)*3000;   run_it $m $n &
done
run_it 45001 46421 &
wait

mkdir -p ${db}_genes_attr
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p ${db}_genes_attr/$chr
  cat ${db}_genes/$chr/$c | head -1 | cut -f 9 | sed 's/; /\n/g' | sed 's/ /\t/g' > ${db}_genes_attr/$chr/$c
done
}

# loop it
for i in {0..14}; do
   let m=($i*3000)+1 ; let n=($i+1)*3000;   run_it $m $n &
done
run_it 45001 46421 &
wait


# chr22	FANTOM	transcript	16601866	16630974	.	+	.	gene_id "ENSG00000100181.17"; transcript_id "ENCT00000276244.1"; 
# chr22_16601866_16630974	gene_id "ENSG00000100181.17"	ENCT00000276244.1 
mkdir -p ${db}_trans
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p ${db}_trans/$chr
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR }' > tmp11_$c
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR-1 }' | sed '1,1d' > tmp22_$c
   cat ${db}_genes/$chr/$c | wc -l >> tmp22_$c
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print $1"_"$4"_"$5"\t"$9 }' | cut -f -2 -d';' | sed 's/; transcript_id "/\t/g'  | sed 's/"//g' | cut -f 1,3 |  sed 's/\t/_/g'  > tmp_tname_$c
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
   let m=($i*3000)+1 ; let n=($i+1)*3000;   run_it $m $n &
done
run_it 45001 46421 &
wait

# gene_id "ENSG00000100181.17"; transcript_id "ENCT00000276244.1"; transcript_type "__na"; transcript_name "__na"; coding_status "nonCoding"; cumulative_support "ENCODE:FANTOM"; TIEScore "82.40089";

mkdir -p ${db}_trans_attr
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p ${db}_trans_attr/$chr
   cat ${db}_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 > tmp_trans_$c
   mkdir -p  ${db}_trans_attr/$chr/$c
   cat tmp_trans_$c | while read line; do 
      ID=`echo -e "$line" | sed 's/;/\n/g' | sed 's/ /\t/g' | sed 's/"//g'   | awk '{print $1"\t"$2}' | grep -w ^transcript_id | cut -f 2 | head -1 `
      echo -e "$line" | sed 's/;/\n/g' | sed 's/ /\t/g' | sed 's/"//g'   | awk '{print $1"\t"$2}'  > ${db}_trans_attr/$chr/$c/${c}_$ID
  done
   rm -f tmp_trans_$c
done
}

# loop it
for i in {0..14}; do
   let m=($i*3000)+1 ; let n=($i+1)*3000;   run_it $m $n &
done
run_it 45001 46421 &
wait


############## 

 cat FANTOM_CAT.lv3_robust.only_lncRNA.hg38.gtf | awk '{if($3=="gene") print $1"\t"$4"\t"$5"\t"$7"\t"$10}' \
   | sed 's/;/\t/g' | sed 's/"//g' | cut -f -5 > chr_bp_gene_dedup





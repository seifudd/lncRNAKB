
ml bedtools

mkdir -p 1_0

F0="../00-CHESS2.1/chr_bp_gene_dedup"
F1="../01-FANTOM5.0.v3/chr_bp_gene_dedup"

# 1) if strand is the same , compare transcript
run_pairwise(){
a=$1 ; b=$2 ; c=$3 ; d=$4 
bedtools intersect -nonamecheck    -a $a -b $b -wa -wb |awk '{if($4==$9) print}' > ${c}_$d/${c}_${d}.overlap
bedtools subtract  -nonamecheck -A -a $a -b $b         > ${c}_$d/${c}_${d}.subtract
}
run_pairwise $F1 $F0 1 0

# chr10	69066079	69066488	+	CATG00000115504.1	chr10	69056418	69106491	+	CHS.6371
# chr10	69086819	69089052	-	CATG00000000008.1	chr10	69088137	69089052	-	CHS.6372
# chr22	29509307	29510167	-	CATG00000058414.1	chr22	29508167	29554254	-	CHS.35544
# chr22	29538361	29538799	-	CATG00000058415.1	chr22	29508167	29554254	-	CHS.35544

#chr=chr10 
#st=69086819
#en=69089052
#strand="-"
#g1=CATG00000000008.1
#chr2=chr10 
#st2=69088137
#en2=69089052
#strand2="-"
#g2=CHS.6372

chr2=chr22
st2=29508167
en2=29554254

run_loop(){

chr_loop=$1

c=1 ; d=0 
cat  ${c}_$d/${c}_${d}.overlap | grep -P "^$chr_loop\t"  | while read chr st en strand g1 chr2 st2 en2 strand2 g2; do

# 2) get the gene_id
gene_id=`cat chess_genes_GTF/$chr2/${chr2}_${st2}_${en2} | head -1 | cut -f 9 | sed 's/;/\n/g' | grep gene_id | cut -f 2 -d' ' | sed 's/"//g'`

# 3) get the last transcript num (id) and add 1 # new transcript is $gene_id.$a
last_trans_num=`cat chess_genes_GTF/$chr2/${chr2}_${st2}_${en2}    1_0/transcript/$chr2/${chr2}_${st2}_${en2}/*  2>/dev/null \
                   | awk '{if($3=="transcript") print}' | cut -f 9 \
                   | sed 's/;/\n/g' | sed "s/^[ \t]*//" | grep -w ^transcript_id  | cut -f 2 -d' ' | sed 's/"//g' | cut -f 3 -d'.' | sort -n | tail -1`
let last_trans_num_plus=$last_trans_num+1
echo "$last_trans_num_plus"  > 1_${chr2}_${st2}_${en2}_last_trans_num_plus

# 4) get transcript bed
cat chess_genes_GTF/$chr2/${chr2}_${st2}_${en2} | awk '{if($3=="transcript") print}' | cut -f 1,4,5  > 1_${chr2}_${st2}_${en2} 

# 5) for each transcript in ../01-FANTOM5.0.v3/fantom_trans/$chr/${chr}_${st}_${en}/* check if it is in boundry 
# 6) for each transcript in  if F1 record doesn't match append it 

mkdir -p 1_0/transcript/$chr2/${chr2}_${st2}_${en2}/
cp chess_transcript_GTF/${chr2}/${chr2}_${st2}_${en2}/*   1_0/transcript/$chr2/${chr2}_${st2}_${en2}/
for i in `ls ../01-FANTOM5.0.v3/fantom_trans/$chr/${chr}_${st}_${en}/*`; do
  cat $i | head -1 | cut -f 1,4,5 | while read a b c; do
     last_trans_num_plus=`cat 1_${chr2}_${st2}_${en2}_last_trans_num_plus`
     #echo "$a $b $c $st2 $en2"
     if [ $b -ge $st2 ] && [ $c -le $en2 ]; then
     m1=`cat 1_${chr2}_${st2}_${en2}  |grep -P "${a}\t${b}\t${c}" | wc -l`
     if [ $m1 -lt 1 ]; then
        cat $i | sed 's/gene_id /gene_id_alias /g' | sed 's/transcript_id /transcript_id_alias /g' \
           | awk -v FS='\t' -v OFS='\t' -v gene_id=$gene_id -v last_trans_num_plus=$last_trans_num_plus '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tgene_id \""gene_id"\"; transcript_id \""gene_id"."last_trans_num_plus"\"; "$9}' \
           >  1_0/transcript/$chr2/${chr2}_${st2}_${en2}/${a}_${b}_${c}_${gene_id}.$last_trans_num_plus 
     let last_trans_num_plus=$last_trans_num_plus+1
     echo "$last_trans_num_plus" > 1_${chr2}_${st2}_${en2}_last_trans_num_plus
     fi; fi
  done
done
rm -f 1_${chr2}_${st2}_${en2} 
rm -f 1_${chr2}_${st2}_${en2}_last_trans_num_plus
done
}


for i in `cat 1_0/1_0.overlap | cut -f 1  | sort | uniq `; do
 run_loop $i &
done
wait

#test:
cd /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/09-merge-DBs/1_0/transcript/chr1
ls * | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq | wc
ls * | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq -c | sort  |tail
cd -

#  gene
ls 1_0/transcript/chr*  | grep -v ':' | awk 'NF' | wc -l
ls 1_0/transcript/chr*  | grep -v ':' | awk 'NF' > 1_0/genes_modified #14300
m=`cat 1_0/genes_modified | wc -l`; n=1
for i in `cat 1_0/genes_modified`; do
 chr=`echo "$i"  | rev | cut -f 3- -d'_' | rev`
 mkdir -p 1_0/gene/$chr
 cat chess_genes_GTF/$chr/$i | head -1 > 1_0/gene/$chr/$i
 cat 1_0/transcript/$chr/$i/* >> 1_0/gene/$chr/$i
 echo -en " of $m ";  echo -en "\r .. $n ";  let n=$n+1
done
echo -en "\r $m genes processed.\n"


########################## subtract add
 

# 1) get the last gene_id 
cat ../00-CHESS2.1/chr_bp_gene_dedup | cut -f 5 | cut -f 2 -d'.' | sort -n  | tail -1 > 1_0/last_gene_id_num  #59970

# 2) for each transcript add .1 .2 .3 etc 
# per gene
m=`cat  ${c}_$d/${c}_${d}.subtract | wc -l`
c=1 ; d=0 ; n=0 ; cat  ${c}_$d/${c}_${d}.subtract |   while read chr st en strand g1 ; do
  last_gene_id_num=`cat 1_0/last_gene_id_num`
  let last_gene_id_num_plus=$last_gene_id_num+1   
  trans_num=1
  for i in `ls ../01-FANTOM5.0.v3/fantom_trans/$chr/${chr}_${st}_${en}/*`; do
    trans_chr_bp=`echo  $i | rev | cut -f 1 -d'/' | cut -f 2- -d'_' | rev`  #chr1_46617031_46617973
    trans_new_id="$last_gene_id_num_plus.$trans_num"  #59971.1  (need CHS.59971.1)
    mkdir -p 1_0/trans_subtract/$chr/${chr}_${st}_${en}
    cat $i | sed 's/gene_id /gene_id_alias /g' | sed 's/transcript_id /transcript_id_alias /g' \
           | awk -v FS='\t' -v OFS='\t' -v gene_id=$last_gene_id_num_plus -v trans_num=$trans_new_id  \
             '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tgene_id \"CHS."gene_id"\"; transcript_id \"CHS."trans_num"\"; "$9}' \
             >  1_0/trans_subtract/$chr/${chr}_${st}_${en}/${trans_chr_bp}_CHS.$trans_new_id
    let trans_num=$trans_num+1  #2
  done
  echo "$last_gene_id_num_plus" > 1_0/last_gene_id_num
  echo -en " of $m ";  echo -en "\r .. $n ";  let n=$n+1
done
echo -en "\r $m genes processed.\n"

# test and count
cd 1_0/trans_subtract
ls * | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq | wc # 7157  gene
ls * | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq -c | sort  |tail
ls chr*/* | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq  | wc  #11583 trans
ls chr*/* | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq -c | sort  |tail
cd -










































wdir=6_54310
wdir_prev1=1_0
wdir_prev2=3_10
wdir_prev3=4_310
wdir_prev4=5_4310

mkdir -p $wdir

# F4310 = 4_310/310.chr_bp_gene_dedup   +   4_310/gene_subtract.chr_bp_gene_dedup  (((46421+7157)53578  +  10506 )  + 20700) + 15164  = 99948
for i in `ls $wdir_prev4/gene_subtract/chr*/chr* `; do 
   cat $i | head -1 | cut -f 1,4,5,7,9 | cut -f 1 -d';' | sed 's/gene_id "//g' | sed 's/"//g' >> $wdir_prev4/gene_subtract.chr_bp_gene_dedup
done
cat $wdir_prev4/gene_subtract.chr_bp_gene_dedup | wc -l  #15164 
cat  $wdir_prev4/4310.chr_bp_gene_dedup   $wdir_prev4/gene_subtract.chr_bp_gene_dedup   > $wdir/54310.chr_bp_gene_dedup  # 99948
F54310="$wdir/54310.chr_bp_gene_dedup"
F6="../06-BIGTRANSCRIPTOMEv1/chr_bp_gene_dedup"

# 1) if strand is the same , compare transcript
run_pairwise(){
ml bedtools
a=$1 ; b=$2 ; c=$3 ; d=$4 
bedtools intersect -nonamecheck    -a $a -b $b -wa -wb |awk '{if($4==$9) print}' > $wdir/$wdir.overlap
bedtools subtract  -nonamecheck -A -a $a -b $b         > $wdir/$wdir.subtract
}
run_pairwise $F6 $F54310 6 54310


run_loop(){

chr_loop=$1

wdir=6_54310
wdir_prev1=1_0
wdir_prev2=3_10
wdir_prev3=4_310
wdir_prev4=5_4310


cat  $wdir/$wdir.overlap | grep -P "^$chr_loop\t"  | while read chr st en strand g1 chr2 st2 en2 strand2 g2; do

# 2) get the gene_id
gene_id=$g2

# 3) get the last transcript num (id) and add 1 # new transcript is $gene_id.$a
last_trans_num=`cat chess_genes_GTF/$chr2/${chr2}_${st2}_${en2}   \
                    $wdir_prev1/gene/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir_prev1/gene_subtract/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir_prev2/gene/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir_prev2/gene_subtract/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir_prev3/gene/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir_prev3/gene_subtract/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir_prev4/gene/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir_prev4/gene_subtract/$chr2/${chr2}_${st2}_${en2}  \
                    $wdir/transcript/$chr2/${chr2}_${st2}_${en2}/*  2>/dev/null \
                   | awk '{if($3=="transcript") print}' | cut -f 9 \
                   | sed 's/;/\n/g' | sed "s/^[ \t]*//" | grep -w ^transcript_id  | cut -f 2 -d' ' | sed 's/"//g' | cut -f 3 -d'.' | sort -n | tail -1`
let last_trans_num_plus=$last_trans_num+1
echo "$last_trans_num_plus" > 1_${chr2}_${st2}_${en2}_last_trans_num_plus
# 4) get transcript bed
if [ -e $wdir_prev4/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev4/gene_subtract/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev4/trans_subtract/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev4/gene/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev4/gene/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev4/transcript/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev3/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev3/gene_subtract/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev3/trans_subtract/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev3/gene/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev3/gene/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev3/transcript/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev2/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev2/gene_subtract/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev2/trans_subtract/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev2/gene/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev2/gene/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev2/transcript/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev1/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev1/gene_subtract/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev1/trans_subtract/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev1/gene/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev1/gene/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev1/transcript/$chr2/${chr2}_${st2}_${en2}"
else
    gene_file="chess_genes_GTF/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="chess_transcript_GTF/$chr2/${chr2}_${st2}_${en2}"
fi
cat  $gene_file | awk '{if($3=="transcript") print}' | cut -f 1,4,5  > 1_${chr2}_${st2}_${en2} 

# 5) for each transcript in this gene (../06-BIGTRANSCRIPTOMEv1/BIGTranscriptome_trans/$chr/${chr}_${st}_${en}/*) check if it is in boundry 
# 6) for each transcript in  if F3 record doesn't match append it 

mkdir -p $wdir/transcript/$chr2/${chr2}_${st2}_${en2}/
cp $trans_dir/*   $wdir/transcript/$chr2/${chr2}_${st2}_${en2}/
for i in `ls ../06-BIGTRANSCRIPTOMEv1/BIGTranscriptome_trans/$chr/${chr}_${st}_${en}/*`; do
  cat $i | head -1 | cut -f 1,4,5 | while read a b c; do
     last_trans_num_plus=`cat 1_${chr2}_${st2}_${en2}_last_trans_num_plus`
     #echo "$a $b $c $st2 $en2  $last_trans_num_plus"
     if [ $b -ge $st2 ] && [ $c -le $en2 ]; then
     m1=`cat 1_${chr2}_${st2}_${en2}  |grep -P "${a}\t${b}\t${c}" | wc -l`
     if [ $m1 -lt 1 ]; then
        cat $i | sed 's/gene_id /gene_id_alias /g' | sed 's/transcript_id /transcript_id_alias /g' \
           | awk -v FS='\t' -v OFS='\t' -v gene_id=$gene_id -v last_trans_num_plus=$last_trans_num_plus '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tgene_id \""gene_id"\"; transcript_id \""gene_id"."last_trans_num_plus"\"; "$9}' \
            >  $wdir/transcript/$chr2/${chr2}_${st2}_${en2}/${a}_${b}_${c}_${gene_id}.$last_trans_num_plus 
     let last_trans_num_plus=$last_trans_num_plus+1
     echo "$last_trans_num_plus" > 1_${chr2}_${st2}_${en2}_last_trans_num_plus
     fi; fi
  done
done
rm -f 1_${chr2}_${st2}_${en2} 
rm -f 1_${chr2}_${st2}_${en2}_last_trans_num_plus
done
}

for i in `cat $wdir/$wdir.overlap | cut -f 1 | grep -v chr22  | sort | uniq `; do
 run_loop $i &
done
wait

#test:
cd 6_54310/transcript/chr22
ls * | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq | wc
ls * | grep -v ':' | awk 'NF' | cut -f 2- -d'.' | sort | uniq -c | sort  |tail
cd -


#  gene
ls $wdir/transcript/chr*  | grep -v ':' | awk 'NF' > $wdir/genes_modified # 13192
m=`cat $wdir/genes_modified | wc -l`; n=1
for i in `cat $wdir/genes_modified`; do
 chr=`echo "$i"  | rev | cut -f 3- -d'_' | rev`
 mkdir -p $wdir/gene/$chr
if [ -e $wdir_prev4/gene_subtract/$chr/$i ]
then
    gene_file="$wdir_prev4/gene_subtract/$chr/$i"
    trans_dir="$wdir_prev4/trans_subtract/$chr/$i"
elif [ -e $wdir_prev4/gene/$chr/$i ]
then
    gene_file="$wdir_prev4/gene/$chr/$i"
    trans_dir="$wdir_prev4/transcript/$chr/$i"
elif [ -e $wdir_prev3/gene_subtract/$chr/$i ]
then
    gene_file="$wdir_prev3/gene_subtract/$chr/$i"
    trans_dir="$wdir_prev3/trans_subtract/$chr/$i"
elif [ -e $wdir_prev3/gene/$chr/$i ]
then
    gene_file="$wdir_prev3/gene/$chr/$i"
    trans_dir="$wdir_prev3/transcript/$chr/$i"
elif [ -e $wdir_prev2/gene_subtract/$chr/$i ]
then
    gene_file="$wdir_prev2/gene_subtract/$chr/$i"
    trans_dir="$wdir_prev2/trans_subtract/$chr/$i"
elif [ -e $wdir_prev2/gene/$chr/$i ]
then
    gene_file="$wdir_prev2/gene/$chr/$i"
    trans_dir="$wdir_prev2/transcript/$chr/$i"
elif [ -e $wdir_prev1/gene_subtract/$chr/$i ]
then
    gene_file="$wdir_prev1/gene_subtract/$chr/$i"
    trans_dir="$wdir_prev1/trans_subtract/$chr/$i"
elif [ -e $wdir_prev1/gene/$chr/$i ]
then
    gene_file="$wdir_prev1/gene/$chr/$i"
    trans_dir="$wdir_prev1/transcript/$chr/$i"
else
    gene_file="chess_genes_GTF/$chr/$i"
    trans_dir="chess_transcript_GTF/$chr/$i"
fi
 cat $gene_file | head -1 > $wdir/gene/$chr/$i
 cat $wdir/transcript/$chr/$i/* >> $wdir/gene/$chr/$i
 echo -en " of $m ";  echo -en "\r .. $n ";  let n=$n+1
done
echo -en "\r $m genes processed.\n"


########################## subtract add # 20700 

# 1) get the last gene_id 
cat $wdir/54310.chr_bp_gene_dedup | cut -f 5 | cut -f 2 -d'.' | sort -n  | tail -1
cat $wdir/54310.chr_bp_gene_dedup | cut -f 5 | cut -f 2 -d'.' | sort -n  | tail -1 > $wdir/last_gene_id_num  #113497 last number

# 2) for each transcript add .1 .2 .3 etc 
# per gene
m=`cat  $wdir/$wdir.subtract | wc -l`  # 333

c=1 ; d=0 ; n=0 ; cat  $wdir/$wdir.subtract  |   while read chr st en strand g1 ; do
  last_gene_id_num=`cat $wdir/last_gene_id_num`
  let last_gene_id_num_plus=$last_gene_id_num+1   
  trans_num=1
  for i in `ls ../06-BIGTRANSCRIPTOMEv1/BIGTranscriptome_trans/$chr/${chr}_${st}_${en}/*`; do
    trans_chr_bp=`echo  $i | rev | cut -f 1 -d'/' | rev | cut -f -3 -d'_'`  #chr1_46617031_46617973
    trans_new_id="$last_gene_id_num_plus.$trans_num"  #67127.1  (need CHS.59971.1)
    mkdir -p $wdir/trans_subtract/$chr/${chr}_${st}_${en}
    cat $i | sed 's/gene_id /gene_id_alias /g' | sed 's/transcript_id /transcript_id_alias /g' \
           | awk -v FS='\t' -v OFS='\t' -v gene_id=$last_gene_id_num_plus -v trans_num=$trans_new_id  \
             '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tgene_id \"CHS."gene_id"\"; transcript_id \"CHS."trans_num"\"; "$9}' \
              >  $wdir/trans_subtract/$chr/${chr}_${st}_${en}/${trans_chr_bp}_CHS.$trans_new_id
    let trans_num=$trans_num+1  #2
  done
  mkdir -p $wdir/gene_subtract/$chr
  cat ../06-BIGTRANSCRIPTOMEv1/BIGTranscriptome_genes/$chr/${chr}_${st}_${en} | head -1 | sed 's/GENE_NAME=/gene_alias "/g'  \
           | sed 's/;gene_alias=/;gene_alias "/g' | sed 's/;/"; /g' | sed 's/""/"/g'  \
           | awk -v FS='\t' -v OFS='\t' -v gene_id=$last_gene_id_num_plus  \
             '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tgene_id \"CHS."gene_id"\"; "$9"\";"}' \
             >  $wdir/gene_subtract/$chr/${chr}_${st}_${en}
  cat $wdir/trans_subtract/$chr/${chr}_${st}_${en}/* >> $wdir/gene_subtract/$chr/${chr}_${st}_${en}
  echo "$last_gene_id_num_plus" > $wdir/last_gene_id_num
  echo -en " of $m ";  echo -en "\r .. $n ";  let n=$n+1
done
echo -en "\r $m genes processed.\n"


############### clean up the gene name
for i in `ls 6_54310/gen*/*/*`; do
 cat $i | head -1 | sed 's/; gene_id/; gene_id_alis/g' | sed 's/"; ";/";/g' > tmp1
 cat $i | sed '1,1d' >> tmp1
 cat tmp1 > $i
done






































































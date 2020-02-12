

wdir=7_654310
wdir_prev1=1_0
wdir_prev2=3_10
wdir_prev3=4_310
wdir_prev4=5_4310
wdir_prev5=6_54310

mkdir -p $wdir
# F4310 = 4_310/310.chr_bp_gene_dedup   +   4_310/gene_subtract.chr_bp_gene_dedup  (((46421+7157)53578  +  10506 )  + 20700) + 15164 + 333 = 100281
for i in `ls $wdir_prev5/gene_subtract/chr*/chr* `; do 
   cat $i | head -1 | cut -f 1,4,5,7,9 | cut -f 1 -d';' | sed 's/gene_id "//g' | sed 's/"//g' >> $wdir_prev5/gene_subtract.chr_bp_gene_dedup
done
cat $wdir_prev5/gene_subtract.chr_bp_gene_dedup | wc -l  #333 
cat  $wdir_prev5/54310.chr_bp_gene_dedup   $wdir_prev5/gene_subtract.chr_bp_gene_dedup   > $wdir/654310.chr_bp_gene_dedup  #100281

# 3 dupicates coming from chess  
chr6_GL000250v2_alt	2862768	2871733
chr6_GL000253v2_alt	3044995	3052474
chr6_GL000255v2_alt	2995767	3017449


echo -e "##description: evidence-based annotation of the human genome (GRCh38), version 1.0 
##provider: lncRNAKB
##contact: fayaz.seifuddin@nih.gov
##format: gtf
##date: 2018-12-31" > lncRNAKB_hg38.gtf



cat 7_654310/654310.chr_bp_gene_dedup | while read chr2 st2 en2 strand gene; do

if [ -e $wdir_prev5/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev5/gene_subtract/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev5/trans_subtract/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev5/gene/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev5/gene/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev5/transcript/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev4/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
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
cat  $gene_file >> lncRNAKB_hg38.gtf


done 2>error.log


cat  lncRNAKB_hg38.gtf | awk '{if($3=="gene") print $2}'  | wc
cat  lncRNAKB_hg38.gtf | awk '{if($3=="transcript") print $2}'  | wc
cat  lncRNAKB_hg38.gtf | awk '{if($3=="transcript") print $2}'  | sort | uniq 
cat  lncRNAKB_hg38.gtf | awk '{if($3=="transcript") print $2}'  | grep -w FANTOM  | wc -l
cat  lncRNAKB_hg38.gtf | awk '{if($3=="transcript") print $2}'  | grep -w lncipedia.org  | wc -l
cat  lncRNAKB_hg38.gtf | awk '{if($3=="transcript") print $2}'  | grep -w mitranscriptome  | wc -l
cat  lncRNAKB_hg38.gtf | awk '{if($3=="transcript") print $2}'  | grep -w CAFE  | wc -l

Total:
Gene numbers: 100281 
Transcripts: 532420 

Some DB counts:
Fantom           : gene 7157  | transcript 76400 
lncipedia        : gene 10506 | transcript 57061
mitranscriptome  : gene 15164 | transcript 49489 
BIGtranscriptome : gene 333   | transcript 6780 
NONCODE          : gene 20701 | transcript 93453

cat  lncRNAKB_hg38.gtf | awk '{if($3=="transcript") print}'  | grep  '"NON'  | wc -l
93453
cat  lncRNAKB_hg38.gtf | awk '{if($3=="gene") print}'  | grep  '"NON'  | wc -l
20701

################### Gene : clean up the bene name records and get CHS ids for each DBs

run_it(){
d="$1"
ls $d/gene/*/* > tmp$d.list ; ls $d/gene_subtract/*/* >> tmp$d.list
for i in `cat tmp$d.list`; do cat $i | head -1 | sed 's/gene_id "/\n/g' | sed '1,1d' | cut -f 1 -d'"' >> tmp2$d.list ; done 
cat tmp2$d.list | sort | uniq > CHS.$d 
rm -f tmp$d.list tmp2$d.list
}

# 00-CHESS2.1       #46421
cat ../00-CHESS2.1/chr_bp_gene_dedup | cut -f 5 > CHS.CHESS 

# 01-FANTOM5.0.v3  # 21457
run_it 1_0 &
# 03-Lncipedia-5.2 # 40188
run_it 3_10 &
# 04-NONCODEv5    # redo
run_it 4_310 &
# 05-MITRANSCRIPTOMEv2
run_it 5_4310 &
# 06-BIGTRANSCRIPTOMEv1
run_it 6_54310 &

ls $d/gene/*/* > tmp$d.list ; ls $d/gene_subtract/*/* >> tmp$d.list
for i in `cat tmp$d.list`; do cat $i | head -1 | sed 's/gene_id "/\n/g' | sed '1,1d' | cut -f 1 -d'"' >> tmp2$d.list ; done 
cat tmp2$d.list | sort | uniq > CHS.$d 

for i in `ls 6_54310/gen*/*/*`; do
 cat $i | head -1 | sed 's/; gene_id/; gene_id_alis/g' | sed 's/"; ";/";/g' > tmp1
 cat $i | sed '1,1d' >> tmp1
 cat tmp1 > $i
done

ls 5_4310/gene_subtract/*/* > tmp20
ls 5_4310/gene/*/* >> tmp20
m=`cat tmp20 | wc -l`; n=0
for i in `cat tmp20`; do
 cat $i | head -1 | sed 's/; gene_id/; gene_id_alis/g' | sed 's/"; ";/";/g' > tmp2
 cat $i | sed '1,1d' >> tmp2
 cat tmp2 > $i
  echo -en " of $m ";  echo -en "\r .. $n ";  let n=$n+1
done
echo -en "\r $m genes processed.\n"


ls 4_310/gene_subtract/*/* > tmp30
ls 4_310/gene/*/* >> tmp30
m=`cat tmp30 | wc -l`; n=0
for i in `cat tmp30`; do
 cat $i | head -1 | sed 's/; gene_id/; gene_id_alis/g' | sed 's/"; ";/";/g' > tmp3
 cat $i | sed '1,1d' >> tmp3
 cat tmp3 > $i
  echo -en " of $m ";  echo -en "\r .. $n ";  let n=$n+1
done
echo -en "\r $m genes processed.\n"

m=`cat tmp30 | wc -l`; n=0
for i in `cat tmp30`; do
  cat $i | head -1 | sed 's/gene_id "/\n/g' | sed '1,1d' | cut -f 1 -d'"' >> tmp230.list 
  echo -en " of $m ";  echo -en "\r .. $n ";  let n=$n+1
done
echo -en "\r $m genes processed.\n"
  cat tmp230.list | sort | uniq > CHS.4_310 

#######################

wc -l CHS.*
  46421 CHS.0
  21457 CHS.1_0
  40188 CHS.3_10
  63055 CHS.4_310
  45282 CHS.5_4310
  13525 CHS.6_54310
 229928 total






chr2=chr9
st2=97805935
en2=97810008

if [ -e $wdir_prev5/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev5/gene_subtract/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev5/trans_subtract/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev5/gene/$chr2/${chr2}_${st2}_${en2} ]
then
    gene_file="$wdir_prev5/gene/$chr2/${chr2}_${st2}_${en2}"
    trans_dir="$wdir_prev5/transcript/$chr2/${chr2}_${st2}_${en2}"
elif [ -e $wdir_prev4/gene_subtract/$chr2/${chr2}_${st2}_${en2} ]
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

echo  "$gene_file" 



# more clean up some left Dbxref= 
cat ../../lncRNAKB_hg38.gtf > ../../lncRNAKB_hg38_pre01.gtf
cat ../../lncRNAKB_hg38_pre01.gtf | sed 's/transcript_id=/transcript_id_alias "/g'  \
    | sed  's/;Dbxref=/"; Dbxref "/g' | sed 's/;/"; /g' | sed 's/=/ "/g' | sed 's/""/"/g' \
     | sed 's/;  [012.] transcript_id /;  transcript_id /g' > ../../lncRNAKB_hg38.gtf

cat lncRNAKB_hg38_pre01.gtf | sed 's/transcript_id=/transcript_id_alias "/g'  \
    | sed  's/;Dbxref=/"; Dbxref "/g' | sed 's/;/"; /g' | sed 's/=/ "/g' | sed 's/""/"/g' \
    | sed 's/;  [012.] transcript_id /;  transcript_id /g' \
    | sed 's/gene_id\t/.\tgene_id /g'  > lncRNAKB_hg38.gtf


Statistics:
  Total: 100279 genes with 532406 transcripts containing 1487422 cds.
  46421 CHS.0	    # CHESS
  21457 CHS.1_0	    # Fantom
  40188 CHS.3_10    # lincipidia
  63055 CHS.4_310   # NONCODE
  45282 CHS.5_4310  # MITranscriptome
  13525 CHS.6_54310 # BIGtranscriptome

http://10.137.19.76/geneset_venn/index2.php





























































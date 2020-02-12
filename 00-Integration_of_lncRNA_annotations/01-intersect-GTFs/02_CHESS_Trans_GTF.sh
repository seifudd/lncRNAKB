
# GTF trans

# chr_bp_gene_dedup
chr1	134773	140566	-	CHS.8
chr1	181049	184258	+	CHS.9
chr1	184878	199860	-	CHS.10

# /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chr_bp_gene_all
4	8	chr1_11874_14409
9	57	chr1_14362_29370
58	65	chr1_29926_31295
958239	958242	chr11_KI270721v1_random_46567_48854
958243	958262	chr11_KI270721v1_random_50358_53012

# Curated\tGenomic test
chr=chr15
b=24823637
c1=24980459

mkdir -p chess_transcript_GTF
run_it(){
m=$1
n=$2
cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chr_bp_gene_dedup | sed -n "$m,${n}p"   | while read a b c1 d1 e1 ; do
   chr=$a
   c="${a}_${b}_${c1}"
   cat chess_genes_GTF/$chr/$c  | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR }' > tmp11_$c
   cat chess_genes_GTF/$chr/$c  | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR-1 }' | sed '1,1d' > tmp22_$c
   cat chess_genes_GTF/$chr/$c | wc -l >> tmp22_$c
   cat chess_genes_GTF/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print $1"_"$4"_"$5"_"$9 }' | cut -f 1 -d';' | sed 's/transcript_id "//g'  | sed 's/"//g' > tmp_tname_$c
   paste tmp11_$c tmp22_$c tmp_tname_$c > tmp_tall_$c
   mkdir -p  chess_transcript_GTF/$chr/$c
   cat tmp_tall_$c | while read d e f; do
     cat chess_genes_GTF/$chr/$c  | sed -n "$d,${e}p" > chess_transcript_GTF/$chr/$c/$f
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









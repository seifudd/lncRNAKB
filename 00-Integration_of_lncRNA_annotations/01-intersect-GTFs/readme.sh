
wc -l  ../0*/chr_bp_gene_dedup  
   46421 ../00-CHESS2.1/chr_bp_gene_dedup
   27871 ../01-FANTOM5.0.v3/chr_bp_gene_dedup
   56224 ../03-Lncipedia-5.2/chr_bp_gene_dedup
   95580 ../04-NONCODEv5/chr_bp_gene_dedup
   63473 ../05-MITRANSCRIPTOMEv2/chr_bp_gene_dedup
   14087 ../06-BIGTRANSCRIPTOMEv1/chr_bp_gene_dedup
  303656 total

"compare gene coordinate from chess with fantom
any overlap:
	go to trans : compare fantom and chess and if any new from fantom whitin boundry bring it in
	get the last transcrip id , add 1 and assing to the new record
	use the same gene id
	the old geneid will be saved but flag to redirect to new gene id
no overlap : 
	bring the whole record
	creat a new lncrnskb id
"
cat ../00-CHESS2.1/chr_bp_gene_dedup 
chr1	11874	14409	+	CHS.1

echo `head -n 1 ../01-FANTOM5.0.v3/chr_bp_gene_dedup` | sed 's/ /\t/g' > tmp1

ml bedtools

wc -l *_overlap
wc -l *_subtract


F0="../00-CHESS2.1/chr_bp_gene_dedup"
F1="../01-FANTOM5.0.v3/chr_bp_gene_dedup"
F3="../03-Lncipedia-5.2/chr_bp_gene_dedup"
F4="../04-NONCODEv5/chr_bp_gene_dedup"
F5="../05-MITRANSCRIPTOMEv2/chr_bp_gene_dedup"
F6="../06-BIGTRANSCRIPTOMEv1/chr_bp_gene_dedup"

# 1. get overlap 2. then exclude multiple hits and get uniq 3. get subtract (novel)

run_pairwise(){
a=$1 ; b=$2 ; c=$3 ; d=$4 
bedtools intersect -nonamecheck    -a $a -b $b -wa -wb > comparisons/${c}_${d}.multi.overlap
cat comparisons/${c}_${d}.multi.overlap | cut -f -5 | sort | uniq -c | awk '{if($1<2) print}' | cut -c 9- > comparisons/${c}_${d}.one.overlap
bedtools subtract  -nonamecheck -A -a $a -b $b                                        > comparisons/${c}_${d}.subtract
}
run_pairwise $F1 $F0 1 0
run_pairwise $F3 $F0 3 0
run_pairwise $F4 $F0 4 0
run_pairwise $F5 $F0 5 0
run_pairwise $F6 $F0 6 0
run_pairwise $F3 $F1 3 1
run_pairwise $F4 $F1 4 1
run_pairwise $F5 $F1 5 1
run_pairwise $F6 $F1 6 1
run_pairwise $F4 $F3 4 3
run_pairwise $F5 $F3 5 3
run_pairwise $F6 $F3 6 3
run_pairwise $F5 $F4 5 4
run_pairwise $F6 $F4 6 4
run_pairwise $F6 $F5 6 5


run_3wise(){
a=$1 ; b1=$2 ; b2=$3 ; c=$4 ; d1=$5 ; d2=$6 
bedtools intersect -nonamecheck    -a $a -b $b1 -b $b2 -wa -wb > comparisons/${c}_${d1}_${d2}.multi.overlap
cat comparisons/${c}_${d1}_${d2}.multi.overlap | cut -f -5 | sort | uniq -c | awk '{if($1<2) print}' | cut -c 9- > comparisons/${c}_${d1}_${d2}.one.overlap
bedtools subtract  -nonamecheck -A -a $a -b $b1 -b $b2         > comparisons/${c}_${d1}_${d2}.subtract
}
run_3wise $F1 $F0 $F3 1 0 3
run_3wise $F1 $F0 $F4 1 0 4
run_3wise $F1 $F0 $F5 1 0 5
run_3wise $F1 $F0 $F6 1 0 6
run_3wise $F3 $F0 $F1 3 0 1
run_3wise $F3 $F0 $F4 3 0 4
run_3wise $F3 $F0 $F5 3 0 5
run_3wise $F3 $F0 $F6 3 0 6
run_3wise $F4 $F0 $F3 4 0 3
run_3wise $F4 $F0 $F1 4 0 1
run_3wise $F4 $F0 $F5 4 0 5
run_3wise $F4 $F0 $F6 4 0 6
run_3wise $F5 $F0 $F3 5 0 3
run_3wise $F5 $F0 $F4 5 0 4
run_3wise $F5 $F0 $F1 5 0 1
run_3wise $F5 $F0 $F6 5 0 6
run_3wise $F6 $F0 $F3 6 0 3
run_3wise $F6 $F0 $F4 6 0 4
run_3wise $F6 $F0 $F5 6 0 5
run_3wise $F6 $F0 $F1 6 0 1

run_4wise(){
a=$1 ; b1=$2 ; b2=$3 ; b3=$4 ; c=$5 ; d1=$6 ; d2=$7 ; d3=$8 
bedtools intersect -nonamecheck    -a $a -b $b1 -b $b2 -b $b3 -wa -wb > comparisons/${c}_${d1}_${d2}_${d3}.multi.overlap
cat comparisons/${c}_${d1}_${d2}_${d3}.multi.overlap | cut -f -5 | sort | uniq -c | awk '{if($1<2) print}' \
  | cut -c 9- > comparisons/${c}_${d1}_${d2}_${d3}.one.overlap
bedtools subtract  -nonamecheck -A -a $a -b $b1 -b $b2 -b $b3      > comparisons/${c}_${d1}_${d2}_${d3}.subtract
}
run_4wise $F1 $F0 $F3 $F4 1 0 3 4
run_4wise $F1 $F0 $F3 $F5 1 0 3 5
run_4wise $F1 $F0 $F3 $F6 1 0 3 6
run_4wise $F1 $F0 $F4 $F5 1 0 4 5
run_4wise $F1 $F0 $F4 $F6 1 0 4 6
run_4wise $F1 $F0 $F5 $F6 1 0 5 6

run_4wise $F3 $F0 $F1 $F4 3 0 1 4
run_4wise $F3 $F0 $F1 $F5 3 0 1 5
run_4wise $F3 $F0 $F1 $F6 3 0 1 6
run_4wise $F3 $F0 $F4 $F5 3 0 4 5
run_4wise $F3 $F0 $F4 $F6 3 0 4 6
run_4wise $F3 $F0 $F5 $F6 3 0 5 6

run_4wise $F4 $F0 $F1 $F3 4 0 1 3
run_4wise $F4 $F0 $F1 $F5 4 0 1 5
run_4wise $F4 $F0 $F1 $F6 4 0 1 6
run_4wise $F4 $F0 $F3 $F5 4 0 3 5
run_4wise $F4 $F0 $F3 $F6 4 0 3 6
run_4wise $F4 $F0 $F5 $F6 4 0 5 6

run_4wise $F5 $F0 $F1 $F3 5 0 1 3
run_4wise $F5 $F0 $F1 $F4 5 0 1 4
run_4wise $F5 $F0 $F1 $F6 5 0 1 6
run_4wise $F5 $F0 $F3 $F4 5 0 3 4
run_4wise $F5 $F0 $F3 $F6 5 0 3 6
run_4wise $F5 $F0 $F4 $F6 5 0 4 6

run_4wise $F6 $F0 $F1 $F3 6 0 1 3
run_4wise $F6 $F0 $F1 $F4 6 0 1 4
run_4wise $F6 $F0 $F1 $F5 6 0 1 5
run_4wise $F6 $F0 $F3 $F4 6 0 3 4
run_4wise $F6 $F0 $F3 $F5 6 0 3 5
run_4wise $F6 $F0 $F4 $F5 6 0 4 5


run_5wise(){
a=$1 ; b1=$F0 ; b2=$2 ; b3=$3 ; b4=$4 ; c=$5 ; d1=0 ; d2=$6 ; d3=$7 ; d3=$8 
bedtools intersect -nonamecheck    -a $a -b $b1 -b $b2 -b $b3 -b $b4 -wa -wb > comparisons/${c}_${d1}_${d2}_${d3}_${d4}.multi.overlap
cat comparisons/${c}_${d1}_${d2}_${d3}_${d4}.multi.overlap | cut -f -5 | sort | uniq -c | awk '{if($1<2) print}' \
  | cut -c 9- > comparisons/${c}_${d1}_${d2}_${d3}_${d4}.one.overlap
bedtools subtract  -nonamecheck -A -a $a -b $b1 -b $b2 -b $b3      > comparisons/${c}_${d1}_${d2}_${d3}_${d4}.subtract
}
run_5wise $F1 $F3 $F4 $F5 1 3 4 5
run_5wise $F1 $F3 $F4 $F6 1 3 4 6
run_5wise $F1 $F3 $F5 $F6 1 3 5 6
run_5wise $F1 $F4 $F5 $F6 1 4 5 6

run_5wise $F3 $F1 $F4 $F5 3 1 4 5
run_5wise $F3 $F1 $F4 $F6 3 1 4 6
run_5wise $F3 $F1 $F5 $F6 3 1 5 6
run_5wise $F3 $F4 $F5 $F6 3 4 5 6

run_5wise $F4 $F1 $F3 $F5 4 1 3 5
run_5wise $F4 $F1 $F3 $F6 4 1 3 6
run_5wise $F4 $F1 $F5 $F6 4 1 5 6
run_5wise $F4 $F3 $F5 $F6 4 3 5 6

run_5wise $F5 $F1 $F3 $F4 5 1 3 4
run_5wise $F5 $F1 $F3 $F6 5 1 3 6
run_5wise $F5 $F1 $F4 $F6 5 1 4 6
run_5wise $F5 $F3 $F4 $F6 5 3 4 6

run_5wise $F6 $F1 $F3 $F4 6 1 3 4
run_5wise $F6 $F1 $F3 $F5 6 1 3 5
run_5wise $F6 $F1 $F4 $F5 6 1 4 5
run_5wise $F6 $F3 $F4 $F5 6 3 4 5

wc -l comparisons/*

#

sed -i ':a;N;$!ba;s/\n/,/g'







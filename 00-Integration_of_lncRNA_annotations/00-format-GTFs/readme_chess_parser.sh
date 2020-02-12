
cat chess2.1.gff | awk -v FS='\t' -v OFS='\t'  '{if($3=="gene") print NR }' > tmp1
cat chess2.1.gff | awk -v FS='\t' -v OFS='\t'  '{if($3=="gene") print NR-1 }' | sed '1,1d' > tmp2
cat chess2.1.gff | wc -l >> tmp2
cat chess2.1.gff | awk -v FS='\t' -v OFS='\t'  '{if($3=="gene") print $1"_"$4"_"$5 }' > tmp_name
paste tmp* > chr_bp_gene_all

mkdir -p chess_genes
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p chess_genes/$chr
  cat chess2.1.gff | sed -n "$a,${b}p" > chess_genes/$chr/$c
done
}

# loop it
run_it 1 3000 &
run_it 3001 6000 &
run_it 6001 9000 &
run_it 9001 12000 &
run_it 12001 15000 &
run_it 15001 18000 &
run_it 18001 21000 &
run_it 21001 24000 &
run_it 24001 27000 &
run_it 27001 30000 &
run_it 30001 33000 &
run_it 33001 36000 &
run_it 36001 39000 &
run_it 39001 42000 &
run_it 42001 45000 &
run_it 45001 46421 &
wait

mkdir -p chess_genes_attr
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p chess_genes_attr/$chr
  cat chess_genes/$chr/$c | head -1 | cut -f 9 | sed 's/;/\n/g' | sed 's/=/\t/g' > chess_genes_attr/$chr/$c
done
}

# loop it
for i in {0..14}; do
   let m=($i*3000)+1 ; let n=($i+1)*3000;   run_it $m $n &
done
run_it 45001 46421 &
wait



mkdir -p chess_trans
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p chess_trans/$chr
   cat chess_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR }' > tmp11_$c
   cat chess_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print NR-1 }' | sed '1,1d' > tmp22_$c
   cat chess_genes/$chr/$c | wc -l >> tmp22_$c
   cat chess_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print $1"_"$4"_"$5"_"$9 }' | cut -f 1 -d';' | sed 's/ID=//g' > tmp_tname_$c
   paste tmp11_$c tmp22_$c tmp_tname_$c > tmp_tall_$c
   mkdir -p  chess_trans/$chr/$c
   cat tmp_tall_$c | while read d e f; do
     cat chess_genes/$chr/$c | sed -n "$d,${e}p" > chess_trans/$chr/$c/$f
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


mkdir -p chess_trans_attr
run_it(){
m=$1
n=$2
cat chr_bp_gene_all | sed -n "$m,${n}p"   | while read a b c; do
  chr=`echo $c | cut -f 1 -d'_'`
  mkdir -p chess_trans_attr/$chr
   cat chess_genes/$chr/$c | awk -v FS='\t' -v OFS='\t'  '{if($3=="transcript") print }' | cut -f 9 > tmp_trans_$c
   mkdir -p  chess_trans_attr/$chr/$c
   cat tmp_trans_$c | while read line; do 
      ID=`echo -e "$line" | sed 's/;/\n/g' | sed 's/=/\t/g' | grep -w ^ID | cut -f 2 | head -1 `
      echo -e "$line" | sed 's/;/\n/g' | sed 's/=/\t/g' > chess_trans_attr/$chr/$c/${c}_$ID
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

cat chess2.1.gff | awk '{if($3=="gene") print $1"\t"$4"\t"$5"\t"$7"\t"$9}' | sed 's/;/\t/g' | sed 's/ID=//g' | cut -f -5 > chr_bp_gene_dedup



############## 

compare gene coordinate from chess with fantom
any overlap:
	go to trans : compare fantom and chess and if any new from fantom whitin boundry bring it in
	get the last transcrip id , add 1 and assing to the new record
	use the same gene id
	the old geneid will be saved but flag to redirect to new gene id
no overlap : 
	bring the whole record
	creat a new lncrnskb id




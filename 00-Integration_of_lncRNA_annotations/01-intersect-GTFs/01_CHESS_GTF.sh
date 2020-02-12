

# cp and rename chess_gene and chess_gene_att : gene and gene_att

# cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chess2.1.gff | grep chr21_GL383579v2_alt

run_it(){
m=$1
n=$2
cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chr_bp_gene_dedup | sed -n "$m,${n}p"  | while read a b c d e ; do
chess_trans="/data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chess_trans"
mkdir -p chess_genes_GTF/$a/
cat chess_genes/$a/${a}_${b}_${c}  | head -1  | awk '{if ($3=="gene") { print } }' | sed 's/ID=/gene_id "/g' | sed 's/;/"; /g'  \
    | sed 's/=/ "/g' | sed 's/GENCODE_GENE_NAME/gene_name/g' | sed 's/GENE_NAME/gene_name/g' | awk '{print $0"\";"}' \
    | sed "s/Curated\tGenomic/Curated_Genomic/g" | sed "s/Curated Genomic/Curated_Genomic/g" | sed 's/ /\t/8' > chess_genes_GTF/$a/${a}_${b}_${c}

gene_id=`cat chess_genes/$a/${a}_${b}_${c}  | head -1  | awk '{if ($3=="gene") { print $9} }' | sed 's/ID=//g' | cut -f 1 -d';'`

for i in `ls $chess_trans/$a/${a}_${b}_${c}`; do
cat $chess_trans/$a/${a}_${b}_${c}/$i | head -1 \
    | sed 's/Parent=/gene_id "/g' | sed 's/;transcript_id/;transcript_id_alias/g' | sed 's/ID=/transcript_id "/g' | sed 's/;/"; /g' | sed 's/=/ "/g' \
    | awk '{print $0"\";"}' | sed 's/\t/ /g' | sed "s/Curated\tGenomic/Curated_Genomic/g" | sed "s/Curated Genomic/Curated_Genomic/g" \
    | sed 's/ /\t/8' >> chess_genes_GTF/$a/${a}_${b}_${c}

cat $chess_trans/$a/${a}_${b}_${c}/$i | sed '1,1d' | grep -v ^# \
    | awk -v gene_id=$gene_id '{$9="gene_id \""gene_id"\"; "$9; print}' | sed 's/Parent=/transcript_id "/g'  \
    | awk '{print $0"\";"}' | sed 's/\t/ /g' | sed "s/Curated\tGenomic/Curated_Genomic/g" | sed "s/Curated Genomic/Curated_Genomic/g" \
    | sed 's/ /\t/8' >> chess_genes_GTF/$a/${a}_${b}_${c}
done

cat chess_genes_GTF/$a/${a}_${b}_${c} | sed 's/\t/ /g' | cut -f -8 -d' ' | sed 's/ /\t/g' > tmp1_${a}_${b}_${c}
cat chess_genes_GTF/$a/${a}_${b}_${c} | sed 's/\t/ /g' | cut -f 9- -d' '  > tmp2_${a}_${b}_${c}
paste tmp1_${a}_${b}_${c} tmp2_${a}_${b}_${c} >  chess_genes_GTF/$a/${a}_${b}_${c}
rm -f tmp1_${a}_${b}_${c} tmp2_${a}_${b}_${c} 
done
}

for i in {0..14}; do
   let m=($i*3000)+1 ; let n=($i+1)*3000;   run_it $m $n &
done
run_it 45001 46421 &
wait

##### correct for chrUn , alt and random # they are in wrong folder: in 00-CHESS2.1/
# chr13_KI270842v1_alt/chr13_KI270842v1_alt_22003_24579 is in chr13/chr13_KI270842v1_alt_22003_24579

cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chr_bp_gene_dedup | grep _alt     > CHESS2.1_chr_bp_gene_dedup_alt_random
cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chr_bp_gene_dedup | grep _random >> CHESS2.1_chr_bp_gene_dedup_alt_random
cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chr_bp_gene_dedup | grep chrUn   >> CHESS2.1_chr_bp_gene_dedup_alt_random

cat CHESS2.1_chr_bp_gene_dedup_alt_random  | while read a b c d e ; do
chess_trans="/data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chess_trans"
mkdir -p chess_genes_GTF/$a/
chr=`echo $a | cut -f 1 -d'_'`
cat chess_genes/$chr/${a}_${b}_${c}  | head -1  | awk '{if ($3=="gene") { print } }' | sed 's/ID=/gene_id "/g' | sed 's/;/"; /g'  \
    | sed 's/=/ "/g' | sed 's/GENCODE_GENE_NAME/gene_name/g' | sed 's/GENE_NAME/gene_name/g' | awk '{print $0"\";"}' \
    | sed "s/Curated\tGenomic/Curated_Genomic/g" | sed "s/Curated Genomic/Curated_Genomic/g" \
    | sed 's/ /\t/8' > chess_genes_GTF/$a/${a}_${b}_${c}

gene_id=`cat chess_genes/$chr/${a}_${b}_${c}  | head -1  | awk '{if ($3=="gene") { print $9} }' | sed 's/ID=//g' | cut -f 1 -d';'`

for i in `ls $chess_trans/$chr/${a}_${b}_${c}`; do
cat $chess_trans/$chr/${a}_${b}_${c}/$i | head -1 \
    | sed 's/Parent=/gene_id "/g' | sed 's/;transcript_id/;transcript_id_alias/g' | sed 's/ID=/transcript_id "/g' | sed 's/;/"; /g' | sed 's/=/ "/g' \
    | awk '{print $0"\";"}' | sed 's/\t/ /g' | sed "s/Curated\tGenomic/Curated_Genomic/g" | sed "s/Curated Genomic/Curated_Genomic/g" \
    | sed 's/ /\t/8' >> chess_genes_GTF/$a/${a}_${b}_${c}

cat $chess_trans/$chr/${a}_${b}_${c}/$i | sed '1,1d' | grep -v ^# \
    | awk -v gene_id=$gene_id '{$9="gene_id \""gene_id"\"; "$9; print}' | sed 's/Parent=/transcript_id "/g'  \
    | awk '{print $0"\";"}' | sed 's/\t/ /g' | sed "s/Curated\tGenomic/Curated_Genomic/g" | sed "s/Curated Genomic/Curated_Genomic/g" \
    | sed 's/ /\t/8' >> chess_genes_GTF/$a/${a}_${b}_${c}
done

cat chess_genes_GTF/$a/${a}_${b}_${c} | sed 's/\t/ /g' | cut -f -8 -d' ' | sed 's/ /\t/g' > tmp1_${a}_${b}_${c}
cat chess_genes_GTF/$a/${a}_${b}_${c} | sed 's/\t/ /g' | cut -f 9- -d' '  > tmp2_${a}_${b}_${c}
paste tmp1_${a}_${b}_${c} tmp2_${a}_${b}_${c} >  chess_genes_GTF/$a/${a}_${b}_${c}
rm -f tmp1_${a}_${b}_${c} tmp2_${a}_${b}_${c} 
done






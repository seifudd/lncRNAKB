#!/bin/bash



# Download lncRNAKB GTF file from: http://psychiatry.som.jhmi.edu/lncrnakb/download.php



awk -F'[; \t"]' 'NR==FNR{a[$0];next} $11 in a' lncRNAKB_gene_ids_lncRNA.txt /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v6.gtf > lncRNAKB_full_record_lncRNA.txt



awk -F'[; \t"]' '$3=="exon" {print $1"\t"$4"\t"$5"\t"$11}' lncRNAKB_full_record_lncRNA.txt > lncRNAKB_exons_lncRNA.bed


sort -k1,1 -k2,2n lncRNAKB_exons_lncRNA.bed > lncRNAKB_exons_lncRNA_sorted.bed

bedtools merge -i lncRNAKB_exons_lncRNA_sorted.bed -c 4 -o distinct > lncRNAKB_exons_lncRNA_merged.bed

awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$4}' lncRNAKB_exons_lncRNA_merged.bed > lncRNAKB_exons_lncRNA_merged_input.bed

./bigWigAverageOverBed hg38.phastCons30way.bw lncRNAKB_exons_lncRNA_merged_input.bed output_lncRNA_exons.tab




awk -F'[; \t"]' 'NR==FNR{a[$0];next} $11 in a' lncRNAKB_gene_ids_antisense_RNA.txt /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v6.gtf > lncRNAKB_full_record_antisense_RNA.txt



awk -F'[; \t"]' '$3=="exon" {print $1"\t"$4"\t"$5"\t"$11}' lncRNAKB_full_record_antisense_RNA.txt > lncRNAKB_exons_antisense_RNA.bed


sort -k1,1 -k2,2n lncRNAKB_exons_antisense_RNA.bed > lncRNAKB_exons_antisense_RNA_sorted.bed

bedtools merge -i lncRNAKB_exons_antisense_RNA_sorted.bed -c 4 -o distinct > lncRNAKB_exons_antisense_RNA_merged.bed

awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$4}' lncRNAKB_exons_antisense_RNA_merged.bed > lncRNAKB_exons_antisense_RNA_merged_input.bed

./bigWigAverageOverBed hg38.phastCons30way.bw lncRNAKB_exons_antisense_RNA_merged_input.bed output_antisense_RNA_exons.tab

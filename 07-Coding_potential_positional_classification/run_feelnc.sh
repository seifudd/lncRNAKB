#!/bin/bash
set -e



source /data/NHLBI_BCB/Abhilash/06_FeelNC/conda/etc/profile.d/conda.sh
conda activate base
which python
#conda update conda
#conda clean --all --yes

conda activate ./feelnc_install_dir

#FEELnc_filter.pl -i "/data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v5_GB.gtf" -a "/data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/00-CHESS2.1/chess2.1.gtf" -b transcript_biotype=protein_coding -p 30 --monoex=1 -f 0.75 > "./candidate_lncRNA_v2.gtf"

#####Run this
#FEELnc_filter.pl -i "/data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v6.gtf" -a "/data/NHLBI_BCB/Abhilash/06_FeelNC/temp/gencode.v29.gtf" -b transcript_biotype=protein_coding -p 30 --monoex=1 -f 0.75 > "./feelNC_filter/candidate_lncRNA_v5.gtf"
#####

# some stats
# Threshold:0.9
# 355,720 unique transcripts in the candidate lncRNA file
# 97,604 unique genes in candidate lncRNA file 

# 530,960 unique transcripts in lncRNAKB v6 file
# 99,720 unique genes in lncRNAKB v6 file

# 206,694 unique transcripts in gencode.v29 annotation file
# 58,721 unique genes in gencode.v29 annotation file
# 83,129 protein coding transcripts in gencode.v29 annotation file

# 152,244 transcripts removed due to overlap with 31,205 protein coding transcripts in gencode.v29 annotation file
# 158,987 total transcripts removed


#####Run this
#FEELnc_codpot.pl -i "./feelNC_filter/candidate_lncRNA_v5.gtf" -a "/data/NHLBI_BCB/Abhilash/06_FeelNC/temp/gencode.v29.gtf" -g "./hg38.fa" -b transcript_biotype=protein_coding -b transcript_status=KNOWN -mode=intergenic --outdir="./codpot_out_v5"
#####

#####Run this
FEELnc_classifier.pl -i "./codpot_out_v5/candidate_lncRNA_v5.gtf.lncRNA.gtf" -a  "/data/NHLBI_BCB/Abhilash/06_FeelNC/temp/gencode.v29.gtf" > "./codpot_out_v5/lncRNA_classes_v5.txt"
#####


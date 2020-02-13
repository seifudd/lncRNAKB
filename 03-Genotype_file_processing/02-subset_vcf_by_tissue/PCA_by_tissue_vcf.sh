#!/bin/bash

tissue=$1
tissuesubsetfilename=$2
vcfbytissuedir=$3
covbytissueoutdir=$4

#### do convert vcf by tissue to bed, bim and fam
# module load plink

# plink --vcf $vcfbytissuedir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.vcf --out $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL

#### do convert vcf by tissue to ped and map
# module load plink

# plink --vcf $vcfbytissuedir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.vcf --recode12 --out $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL

#### do perform pca by tissue

module load gnuplot
module load eigensoft

smartpca.perl -i $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.sorted.FINAL.ped \
	      -a $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.sorted.FINAL.map \
              -b $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.sorted.FINAL.ped \
	      -k 10 \
              -m 0 \
              -o $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.pca \
              -p $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.plot \
              -e $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.eval \
              -l $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.pca.log 


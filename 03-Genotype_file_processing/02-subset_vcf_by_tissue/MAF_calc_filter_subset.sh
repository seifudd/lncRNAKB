#!/bin/bash

set -e #stop scrip if a command fails
function fail {
echo “FAIL: $@“ >&2
exit 1 # signal failure
}

tissue=$1
tissuebysubsetfile=$2
vcfsubsetdir=$3
covbytissueoutdir=$4

module load bcftools
module load plink

function do_MAF_calc_filter_subset {
	tissuesubsetfilename=`echo $tissuebysubsetfile | sed 's/.txt//'`

#	bcftools +fill-tags $vcfsubsetdir/$tissue/$tissue.GTEx.subset.vcf -o $vcfsubsetdir/$tissue/$tissue.GTEx.subset.MAF.vcf -- -t MAF
#	bcftools view -i 'MAF[0]>0.05' $vcfsubsetdir/$tissue/$tissue.GTEx.subset.MAF.vcf > $vcfsubsetdir/$tissue/$tissue.GTEx.subset.MAF0.05.FINAL.vcf
#	bcftools +fill-tags $vcfsubsetdir/$tissue/$tissue.GTEx.subset.vcf -- -t MAF | bcftools view -i 'MAF[0]>0.05' - > $vcfsubsetdir/$tissue/$tissue_subset_file_name.GTEx.subset.MAF0.05.FINAL.vcf

	bcftools +fill-tags $vcfsubsetdir/$tissue/$tissue.GTEx.subset.vcf -- -t MAF | bcftools view -i 'MAF[0]>0.05' - | grep -v contig | awk '/^#CHROM/ { printf("##contig=<ID=1,length=248956422>\n##contig=<ID=10,length=133797422>\n##contig=<ID=11,length=135086622>\n##contig=<ID=12,length=133275309>\n##contig=<ID=13,length=114364328>\n##contig=<ID=14,length=107043718>\n##contig=<ID=15,length=101991189>\n##contig=<ID=16,length=90338345>\n##contig=<ID=17,length=83257441>\n##contig=<ID=18,length=80373285>\n##contig=<ID=19,length=58617616>\n##contig=<ID=2,length=242193529>\n##contig=<ID=20,length=64444167>\n##contig=<ID=21,length=46709983>\n##contig=<ID=22,length=50818468>\n##contig=<ID=3,length=198295559>\n##contig=<ID=4,length=190214555>\n##contig=<ID=5,length=181538259>\n##contig=<ID=6,length=170805979>\n##contig=<ID=7,length=159345973>\n##contig=<ID=8,length=145138636>\n##contig=<ID=9,length=138394717>\n##contig=<ID=X,length=156040895>\n");} {print;}' | awk '/#/{print;next}{print $0 | "sort -k1,1 -k2,2n"}' | sed 's/##contig=<ID=/##contig=<ID=chr/' | awk '/#/{print;next}{print $0="chr"$0}' > $vcfsubsetdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.sorted.vcf
	
	bgzip $vcfsubsetdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.sorted.vcf && tabix -f -p vcf $vcfsubsetdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.sorted.vcf.gz

	#### do convert vcf by tissue to bed, bim and fam
	plink --vcf $vcfsubsetdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.sorted.vcf.gz --out $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.sorted.FINAL

	#### do convert vcf by tissue to ped and map
	plink --vcf $vcfsubsetdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.FINAL.sorted.vcf.gz --recode12 --out $covbytissueoutdir/$tissue/$tissuesubsetfilename.GTEx.subset.MAF0.05.sorted.FINAL
}

do_MAF_calc_filter_subset


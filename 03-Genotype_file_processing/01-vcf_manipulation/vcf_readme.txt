Websites: http://vcftools.sourceforge.net/perl_module.html#vcf-indel-stats
cat hg38.fa | sed 's/>chr/>/g' > hg38.nochr.fa
 
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=InbreedingCoeff,Description="InbreedingCoeff < -0.3">
##FILTER=<ID=InbreedingCoeff_Filter,Description="InbreedingCoeff <= -0.3">
##FILTER=<ID=LCR,Description="Overlaps a user-input mask">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=VQSRTrancheINDEL99.90to99.95,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -6.1567 <= x < -3.9332">
##FILTER=<ID=VQSRTrancheINDEL99.95to100.00+,Description="Truth sensitivity tranche level for INDEL model at VQS Lod < -40052.6727">
##FILTER=<ID=VQSRTrancheINDEL99.95to100.00,Description="Truth sensitivity tranche level for INDEL model at VQS Lod: -40052.6727 <= x < -6.1567">
##FILTER=<ID=VQSRTrancheSNP99.80to99.90,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -13.0817 <= x < -4.0054">
##FILTER=<ID=VQSRTrancheSNP99.90to99.95,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -24.4173 <= x < -13.0817">
##FILTER=<ID=VQSRTrancheSNP99.95to100.00+,Description="Truth sensitivity tranche level for SNP model at VQS Lod < -96563.5263">
##FILTER=<ID=VQSRTrancheSNP99.95to100.00,Description="Truth sensitivity tranche level for SNP model at VQS Lod: -96563.5263 <= x < -24.4173">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">




#columns in vcd
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	GTEX-T5JW	GTEX-SUCS	GTEX-U3ZM	GTEX-U4B1	GTEX-S4Q7	GTEX-TMMY	GTEX-SE5C	GTEX-S4UYGTEX-S32W	GTEX-S3XE	GTEX-U3ZN

#1-9 are common information
###########
grep -v "##" Bladderv2.vcf |  awk -F "\t" '{print $7}' | head

FILTER
LCR;VQSRTrancheINDEL99.90to99.95
LCR;VQSRTrancheINDEL99.95to100.00
LCR;VQSRTrancheINDEL99.90to99.95
LCR
LCR
LCR;VQSRTrancheINDEL99.95to100.00
LCR;VQSRTrancheINDEL99.95to100.00
LCR;VQSRTrancheINDEL99.95to100.00
LCR;VQSRTrancheINDEL99.95to100.00

grep -v "##" Bladderv2.vcf |  awk -F "\t" '{print $7}' | wc -l
11959664

grep -v "##" Bladderv2.vcf |  awk -F "\t" '{print $7}' | grep "PASS" | wc -l
8750445
grep -v "##" Bladderv2.vcf |  awk -F "\t" '{print $7}' | grep “LCR” | wc -l
2237459
grep -v "##" Bladderv2.vcf |  awk -F "\t" '{print $7}' | grep “LCR;VSQRT” | wc -l

grep -v "##" Bladderv2PASS.vcf | wc -l
8750446

grep -v "##" Bladderv2PASS.vcf | awk -F "\t" '(substr($5,1,2)=="*,")' | wc -l
70649

grep -v "##" Bladderv2PASS.vcf | awk -F "\t" '(substr($5,2,3)==",*")' | wc -l
99667
70649+99667=170316
grep -v "##" Bladderv2PASS.vcf | grep "*" | wc -l
224823

[singhk4@biowulf vcf_by_tissue]$ grep -v "##" Bladderv2PASS.vcf | grep "*" | awk 'length($5)>3 { print; }' | awk -F "\t" '{print $4,$5}' | wc -l
72768
[singhk4@biowulf vcf_by_tissue]$ grep -v "##" Bladderv2PASS.vcf | grep "*" | awk 'length($5)==3 { print; }' | awk -F "\t" '{print $4,$5}' | wc -l
152055


grep -v "##" Bladderv2PASS.vcf | grep "*" | awk -F "\t" '(substr($5,2,3)!=",*")' | wc -l
125156
grep -v "##" Bladderv2PASS.vcf | grep "*" | awk -F "\t" '(substr($5,1,2)!="*,")' | wc -l
154174
125156+154174=279330
279330-224843=54507
grep -v "##" Bladderv2PASS.vcf | grep "*" | awk -F "\t" '(substr($5,2,3)!=",*")' | awk -F "\t" '(substr($5,1,2)!="*,")' | wc -l
54507

What filter tex used
maf and pass

liftover vcf 


#####gtex vcd gawk.haplotypecaller information


bcftools index some.vcf.gz # create the tbi
bcftools index --nrecords some.vcf.gz # get total variant count
bcftools index --stats some.vcf.gz # get variant count per chromsome

vcf-stats file.vcf.gz #to get sap stats


bcftools filter -e'%TYPE="snp"' in.vcf > indels.vcf
vcftools --vcf DP3g95p5maf05.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.DP3g95p5maf05

zcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz | grep -v "##" | wc -l

#to download taxi
Download the latest tabix jar at http://sourceforge.net/projects/samtools/files/tabix/.
Install C compiler and the ‘make’ utility with the package management system of your specific Linux distribution. In the case of Ubuntu, the following command should suffice:
	$ sudo apt-get install build-essential
	$ tar -xjf tabix-0.2.6.tar.bz2
	$ cd tabix-0.2.6
	$ make
#to make vcd index file using tabix
$ /path/to/tabix/bgzip myvcf.vcf
	$ /path/to/tabix/tabix -p vcf myvcf.vcf.gz

newwall --jobid 15497649 --time 4:00:00:00

https://www.biostars.org/p/59492/
Question: Generate vcf.gz file and its index file vcf.gz.tbi
bgzip -c file.vcf > file.vcf.gz

MOST of the variants in the vcd file are private variants

tabix -p vcf file.vcf.gz

vcfstats for original file:
GTEX-16NFA:
'missing' => 608702,
                                           'indel_count' => 1004475,
                                           'snp_count' => 4815639, (but totaling SNPs gives:5117515 (301876 diff)
                                           'ref' => 48347981,
                                           'hom_RR_count' => 44338346,
                                           'count' => 50253761,
                                           'unphased' => 50253761, (Phased data are ordered along one chromosome and so from these data you know the haplotype. Unphased data are simply the genotypes without regard to which one of the pair of chromosomes holds that allele.)
                                           'snp' => {
                                                      'T>G' => 192562,
                                                      'C>A' => 204746,
                                                      'G>*' => 34735,
                                                      'A>T' => 170254,
                                                      'C>T' => 799833,
                                                      'T>*' => 58067,
                                                      'G>A' => 795612,
                                                      'C>*' => 32227,
                                                      'T>C' => 749648,
                                                      'A>G' => 747015,
                                                      'G>T' => 211480,
                                                      'T>A' => 168328,
                                                      'A>C' => 196829,
                                                      'C>G' => 202641,
                                                      'G>C' => 205198,
                                                      'A>*' => 58340
                                                    },


zcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz | grep -v "##" | awk -F "\t" '{print $1}' | wc -l
50862464
zcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf.gz | grep -v "##" | awk -F "\t" '{print $1}' | wc -l
43752690
INDEL number: 7109835

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf |grep -v "*" | wc -l
42529769

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf |grep  "*" | awk -F "\t" '{print $2"\t"$5}' | awk -F '[\t,]' '{print $1,NF-1}' | head
10169 2
10180 2
10181 2
10250 2
10257 2
10327 2
10423 2
10469 2
10470 2
10473 2
grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf |grep  "*" | awk -F "\t" '{print $1"\t"$2 "\t"$3"\t"$4"\t"$5}' | tail
X	155259777	X_155259777_G_*_b37;X_155259777_G_A_b37	G	*,A
X	155259791	X_155259791_A_*_b37;X_155259791_A_C_b37	A	*,C
X	155259833	X_155259833_G_T_b37;X_155259833_G_*_b37	G	T,*
X	155259835	X_155259835_G_*_b37;X_155259835_G_A_b37	G	*,A
X	155259839	X_155259839_G_*_b37;X_155259839_G_T_b37	G	*,T
X	155259869	X_155259869_T_*_b37;X_155259869_T_G_b37	T	*,G
X	155259905	X_155259905_G_*_b37;X_155259905_G_T_b37	G	*,T
X	155259908	X_155259908_A_*_b37;X_155259908_A_G_b37	A	*,G
X	155260412	X_155260412_G_*_b37;X_155260412_G_T_b37	G	*,T
X	155260419	X_155260419_G_*_b37;X_155260419_G_T_b37	G	*,T
###ALLL the Y chr has dropped off!!!!

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep  "*" | awk -F "\t" '{print $2 "\t"$3"_"$4"_"$5}' | awk -F "_" '{print $3"_"$4 "\t" $7"_"$8"\t"$10"_"$11}' | sed 's|[,*]||g' | head
T_G	T_	T_G
T_C	T_	T_C
A_T	A_	A_T
A_C	A_	A_C
A_C	A_	A_C
T_C	T_	T_C
C_	C_G	C_G
C_	C_G	C_G
G_	G_C	G_C
G_	G_A	G_A
grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep  "*" | awk -F "\t" '{print $2 "\t"$3"_"$4"_"$5}' | awk -F "_" '{print $3"_"$4 "\t" $7"_"$8"\t"$10"_"$11}' | sed 's|[,*]||g' | awk -F "\t" '$1==$3 {print $1,$3}' | head
T_G T_G
T_C T_C
A_T A_T
A_C A_C
A_C A_C
T_C T_C
T_G T_G
T_C T_C
G_C G_C
G_T G_T
grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep  "*" | awk -F "\t" '{print $2 "\t"$3"_"$4"_"$5}' | awk -F "_" '{print $3"_"$4 "\t" $7"_"$8"\t"$10"_"$11}' | sed 's|[,*]||g' | awk -F "\t" '$1==$3 {print $1,$3}' | wc -l
622269

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep  "*" | awk -F "\t" '{print $2 "\t"$3"_"$4"_"$5}' | awk -F "_" '{print $3"_"$4 "\t" $7"_"$8"\t"$10"_"$11}' | sed 's|[,*]||g' | awk -F "\t" '$1!=$3 {print $1,$3}' | wc -l
600652


grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep  "*" | awk -F "\t" '{print $2 "\t"$3"_"$4"_"$5}' | awk -F "_" '{print $3"_"$4 "\t" $7"_"$8"\t"$10"_"$11}' | sed 's|[*,]||g' | awk -F "\t" '$2==$3 {print $2,$3}' | wc -l
581372
 total matches: 622269+ 581372=1203641

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep  "*" | awk -F "\t" '{print $2 "\t"$3"_"$4"_"$5}' | awk -F "_" '{print $3"_"$4 "\t" $7"_"$8"\t"$10"_"$11}' | sed 's|[*,]||g' | awk -F "\t" '$2!=$3 {print $2,$3}' | wc -l
641549
600652+641549=1242201
38506 * entries: check out what they are?????

or possibly this?:
641549-622269=19280
600652-581372=19280


grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep  "*" | wc -l
1222921
1222921-1203641=19280


grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep "*" | awk -F "\t" '(substr($5,1,2)!="*,")' | awk -F "\t" '(substr($5,2,3)!=",*")' | awk -F "\t" '{print $4,$5}' | head
A G,T,*
G A,*,C
G A,*,C
G C,T,*
G T,A,*
C G,*,A
G C,*,T
A G,*,T
A C,*,T
T G,A,C,*

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep "*" | awk -F "\t" '(substr($5,2,3)!=",*")' | awk -F "\t" '(substr($5,1,2)!="*,")' | awk -F "\t" '{print $4,$5}' | wc -l
15554


grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep "*" | awk -F "\t" '(substr($5,1,2)!="*,")' | awk -F "\t" '{print $4,$5}' | head
T G,*
T C,*
A T,*
A C,*
A C,*
T C,*
T G,*
T C,*
G C,*
G T,*

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep "*" | awk -F "\t" '(substr($5,1,2)=="*,")' | awk -F "\t" '{print $4,$5}' | wc -l
585098
grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep "*" | awk -F "\t" '(substr($5,2,3)==",*")' | awk -F "\t" '{print $4,$5}' | wc -l
622269
585098+622269=1207367

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep "*" | awk 'length($5)==3 { print; }' | wc -l
1203641
grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.vcf.gz.recode.vcf | grep "*" | awk 'length($5)>3 { print; }' | wc -l
19280

1222921 total *

wc -l GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissing.vcf 
42529769 

grep "PASS" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissing.vcf | wc -l
38532058
 3997711: non PASS entries
grep -v "PASS" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissing.vcf | wc -l
3997711


#filter out entries without PASS
bcftools view -f PASS

bcftools view  -i  'MIN(FMT/DP)>10 & MIN(FMT/GQ)>15'   my.vcf.gz

to filter maf
To filter on minor allele frequency you need to add :minor after your float like so:

bcftools view -q 0.05:minor chr.all.vcf.gz
SEE: -q, --min-af FLOAT[:nref|:alt1|:minor|:major|:nonmajor]
minimum allele frequency (INFO/AC / INFO/AN) of sites to be printed. Specifying the type of allele is optional and can be set to non-reference (nref, the default), 1st alternate (alt1), the least frequent (minor), the most frequent (major) or sum of all but the most frequent (nonmajor) alleles.


You can add the AF and MAF via the BCFtools +fill-tags plugin.
bcftools +fill-tags test.vcf
#########
AN .. Total number of alleles in called genotypes

AC .. Allele count in genotypes

NS .. Number of samples with data

AC_Hom .. Allele counts in homozygous genotypes

AC_Het .. Allele counts in heterozygous genotypes

AC_Hemi .. Allele counts in hemizygous genotypes (haploid)

AF .. Allele frequency

MAF .. Minor Allele frequency

HWE .. HWE test (PMID:15789306)
About: Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, HWE, MAF, NS.
Usage: bcftools +fill-tags [General Options] -- [Plugin Options]
Options:
   run "bcftools plugin" for a list of common options

Plugin options:
   -d, --drop-missing          do not count half-missing genotypes "./1" as hemizygous
   -t, --tags LIST             list of output tags. By default, all tags are filled.
   -S, --samples-file FILE     list of samples (first column) and comma-separated list of populations (second column)

Example:
   bcftools +fill-tags in.bcf -Ob -o out.bcf
   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t AN,AC
   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -d
   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -S sample-group.txt -t HWE

bcftools +fill-tags Bladderv2PASS.vcf -o Bladderv2PASScalcMAF.vcf -- -t MAF

bcftools view -i 'MAF[0]<0.01' in.vcf > out.vcf

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF.vcf | grep -v "#" | awk -F "\t" '{print $8}' | awk -F ";" '{print $2 "\t" $NF}' | sed 's/\AF=//g' | sed 's/\M//g' | head
0.003263	0.00326264
0.001603	0.00160256
0.0007987	0.000798722
0.0007788	0.000778816
0.093	0.092535
0.001546	0.0015456
0.0007704	0.000770416
0.0007704	0.000770416
0.0007704	0.000770416
0.001543	0.00154321


###find the difference between AF and MAF
grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF.vcf | grep -v "#" | awk -F "\t" '{print $8}' | awk -F ";" '{print $2 "\t" $NF}' | sed 's/\AF=//g' | sed 's/\M//g' | awk 'NR>1 {{a=$1-$2}{print a}}' | sort | tail
9e-09
9e-09
9e-09
9e-09
9e-09
9e-09
9e-09
9e-09
9e-09
9e-09


grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASS.vcf | wc -l
38532059

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.vcf | wc -l
6130330

grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.vcf | awk -F "\t" '{print $8}' | awk -F ";" '{print $NF}' | sed 's/\MAF=//g' | sort | head
0.050077
0.050077
0.050077
0.050077
0.050077
0.050077
0.050077
0.050077
0.050077
0.050077

# after cross mapping to hrch38
grep -v "##" GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38n.vcf | cut -f 2,3 | wc -l
6096726
6130330-6096726=33604


####OVERLOOKED!!!!!
when checking
grep -v ^# GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38n.vcf | awk '{print $3}' | grep ";" | wc -l
78677

####one last step to remove multi allelic from passmaf0.05.hg38n.vcf file
awk '/#/{print;next}{if($5 !~ /,/ && length($5)==1 && length($4)==1){print}}' file.vcf

grep -v ^# GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.vcf | wc -l
6018048

grep -v ^# GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38.FINAL.vcf | awk -F "\t" '{print $1,$2,$3}' > gtex.vcf.b37snpID.positions.txt


wc -l gtex.vcf.b37snpID.positions.txt 
6018048 gtex.vcf.b37snpID.positions.txt


crossmap.std.err
[+] Loading crossmap  0.2.7 
@ 2019-01-08 07:00:05: Read chain_file:  GRCh37_to_GRCh38.chain.gz
@ 2019-01-08 07:08:46: Total entries: 6130329
@ 2019-01-08 07:08:46: Failed to map: 33604

grep -v ^# GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCallerSNPonly.NoMulti.noMissingV2.PASSMAF0.05.hg38n.vcf | awk '{print $3,$2}' | awk -F"_" '{print $2, $NF}' | awk '$1!=$3 {print $1,$3}' | wc -l
6041376

##########

DEALING with multiallelic sips
 split these multi-allelic calls into multiple lines.

bcftools norm -Ov -m-any MyFile.vcf > MyFile.Split.vcf
Also get into the habit of checking your ref alleles against an existing reference genome (this also left-aligns indels):

bcftools norm -m-any MyFile.vcf | bcftools norm -Ov --check-ref w -f /ReferenceMaterial/1000Genomes/human_g1k_v37.fasta > MyFile.Split.RefCheck.vcf

OR------
https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_LeftAlignAndTrimVariants.php
-------------------------
https://samtools.github.io/bcftools/howtos/query.html (IMPOrtANT!!!!)

###########################################################################################
###########################details of bed bim and fam files
PLINK Binary files (BED/BIM/FAM)
PLINK is a very widely used application for analyzing genotypic data. It can be considered the “de-facto” standard of the field, although newer formats are starting to be widespread as well. The binary PLINK format description can be accessed  at their site

The binary PLINK fromat contains the same information as the flat file PLINK format but in a compressed and signifficantly more efficient form. You may use the FAM file to deliver sex and affection data for GWASpi to perform GWS studies. Optionally, if you want to provide the additional data fields that GWASpi supports, you may chose to use the Sample Info & Phenotype format.

BED files
The BED files are encoded in binary format. A description of it’s content con be found here

BIM files
The fields in a MAP file are:

Chromosome
Marker ID
Genetic distance
Physical position
Allele 1
Allele 2
Example of a BIM file of the binary PLINK format:
21	rs11511647	0	26765	A	T
X	rs3883674	0	32380	C	G
X	rs12218882	0	48172	T	T
9	rs10904045	0	48426	A	T
9	rs10751931	0	49949	C	T
8	rs11252127	0	52087	A	C
10	rs12775203	0	52277	A	A
8	rs12255619	0	52481	G	T
FAM files
The fields in a FAM file are

Family ID
Sample ID
Paternal ID
Maternal ID
Sex (1=male; 2=female; other=unknown)
Affection (0=unknown; 1=unaffected; 2=affected)
Example of a FAM file of the binary PLINK format:
FAM1	NA06985	0	0	1	1
FAM1	NA06991	0	0	1	1
0	NA06993	0	0	1	1
0	NA06994	0	0	1	1
0	NA07000	0	0	2	1
0	NA07019	0	0	1	1
0	NA07022	0	0	2	1
0	NA07029	0	0	1	1
FAM2	NA07056	0	0	0	2
FAM2	NA07345	0	0	1	1
###########################details of bed bim and fam files################################
###########################################################################################

#############################################################################
#############################################################################
#############################################################################
##after tissue subsetting of vcf file:

wc -l blood.GTEx.subset.MAF0.05.FINAL.vcf 
5948231 blood.GTEx.subset.MAF0.05.FINAL.vcf


cat blood.bed | sed 1d | sort -k1,1 -k 2,2n -k 3,3n | bgzip > blood.bed.gz

 tabix -f -p bed blood.bed.gz > blood.bed.gz.tbi

bgzip blood.GTEx.subset.MAF0.05.FINAL.vcf

################################################
#############to vcf sort
vcf-sort "$script_path""/""$input" >
> "$script_path""/""$input"".sorted"; 
###but this sorts the first column only
########instead use this!!!!!
cat in.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > out_sorted.vcf
cat in.vcf | awk '/#/{print;next}{print $0 | "sort -k1,1 -k2,2n"}' >

######sort bed file
cat bed | awk 'NR == 1; NR > 1 {print $0 | "sort -k1,1 -k2,2n -k3,3n"}'

################################################
##########Use awk to insert the contig before the #CHROM line. Something like
awk '/^#CHROM/ { printf("##contig=<ID=1,length=195471971>\n##contig=<ID=2,length=182113224>\n");} {print;}' in.vcf > out.vcf
##########you can generate the printf line above from the ref.fai file
awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ref.fai

#######some nuances
 awk '$1 ~ /a/' cols.txt ##a is in the first column
 awk '$1 ~ /^a/' cols.txt ##a is in the start of the first column
 awk '/a/' cols.txt ##a is in any column



################################
###################
######VCF file and I want to change part of the header:
	##contig=<ID=1,length=195471971>
	##contig=<ID=10,length=130694993>
to
	##contig=<ID=chr1,length=195471971>
	##contig=<ID=chr10,length=130694993>

########## 1) change the header(s)
bcftools view --header-only $INPUT_FILE | sed 's/##contig=<ID=/##contig=<ID=chr/' | sed 's/##contig=<ID=chrMT/##contig=<ID=chrM/' > $OUTPUT_FILE
############and 2) change the data field (as questioned by @charlesberkn).

bcftools view --no-header $INPUT_FILE | sed 's/^/chr/' | sed 's/^chrMT/chrM/' >> $OUTPUT_FILE


###################################
############to update the contiguous info
awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/crossmap_files/hg38.nochr.fa.fai


##### to only get chr1-22 and chrX do this
cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/crossmap_files/hg38.nochr.fa.fai | grep -v "_" | grep -v M | grep -v Y | awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}'

####you will get this:

##contig=<ID=1,length=248956422>\n##contig=<ID=10,length=133797422>\n##contig=<ID=11,length=135086622>\n##contig=<ID=12,length=133275309>\n##contig=<ID=13,length=114364328>\n##contig=<ID=14,length=107043718>\n##contig=<ID=15,length=101991189>\n##contig=<ID=16,length=90338345>\n##contig=<ID=17,length=83257441>\n##contig=<ID=18,length=80373285>\n##contig=<ID=19,length=58617616>\n##contig=<ID=2,length=242193529>\n##contig=<ID=20,length=64444167>\n##contig=<ID=21,length=46709983>\n##contig=<ID=22,length=50818468>\n##contig=<ID=3,length=198295559>\n##contig=<ID=4,length=190214555>\n##contig=<ID=5,length=181538259>\n##contig=<ID=6,length=170805979>\n##contig=<ID=7,length=159345973>\n##contig=<ID=8,length=145138636>\n##contig=<ID=9,length=138394717>\n##contig=<ID=X,length=156040895>\n[singhk4@cn0629 blood]$ 


#######################################
########### bgzip and tabix
bgzip genotypes.vcf && tabix -p vcf genotypes.vcf.gz
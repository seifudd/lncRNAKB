

<b>Scripts to proces the whole genome sequece (WGS) genotype data from GTEx Release v7 project</b>

VCF files were processed using the following steps with a combination of PLINK v1.9 vcftools v0.1.15 and bcftools  v1.9

<ol type="1">
<li>Remove indels.</li>
  <li>Exclude missing and multi-allelic variants.</li>
  <li>Selected "FILTER == 'PASS'" variants.</li>
  <li>Exclude variants with minor allele frequency (MAF) <5%.</li>
  <li>Update the coordinates of single nucleotide polymorphisms (SNPs) using the UCSC liftOver tool from hg19 to hg38 (latest genome build).</li> 
  <li>Change the SNPs IDs to dbSNP rsID using dbSNP Build 150.</li>
  <li>Convert to bed, bim and fam format.</li>
  </ol>

The VCF file was subset by tissue and the MAF recalculated to exclude variants with MAF <5%. 

After converting to ped and map format, we ran principal component analysis (PCA) on each tissue to get a set of genotype
covariates using eigensoft v6.1.4. 

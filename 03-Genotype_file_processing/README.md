

<b>Scripts to proces the whole genome sequece (WGS) genotype data from GTEx Release v7 project</b>

VCF files were processed using the following steps with a combination of PLINK v1.9 vcftools v0.1.15 and bcftools  v1.9

<ol type="1">
<li>Remove indels.</li>
  <li>Exclude missing and multi-allelic variants.</li>
  <li>Selected "FILTER == 'PASS'" variants.</li>
  <li>Exclude variants with minor allele frequency (MAF) <5%.</li>
  <li>Update the coordinates of single nucleotide polymorphisms (SNPs) using the UCSC liftOver tool from hg19 to hg38 (latest genome build).</li> 
  <li>Change the SNPs IDs to dbSNP50 rsID using dbSNP Build 151.</li>
  <li>Convert to bed, bim and fam format.</li>
  </ol>

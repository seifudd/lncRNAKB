

<b>Scripts to to perform expression quantitative trait loci (eQTL) analysis</b>

<ol type="1">
<li>Genes were first filtered based on TPM to include genes with TPM >0.50 in at least 20% of the samples.</li> 
<li>Genes were filtered based on raw counts to include protein-coding genes and non-coding genes with counts >2 and >1 in at least 20% of the samples, respectively.</li>
<li>The edgeR and limma-voom package in R were used to process the filtered read counts into log2 counts per million (log2CPM) that were normalized using trimmed mean of M-values (TMM).</li>
<li>The expression files were then sorted by gene start and stop, compressed with BGZIP and indexed with TABIX.</li> 
<li>Only tissues with >80 samples were included in the cis-eQTL analysis.</li> 
<li>For eQTL analysis, the first five principal components (PCs) (see Genotype file processing), sex, genotype platform and 15
  PEER factors were included as covariates.</li> 
<li>Within each tissue, cis-eQTLs were identified by linear regression, as implemented in FastQTLv2.0 (threaded option), adjusting for all the covariates.</li> 
<li>We restricted our search to variants within 1 megabase (Mb) of the transcription start site (TSS) of each gene.</li> 
<li>We used the Benjamini and Hochberg correction method to calculate the false discovery rate (FDR) in R statistical programming language (R) across all SNP-gene pairs.</li>
<li>In addition, we ran the adaptive permutations option in FastQTL between 1000 and 10000 permutations.</li> 
<li>Once we obtained the permutation p-values for all the genes, we accounted for multiple testing using FDR to determine the top cis-eQTLs.</li> 
<li>For each tissue, all cis-eQTL results were visualized using a Manhattan plot created using the qqman package in R.</li>
</ol>

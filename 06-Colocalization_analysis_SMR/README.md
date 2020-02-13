

<b>Scripts to analyze lncRNA-GWAS trait association using Summary Mendelian Randomization (SMR)</b>

<ol type="1">
  <li>SMR and HEIDI (Heterogeneity in dependent instruments) methods implemented in the SMR package were used to test the association between lncRNA gene expression and traits tested by means of colocalization of summary GWAS and cis-eQTL signals.   </li> 
  <li>Particularly, HEIDI uses multiple SNPs (n = 20) in a cis-eQTL region to distinguish pleiotropy from linkage, and a pHEIDI >0.05 suggests non-heterogeneity, thus colocalized.</li>
  <li>Summary GWAS data for 323 traits with >5,000 cases available in the UK Biobank were downloaded</li>
  <li>Formatted into .ma format as specified on the CNS genomics’ website (http://cnsgenomics.com/software/smr/).</li>
  <li>Results from the eQTL analysis were filtered by FDR ≤0.05 and formatted into BESD format.</li> 
  <li>SMR was then conducted separately using GWAS meta-analyses summary data for each of the 323 traits using a default cis window of 2000 Kb and peQTL set to 5 × 10−4 for selecting top cis-eQTL SNPs in all tissues with eQTL information.</li>
</ol>

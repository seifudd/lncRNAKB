
<b>Scripts to analyze conservation scores between protein coding genes (PCGs) and lncRNAs</b>

<ol type="1">
<li>Conservation of exons between protein-coding genes and lncRNAs in the lncRNAKB annotation database was analyzed using the bigWigAverageOverBed70 and the cons30way (hg38) track71 both downloaded from the UCSC genome browser.</li>
<li> An exon-level BED file was created using the lncRNAKB GTF annotation file separately for protein-coding genes and lncRNAs.</li> 
<li>We merged (using bedtools) overlapping exons within transcripts to avoid counting conservation scores of overlapping base pairs more than once.</li>
<li>For each exon, the bigWigAverageOverBed function calculates the average conservation score across all base pairs.</li>
</ol>

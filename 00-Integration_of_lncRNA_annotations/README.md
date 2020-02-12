<b>Scripts to systematically combine frequently used lncRNA annotation resources (GTFs/GFFs) using a cumulative stepwise intersection method</b>. 

<ol type="1">
<li>The gene transfer format (GTF) or gene feature format (GFF) from all six annotations (links in Table 1 in the manuscript) were downloaded.</li>
<li>To streamline the data integration step, all the GTF or GFF annotations were parsed to the same format using the following steps:</li>
<ul>
<li>Where required, the coordinates of annotation were updated using the UCSC liftOver tool from
hg19 to hg38 (latest genome build).</li>
<li>For each chromosome, the gene and transcript records were split into individual files labelled by chromosome, strand, start and end base pair locations.</li>
<li>Each gene block file contained the transcripts information and the transcript block file contained the exons information </li>
<li>In cases where the transcripts or exons records lacked genes information, a gene entry was manually created using the gene ids in the transcripts or exons records and combined with the base pair locations of the first exon (as gene start), of the last exon (as gene end), and transcript strand to represent the gene.</li>
strand.
<li>All the redundant records within each annotation file were removed in this process.</li>
</ul>
</ol>

Figure 2 illustrates the cumulative stepwise intersection method for two annotations as an example, D1 (CHESS) in blue and D2 (FANTOM-lncRNA only) in green. 

<ol type="1">
<li>For each gene entry in D1 (top blue panel), we kept genes from D2 (green panel) that had full overlap and were also within D1’s gene boundary. The resulting intersection is shown in orange.</li> 
<li>D2’s gene that had partial overlap with D1’s gene (marked using a red X) were discarded as we did not want to re-define gene boundaries in the reference annotation.</li>
<li>For genes that intersected, the transcript records (shown as smaller bars connected by lines to represent exons and introns, respectively) from D1 and D2 were compared.</li>
<li>Similarly, to the gene intersection, transcript entries whose start and end were within the gene boundaries were included.</li> 
<li>Several transcripts (marked using a red X) that fell outside the gene boundary and were probably incorrectly assigned to genes were removed in this process.</li> 
<li>In addition, if a transcript in D2 had partial or no overlap with transcripts in D1, we incorporated that transcript (marked using red checks) including all the exons to the gene record accordingly.</li> 
<li>For genes with no overlap in D1, we added all the transcripts and corresponding exons to the merged annotation as a lncRNA entry (marked using red checks).</li>
</ol>

![Integration_Procdure](/09-Figures/Figure2.png)

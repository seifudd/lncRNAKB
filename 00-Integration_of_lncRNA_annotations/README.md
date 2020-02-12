

<b>Scripts to systematically combine frequently used lncRNA annotation resources using a cumulative stepwise intersection method</b>. 

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

   
![Repo_List](09-Figures/Figure2.png)

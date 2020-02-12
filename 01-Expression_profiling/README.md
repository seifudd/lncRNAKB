

<b>Scripts to analyze the gene expression data from GTEx Release v7 project</b>
<ol type="1">
  <li>Quality control of paired-end reads were assessed using FastQC tools.</li>
  <li>Adapter sequences and low-quality bases were trimmed using Trimmomatic.</li>
  <li>Aligned to the human reference genome (H. sapiens, GRCh38) using HISAT2.</li>
  <li>Utilizing uniquely aligned reads to the human genome, gene-level expression (raw read counts) were generated with the featureCounts software guided by the lncRNAKB GTF annotation.</li>
  <li>Transcript-level expression was quantified using Salmon guided by the lncRNAKB transcripts FASTA file.</li>
  <li>Raw read counts were transformed to Transcripts Per Kilobase Million (TPM)</li>
<ol>

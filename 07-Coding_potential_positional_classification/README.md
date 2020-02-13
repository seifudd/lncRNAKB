

#### Scripts to analyze the coding potential and positional classification (overlap with mRNA) of lncRNAKB GTF annotation

<ol type="1">
  <li>FlExible Extraction of LncRNAs (FEELnc-https://github.com/tderrien/FEELnc) was used to classify/annotate and calculate
    the coding potential of all the gene entries in the lncRNAKB.</li>
  <li>FEELnc annotates lncRNAs based on a machine learning method, Random Forest (RF), trained with general features such as
    multi k-mer frequencies, RNA sequence length and open reading frames (ORFs) size.</li>
  <li>It is comprised of three modules:</li>
    <ul>
      <li>Filter.</li>
        <ul>
          <li>The filter module flags and removes transcripts overlapping (in sense) exons of the reference annotation, specifically the protein-coding exons.</li>
          <li>GENCODEv29 GTF file was used as the reference annotation to get an estimate of the number of transcripts from lncRNAKB overlapping with “protein_coding” transcripts.</li>
          <li>The minimal fraction out of the candidate lncRNAs size to be considered for overlap to be excluded was arbitrarily set as 0.75 (>75%overlap) to retain many lncRNAs transcripts.</li>
          <li>Transcripts <200 base pairs (bp) long were filtered out and monoexonic transcripts were retained.</li>
        </ul>
      <li>Coding potential.</li>
          <ul>
            <li>The filtered GTF annotation output file from the filter module was used to calculate a coding potential score (CPS) for each transcript using the coding potential module.</li>
            <li>Due to the lack of a gold standard/known human lncRNAs data set for training, we used the “intergenic” mode in the module.</li>
            <li>We used the human reference genome FASTA file (hg38) and the GENCODE GTF file as the reference annotation.</li>
            <li>To get the best training set of known mRNA, we used “transcript_biotype=protein_coding” and “transcript_status=KNOWN” for the RF model.</li>
            <li>We used the default values for the k-mer sizes, number of trees and ORF type.</li>
          </ul>
      <li>Classifier.</li>
        <ul>
          <li>We used the final set of lncRNAs transcripts output from the coding potential module and classified them using the GENCODEv29 GTF file as the reference annotation.</li>
          <li>We used a minimum and maximum window size of 10 kilobase (kb) and 100kb respectively.</li>
        </ul>
    </ul>
</ol>

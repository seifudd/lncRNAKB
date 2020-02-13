

<b>Scripts to perform functional characterization of lncRNA using a network-based approach (WGCNA)</b>

<ol type="1">
<li>Weighted gene co-expression network analysis (WGCNA) approach as implemented in the Co-Expression Modules identification Tool (CEMiTool) package in R was used to identify modules of lncRNA-mRNA clusters that are co-expressed and therefore likely work in concert to carry out various biological functions.</li>
  <ul>
    <li>Gene expression data by tissue (check folder: 04-eQTL_analysis) was filtered by log2CPM >2 in at least 50% of the samples to avoid random correlations between low-expressing genes.</li>
    <li>The default CEMiTool parameters were used with the following exceptions:</li>
      <ul>
        <li>Pearson method was used for calculating the correlation coefficients.</li>
        <li>Network type used was unsigned.</li>
        <li>No additional filter was used for the expression data.</li>
        <li>Applied Variance Stabilizing Transformation (VST).</li>
        <li>Correlation threshold for merging similar modules were set to 0.90.</li>
      </ul>
    <li>All the co-expressed modules were subjected to over-representation analysis (ORA) by module based on the hypergeometric test.</li>  
    <li>We used Gene Ontology (GO) pathways to check for overrepresentation of genes and determined the most significant module functions based on pathways FDR q-value ≤0.0567.</li>    
    <li>The background set used for the pathway enrichment analysis was genes represented across all GO pathways.</li> 
    <li>To visualize the interactions between the genes in each co-expression module, we selected 25 notable pathways for each tissue.</li> 
    <li>The module adjacency matrices for each of these pathways were filtered based on correlations >0.20 across all genes in each pathway.</li>
    <li> A JSON file (one per pathway) was created to produce interactive networks using Cytoscape v3.6.0 JavaScript modules. </li> 
    </ul>
</ol>

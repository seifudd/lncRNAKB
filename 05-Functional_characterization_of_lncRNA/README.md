

<b>Scripts to perform functional characterization of lncRNA using a network-based approach (WGCNA)</b>

<ol type="1">
<li>Weighted gene co-expression network analysis (WGCNA) approach as implemented in the Co-Expression Modules identification Tool (CEMiTool) package in R was used to identify modules of lncRNA-mRNA clusters that are co-expressed and therefore likely work in concert to carry out various biological functions.</li>
  <ul>
    <li>Gene expression data by tissue was filtered by log2CPM >2 in at least 50% of the samples to avoid random
correlations between low-expressing genes.<\li> 
  
  
  The default CEMiTool parameters were used
with the following exceptions: (i) Pearson method was used for calculating the correlation
coefficients, (ii) the network type used was unsigned, (iii) no filter was used for the
expression data, (iv) applied Variance Stabilizing Transformation (VST) and the correlation
threshold for merging similar modules were set to 0.90. All the co-expressed modules were
subjected to over-representation analysis (ORA) by module based on the hypergeometric
test63. We used Gene Ontology (GO) pathways64–66 to check for overrepresentation of genes
and determined the most significant module functions based on pathways FDR q-value
≤0.0567. The background set used for the pathway enrichment analysis was genes
represented across all GO pathways. To visualize the interactions between the genes in each
co-expression module, we selected 25 notable pathways for each tissue. The module
adjacency matrices for each of these pathways were filtered based on correlations >0.20
across all genes in each pathway. A JSON file (one per pathway) was created to produce
interactive networks using Cytoscape v3.6.0 JavaScript modules68. The network files and the
module adjacency/correlation matrix files are available for downloading on lncRNAKB.
</ul>
</ol>

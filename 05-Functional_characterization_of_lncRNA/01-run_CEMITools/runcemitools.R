
library("CEMiTool")
library("AnnotationDbi")
library("org.Hs.eg.db")
library(edgeR)

args <- commandArgs()
tissue <- args[6]
cpmexpressionfile <- args[7]

setwd(paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/",tissue,sep=""))

/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/Bone_Marrow/Bone_Marrow.91samples.12552pcs.9507lincs.gene.CPM.GTExId.txt

data = read.table(file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6/",tissue,"/",cpmexpressionfile,sep=""), header=T, row.names=1, check.names=F)
dim(data)
head(data[1:10,1:10])
annotation = read.table(file="/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/lncRNAKB_hg38_v6_query_table.genes", header=T, row.names=1, colClasses=c("character"))
datawgenenames = merge(x=data, y=annotation, by="row.names", all.x=T)
missingrownames = which(is.na(datawgenenames$gene_name))
datawgenenames$gene_name[missingrownames] = datawgenenames$Row.names[missingrownames]

datadf = data.frame(datawgenenames[,7:(dim(datawgenenames)[2])-1], row.names=make.unique(datawgenenames$gene_name))
# lncrnas = read.table(file="/data/NHLBI_BCB/Cao_Lab/08-CEMITOOLS/03-data-from-Yi-redo-cutoffs/sig_genes.tsv", header=T)
## filtering
#Keep genes with at least CPM>2 in at least 20%/50% replicates/samples?
SampleSize = dim(datadf)[2]
sampleCutoff = round(0.5*SampleSize)
isexpr = rowSums(datadf>1) >= sampleCutoff
table(isexpr)
datadf = datadf[isexpr,]
dim(datadf)

# to get lnckbids only
datawgenenamesfiltered = datawgenenames[datawgenenames$gene_name %in% rownames(datadf),]
dim(datawgenenamesfiltered)
datawgenenamesfilteredlncrnas = datawgenenamesfiltered[datawgenenamesfiltered$gene_type %in% c("lncRNA"),]
datalncrnas = datawgenenamesfilteredlncrnas$gene_name
datalncrnaslnckbidonly = datawgenenamesfilteredlncrnas$Row.names
length(datalncrnas)
length(datalncrnaslnckbidonly)

gmt.cemitool = read_gmt("/data/NHLBI_BCB/Cao_Lab/08-CEMITOOLS/c5.GO.all.v6.1.symbols.gmt")
cemres = cemitool(expr=datadf, gmt=gmt.cemitool, verbose=T, apply_vst=T, cor_method="pearson", filter=F, ora_pval=0.05, network_type="unsigned", diss_thresh=0.90)

module_names = mod_names(cemres)
moduleinfo = matrix(nrow=length(module_names), ncol=8, byrow=T)
for(i in 1:length(module_names)) {
	mgenes = module_genes(cemres, module=module_names[i])
	mgenesvec = mgenes$genes

	c = which(datawgenenames$gene_name %in% mgenesvec)
	d = datawgenenames$Row.names[c]
	if(length(d) > 0){e = paste(d, collapse = "/")} else{ y="NA" }

	# which lncRNAs are in the modules (gene_name)
	a = which(datalncrnas %in% mgenesvec)
	b = datalncrnas[a]
	if(length(b) > 0){z = paste(b, collapse = "/")} else{ y="NA" }

	x = datalncrnaslnckbidonly[a]
	if(length(x) > 0){y = paste(x, collapse = "/")} else{ y="NA" }

	moduleinfo[i,1] = tissue # tissue name
	moduleinfo[i,2] = module_names[i] # module names
	moduleinfo[i,3] = length(mgenes$genes) # number of genes in the module
	moduleinfo[i,4] = length(x) # number of lncRNAs in the module
	moduleinfo[i,5] = y # list of lncRNAs in the module (lnckbid)
	moduleinfo[i,6] = z # list of lncRNAs in the module (gene_name)
	moduleinfo[i,7] = e # list of genes in the module (lnckbid)
}
hubs = get_hubs(cem=cemres, n=20) # Returns n genes in each module with high connectivity
for(i in 1:length(module_names)) {
	moduleinfo[i,8] = paste(names(hubs[[module_names[i]]]), collapse="/") #top 20 most connected genes in each module
}
moduleinfo = as.data.frame(moduleinfo)
colnames(moduleinfo) = c("Tissue","ModuleName","number_genes_in_module","number_lncrnas_in_module","lncrnas_in_module","lncrnas_in_module_v2","genes_in_module","top20_hub_genes")
write.csv(moduleinfo, file="Modules_Identified.csv", row.names=F)

orares = ora_data(cem=cemres) # Retrieve over representation analysis (ORA) results
moduleenrichmentinfo = data.frame(head(orares[orares$Module %in% module_names[1],], n=3))
for(i in 2:length(module_names)) {
	module_name = module_names[i]
	moduleenrichmentinfo = rbind(moduleenrichmentinfo, head(orares[orares$Module %in% module_name,], n=3))
}
colnames(moduleenrichmentinfo) = c("Module","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")
moduleenrichmentinfo = as.data.frame(moduleenrichmentinfo, colClasses=c("character","character","character","numeric","numeric","numeric","numeric","numeric","character","numeric"), stringsAsFactors=FALSE, make.names=T)
moduleenrichmentinfo_moduleinfo_merged = merge(x=moduleenrichmentinfo, y=moduleinfo, by.x="Module", by.y="ModuleName", all.x=T)
write.csv(moduleenrichmentinfo, file="Modules_Identified_Enrichment_top3_pathways.csv", row.names=F)
write.csv(moduleenrichmentinfo_moduleinfo_merged, file="Modules_Identified_Enrichment_top3_pathways_with_module_info.csv", row.names=F)
write.csv(orares, file="Modules_Identified_Enrichment_all_pathways.csv")

orares_moduleinfo_merged = merge(x=orares, y=moduleinfo, by.x="Module", by.y="ModuleName", all.x=T)
write.csv(orares_moduleinfo_merged, file="Modules_Identified_Enrichment_all_pathways_with_module_info.csv")

# get adjacency matrix for these modules
# modules = c("M3","M11","M18","M20")
adj = adj_data(cemres)
for(i in 1:length(module_names)) {
	mgenes = module_genes(cemres, module=module_names[i])
	adjmod = adj[rownames(adj) %in% mgenes$genes, colnames(adj) %in% mgenes$genes]
	adjmoddf = data.frame(row=rownames(adjmod)[row(adjmod)[upper.tri(adjmod)]], col=colnames(adjmod)[col(adjmod)[upper.tri(adjmod)]], corr=adjmod[upper.tri(adjmod)])
	adjmoddf_corr = adjmoddf[(adjmoddf$row %in% mgenes$genes & adjmoddf$col %in% mgenes$genes),]
	write.table(adjmoddf_corr, file=paste("adjacency_matrix_corr_module","_",module_names[i],".txt",sep=""), row.names=F)
}


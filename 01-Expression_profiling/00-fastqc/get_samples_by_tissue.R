
setwd("/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-by-tissue-files-and-ids-rnaseq-only")

gtexdownloaded = read.table(file="/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-fastqc-gtex/gtexfastqdownloaded.txt", header=T, row.names=1)
sraRunTable = read.table(file="/data/NHLBI_BCB/Fayaz/21-GTEx-data/SraRunTable.txt", header=T, row.names=17, sep="\t")
gtexdownloaded_sraRunTable = merge(x=gtexdownloaded, y=sraRunTable, by="row.names")
table(gtexdownloaded_sraRunTable$body_site_s)
data.frame(table(gtexdownloaded_sraRunTable$body_site_s))
table(gtexdownloaded_sraRunTable$histological_type_s)
data.frame(table(gtexdownloaded_sraRunTable$histological_type_s))
gtexdownloaded_sraRunTable = gtexdownloaded_sraRunTable[gtexdownloaded_sraRunTable$Assay_Type_s=="RNA-Seq",]
# tissue = subset(gtexdownloaded_sraRunTable, gtexdownloaded_sraRunTable$histological_type_s=="Heart")
# tissue_ID_R1_R2 = cbind(tissue$Row.names, as.character(tissue$gtex_d_R1), as.character(tissue$gtex_d_R2))
# write.table(tissue_ID_R1_R2, file="tissue_ID_R1_R2.txt", row.names=F, col.names=F, quote=F)
tissuetypes = as.character(data.frame(table(gtexdownloaded_sraRunTable$histological_type_s))[,1])
for(i in 1:length(tissuetypes)){
	tissue = subset(gtexdownloaded_sraRunTable, gtexdownloaded_sraRunTable$histological_type_s==tissuetypes[i])
	tissue_ID_R1_R2 = cbind(tissue$Row.names, as.character(tissue$gtex_d_R1), as.character(tissue$gtex_d_R2))
	colnames(tissue_ID_R1_R2) = c("Run_s","R1","R2")
#	write.table(tissue_ID_R1_R2, file=paste(tissuetypes[i],"_ID_R1_R2.txt",sep=""), row.names=F, col.names=F, quote=F)
	write.table(tissue_ID_R1_R2, file=paste(tissuetypes[i],"_ID_R1_R2.txt",sep=""), row.names=F, col.names=T, quote=F)
}



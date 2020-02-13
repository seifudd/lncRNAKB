
library(tidyr)

setwd("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/")

args <- commandArgs()
tissue <- args[6]
tissueSRA1 <- args[7]
tissueSRA2 <- args[8]

if (is.na(tissueSRA2)) {
	tissueSRA <- tissueSRA1
} else {
	tissueSRA <- paste(tissueSRA1,tissueSRA2,sep=" ")	
}

print(tissue)
print(tissueSRA)

sampleInfo = read.delim("/data/NHLBI_BCB/Fayaz/21-GTEx-data/SraRunTable.txt")
sampleInfo2 = sampleInfo[,c(1,33,32,24,17)]
RNAseq = sampleInfo2[sampleInfo2$Assay_Type_s=='RNA-Seq',]
TissueSRA = tissueSRA
Tissue = tissue
TissueSpecific = RNAseq[RNAseq$histological_type_s==TissueSRA,]
head(TissueSpecific)
dim(TissueSpecific)

###get covariate file

if (tissue %in% c("Pancreas","Prostate")) {
	covfilepath <- Sys.glob(file.path(paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6/",tissue,"/",tissue,".*.subset.GTEx.subset.MAF0.05.FINAL.pca.evec",sep="")))
} else {
	covfilepath <- Sys.glob(file.path(paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files/",tissue,"/",tissue,".*.subset.GTEx.subset.MAF0.05.FINAL.pca.evec",sep="")))
}
covariate = read.table(file=covfilepath, header=F, skip=1, check.names=F)
covariate = covariate[,1:6]
covariate = separate(data = covariate, col = V1, into = c("samples"), sep = "\\:")
names(covariate)<-c("samples","PC1","PC2","PC3","PC4","PC5")
head(covariate)
dim(covariate)
sampleInfo<-read.delim("/data/NHLBI_BCB/Fayaz/21-GTEx-data/SraRunTable.txt")
sampleInfo2<-sampleInfo[,c(1,8,33,32,24,17)]
WGSeq<-sampleInfo2[sampleInfo2$Assay_Type_s!='RNA-Seq',]
uniqueGTEx<-WGSeq[!duplicated(WGSeq$submitted_subject_id_s),]
values<-c("male","female")
index<-c("1","2")
uniqueGTEx$genderCode<-index[match(uniqueGTEx$sex_s,values)]
values<-c("HiSeq X Ten","Illumina HiSeq 2000","Illumina HiSeq 2500","Illumina MiSeq")
index<-c("a","b","c","d")
uniqueGTEx$seqplatform<-index[match(uniqueGTEx$Instrument_s,values)]
names(uniqueGTEx)<-c("assay","Instrument_s","samples","gender","tissue","SRR","sex","sequencingplatform")
head(uniqueGTEx)
dim(uniqueGTEx)
mergeFile<-merge(covariate, uniqueGTEx, by = "samples", all.x=TRUE)
head(mergeFile)
dim(mergeFile)
covariateFile<-mergeFile[,c(1:6,12,13)]
head(covariateFile)
dim(covariateFile)

TPMData = read.delim(paste0("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb-v6/",Tissue,"/",paste(Tissue),".gene.TPM.with.annotation.filtered.txt"),sep = " ")
annot = TPMData[,c(1:4)]
cMatrix = TPMData[,c(5:ncol(TPMData))]
TrueMatch = which(names(cMatrix)%in%TissueSpecific$Run_s, useNames = TRUE)
cMatrix2 = cMatrix[,TrueMatch]
head(cMatrix2)[1:5,1:5]
names(cMatrix2) = TissueSpecific$submitted_subject_id_s[match(names(cMatrix2),TissueSpecific$Run_s)]
head(cMatrix2)[1:5,1:5]

###find intersect between all three
ans = intersect(intersect(TissueSpecific$submitted_subject_id_s,names(cMatrix2)),covariateFile$samples)
SampleSize = length(ans)
ans
# write.table(ans, paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/",tissue,"/", Tissue,".",SampleSize,".subset.txt", sep=""),row.names = FALSE,col.names = FALSE, quote = FALSE)

###subset covariate based on samples
subCovariate = covariateFile[covariateFile$samples%in%ans,]
head(subCovariate)
subCovariate = subCovariate[match(ans,subCovariate$sample),]
head(subCovariate)
FinalCovariate = as.data.frame(t(subCovariate))
head(FinalCovariate)
write.table(FinalCovariate, paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/",tissue,"/",Tissue,".",SampleSize,"samples.covariate.sorted.txt", sep=""),col.names = FALSE, quote = FALSE, sep = "\t")

###get gene expression cutoff based on 20% sample number
sampleCutoff = round(0.2*SampleSize)
subTPM = cMatrix2[,match(ans,names(cMatrix2))]
subTPM$rowsum = rowSums(subTPM[,c(1:ncol(subTPM))]>0.5)
FinalTPM = data.frame(annot,subTPM)
FinalTPM2 = subset(FinalTPM,FinalTPM$rowsum>=sampleCutoff)
totalLincs = dim(subset(FinalTPM2,FinalTPM2$gene_type!="protein_coding"))[1]
totalGenes = dim(subset(FinalTPM2,FinalTPM2$gene_type=="protein_coding"))[1]
table(FinalTPM2$gene_type)
FinalTPM3 = FinalTPM2[,-c(ncol(FinalTPM2))]
# FinalTPM3 = format(FinalTPM3, digits=1, scientific=FALSE, trim=TRUE)
# write.table(FinalTPM3,file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/",tissue,"/",Tissue,".",SampleSize,"samples.",totalGenes,"pcs.",totalLincs,"lincs",".gene.TPM.GTExId.txt", sep=""),sep = "\t", col.names = NA, quote = FALSE)
GeneToInclude = rownames(FinalTPM3)
ColumnsToInclude = colnames(FinalTPM3)

#############################################################################################
###bring count data and import GTEx IDs and match to get GTEx IDs and subset with TPM file
#############################################################################################
CountData = read.delim(paste0("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb-v6/",Tissue,"/",paste(Tissue),".gene.featureCounts.with.annotation.filtered.txt"),sep = " ")
annotC = CountData[,c(1:4)]
cMatrixC = CountData[,c(5:ncol(CountData))]
# matches<-names(cMatrix)%in%TissueSpecific$Run_s
TrueMatchC = which(names(cMatrixC)%in%TissueSpecific$Run_s, useNames = TRUE)
cMatrix2C = cMatrixC[,TrueMatchC]
head(cMatrix2C)[1:5,1:5]
names(cMatrix2C) = TissueSpecific$submitted_subject_id_s[match(names(cMatrix2C),TissueSpecific$Run_s)]
head(cMatrix2C)[1:5,1:5]
FinalCount = data.frame(annotC,cMatrix2C)
FinalCount1 = FinalCount[,colnames(FinalCount)%in%ColumnsToInclude]
FinalCount2 = FinalCount1[rownames(FinalCount1)%in%GeneToInclude,]
FinalCount2 = FinalCount2[,ColumnsToInclude]

# filter by count
# filter, protein-coding genes > 2 in at least 50% of samples
# filter, lncRNAs > 1 in 50% samples
sampleCutoff = round(0.5*SampleSize) # change to at least 50% of samples from at least 20% of samples

FinalCount2pcgs = FinalCount2[FinalCount2$gene_type=="protein_coding",]
dim(FinalCount2pcgs)
exprpcgs = rowSums(FinalCount2pcgs[,-c(1:4)]>2) >= sampleCutoff
table(exprpcgs)
counts_datapcgs = FinalCount2pcgs[exprpcgs,]
dim(counts_datapcgs)

FinalCount2lncRNAs = FinalCount2[FinalCount2$gene_type!="protein_coding",]
dim(FinalCount2lncRNAs)
exprlncRNAs = rowSums(FinalCount2lncRNAs[,-c(1:4)]>1) >= sampleCutoff
table(exprlncRNAs)
counts_datalncRNAs = FinalCount2lncRNAs[exprlncRNAs,]
dim(counts_datalncRNAs)

FinalCount3 = rbind(counts_datapcgs, counts_datalncRNAs)
dim(FinalCount3)
table(FinalCount3$gene_type)
lincsCount = dim(subset(FinalCount3,FinalCount3$gene_type!="protein_coding"))[1]
pcCount = dim(subset(FinalCount3,FinalCount3$gene_type=="protein_coding"))[1]
# write.table(FinalCount3,file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/",tissue,"/",Tissue,".",SampleSize,"samples.",pcCount,"pcs.",lincsCount,"lincs",".gene.featureCounts.GTExId.txt", sep=""),sep = "\t", col.names = NA, quote = FALSE)

################################################################################################
###calculate CPM
################################################################################################
counts_data = FinalCount3
library(edgeR)
library(limma)
dge = DGEList(counts = counts_data[,-c(1:4)])
dge = calcNormFactors(dge, method="TMM")
CPMbyTMM<-cpm(dge)
dim(CPMbyTMM)
head(CPMbyTMM)[1:5,1:7]
annotCPM<-counts_data[,c(1:4)]
FinalCPM<-data.frame(annotCPM,CPMbyTMM)
FinalCPM = format(FinalCPM, digits=1, scientific=FALSE, trim=TRUE)
CPMgenes<-dim(subset(FinalCPM,FinalCPM$gene_type=="protein_coding"))[1]
CPMlincs<-dim(subset(FinalCPM,FinalCPM$gene_type!="protein_coding"))[1]
write.table(FinalCPM,file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/covariate_files-v6-FS/",tissue,"/",Tissue,".",SampleSize,"samples.",CPMgenes,"pcs.",CPMlincs,"lincs",".gene.CPM.GTExId.txt", sep=""),sep = "\t", col.names = NA, quote = FALSE)


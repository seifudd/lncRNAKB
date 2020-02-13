
args <- commandArgs()
tissue <- args[6]

gene_list = read.table(file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/00-SMR-lncRNAKB-UKBiobank/",tissue,"/",tissue,"_genelist_threaded_fdrbh_less_than_0.05.txt",sep=""), header=F)
eQTL_for_esd = read.table(file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/00-SMR-lncRNAKB-UKBiobank/",tissue,"/",tissue,"_eQTL_for.esd", sep=""), header=F)

for (i in 1:dim(gene_list)[1]) {
	gene = as.character(gene_list[i,1])
	print(gene)
	gene_eqtl_table = eQTL_for_esd[eQTL_for_esd$V1 %in% gene,]
	gene_eqtl_table = data.frame(Chr=substr(gene_eqtl_table$V2,4,nchar(as.character(gene_eqtl_table$V2))), SNP=gene_eqtl_table$V3, Bp=gene_eqtl_table$V4, A1=gene_eqtl_table$V5, A2=gene_eqtl_table$V6, Freq=gene_eqtl_table$V7, Beta=gene_eqtl_table$V8, se=gene_eqtl_table$V9, p=gene_eqtl_table$V10)
	write.table(gene_eqtl_table, paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/00-SMR-lncRNAKB-UKBiobank/",tissue,"/ESD/",gene,".esd",sep=""), row.names=F, quote=F)
}


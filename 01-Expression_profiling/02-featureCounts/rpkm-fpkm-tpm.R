
setwd("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb/")

args <- commandArgs()
tissue <- args[6]

counts_data = read.table(file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb/",tissue,"/gene.featureCounts.txt",sep=""), header=T, row.names=1, check.names=F, sep="\t")
dim(counts_data)
subjects_to_include = read.table(file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-by-tissue-files-and-ids-rnaseq-only/",tissue,"_ID_R1_R2.txt",sep=""), header=F, row.names=1)
subjects_to_include = rownames(subjects_to_include)
print("SUBJECTS_TO_INCLUDE_RNASEQ_ONLY")
length(subjects_to_include)
counts_data = counts_data[,(colnames(counts_data) %in% subjects_to_include)]
dim(counts_data)

info = file.info(paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb/",tissue,"/",tissue,"_subjects_to_exclude.txt",sep=""))
if(info$size > 0){
	print("SUBJECTS_TO_EXCLUDE")
	subjects_to_exclude = read.table(file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb/",tissue,"/",tissue,"_subjects_to_exclude.txt",sep=""), header=F, row.names=1)
	subjects_to_exclude = rownames(subjects_to_exclude)
	print(length(subjects_to_exclude))
	counts_data = counts_data[,!(colnames(counts_data) %in% subjects_to_exclude)]
	print(dim(counts_data))
}

annotation = read.table(file="/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb/lncRNAKB_hg38_v1_query_table.txt", header=T, check.names=F)
dim(annotation)
counts_data_w_annotation = merge(x=annotation, y=counts_data, by.x="chessID", by.y="row.names")
counts_data_w_annotation = data.frame(counts_data_w_annotation[,2:dim(counts_data_w_annotation)[2]], row.names=counts_data_w_annotation$chessID)
dim(counts_data_w_annotation)
write.table(counts_data_w_annotation, file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb/",tissue,"/",tissue,".gene.featureCounts.with.annotation.filtered.txt",sep=""))

counts_data = counts_data_w_annotation[,5:dim(counts_data_w_annotation)[2]]
gene_lenghts_kb = (counts_data_w_annotation$gend - counts_data_w_annotation$gstart)/1000

# FPKM
# lib_sizes = colSums(counts_data)
# scaling_factor = lib_sizes/1000000
# counts_data_scaled_by_scaling_factor = sapply(1:ncol(counts_data), function(i) counts_data[,i]/scaling_factor[i])
# counts_data_scaled_by_scaling_factor_and_gene_length = t(sapply(1:nrow(counts_data_scaled_by_scaling_factor), function(i) counts_data_scaled_by_scaling_factor[i,]/gene_lenghts_kb[i]))
# rownames(counts_data_scaled_by_scaling_factor_and_gene_length) = rownames(counts_data)
# colnames(counts_data_scaled_by_scaling_factor_and_gene_length) = colnames(counts_data)
# fpkm = counts_data_scaled_by_scaling_factor_and_gene_length
# write.table(fpkm, file="")

# TPM
counts_data_scaled_by_gene_length = sapply(1:nrow(counts_data), function(i) counts_data[i,]/gene_lenghts_kb[i])
counts_data_scaled_by_gene_length_df = data.frame(matrix(unlist(counts_data_scaled_by_gene_length), nrow=dim(counts_data)[2], byrow=F), stringsAsFactors=FALSE)
counts_data_scaled_by_gene_length_df = t(counts_data_scaled_by_gene_length_df)
dim(counts_data_scaled_by_gene_length_df)
lib_sizes = colSums(counts_data_scaled_by_gene_length_df)
scaling_factor = lib_sizes/1000000
counts_data_scaled_by_gene_length_and_scaling_factor_df = sapply(1:ncol(counts_data_scaled_by_gene_length_df), function(i) counts_data_scaled_by_gene_length_df[,i]/scaling_factor[i])
counts_data_scaled_by_gene_length_and_scaling_factor_df = data.frame(counts_data_w_annotation[,1:4], counts_data_scaled_by_gene_length_and_scaling_factor_df, row.names=rownames(counts_data))
colnames(counts_data_scaled_by_gene_length_and_scaling_factor_df) = colnames(counts_data_w_annotation)
dim(counts_data_scaled_by_gene_length_and_scaling_factor_df)
tpm = counts_data_scaled_by_gene_length_and_scaling_factor_df
write.table(tpm, file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-featurecounts-lncRNAKb/",tissue,"/",tissue,".gene.TPM.with.annotation.filtered.txt",sep=""))


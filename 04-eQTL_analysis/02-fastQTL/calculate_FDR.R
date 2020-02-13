
args <- commandArgs()
tissue <- args[6]
tissue_samplesize <- args[7]

eqtl_results = read.table(gzfile(paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS/",tissue,"/",tissue_samplesize,".subset.GTEx.subset.MAF0.05.FINAL.win.500000.bp.perm.0.cis.nominal.threaded.eqtl.allpairs.txt.gz",sep="")), header=T)

write.table(data.frame(eqtl_results, fdrbh=p.adjust(eqtl_results$pval_nominal,method="BH")), file=paste("/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS/",tissue,"/",tissue_samplesize,".subset.GTEx.subset.MAF0.05.FINAL.win.500000.bp.perm.0.cis.nominal.threaded.eqtl.allpairs.fdrbh.txt", sep=""), row.names=F)


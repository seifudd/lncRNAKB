#!/bin/bash

ml SMR
ml plink
ml crossmap
ml python/3.7

script_dir=$1
tissue=$2
vcf_output_filepath=$3
eqtl_output_filepath=$4
gwas_summary_filename=$5
gwas_summary_folder=$6
smr_dir=$7
threads=$8

# for testing purposes only
# vcf_filename="Liver.118.subset.GTEx.subset.MAF0.05.FINAL.sorted.vcf.gz"
# eqtl_results_filename="Heart.251samples.subset.GTEx.subset.MAF0.05.FINAL.win.500000.bp.perm.0.cis.nominal.threaded.eqtl.allpairs.fdrbh.txt"
# gwas_summary_filename="20002_1065.gwas.imputed_v3.both_sexes.tsv"
# tissue="Liver"
# eqtl_output_filepath="/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS/$tissue/$eqtl_results_filename"
# vcf_output_filepath="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_FS/$tissue/$vcf_filename"
# smr_dir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/00-SMR-lncRNAKB-UKBiobank"
# echo "VCF: $vcf_output_filepath"
# echo "eQTL: $eqtl_output_filepath"

do_preprocessing_eQTL_results_for_SMR_by_tissue (){
	date

#	cat $eqtl_output_filepath | sed 's/["]//g' > /data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS/$tissue/temp
#	mv -f /data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS/$tissue/temp $eqtl_output_filepath

#	cat $eqtl_output_filepath | awk '{if($10<0.05){print $1}}' | sed 's/["]//g' | sed '1,1d' | sort -k1 | uniq > $smr_dir/$tissue/${tissue}_genelist_threaded_fdrbh_less_than_0.05.txt
	
	#Step A-get gene annotation information from lncRNAKBv6 GTF file
#	cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v6.gtf | grep -v "##" | awk -F"\t" '{if ($3=="gene") print $1,$4,$7,$9}' | awk -F";" '{print $1}' | sed 's/gene_id//' | sed 's/\"//g' | awk -F" " '{print $4,$1,$2,$3}' | sort | uniq | awk -F" " -v OFS="\t" '{print $2,$1,"0",$3,$1,$4}' > $smr_dir/$tissue/lncRNAKB_v6_GeneID_position_strand_for_epi.txt

	#Step B-subset lncRNAKB_v6_GeneID_position_strand_for_epi.txt for genes with eQTL fdrbh < 0.05 (see above)
#	awk 'NR==FNR{a[$0];next}NR!=FNR{c=$2; if (c in a) {print $0 "\t" a[c]}}' $smr_dir/$tissue/${tissue}_genelist_threaded_fdrbh_less_than_0.05.txt $smr_dir/$tissue/lncRNAKB_v6_GeneID_position_strand_for_epi.txt > $smr_dir/$tissue/${tissue}_genelist_threaded_fdrbh_less_than_0.05.epi

	#Step C-create flist file, remove "chr" from chromosome label
#	cat $smr_dir/$tissue/${tissue}_genelist_threaded_fdrbh_less_than_0.05.epi | sed 's/chr//' | awk -v smr_dir="$smr_dir" -v tissue="$tissue" -v OFS="\t" '{print $0, smr_dir"/"tissue"/ESD/"$2".esd"}' | awk 'BEGIN{print "Chr" "\t"	"ProbeID" "\t"	"GeneticDistance" "\t" "ProbeBp" "\t" "Gene" "\t" "Orientation" "\t" "PathOfEsd"};{print}' > $smr_dir/$tissue/${tissue}_genelist_threaded_fdrbh_less_than_0.05_epi_nochr.flist

	#Step D-create a snp info txt file containing: chr, position, snpID, ref, alt from the VCF file
#	zcat $vcf_output_filepath | grep -v ^# |  awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$4,$5}' >  $smr_dir/$tissue/${tissue}_snps_chr_position_ref_alt.txt

	#Step E-merge snp info txt file containing: chr, position, snpID, ref, alt from the VCF file with eQTL results
#	awk -v OFS="\t" 'NR==FNR{a[$3]=$1 "\t" $2 "\t" $4 "\t" $5; next}NR!=FNR{c=$2; if (c in a) {print $0 "\t" a[c]}}' $smr_dir/$tissue/${tissue}_snps_chr_position_ref_alt.txt $eqtl_output_filepath > $smr_dir/$tissue/${tissue}_threaded_eQTL_results_with_chr_position_ref_alt.txt

	#Step F-create ESD file for SMR; pre-step for creating BESD files for SMR	
#	cat $smr_dir/$tissue/${tissue}_threaded_eQTL_results_with_chr_position_ref_alt.txt | awk -v OFS="\t" '{print $1,$11,$2,$12,$13,$14,$6,$8,$9,$7}' > $smr_dir/$tissue/${tissue}_eQTL_for.esd

	#Step G-create ESD file for SMR by gene
#	mkdir -p $smr_dir/$tissue/ESD/
#	for i in `cat $smr_dir/$tissue/${tissue}_genelist_threaded_fdrbh_less_than_0.05.txt`; do echo $i; grep -w "$i" $smr_dir/$tissue/${tissue}_eQTL_for.esd | cut -f2- |sed 's/chr//' | awk 'BEGIN{print "Chr" "\t" "SNP"  "\t" "Bp" "\t" "A1" "\t" "A2" "\t" "Freq" "\t" "Beta" "\t" "se" "\t" "p"};{print}' > $smr_dir/$tissue/ESD/$i.esd; done

#	echo "Creating ESD files per gene..."
#	echo "..."
#	echo "..."

	#Step G-create ESD file for SMR by gene, faster than above
#	Rscript /data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/SMR_make_esd_file.R $tissue

#	echo "Creating BESD file..."
#	echo "..."
#	echo "..."

	#create besd files
#	smr --eqtl-flist $smr_dir/$tissue/${tissue}_genelist_threaded_fdrbh_less_than_0.05_epi_nochr.flist --make-besd --peqtl-other 5.0e-4 --out $smr_dir/$tissue/${tissue}_threaded_subset_fdrbh_less_than_0.05_eQTL.nochr

#	echo "Creating BED, BIM and FAM file from VCF..."
#	echo "..."
#	echo "..."

	##create bed bim fam from VCF file
	##create no header VCF file
	vcf_filename_mod=`echo $vcf_output_filepath | awk -F"/" '{print $NF}' | sed 's/.vcf.gz//g'`
	zcat $vcf_output_filepath | grep -v "#" | cut -c4- > $smr_dir/$tissue/$vcf_filename_mod.no.header.vcf
	##get header from original VCF file
	zcat $vcf_output_filepath | grep ^# > $smr_dir/$tissue/$vcf_filename_mod.header.vcf
	##append header from VCF file to no header VCF file
	cat $smr_dir/$tissue/$vcf_filename_mod.no.header.vcf >> $smr_dir/$tissue/$vcf_filename_mod.header.vcf

	mv -f $smr_dir/$tissue/$vcf_filename_mod.header.vcf $smr_dir/$tissue/$vcf_filename_mod.vcf
	rm -f $smr_dir/$tissue/$vcf_filename_mod.no.header.vcf

	plink --vcf $smr_dir/$tissue/$vcf_filename_mod.vcf --out $smr_dir/$tissue/${tissue}_vcf_plink_format
	
	date
}

do_preprocessing_gwas_results_for_SMR (){
	date
	
	gwas_summary_filename_mod=`echo $gwas_summary_filename | sed s/.tsv//g`

	#Step A - liftOver UK Biobank summary GWAS coordinates from hg19 to hg38
	cat $gwas_summary_folder/$gwas_summary_filename | cut -f 1 | sed 's/:/\t/g' | awk '{print "chr"$1"\t"$2"\t"$2}' | sed '1,1d' > $gwas_summary_folder/${gwas_summary_filename_mod}.only.coordinates.hg37.bed
	crossmap bed $script_dir/GRCh37_to_GRCh38.chain.gz $gwas_summary_folder/${gwas_summary_filename_mod}.only.coordinates.hg37.bed > $gwas_summary_folder/${gwas_summary_filename_mod}.only.coordinates.hg38.bed
	
	#Step B - do some reformatting and add the hg38 coordinates to the UK Biobank summary GWAS data
	cat $gwas_summary_folder/${gwas_summary_filename_mod}.only.coordinates.hg38.bed | awk '{print $1"_"$2"\t"$5"_"$6}' | sed 's/chr//g' > $gwas_summary_folder/tmp.txt
	cat $gwas_summary_folder/tmp.txt | awk 'BEGIN{print "variant_hg37""\t""variant_hg38"}; {print $0}' > $gwas_summary_folder/tmp1.txt
	cat $gwas_summary_folder/$gwas_summary_filename | cut -f 1 | sed 's/:/\t/g' | awk '{print $1"_"$2}' > $gwas_summary_folder/tmp2.txt
	paste $gwas_summary_folder/tmp2.txt $gwas_summary_folder/$gwas_summary_filename > $gwas_summary_folder/tmp3.txt

	#Step C - take the dbSNP150_hg38 file from Annovar with some reformatting already done (ask Komudi) and reformat is again so it could be matched with the UK Biobank summary GWAS data from Step B above, run only once
#	cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/hg38_avsnp150_SNPonlyV3.txt | sed 's/_/\t/g' | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5}' | sed 's/chr//g' > /data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/hg38_avsnp150_SNPonlyV4.txt

	echo "Merging with dbSNP150 to retrieve rsIDs..."
	echo "..."
	echo "..."
	#Step D - merge with reformatted dbSNP150_hg38 (from #Step C above) to finally create the GWAS .ma file format for SMR analysis [https://cnsgenomics.com/software/smr/#SMR&HEIDIanalysis]
	python $script_dir/SMR_update_GWAS_from_hg19_to_hg38.py "$gwas_summary_folder/tmp3.txt" "$gwas_summary_folder/tmp1.txt" "$script_dir/hg38_avsnp150_SNPonlyV4.txt" $gwas_summary_folder $gwas_summary_filename_mod

	cat $gwas_summary_folder/${gwas_summary_filename_mod}_hg38.tsv | awk '{print $NF"\t"$2"\t"$4"\t"$10"\t"$11"\t"$13"\t""NA"}' | sed 's/:/\t/g' | cut -f1,4-5,6- | sed '1,1d' | sort -k1 | uniq | awk 'BEGIN{print "SNP\tA1\tA2\tfreq\tb\tse\tp\tn"}; {print $0}' > $gwas_summary_folder/${gwas_summary_filename_mod}_hg38.ma
	cat $gwas_summary_folder/${gwas_summary_filename_mod}_hg38.ma | awk '{if(NF<8){}else{print $0}}' > $gwas_summary_folder/${gwas_summary_filename_mod}_hg38_temp.ma
	mv -f $gwas_summary_folder/${gwas_summary_filename_mod}_hg38_temp.ma $gwas_summary_folder/${gwas_summary_filename_mod}_hg38.ma

	rm -f $gwas_summary_folder/tmp.txt
	rm -f $gwas_summary_folder/tmp1.txt
	rm -f $gwas_summary_folder/tmp2.txt
	rm -f $gwas_summary_folder/tmp3.txt

	date
}

do_SMR (){
	date

	gwas_summary_filename_mod=`echo $gwas_summary_filename | sed s/.tsv//g`
	gwas_summary_filename_ukbiobank_code=`echo $gwas_summary_filename | awk -F"." '{print $1}'`
	smr_output_filename=`echo $gwas_summary_folder | awk -F"/" '{print $(NF-1)}'`

	smr --bfile $smr_dir/$tissue/${tissue}_vcf_plink_format \
	--gwas-summary $gwas_summary_folder/${gwas_summary_filename_mod}_hg38.ma \
	--beqtl-summary $smr_dir/$tissue/${tissue}_threaded_subset_fdrbh_less_than_0.05_eQTL.nochr \
	--peqtl-smr 5.0e-4 \
	--thread-num $threads \
	--out $smr_dir/$tissue/00-smr_out/${gwas_summary_filename_ukbiobank_code}_${smr_output_filename}_hg38.smr.out.txt

	date
}

# do_preprocessing_eQTL_results_for_SMR_by_tissue
# do_preprocessing_gwas_results_for_SMR
do_SMR


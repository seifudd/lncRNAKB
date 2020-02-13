#!/bin/bash

function do_smr (){
    vcf_filename=$1
    eqtl_results_filename=$2
    tissue=$3
    gwas_summary_filename=$4
    gwas_summary_folder=$5
    eqtl_output_filepath="/data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded-FS/$tissue/$eqtl_results_filename"
    vcf_output_filepath="/data/NHLBI_BCB/Fayaz/21-GTEx-data/dbGaP-13516/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/vcf_tissue_subset_FS/$tissue/$vcf_filename"
    smr_dir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/00-SMR-lncRNAKB-UKBiobank"
    script_dir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK"
    threads=8

    mkdir -p $smr_dir/$tissue
    mkdir -p $smr_dir/$tissue/00-smr_out

    echo "TISSUE=$tissue,$gwas_summary_filename"
    error_filename=`echo $gwas_summary_folder | awk -F"/" '{print $(NF-1)}'`
    sbatch --job-name="${tissue}.${gwas_summary_filename}.smr" \
    	   --partition=norm \
           --time=5:00:00 \
           --mem=100g \
           --gres=lscratch:100 \
           --cpus-per-task=$threads \
           --error="$smr_dir/$tissue/${error_filename}.stderr" \
           --output="$smr_dir/$tissue/${error_filename}.stdout" \
	$script_dir/run_smr.sh \
        $script_dir \
        $tissue \
        $vcf_output_filepath \
        $eqtl_output_filepath \
        $gwas_summary_filename \
        $gwas_summary_folder \
        $smr_dir \
        $threads
}

let "n = 0"
while read vcf_filename eqtl_results_filename tissue
do
	while read gwas_summary_filename gwas_summary_folder
	do
#		echo $tissue,$gwas_summary_filename,$gwas_summary_folder
		let "n++"
#		echo "initial_value",$n
     		do_smr $vcf_filename $eqtl_results_filename $tissue $gwas_summary_filename $gwas_summary_folder
		#do_smr $vcf_filename $eqtl_results_filename $tissue
		#do_smr $gwas_summary_filename $gwas_summary_folder
		#do_smr $vcf_filename $eqtl_results_filename $tissue $gwas_summary_filename $gwas_summary_folder
		#do_smr "20117_0.gwas.imputed_v3.both_sexes.tsv" "/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/Alcohol_drinker_status_Never/"
		#do_smr
		if [ $n -gt 1500 ]
		then
#			echo "waiting_value",$n
			echo "Waiting for 60 min..."
			echo "..."
			echo "..."
			echo "..."
			sleep 60m
			let "n = 0"
		fi
#		echo "reset_value",$n
	done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/list_of_gwas_filenames_and_folderpaths_for_smr.txt"
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/list_of_vcf_eqtl_filenames_for_smr.txt"

######
# vcf_filename eqtl_results_filename tissue
# /data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/list_of_vcf_eqtl_filenames_for_smr.txt
# gwas_summary_filename gwas_summary_folder
# /data/NHLBI_BCB/Fayaz/21-GTEx-data/12-SMR-UKBIOBANK/list_of_gwas_filenames_and_folderpaths_for_smr.txt


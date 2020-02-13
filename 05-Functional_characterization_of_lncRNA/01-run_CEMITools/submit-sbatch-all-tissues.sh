#!/bin/bash
#$ -cwd

function do_cemitools () {
	scriptdir="/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo"
	tissue=$1
	cpmexpressionfile=$2
	numcpus=8
	echo "TISSUE=$tissue"
	mkdir $scriptdir/$tissue
	sbatch  --job-name="$tissue" \
		--partition=norm \
		--time=8:00:00 \
		--mem=200g \
		--gres=lscratch:150 \
		--cpus-per-task=$numcpus \
		--error="$scriptdir/$tissue.slurm.err.txt" \
		--output="$scriptdir/$tissue.slurm.out.txt" \
		$scriptdir/sbatch-runcemitools.sh $scriptdir $tissue $cpmexpressionfile
}

# do_cemitools Adipose_Tissue	Adipose_Tissue.363samples.15175pcs.11854lincs.gene.CPM.GTExId.txt

while read tissue cpmexpression; do
   	do_cemitools ${tissue} $cpmexpression
done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/cpm_filteredbycount_filteredbytpm_file_list_for_cemitools.txt"

# while read tissue cpmexpression; do
#	k=`cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/$tissue/Modules_Identified.csv | grep -v "ModuleName" | grep -v "Not.Correlated" | wc -l`
#	echo $tissue,$k
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/cpm_filteredbycount_filteredbytpm_file_list_for_cemitools_copy.txt"

# while read tissue cpmexpression; do
#	k=`cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/$tissue/Modules_Identified.csv | grep -v "ModuleName" | grep -v "Not.Correlated" | awk -F"," '{print $3}' | sed 's/"//g' | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }'`
#	echo $tissue,$k
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/cpm_filteredbycount_filteredbytpm_file_list_for_cemitools_copy.txt"

# while read tissue cpmexpression; do
#	k=`cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/$tissue/Modules_Identified.csv | grep -v "ModuleName" | grep -v "Not.Correlated" | awk -F"," '{print $2}' | sed 's/"//g' | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }'`
#	echo $tissue,$k
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/cpm_filteredbycount_filteredbytpm_file_list_for_cemitools_copy.txt"

# while read tissue cpmexpression; do
#	k=`cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/$tissue/Modules_Identified_Enrichment_all_pathways_with_module_info.csv | grep -v "ModuleName" | grep -v "Not.Correlated" | wc -l`
#	echo $tissue,$k
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/cpm_filteredbycount_filteredbytpm_file_list_for_cemitools_copy.txt"

# while read tissue cpmexpression; do
#	k=`cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/$tissue/Modules_Identified_Enrichment_all_pathways_with_module_info.csv | grep -v "ModuleName" | grep -v "Not.Correlated" | awk -F"," '{if($7<=0.05)print $0}' | wc -l`
#	echo $tissue,$k
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/cpm_filteredbycount_filteredbytpm_file_list_for_cemitools_copy.txt"

# while read tissue cpmexpression; do
#	k=`cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/$tissue/Modules_Identified_Enrichment_all_pathways_with_module_info.csv | grep -v "ModuleName" | grep -v "Not.Correlated" | awk -F"," '{if($9<=0.05)print $0}' | wc -l`
#	echo $tissue,$k
# done < "/data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools/cpm_filteredbycount_filteredbytpm_file_list_for_cemitools_copy.txt"



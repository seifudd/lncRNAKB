#!/bin/bash

set -e
function fail {
echo “FAIL: $@“ >&2
exit 1 # signal failure
}


tissue=$1
dir=$2

function do_merge {
		for i in `ls /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/$tissue | grep "info" | grep "top3"`; do
	echo $i
	cat /data/NHLBI_BCB/Fayaz/21-GTEx-data/10-CEMITools-redo/$tissue/$i | sed '1d' | grep -v "Not.Correlated" | awk -F"," -v OFS="\t" '{print $0}' > $dir/$tissue.tmp.txt
	
	done 
#cat *.tmp >> $dir/Top3_pathways_all_tissues-2.txt
#rm *.tmp
date
}
do_merge
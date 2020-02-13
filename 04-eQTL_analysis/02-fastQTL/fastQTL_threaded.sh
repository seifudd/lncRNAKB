#!/bin/bash

ml singularity
ml python

export SINGULARITY_BINDPATH="/data/NHLBI_BCB"

inputdir=$1
vcfinputdir=$2
expressiondatafile=$3
tissue=$4
window=$5
permutations=$6
chunks=$7
threads=$8
outputdir=$9

expressiondatafilename=`echo $expressiondatafile | sed 's/.txt//g'`
covariatedatafilename=`echo $expressiondatafile | cut -d"." -f1,2`
vcffilename=`echo $covariatedatafilename | sed 's/samples//'`

singularity exec /data/NHLBI_BCB/Fayaz/21-GTEx-data/09-fastQTL_threaded/gtex_eqtl_V8.sif /opt/fastqtl/python/run_FastQTL_threaded.py --covariates  $inputdir/$tissue/$covariatedatafilename.covariate.sorted_v2.txt \
					--window $window \
					--chunks $chunks \
					--threads $threads \
					-o $outputdir/$tissue/ \
$vcfinputdir/$tissue/$vcffilename.subset.GTEx.subset.MAF0.05.FINAL.sorted.vcf.gz $inputdir/$tissue/$expressiondatafilename.sorted.bed.gz $covariatedatafilename.subset.GTEx.subset.MAF0.05.FINAL.win.$window.bp.perm.$permutations.cis.nominal.threaded.eqtl

echo -e "
	usage: run_FastQTL_threaded.py [-h] [--covariates COVARIATES]
                               [--phenotype_groups PHENOTYPE_GROUPS]
                               [--chunks CHUNKS]
                               [--permute PERMUTE [PERMUTE ...]]
                               [--interaction INTERACTION]
                               [--best_variant_only] [--window WINDOW]
                               [--threshold THRESHOLD]
                               [--maf_threshold MAF_THRESHOLD]
                               [--ma_sample_threshold MA_SAMPLE_THRESHOLD]
                               [--interaction_maf_threshold INTERACTION_MAF_THRESHOLD]
                               [--fdr FDR] [--seed SEED]
                               [--exclude_samples EXCLUDE_SAMPLES]
                               [--exclude_sites EXCLUDE_SITES]
                               [--qvalue_lambda QVALUE_LAMBDA] [-t THREADS]
                               [-o OUTPUT_DIR]
                               vcf bed prefix

	Run FastQTL

	positional arguments:
	  vcf                   Genotypes in VCF 4.1 format
	  bed                   Phenotypes in UCSC BED extended format
	  prefix                Prefix for output file name

	optional arguments:
	  -h, --help            show this help message and exit
	  --covariates COVARIATES
		                Covariates
	  --phenotype_groups PHENOTYPE_GROUPS
		                File with mapping of phenotype_id to group_id
		                (gene_id)
	  --chunks CHUNKS       Number of chunks, minimum: #chromosomes
	  --permute PERMUTE [PERMUTE ...]
		                Number of permutations, e.g. [1000, 10000] (adaptive).
		                Default: None (run nominal pass)
	  --interaction INTERACTION
		                Interaction term
	  --best_variant_only
	  --window WINDOW       Cis-window size. Default values is 1Mb (1e6).
	  --threshold THRESHOLD
		                Output only significant phenotype-variant pairs with a
		                p-value below threshold (default 1)
	  --maf_threshold MAF_THRESHOLD
		                Include only genotypes with minor allele frequency
		                >=maf_threshold (default 0)
	  --ma_sample_threshold MA_SAMPLE_THRESHOLD
		                Include only genotypes with >=ma_sample_threshold
		                samples carrying the minor allele (default 0)
	  --interaction_maf_threshold INTERACTION_MAF_THRESHOLD
		                MAF threshold for interactions, applied to lower and
		                upper half of samples
	  --fdr FDR
	  --seed SEED           Random number generator seed
	  --exclude_samples EXCLUDE_SAMPLES
	  --exclude_sites EXCLUDE_SITES
	  --qvalue_lambda QVALUE_LAMBDA
		                lambda parameter for pi0est in qvalue.
	  -t THREADS, --threads THREADS
		                Number of threads
	  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
		                Output directory

" > /dev/null


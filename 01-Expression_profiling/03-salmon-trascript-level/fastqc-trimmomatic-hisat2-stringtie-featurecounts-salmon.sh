#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

SAMPLE=$1
DATAPATH=$2
READ1=$3
READ2=$4
REF=$5
outdir=$6
numcpus=$7
transcript_ref=$8
tissue=$9
tissue_R1_R2_filelist=${10}

module load salmon
# module load Software/hisat2/2.0.5
# hisat2="/home/seifuddinft/miniconda3/bin/hisat2" #version 2.1.0
# module load Software/samtools/0.1.18
# samtools="/home/seifuddinft/miniconda3/bin/samtools" #Version: 1.9 (using htslib 1.9)
# module load Software/Trimmomatic/0.36
# trimmomatic="/home/seifuddinft/miniconda3/bin/trimmomatic" #0.38
# module load Software/FastQC/0.11.5
# fastqc="/home/seifuddinft/miniconda3/bin/fastqc" #FastQC v0.11.8
# module load stringtie
# stringtie="/home/seifuddinft/miniconda3/bin/stringtie" #stringtie-1.3.4
# module load subread
# featureCounts="/home/seifuddinft/miniconda3/bin/featureCounts" #featureCounts v1.6.3
# salmon="/home/seifuddinft/miniconda3/bin/salmon" #salmon v0.12.0

function do_fastqc () {
	date
	########################################################################################################################
	# if trimming, change DATAPATH
#	DATAPATH="/lscratch/${SLURM_JOBID}"

	# if trimming, change $READ1 and $READ2
#	READ1="${SAMPLE}_1P.fastq.gz"
#	READ2="${SAMPLE}_2P.fastq.gz"

	# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
	out_dir="fastqc"
#	out_dir="fastqc_post_trimming"

#	fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

	mkdir -p $outdir/$SAMPLE/$out_dir

	$fastqc -o "$outdir/$SAMPLE/$out_dir"  \
	--nogroup \
	"$DATAPATH/$READ1"  \
	"$DATAPATH/$READ2"  \
	|| fail "fastqc failed"

	echo "fastqc done"
	########################################################################################################################
	date
}

function do_trimmomatic () {
	date
	########################################################################################################################

	out_dir="trimmedfastqs"
	mkdir -p $outdir/$SAMPLE/$out_dir

	java -jar $trimmomatic PE \
	            -threads $numcpus \
	            "$DATAPATH/$READ1" \
	            "$DATAPATH/$READ2" \
	 	    -baseout "$outdir/$SAMPLE/$out_dir/${SAMPLE}.fastq.gz" \
	           ILLUMINACLIP:"/home/seifuddinft/miniconda3/bin/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa":2:30:10 \
	           MINLEN:50

	echo "trimmomatic done"
	########################################################################################################################
	date
}

function do_hisat2 () {
	date
	########################################################################################################################
	# if trimming, please change DATAPATH
#	DATAPATH="/lscratch/${SLURM_JOBID}"

	# if trimming, please change $READ1 and $READ2
#	READ1="${SAMPLE}_1P.fastq.gz"
#	READ2="${SAMPLE}_2P.fastq.gz"

	out_dir="hisat2"
	mkdir -p $outdir/$SAMPLE/$out_dir

	$hisat2 	-p $numcpus \
			-x $REF/genome_tran \
			--downstream-transcriptome-assembly \
			-1 "$DATAPATH/$READ1" \
			-2 "$DATAPATH/$READ2" \
			--rg-id $SAMPLE --rg SM:$SAMPLE \
	| $samtools view -h -f 3 -O SAM - \
	| perl -nle  'print if m/^@(?:[A-Z]{2})\s|\bNH:i:1\b/' \
	| $samtools sort -@ $numcpus \
		-o "$outdir/$SAMPLE/$out_dir/$SAMPLE.unique.bam"
#		-T /lscratch/${SLURM_JOB_ID}/${SAMPLE}_chunk -
	   
	$samtools index "$outdir/$SAMPLE/$out_dir/$SAMPLE.unique.bam"
	echo "hisat2 done"
	########################################################################################################################
	date
}

function do_stringtie () {
	date
	########################################################################################################################
	GTF="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/GENCODE_human_v28/gencode.v28.annotation.gtf"
	bamfile="$outdir/$SAMPLE/hisat2/$SAMPLE.unique.bam"
	out_dir="stringtie"
	mkdir -p $outdir/$SAMPLE/$out_dir

	# only reference
	stringtie -p 4 \
	  -o $outdir/$SAMPLE/$out_dir/$SAMPLE.gtf \
	  -e -G $GTF \
	  -l $SAMPLE \
	  -v -B -c 10 -j 5  -f 0.1 \
	  $bamfile

	# novel
	# stringtie -p 4 \
	#  -o $outdir/$SAMPLE/$out_dir/$SAMPLE.gtf \
	#  -G $GTF \
	#  -l $SAMPLE \
	#  -v -B -c 10 -j 5  -f 0.1 \
	#  $bamfile

	# https://github.com/gpertea/stringtie
	echo -e "
	--version : print current version at stdout
	 -h print this usage message
	 -G reference annotation to use for guiding the assembly process (GTF/GFF3)
	 -l name prefix for output transcripts (default: STRG)
	 -f minimum isoform fraction (default: 0.1)
	 -m minimum assembled transcript length to report (default 100bp)
	 -o output path/file name for the assembled transcripts GTF (default: stdout)
	 -a minimum anchor length for junctions (default: 10)
	 -j minimum junction coverage (default: 1)
	 -t disable trimming of predicted transcripts based on coverage
	    (default: coverage trimming is enabled)
	 -c minimum reads per bp coverage to consider for transcript assembly (default: 2.5)
	 -v verbose (log bundle processing details)
	 -g gap between read mappings triggering a new bundle (default: 50)
	 -C output file with reference transcripts that are covered by reads
	 -M fraction of bundle allowed to be covered by multi-hit reads (default:0.95)
	 -p number of threads (CPUs) to use (default: 1)
	 -A gene abundance estimation output file name
	 -B enable output of Ballgown table files which will be created in the
	    same directory as the output GTF (requires -G, -o recommended)
	 -b enable output of Ballgown table files but these files will be 
	    created under the directory path given as <dir_path>
	 -e only estimates the abundance of given reference transcripts (requires -G)
	 -x do not assemble any transcripts on the given reference sequence(s)

	Transcript merge usage mode:

	 stringtie --merge [Options] { gtf_list | strg1.gtf ...}
	With this option StringTie will assemble transcripts from multiple
	input files generating a unified non-redundant set of isoforms. In this
	usage mode the following options are available:
	  -G <guide_gff>   reference annotation to include in the merging (GTF/GFF3)
	  -o <out_gtf>     output file name for the merged transcripts GTF
		            (default: stdout)
	  -m <min_len>     minimum input transcript length to include in the merge
		            (default: 50)
	  -c <min_cov>     minimum input transcript coverage to include in the merge
		            (default: 0)
	  -F <min_fpkm>    minimum input transcript FPKM to include in the merge
		            (default: 1.0)
	  -T <min_tpm>     minimum input transcript TPM to include in the merge
		            (default: 1.0)
	  -f <min_iso>     minimum isoform fraction (default: 0.01)
	  -g <gap_len>     gap between transcripts to merge together (default: 250)
	  -i               keep merged transcripts with retained introns; by default
		           these are not kept unless there is strong evidence for them
	  -l <label>       name prefix for output transcripts (default: MSTRG)
	" > /dev/null
	echo "stringtie done"
	########################################################################################################################
	date
}

function do_featurecounts () {
	date
	########################################################################################################################
	GTF="/home/seifuddinft/gencode.v29.annotation.gtf"
	bamfile="$outdir/$SAMPLE/hisat2/$SAMPLE.unique.bam" # only using uniquely mapped BAM to count, if using nonunique BAM please check -M option for featureCounts

#	out_dir="featurecounts"
	out_dir="featurecounts_exon"
	mkdir -p $outdir/$SAMPLE/$out_dir

	s=2  # -s strand-specific : 0 (unstranded), 1 (stranded) 2 (reversely stranded). 0 default.
		
	end_num="" # if single end sequencing, please change to "single"
	feature="" # if not gene, please leave "blank"
	
	if [[ $feature == "gene" ]]
		then
			if [[ $end_num == "single" ]]
			then

				$featureCounts -T $numcpus \
					-t exon \
					-g gene_id \
					-a $GTF \
					-s $s \
					-o $outdir/$SAMPLE/$out_dir/$SAMPLE.genefeatureCounts.txt \
					$bamfile
			else

				$featureCounts -T $numcpus \
					-t exon \
					-g gene_id \
					-p \
					-a $GTF \
					-s $s \
					-O \
					-o $outdir/$SAMPLE/$out_dir/$SAMPLE.genefeatureCounts.txt \
					$bamfile
			fi
		
	else
			if [[ $end_num == "single" ]]
			then

				$featureCounts -T $numcpus \
					-t exon \
					-g exon_id \
					-f \
					-O \
					-a $GTF \
					-s $s \
					-o $outdir/$SAMPLE/$out_dir/$SAMPLE.exonfeatureCounts.txt \
					$bamfile

			else
				$featureCounts -T $numcpus \
					-t exon \
					-g exon_id \
					-f \
					-a $GTF \
					-s $s \
					-p \
					-O \
					-J \
					-o $outdir/$SAMPLE/$out_dir/$SAMPLE.exonfeatureCounts.txt \
					$bamfile
			fi
	fi

	echo -e "
	Version 1.6.3

	Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 

	## Mandatory arguments:

	  -a <string>         Name of an annotation file. GTF/GFF format by default. See
		              -F option for more format information. Inbuilt annotations
		              (SAF format) is available in 'annotation' directory of the
		              package. Gzipped file is also accepted.

	  -o <string>         Name of the output file including read counts. A separate
		              file including summary statistics of counting results is
		              also included in the output ('<string>.summary')

	  input_file1 [input_file2] ...   A list of SAM or BAM format files. They can be
		              either name or location sorted. If no files provided,
		              <stdin> input is expected. Location-sorted paired-end reads
		              are automatically sorted by read names.

	## Optional arguments:
	# Annotation

	  -F <string>         Specify format of the provided annotation file. Acceptable
		              formats include 'GTF' (or compatible GFF format) and
		              'SAF'. 'GTF' by default.  For SAF format, please refer to
		              Users Guide.

	  -t <string>         Specify feature type in GTF annotation. 'exon' by 
		              default. Features used for read counting will be 
		              extracted from annotation using the provided value.

	  -g <string>         Specify attribute type in GTF annotation. 'gene_id' by 
		              default. Meta-features used for read counting will be 
		              extracted from annotation using the provided value.

	  --extraAttributes   Extract extra attribute types from the provided GTF
		              annotation and include them in the counting output. These
		              attribute types will not be used to group features. If
		              more than one attribute type is provided they should be
		              separated by comma.

	  -A <string>         Provide a chromosome name alias file to match chr names in
		              annotation with those in the reads. This should be a two-
		              column comma-delimited text file. Its first column should
		              include chr names in the annotation and its second column
		              should include chr names in the reads. Chr names are case
		              sensitive. No column header should be included in the
		              file.

	# Level of summarization

	  -f                  Perform read counting at feature level (eg. counting 
		              reads for exons rather than genes).

	# Overlap between reads and features

	  -O                  Assign reads to all their overlapping meta-features (or 
		              features if -f is specified).

	  --minOverlap <int>  Minimum number of overlapping bases in a read that is
		              required for read assignment. 1 by default. Number of
		              overlapping bases is counted from both reads if paired
		              end. If a negative value is provided, then a gap of up
		              to specified size will be allowed between read and the
		              feature that the read is assigned to.

	  --fracOverlap <float> Minimum fraction of overlapping bases in a read that is
		              required for read assignment. Value should be within range
		              [0,1]. 0 by default. Number of overlapping bases is
		              counted from both reads if paired end. Both this option
		              and '--minOverlap' option need to be satisfied for read
		              assignment.

	  --fracOverlapFeature <float> Minimum fraction of overlapping bases in a
		              feature that is required for read assignment. Value
		              should be within range [0,1]. 0 by default.

	  --largestOverlap    Assign reads to a meta-feature/feature that has the 
		              largest number of overlapping bases.

	  --nonOverlap <int>  Maximum number of non-overlapping bases in a read (or a
		              read pair) that is allowed when being assigned to a
		              feature. No limit is set by default.

	  --nonOverlapFeature <int> Maximum number of non-overlapping bases in a feature
		              that is allowed in read assignment. No limit is set by
		              default.

	  --readExtension5 <int> Reads are extended upstream by <int> bases from their
		              5' end.

	  --readExtension3 <int> Reads are extended upstream by <int> bases from their
		              3' end.

	  --read2pos <5:3>    Reduce reads to their 5' most base or 3' most base. Read
		              counting is then performed based on the single base the 
		              read is reduced to.

	# Multi-mapping reads

	  -M                  Multi-mapping reads will also be counted. For a multi-
		              mapping read, all its reported alignments will be 
		              counted. The 'NH' tag in BAM/SAM input is used to detect 
		              multi-mapping reads.

	# Fractional counting

	  --fraction          Assign fractional counts to features. This option must
		              be used together with '-M' or '-O' or both. When '-M' is
		              specified, each reported alignment from a multi-mapping
		              read (identified via 'NH' tag) will carry a fractional
		              count of 1/x, instead of 1 (one), where x is the total
		              number of alignments reported for the same read. When '-O'
		              is specified, each overlapping feature will receive a
		              fractional count of 1/y, where y is the total number of
		              features overlapping with the read. When both '-M' and
		              '-O' are specified, each alignment will carry a fractional
		              count of 1/(x*y).

	# Read filtering

	  -Q <int>            The minimum mapping quality score a read must satisfy in
		              order to be counted. For paired-end reads, at least one
		              end should satisfy this criteria. 0 by default.

	  --splitOnly         Count split alignments only (ie. alignments with CIGAR
		              string containing 'N'). An example of split alignments is
		              exon-spanning reads in RNA-seq data.

	  --nonSplitOnly      If specified, only non-split alignments (CIGAR strings do
		              not contain letter 'N') will be counted. All the other
		              alignments will be ignored.

	  --primary           Count primary alignments only. Primary alignments are 
		              identified using bit 0x100 in SAM/BAM FLAG field.

	  --ignoreDup         Ignore duplicate reads in read counting. Duplicate reads 
		              are identified using bit Ox400 in BAM/SAM FLAG field. The 
		              whole read pair is ignored if one of the reads is a 
		              duplicate read for paired end data.

	# Strandness

	  -s <int or string>  Perform strand-specific read counting. A single integer
		              value (applied to all input files) or a string of comma-
		              separated values (applied to each corresponding input
		              file) should be provided. Possible values include:
		              0 (unstranded), 1 (stranded) and 2 (reversely stranded).
		              Default value is 0 (ie. unstranded read counting carried
		              out for all input files).

	# Exon-exon junctions

	  -J                  Count number of reads supporting each exon-exon junction.
		              Junctions were identified from those exon-spanning reads
		              in the input (containing 'N' in CIGAR string). Counting
		              results are saved to a file named '<output_file>.jcounts'

	  -G <string>         Provide the name of a FASTA-format file that contains the
		              reference sequences used in read mapping that produced the
		              provided SAM/BAM files. This optional argument can be used
		              with '-J' option to improve read counting for junctions.

	# Parameters specific to paired end reads

	  -p                  If specified, fragments (or templates) will be counted
		              instead of reads. This option is only applicable for
		              paired-end reads.

	  -B                  Only count read pairs that have both ends aligned.

	  -P                  Check validity of paired-end distance when counting read 
		              pairs. Use -d and -D to set thresholds.

	  -d <int>            Minimum fragment/template length, 50 by default.

	  -D <int>            Maximum fragment/template length, 600 by default.

	  -C                  Do not count read pairs that have their two ends mapping 
		              to different chromosomes or mapping to same chromosome 
		              but on different strands.

	  --donotsort         Do not sort reads in BAM/SAM input. Note that reads from 
		              the same pair are required to be located next to each 
		              other in the input.

	# Number of CPU threads

	  -T <int>            Number of the threads. 1 by default.

	# Read groups

	  --byReadGroup       Assign reads by read group. RG tag is required to be
		              present in the input BAM/SAM files.
		              

	# Long reads

	  -L                  Count long reads such as Nanopore and PacBio reads. Long
		              read counting can only run in one thread and only reads
		              (not read-pairs) can be counted. There is no limitation on
		              the number of 'M' operations allowed in a CIGAR string in
		              long read counting.

	# Assignment results for each read

	  -R <format>         Output detailed assignment results for each read or read-
		              pair. Results are saved to a file that is in one of the
		              following formats: CORE, SAM and BAM. See Users Guide for
		              more info about these formats.

	  --Rpath <string>    Specify a directory to save the detailed assignment
		              results. If unspecified, the directory where counting
		              results are saved is used.

	# Miscellaneous

	  --tmpDir <string>   Directory under which intermediate files are saved (later
		              removed). By default, intermediate files will be saved to
		              the directory specified in '-o' argument.

	  --maxMOp <int>      Maximum number of 'M' operations allowed in a CIGAR
		              string. 10 by default. Both 'X' and '=' are treated as 'M'
		              and adjacent 'M' operations are merged in the CIGAR
		              string.

	  --verbose           Output verbose information for debugging, such as un-
		              matched chromosome/contig names.

	  -v                  Output version of the program.
	"> /dev/null
	echo "featurecounts done"
	########################################################################################################################
	date
}

function do_salmon_index () {

	date
	########################################################################################################################

### Get lncRNAKB transcript FASTA file - Run only once

#	cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.gtf | awk '{if($3=="transcript") print $1"\t"$4"\t"$5}' > /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp1
#	cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.gtf | awk '{if($3=="transcript") print $0}' | awk 'match($0,/transcript_id\s+["]lnckb[.][0-9]+[.][0-9]+["]/) {print substr($0,RSTART,RLENGTH)}' | sed 's/transcript_id //g' | sed 's/["]//g' > /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp2
#	cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.gtf | awk '{if($3=="transcript") print $0}' | awk 'match($0,/gene_id\s+["]lnckb[.][0-9]+["]/) {print substr($0,RSTART,RLENGTH)}' | sed 's/gene_id //g' | sed 's/["]//g' > /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp3
#	paste -d "|" /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp2 /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp3 > /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp4
#	paste /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp1 /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp4 > /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.bed

#	rm -f /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp1
#	rm -f /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp2
#	rm -f /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp3
#	rm -f /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/temp4

#	cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.bed | grep -v "_" | cut -c4- > /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.primary_assembly.bed
#	cat /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.bed | grep "_" > /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.patches.bed

#	bedtools getfasta -fi /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.primary_assembly.bed -name -fo /data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.transcript.primary.assembly.fa

### Get lncRNAKB transcript FASTA file - Run only once

	salmon index --type quasi \
		      -p $numcpus \
                      -k 31 \
                      --gencode \
		      -t $transcript_ref \
		      -i /data/NHLBI_BCB/Fayaz/21-GTEx-data/05-salmon/salmon_transcript_index_hg38_lncRNAKB_v7
	
	echo -e "
	Version Info: This is the most recent version of salmon.

	Index
	==========
	Creates a salmon index.

	Command Line Options:
	  -v [ --version ]           print version string
	  -h [ --help ]              produce help message
	  -t [ --transcripts ] arg   Transcript fasta file.
	  -k [ --kmerLen ] arg (=31) The size of k-mers that should be used for the 
		                     quasi index.
	  -i [ --index ] arg         salmon index.
	  --gencode                  This flag will expect the input transcript fasta 
		                     to be in GENCODE format, and will split the 
		                     transcript name at the first '|' character.  These
		                     reduced names will be used in the output and when 
		                     looking for these transcripts in a gene to 
		                     transcript GTF.
	  --keepDuplicates           This flag will disable the default indexing 
		                     behavior of discarding sequence-identical 
		                     duplicate transcripts.  If this flag is passed, 
		                     then duplicate transcripts that appear in the 
		                     input will be retained and quantified separately.
	  -p [ --threads ] arg (=2)  Number of threads to use (only used for computing 
		                     bias features)
	  --perfectHash              [quasi index only] Build the index using a perfect
		                     hash rather than a dense hash.  This will require 
		                     less memory (especially during quantification), 
		                     but will take longer to construct
	  --type arg (=quasi)        The type of index to build; the only option is 
		                     quasi in this version of salmon.
	"> /dev/null
	########################################################################################################################
	date      
}

function do_salmon_quant () {
	GTF="/data/NHLBI_BCB/Fayaz/24-lncRNA-database-curation/10-modify_lncRNAKB_GTF/lncRNAKB_hg38_v7.gtf"

	# if trimming, please change DATAPATH
	# DATAPATH="/lscratch/${SLURM_JOBID}"

	# if trimming, please change $READ1 and $READ2
#	READ1="${SAMPLE}_1P.fastq.gz"
#	READ2="${SAMPLE}_2P.fastq.gz"

# salmon index path=/data/NHLBI_BCB/Fayaz/21-GTEx-data/05-salmon/salmon_transcript_index_hg38_lncRNAKB_v7
	
	out_dir="salmon"
	mkdir -p $outdir/$tissue/$SAMPLE/$out_dir

	salmon quant -i  /scratch/seifuddinft/salmon_transcript_index_hg38_lncRNAKB_v7 \
		      -l A \
		      -1 "$DATAPATH/$READ1" \
		      -2 "$DATAPATH/$READ2" \
		      -p $numcpus \
		      -g $GTF \
		      --validateMappings \
	              -o $outdir/$tissue/$SAMPLE/$out_dir/${SAMPLE}.salmon

	echo -e"
	Fragment Library Types=https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype
	Quant
	==========
	Perform dual-phase, mapping-based estimation of
	transcript abundance from RNA-seq reads

	salmon quant options:


	mapping input options:
	  -l [ --libType ] arg                  Format string describing the library 
		                                type
	  -i [ --index ] arg                    salmon index
	  -r [ --unmatedReads ] arg             List of files containing unmated reads 
		                                of (e.g. single-end reads)
	  -1 [ --mates1 ] arg                   File containing the #1 mates
	  -2 [ --mates2 ] arg                   File containing the #2 mates


	basic options:
	  -v [ --version ]                      print version string
	  -h [ --help ]                         produce help message
	  -o [ --output ] arg                   Output quantification directory.
	  --seqBias                             Perform sequence-specific bias 
		                                correction.
	  --gcBias                              [beta for single-end reads] Perform 
		                                fragment GC bias correction
	  -p [ --threads ] arg (=12)            The number of threads to use 
		                                concurrently.
	  --incompatPrior arg (=0)              This option sets the prior probability 
		                                that an alignment that disagrees with 
		                                the specified library type (--libType) 
		                                results from the true fragment origin. 
		                                Setting this to 0 specifies that 
		                                alignments that disagree with the 
		                                library type should be impossible, 
		                                while setting it to 1 says that 
		                                alignments that disagree with the 
		                                library type are no less likely than 
		                                those that do
	  -g [ --geneMap ] arg                  File containing a mapping of 
		                                transcripts to genes.  If this file is 
		                                provided salmon will output both 
		                                quant.sf and quant.genes.sf files, 
		                                where the latter contains aggregated 
		                                gene-level abundance estimates.  The 
		                                transcript to gene mapping should be 
		                                provided as either a GTF file, or a in 
		                                a simple tab-delimited format where 
		                                each line contains the name of a 
		                                transcript and the gene to which it 
		                                belongs separated by a tab.  The 
		                                extension of the file is used to 
		                                determine how the file should be 
		                                parsed.  Files ending in '.gtf', '.gff'
		                                or '.gff3' are assumed to be in GTF 
		                                format; files with any other extension 
		                                are assumed to be in the simple format.
		                                In GTF / GFF format, the 
		                                transcript_id is assumed to contain 
		                                the transcript identifier and the 
		                                gene_id is assumed to contain the 
		                                corresponding gene identifier.
	  --meta                                If you're using Salmon on a metagenomic
		                                dataset, consider setting this flag to 
		                                disable parts of the abundance 
		                                estimation model that make less sense 
		                                for metagenomic data.


	options specific to mapping mode:
	  --discardOrphansQuasi                 [Quasi-mapping mode only] : Discard 
		                                orphan mappings in quasi-mapping mode. 
		                                If this flag is passed then only paired
		                                mappings will be considered toward 
		                                quantification estimates.  The default 
		                                behavior is to consider orphan mappings
		                                if no valid paired mappings exist.  
		                                This flag is independent of the option 
		                                to write the orphaned mappings to file 
		                                (--writeOrphanLinks).
	  --validateMappings                    [Quasi-mapping mode only] : Validate 
		                                mappings using alignment-based 
		                                verifcation. If this flag is passed, 
		                                quasi-mappings will be validated to 
		                                ensure that they could give rise to a 
		                                reasonable alignment before they are 
		                                further used for quantification.
	  --consensusSlack arg (=0)             [Quasi-mapping mode only] : The amount 
		                                of slack allowed in the quasi-mapping 
		                                consensus mechanism.  Normally, a 
		                                transcript must cover all hits to be 
		                                considered for mapping.  If this is set
		                                to a value, X, greater than 0, then a 
		                                transcript can fail to cover up to X 
		                                hits before it is discounted as a 
		                                mapping candidate.  The default value 
		                                of this option is 1 if 
		                                --validateMappings is given and 0 
		                                otherwise.
	  --minScoreFraction arg                [Quasi-mapping mode (w / mapping 
		                                validation) only] : The fraction of the
		                                optimal possible alignment score that a
		                                mapping must achieve in order to be 
		                                considered valid --- should be in 
		                                (0,1].
		                                Salmon Default 0.65 and Alevin Default 
		                                0.8
	  --maxMMPExtension arg (=7)            [Quasi-mapping mode (w / mapping 
		                                validation) only] : Sets the maximum 
		                                allowable MMP extension when collecting
		                                suffix array intervals to be used in 
		                                chaining.  This prevents MMPs from 
		                                becoming too long, and potentially 
		                                masking intervals that would uncover 
		                                other good quasi-mappings for the read.
		                                  This heuristic mimics the idea of the
		                                maximum mappable safe prefix (MMSP) in 
		                                selective alignment.  Setting a smaller
		                                value will potentially allow for more 
		                                sensitive, but slower, mapping.
	  --ma arg (=2)                         [Quasi-mapping mode (w / mapping 
		                                validation) only] : The value given to 
		                                a match between read and reference 
		                                nucleotides in an alignment.
	  --mp arg (=-4)                        [Quasi-mapping mode (w / mapping 
		                                validation) only] : The value given to 
		                                a mis-match between read and reference 
		                                nucleotides in an alignment.
	  --go arg (=5)                         [Quasi-mapping mode (w / mapping 
		                                validation) only] : The value given to 
		                                a gap opening in an alignment.
	  --ge arg (=3)                         [Quasi-mapping mode (w / mapping 
		                                validation) only] : The value given to 
		                                a gap extension in an alignment.
	  --mimicStrictBT2                      [Quasi-mapping mode (w / mapping 
		                                validation) only] : Set flags to mimic 
		                                the very strict parameters used by 
		                                RSEM+Bowtie2.  This increases 
		                                --minScoreFraction to 0.8, disallows 
		                                dovetailing reads, discards orphans, 
		                                and disallows gaps in alignments.
	  --allowOrphansFMD                     [FMD-mapping mode only] : Consider 
		                                orphaned reads as valid hits when 
		                                performing lightweight-alignment.  This
		                                option will increase sensitivity (allow
		                                more reads to map and more transcripts 
		                                to be detected), but may decrease 
		                                specificity as orphaned alignments are 
		                                more likely to be spurious.
	  -z [ --writeMappings ] [=arg(=-)]     If this option is provided, then the 
		                                quasi-mapping results will be written 
		                                out in SAM-compatible format.  By 
		                                default, output will be directed to 
		                                stdout, but an alternative file name 
		                                can be provided instead.
	  -c [ --consistentHits ]               Force hits gathered during 
		                                quasi-mapping to be consistent (i.e. 
		                                co-linear and approximately the right 
		                                distance apart).
	  --strictIntersect                     Modifies how orphans are assigned.  
		                                When this flag is set, if the 
		                                intersection of the quasi-mappings for 
		                                the left and right is empty, then all 
		                                mappings for the left and all mappings 
		                                for the right read are reported as 
		                                orphaned quasi-mappings
	  --fasterMapping                       [Developer]: Disables some extra checks
		                                during quasi-mapping. This may make 
		                                mapping a little bit faster at the 
		                                potential cost of returning too many 
		                                mappings (i.e. some sub-optimal 
		                                mappings) for certain reads.
	  -x [ --quasiCoverage ] arg (=0)       [Experimental]: The fraction of the 
		                                read that must be covered by MMPs (of 
		                                length >= 31) if this read is to be 
		                                considered as mapped.  This may help 
		                                to avoid spurious mappings. A value 
		                                of 0 (the default) denotes no coverage 
		                                threshold (a single 31-mer can yield a 
		                                mapping).  Since coverage by exact 
		                                matching, large, MMPs is a rather 
		                                strict condition, this value should 
		                                likely be set to something low, if 
		                                used.


	advanced options:
	  --alternativeInitMode                 [Experimental]: Use an alternative 
		                                strategy (rather than simple 
		                                interpolation between) the online and 
		                                uniform abundance estimates to 
		                                initalize the EM / VBEM algorithm.
	  --auxDir arg (=aux_info)              The sub-directory of the quantification
		                                directory where auxiliary information 
		                                e.g. bootstraps, bias parameters, etc. 
		                                will be written.
	  --dumpEq                              Dump the equivalence class counts that 
		                                were computed during quasi-mapping
	  -d [ --dumpEqWeights ]                Includes rich equivlance class 
		                                weights in the output when equivalence 
		                                class information is being dumped to 
		                                file.
	  --minAssignedFrags arg (=10)          The minimum number of fragments that 
		                                must be assigned to the transcriptome 
		                                for quantification to proceed.
	  --reduceGCMemory                      If this option is selected, a more 
		                                memory efficient (but slightly slower) 
		                                representation is used to compute 
		                                fragment GC content. Enabling this will
		                                reduce memory usage, but can also 
		                                reduce speed.  However, the results 
		                                themselves will remain the same.
	  --biasSpeedSamp arg (=5)              The value at which the fragment length 
		                                PMF is down-sampled when evaluating 
		                                sequence-specific & GC fragment bias.  
		                                Larger values speed up effective length
		                                correction, but may decrease the 
		                                fidelity of bias modeling results.
	  --fldMax arg (=1000)                  The maximum fragment length to consider
		                                when building the empirical 
		                                distribution
	  --fldMean arg (=250)                  The mean used in the fragment length 
		                                distribution prior
	  --fldSD arg (=25)                     The standard deviation used in the 
		                                fragment length distribution prior
	  -f [ --forgettingFactor ] arg (=0.65000000000000002)
		                                The forgetting factor used in the 
		                                online learning schedule.  A smaller 
		                                value results in quicker learning, but 
		                                higher variance and may be unstable.  A
		                                larger value results in slower learning
		                                but may be more stable.  Value should 
		                                be in the interval (0.5, 1.0].
	  --initUniform                         initialize the offline inference with 
		                                uniform parameters, rather than seeding
		                                with online parameters.
	  -w [ --maxReadOcc ] arg (=200)        Reads mapping to more than this many 
		                                places won't be considered.
	  --noLengthCorrection                  [experimental] : Entirely disables 
		                                length correction when estimating the 
		                                abundance of transcripts.  This option 
		                                can be used with protocols where one 
		                                expects that fragments derive from 
		                                their underlying targets without regard
		                                to that target's length (e.g. QuantSeq)
	  --noEffectiveLengthCorrection         Disables effective length correction 
		                                when computing the probability that a 
		                                fragment was generated from a 
		                                transcript.  If this flag is passed in,
		                                the fragment length distribution is not
		                                taken into account when computing this 
		                                probability.
	  --noFragLengthDist                    [experimental] : Don't consider 
		                                concordance with the learned fragment 
		                                length distribution when trying to 
		                                determine the probability that a 
		                                fragment has originated from a 
		                                specified location.  Normally, 
		                                Fragments with unlikely lengths will be
		                                assigned a smaller relative probability
		                                than those with more likely lengths.  
		                                When this flag is passed in, the 
		                                observed fragment length has no effect 
		                                on that fragment's a priori 
		                                probability.
	  --noBiasLengthThreshold               [experimental] : If this option is 
		                                enabled, then no (lower) threshold will
		                                be set on how short bias correction can
		                                make effective lengths. This can 
		                                increase the precision of bias 
		                                correction, but harm robustness.  The 
		                                default correction applies a threshold.
	  --numBiasSamples arg (=2000000)       Number of fragment mappings to use when
		                                learning the sequence-specific bias 
		                                model.
	  --numAuxModelSamples arg (=5000000)   The first <numAuxModelSamples> are used
		                                to train the auxiliary model parameters
		                                (e.g. fragment length distribution, 
		                                bias, etc.).  After ther first 
		                                <numAuxModelSamples> observations the 
		                                auxiliary model parameters will be 
		                                assumed to have converged and will be 
		                                fixed.
	  --numPreAuxModelSamples arg (=1000000)
		                                The first <numPreAuxModelSamples> will 
		                                have their assignment likelihoods and 
		                                contributions to the transcript 
		                                abundances computed without applying 
		                                any auxiliary models.  The purpose of 
		                                ignoring the auxiliary models for the 
		                                first <numPreAuxModelSamples> 
		                                observations is to avoid applying these
		                                models before thier parameters have 
		                                been learned sufficiently well.
	  --useEM                               Use the traditional EM algorithm for 
		                                optimization in the batch passes.
	  --useVBOpt                            Use the Variational Bayesian EM 
		                                [default]
	  --rangeFactorizationBins arg (=0)     Factorizes the likelihood used in 
		                                quantification by adopting a new notion
		                                of equivalence classes based on the 
		                                conditional probabilities with which 
		                                fragments are generated from different 
		                                transcripts.  This is a more 
		                                fine-grained factorization than the 
		                                normal rich equivalence classes.  The 
		                                default value (0) corresponds to the 
		                                standard rich equivalence classes, and 
		                                larger values imply a more fine-grained
		                                factorization.  If range factorization 
		                                is enabled, a common value to select 
		                                for this parameter is 4.
	  --numGibbsSamples arg (=0)            Number of Gibbs sampling rounds to 
		                                perform.
	  --noGammaDraw                         This switch will disable drawing 
		                                transcript fractions from a Gamma 
		                                distribution during Gibbs sampling.  In
		                                this case the sampler does not account 
		                                for shot-noise, but only assignment 
		                                ambiguity
	  --numBootstraps arg (=0)              Number of bootstrap samples to 
		                                generate. Note: This is mutually 
		                                exclusive with Gibbs sampling.
	  --bootstrapReproject                  This switch will learn the parameter 
		                                distribution from the bootstrapped 
		                                counts for each sample, but will 
		                                reproject those parameters onto the 
		                                original equivalence class counts.
	  --thinningFactor arg (=16)            Number of steps to discard for every 
		                                sample kept from the Gibbs chain. The 
		                                larger this number, the less chance 
		                                that subsequent samples are 
		                                auto-correlated, but the slower 
		                                sampling becomes.
	  -q [ --quiet ]                        Be quiet while doing quantification 
		                                (don't write informative output to the 
		                                console unless something goes wrong).
	  --perTranscriptPrior                  The prior (either the default or the 
		                                argument provided via --vbPrior) will 
		                                be interpreted as a transcript-level 
		                                prior (i.e. each transcript will be 
		                                given a prior read count of this value)
	  --sigDigits arg (=3)                  The number of significant digits to 
		                                write when outputting the 
		                                EffectiveLength and NumReads columns
	  --vbPrior arg (=1.0000000000000001e-05)
		                                The prior that will be used in the VBEM
		                                algorithm.  This is interpreted as a 
		                                per-nucleotide prior, unless the 
		                                --perTranscriptPrior flag is also 
		                                given, in which case this is used as a 
		                                transcript-level prior
	  --writeOrphanLinks                    Write the transcripts that are linked 
		                                by orphaned reads.
	  --writeUnmappedNames                  Write the names of un-mapped reads to 
		                                the file unmapped_names.txt in the 
		                                auxiliary directory.
	"> /dev/null
}

function do_salmon_quant_alignment_based_mode (){

	echo -e"

		Version Info: Could not resolve upgrade information in the alotted time.
		Check for upgrades manually at https://combine-lab.github.io/salmon

		Quant
		==========
		Perform dual-phase, alignment-based estimation of
		transcript abundance from RNA-seq reads

		salmon quant options:


		alignment input options:
		  --discardOrphans                      [Alignment-based mode only] : Discard 
				                        orphan alignments in the input .  If 
				                        this flag is passed, then only paired 
				                        alignments will be considered toward 
				                        quantification estimates.  The default 
				                        behavior is to consider orphan 
				                        alignments if no valid paired mappings 
				                        exist.
		  -l [ --libType ] arg                  Format string describing the library 
				                        type
		  -a [ --alignments ] arg               input alignment (BAM) file(s).
		  -t [ --targets ] arg                  FASTA format file containing target 
				                        transcripts.


		basic options:
		  -v [ --version ]                      print version string
		  -h [ --help ]                         produce help message
		  -o [ --output ] arg                   Output quantification directory.
		  --seqBias                             Perform sequence-specific bias 
				                        correction.
		  --gcBias                              [beta for single-end reads] Perform 
				                        fragment GC bias correction
		  -p [ --threads ] arg (=8)             The number of threads to use 
				                        concurrently.
		  --incompatPrior arg (=0)              This option sets the prior probability 
				                        that an alignment that disagrees with 
				                        the specified library type (--libType) 
				                        results from the true fragment origin. 
				                        Setting this to 0 specifies that 
				                        alignments that disagree with the 
				                        library type should be impossible, 
				                        while setting it to 1 says that 
				                        alignments that disagree with the 
				                        library type are no less likely than 
				                        those that do
		  -g [ --geneMap ] arg                  File containing a mapping of 
				                        transcripts to genes.  If this file is 
				                        provided salmon will output both 
				                        quant.sf and quant.genes.sf files, 
				                        where the latter contains aggregated 
				                        gene-level abundance estimates.  The 
				                        transcript to gene mapping should be 
				                        provided as either a GTF file, or a in 
				                        a simple tab-delimited format where 
				                        each line contains the name of a 
				                        transcript and the gene to which it 
				                        belongs separated by a tab.  The 
				                        extension of the file is used to 
				                        determine how the file should be 
				                        parsed.  Files ending in '.gtf', '.gff'
				                        or '.gff3' are assumed to be in GTF 
				                        format; files with any other extension 
				                        are assumed to be in the simple format.
				                        In GTF / GFF format, the 
				                        transcript_id is assumed to contain 
				                        the transcript identifier and the 
				                        gene_id is assumed to contain the 
				                        corresponding gene identifier.
		  --meta                                If you're using Salmon on a metagenomic
				                        dataset, consider setting this flag to 
				                        disable parts of the abundance 
				                        estimation model that make less sense 
				                        for metagenomic data.


		alignment-specific options:
		  --noErrorModel                        Turn off the alignment error model, 
				                        which takes into account the the 
				                        observed frequency of different types 
				                        of mismatches / indels when computing 
				                        the likelihood of a given alignment. 
				                        Turning this off can speed up 
				                        alignment-based salmon, but can harm 
				                        quantification accuracy.
		  --numErrorBins arg (=6)               The number of bins into which to divide
				                        each read when learning and applying 
				                        the error model.  For example, a value 
				                        of 10 would mean that effectively, a 
				                        separate error model is leared and 
				                        applied to each 10th of the read, while
				                        a value of 3 would mean that a separate
				                        error model is applied to the read 
				                        beginning (first third), middle (second
				                        third) and end (final third).
		  -s [ --sampleOut ]                    Write a postSample.bam file in the 
				                        output directory that will sample the 
				                        input alignments according to the 
				                        estimated transcript abundances. If 
				                        you're going to perform downstream 
				                        analysis of the alignments with tools 
				                        which don't, themselves, take fragment 
				                        assignment ambiguity into account, you 
				                        should use this output.
		  -u [ --sampleUnaligned ]              In addition to sampling the aligned 
				                        reads, also write the un-aligned reads 
				                        to postSample.bam.
		  --gencode                             This flag will expect the input 
				                        transcript fasta to be in GENCODE 
				                        format, and will split the transcript 
				                        name at the first '|' character.  These
				                        reduced names will be used in the 
				                        output and when looking for these 
				                        transcripts in a gene to transcript 
				                        GTF.
		  --mappingCacheMemoryLimit arg (=2000000)
				                        If the file contained fewer than this 
				                        many mapped reads, then just keep the 
				                        data in memory for subsequent rounds of
				                        inference. Obviously, this value should
				                        not be too large if you wish to keep a 
				                        low memory usage, but setting it large 
				                        enough to accommodate all of the mapped
				                        read can substantially speed up 
				                        inference on small files that contain
				                        only a few million reads.


		advanced options:
		  --alternativeInitMode                 [Experimental]: Use an alternative 
				                        strategy (rather than simple 
				                        interpolation between) the online and 
				                        uniform abundance estimates to 
				                        initalize the EM / VBEM algorithm.
		  --auxDir arg (=aux_info)              The sub-directory of the quantification
				                        directory where auxiliary information 
				                        e.g. bootstraps, bias parameters, etc. 
				                        will be written.
		  --skipQuant                           Skip performing the actual transcript 
				                        quantification (including any Gibbs 
				                        sampling or bootstrapping).
		  --dumpEq                              Dump the equivalence class counts that 
				                        were computed during quasi-mapping
		  -d [ --dumpEqWeights ]                Includes rich equivlance class 
				                        weights in the output when equivalence 
				                        class information is being dumped to 
				                        file.
		  --minAssignedFrags arg (=10)          The minimum number of fragments that 
				                        must be assigned to the transcriptome 
				                        for quantification to proceed.
		  --reduceGCMemory                      If this option is selected, a more 
				                        memory efficient (but slightly slower) 
				                        representation is used to compute 
				                        fragment GC content. Enabling this will
				                        reduce memory usage, but can also 
				                        reduce speed.  However, the results 
				                        themselves will remain the same.
		  --biasSpeedSamp arg (=5)              The value at which the fragment length 
				                        PMF is down-sampled when evaluating 
				                        sequence-specific & GC fragment bias.  
				                        Larger values speed up effective length
				                        correction, but may decrease the 
				                        fidelity of bias modeling results.
		  --fldMax arg (=1000)                  The maximum fragment length to consider
				                        when building the empirical 
				                        distribution
		  --fldMean arg (=250)                  The mean used in the fragment length 
				                        distribution prior
		  --fldSD arg (=25)                     The standard deviation used in the 
				                        fragment length distribution prior
		  -f [ --forgettingFactor ] arg (=0.65000000000000002)
				                        The forgetting factor used in the 
				                        online learning schedule.  A smaller 
				                        value results in quicker learning, but 
				                        higher variance and may be unstable.  A
				                        larger value results in slower learning
				                        but may be more stable.  Value should 
				                        be in the interval (0.5, 1.0].
		  --initUniform                         initialize the offline inference with 
				                        uniform parameters, rather than seeding
				                        with online parameters.
		  -w [ --maxReadOcc ] arg (=200)        Reads mapping to more than this many 
				                        places won't be considered.
		  --noLengthCorrection                  [experimental] : Entirely disables 
				                        length correction when estimating the 
				                        abundance of transcripts.  This option 
				                        can be used with protocols where one 
				                        expects that fragments derive from 
				                        their underlying targets without regard
				                        to that target's length (e.g. QuantSeq)
		  --noEffectiveLengthCorrection         Disables effective length correction 
				                        when computing the probability that a 
				                        fragment was generated from a 
				                        transcript.  If this flag is passed in,
				                        the fragment length distribution is not
				                        taken into account when computing this 
				                        probability.
		  --noFragLengthDist                    [experimental] : Don't consider 
				                        concordance with the learned fragment 
				                        length distribution when trying to 
				                        determine the probability that a 
				                        fragment has originated from a 
				                        specified location.  Normally, 
				                        Fragments with unlikely lengths will be
				                        assigned a smaller relative probability
				                        than those with more likely lengths.  
				                        When this flag is passed in, the 
				                        observed fragment length has no effect 
				                        on that fragment's a priori 
				                        probability.
		  --noBiasLengthThreshold               [experimental] : If this option is 
				                        enabled, then no (lower) threshold will
				                        be set on how short bias correction can
				                        make effective lengths. This can 
				                        increase the precision of bias 
				                        correction, but harm robustness.  The 
				                        default correction applies a threshold.
		  --numBiasSamples arg (=2000000)       Number of fragment mappings to use when
				                        learning the sequence-specific bias 
				                        model.
		  --numAuxModelSamples arg (=5000000)   The first <numAuxModelSamples> are used
				                        to train the auxiliary model parameters
				                        (e.g. fragment length distribution, 
				                        bias, etc.).  After ther first 
				                        <numAuxModelSamples> observations the 
				                        auxiliary model parameters will be 
				                        assumed to have converged and will be 
				                        fixed.
		  --numPreAuxModelSamples arg (=5000)   The first <numPreAuxModelSamples> will 
				                        have their assignment likelihoods and 
				                        contributions to the transcript 
				                        abundances computed without applying 
				                        any auxiliary models.  The purpose of 
				                        ignoring the auxiliary models for the 
				                        first <numPreAuxModelSamples> 
				                        observations is to avoid applying these
				                        models before thier parameters have 
				                        been learned sufficiently well.
		  --useEM                               Use the traditional EM algorithm for 
				                        optimization in the batch passes.
		  --useVBOpt                            Use the Variational Bayesian EM 
				                        [default]
		  --rangeFactorizationBins arg (=0)     Factorizes the likelihood used in 
				                        quantification by adopting a new notion
				                        of equivalence classes based on the 
				                        conditional probabilities with which 
				                        fragments are generated from different 
				                        transcripts.  This is a more 
				                        fine-grained factorization than the 
				                        normal rich equivalence classes.  The 
				                        default value (0) corresponds to the 
				                        standard rich equivalence classes, and 
				                        larger values imply a more fine-grained
				                        factorization.  If range factorization 
				                        is enabled, a common value to select 
				                        for this parameter is 4.
		  --numGibbsSamples arg (=0)            Number of Gibbs sampling rounds to 
				                        perform.
		  --noGammaDraw                         This switch will disable drawing 
				                        transcript fractions from a Gamma 
				                        distribution during Gibbs sampling.  In
				                        this case the sampler does not account 
				                        for shot-noise, but only assignment 
				                        ambiguity
		  --numBootstraps arg (=0)              Number of bootstrap samples to 
				                        generate. Note: This is mutually 
				                        exclusive with Gibbs sampling.
		  --bootstrapReproject                  This switch will learn the parameter 
				                        distribution from the bootstrapped 
				                        counts for each sample, but will 
				                        reproject those parameters onto the 
				                        original equivalence class counts.
		  --thinningFactor arg (=16)            Number of steps to discard for every 
				                        sample kept from the Gibbs chain. The 
				                        larger this number, the less chance 
				                        that subsequent samples are 
				                        auto-correlated, but the slower 
				                        sampling becomes.
		  -q [ --quiet ]                        Be quiet while doing quantification 
				                        (don't write informative output to the 
				                        console unless something goes wrong).
		  --perTranscriptPrior                  The prior (either the default or the 
				                        argument provided via --vbPrior) will 
				                        be interpreted as a transcript-level 
				                        prior (i.e. each transcript will be 
				                        given a prior read count of this value)
		  --sigDigits arg (=3)                  The number of significant digits to 
				                        write when outputting the 
				                        EffectiveLength and NumReads columns
		  --vbPrior arg (=1.0000000000000001e-05)
				                        The prior that will be used in the VBEM
				                        algorithm.  This is interpreted as a 
				                        per-nucleotide prior, unless the 
				                        --perTranscriptPrior flag is also 
				                        given, in which case this is used as a 
				                        transcript-level prior
		  --writeOrphanLinks                    Write the transcripts that are linked 
				                        by orphaned reads.
		  --writeUnmappedNames                  Write the names of un-mapped reads to 
				                        the file unmapped_names.txt in the 
				                        auxiliary directory. 
	">/dev/null

}


function do_dexseq_htseq_counting_reads () {

	date
	########################################################################################################################
#	GTF="/home/seifuddinft/gencode.v29.annotation.gtf"

#	python /home/seifuddinft/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py /home/seifuddinft/gencode.v29.annotation.gtf dexseq.gencode.v29.annotation.gff	
	GTF="/mnt/lab-thein/Fayaz/03-MNS/00-DEXSEq-MNS/dexseq.gencode.v29.annotation.v2.without.chr.gff"
	outdir="/mnt/lab-thein/Fayaz/00-fastqc-trimmomatic-hisat2-stringtie-featurecounts"
	bamfile="$outdir/$SAMPLE/hisat2/$SAMPLE.unique.bam" # only using uniquely mapped BAM to count, if using nonunique BAM please check -M option for featureCounts
		
	python /home/seifuddinft/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -r pos -s reverse -p yes -f bam $GTF $bamfile /mnt/lab-thein/Fayaz/03-MNS/00-DEXSEq-MNS/00-DEXSeq-HTSeq-counts/$SAMPLE/$SAMPLE.dexseq.exon.counts.txt
	
	echo -e"
	There are a number of crucial points to pay attention to when using the python_count.py script:

	Paired-end data: If your data is from a paired-end sequencing run, you need to add the option  -p yes to the command to call the script. (As usual, options have to be placed before the file names, surrounded by spaces.) In addition, the SAM file needs to be sorted, either by read name or by position. Most aligners produce sorted SAM files; if your SAM file is not sorted, use  samtools sort -n to sort by read name (or samtools sort) to sort by position. (See Anders et al. (2013), if you need further explanations on how to sort SAM files.) Use the option -r pos or  -r name to indicate whether your paired-end data is sorted by alignment position or by read name.

	Strandedness: By default, the counting script assumes your library to be strand-specific, i.e., reads are aligned to the same strand as the gene they originate from. If you have used a library preparation protocol that does not preserve strand information (i.e., reads from a given gene can appear equally likely on either strand), you need to inform the script by specifying the option  -s no. If your library preparation protocol reverses the strand (i.e., reads appear on the strand opposite to their gene of origin), use -s reverse. In case of paired-end data, the default (-s yes) means that the read from the first sequence pass is on the same strand as the gene and the read from the second pass on the opposite strand (“forward-reverse” or “fr” order in the parlance of the Bowtie/TopHat manual) and the options -s reverse specifies the opposite case.

	SAM and BAM files: By default, the script expects its input to be in plain-text SAM format. However, it can also read BAM files, i.e., files in the the compressed binary variant of the SAM format. If you wish to do so, use the option -f bam. This works only if you have installed the Python package pysam.

	Alignment quality: The scripts takes a further option, -a to specify the minimum alignment quality (as given in the fifth column of the SAM file). All reads with a lower quality than specified (with default -a 10) are skipped.

	Help pages: Calling either script without arguments displays a help page with an overview of all options and arguments.
	"> /dev/null
	########################################################################################################################
	date

}

function do_dexseq_counts_merge () {
	date
	########################################################################################################################
        sids="/mnt/lab-thein/Fayaz/03-MNS/mns-patients-rbc-system-antigens-history-alloimunization_v2.txt"
	feature="exon"
#	out_dir="featurecounts"
	out_dir="/mnt/lab-thein/Fayaz/03-MNS/00-DEXSEq-MNS/00-DEXSeq-HTSeq-counts"

	for i in `cat $sids | head -1`; do
           SAMPLE=$i
	   echo -e "$feature\t$i" > tmp1
	   cat $out_dir/$SAMPLE/$i.dexseq.$feature.counts.txt | sort -k1,1 >> tmp1
	done

	for i in `cat $sids | sed '1,1d'`; do
           SAMPLE=$i
	   echo -e "$i" > tmp2
	   cat $out_dir/$SAMPLE/$i.dexseq.$feature.counts.txt | sort -k1,1 | cut -f2 >> tmp2
	   paste tmp1 tmp2 > tmp3
	   mv -f tmp3 tmp1
	done

	mv -f tmp1  $feature.dexseq.counts.txt
	rm -f tmp2
	rm -f tmp3
	########################################################################################################################
	date
}

function do_featurecounts_merge () {
	date
	########################################################################################################################
        sids="/mnt/lab-thein/Fayaz/sids.all.txt"
	feature="exon"
#	out_dir="featurecounts"
	out_dir="featurecounts_exon"

	for i in `cat $sids | cut -f1 | head -1`; do
           SAMPLE=$i
	   echo -e "$feature\t$i" > tmp1
	   cat $SAMPLE/$out_dir/$i.${feature}featureCounts.txt | sed '1,2d' | sort -k1,1 | cut -f 1,7 >> tmp1
	done

	for i in `cat $sids | cut -f1 | sed '1,1d' `; do
           SAMPLE=$i
	   echo -e "$i" > tmp2
	   cat $SAMPLE/$out_dir/$i.${feature}featureCounts.txt | sed '1,2d' | sort -k1,1 | cut -f 7 >> tmp2
	   paste tmp1 tmp2 > tmp3
	   mv -f tmp3 tmp1
	done

	mv -f tmp1  $feature.featurecount.txt
	rm -f tmp2
	rm -f tmp3
	########################################################################################################################
	date
}

function do_salmon_merge () {
	date
	########################################################################################################################
        sids="/data/NHLBI_BCB/Fayaz/21-GTEx-data/01-by-tissue-files-and-ids-rnaseq-only/${tissue_R1_R2_filelist}"
	feature="transcript"
	typeofcount="tpm"

	out_dir="$outdir/$tissue"
	
	for i in `cat $sids | cut -d" " -f1 | head -1`; do
           SAMPLE=$i
	   echo -e "$feature\t$i" > $out_dir/tmp1
	   cat $out_dir/$SAMPLE/salmon/$SAMPLE.salmon/quant.sf | sed '1,1d' | sort -k1,1 | cut -f 1,4 >> $out_dir/tmp1
	done

	for i in `cat $sids | cut -d" " -f1 | sed '1,1d'`; do
           SAMPLE=$i
	   echo -e "$i" > $out_dir/tmp2
	   cat $out_dir/$SAMPLE/salmon/$SAMPLE.salmon/quant.sf | sed '1,1d' | sort -k1,1 | cut -f 4 >> $out_dir/tmp2
	   paste $out_dir/tmp1 $out_dir/tmp2 > $out_dir/tmp3
	   mv -f $out_dir/tmp3 $out_dir/tmp1
	done

	mv -f $out_dir/tmp1  $out_dir/$tissue.$feature.salmon.$typeofcount.txt
	rm -f $out_dir/tmp2
	rm -f $out_dir/tmp3
	########################################################################################################################
	date
}

function do_get_hisat2_stats () {
	date
	########################################################################################################################
	# change log files directory to absolute path of logfiles 
#	logfiles="/mnt/lab-thein/Fayaz/00-fastqc-trimmomatic-hisat2-stringtie-featurecounts/01-set1-log-files-hisat2"	
#	logfiles="/mnt/lab-thein/Fayaz/00-fastqc-trimmomatic-hisat2-stringtie-featurecounts/01-set2-log-files-hisat2"
	logfiles="/mnt/lab-thein/Fayaz/00-fastqc-trimmomatic-hisat2-stringtie-featurecounts/01-set3-log-files-hisat2"

	cd $logfiles

	ls *.slurm.err.txt  | cut -f 1 -d'.' > tmp1
	grep 'reads; of these' *.slurm.err.txt | awk '{print $1}' | cut -f 2 -d':' | awk '{ printf("%'"'"'d\n",$1); }' > tmp2
	grep 'overall alignment rate' *.slurm.err.txt | awk '{print $1}' | cut -f 2 -d':' > tmp3
	grep 'aligned concordantly exactly 1 time' *.slurm.err.txt | awk '{print $2}' | awk '{ printf("%'"'"'d\n",$1); }' > tmp4
	grep 'aligned concordantly exactly 1 time' *.slurm.err.txt | awk '{print $3}' | sed 's/(//g' | sed 's/)//g' > tmp5
	echo -e "Sample\tTotal_Reads\tOverall_Alignment_Rate\tUniq_alignment\tUniq_alignment_%" > alignment.stats.txt
	paste tmp1 tmp2 tmp3 tmp4 tmp5 >> alignment.stats.txt
	rm -f tmp*
	a1=`cat alignment.stats.txt | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $2; n++ } END { print int(sum / n); }' | awk '{ printf("%'"'"'d\n",$1); }'`
	a2=`cat alignment.stats.txt | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $3; n++ } END { print sum / n; }' | awk -F. '{print $1"."substr($2,1,2)"%"}'` 
	a3=`cat alignment.stats.txt | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $4; n++ } END { print int(sum / n); }' | awk '{ printf("%'"'"'d\n",$1); }'`
	a4=`cat alignment.stats.txt | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $5; n++ } END { print sum / n; }' | awk -F. '{print $1"."substr($2,1,2)"%"}' `
	echo -e "Average\t$a1\t$a2\t$a3\t$a4" >> alignment.stats.txt
	########################################################################################################################
	date
}

function do_get_featurecount_stats () {
	date
	########################################################################################################################
 	sids="/mnt/lab-thein/Fayaz/sids.all.txt"
	feature="exon"	
	
	# change log files directory to absolute path of logfiles 
	logfiles="/mnt/lab-thein/Fayaz/00-fastqc-trimmomatic-hisat2-stringtie-featurecounts/03-set1-set2-log-files-featurecounts-$feature"	

	cd $logfiles
	
	echo -e "Sample\t${feature}_Count\tTotal_fragments\tSuccessfully_assigned\t%" > $feature.featurecount.summary.txt

	for i in `cat $sids | cut -f1`; do
	    count_of_features=`cat $i.slurm.err.txt | grep 'Meta-features : ' | awk '{print $4}'`
	   total=`cat $i.slurm.err.txt | grep 'Total alignments : ' | awk '{print $5}'`
	       s=`cat $i.slurm.err.txt | grep 'Successfully assigned alignments : ' | awk '{print $6"\t"$7}'`
	  echo -e "$i\t$count_of_features\t$total\t$s" >> $feature.featurecount.summary.txt
	done
	########################################################################################################################
	date
}

# do_fastqc
# do_trimmomatic
# do_fastqc
# do_hisat2
# do_featurecounts
# do_dexseq_htseq_counting_reads
# do_dexseq_counts_merge
# do_featurecounts_merge
do_salmon_merge
# do_stringtie
# do_salmon_index
# do_salmon_quant
# do_get_hisat2_stats
# do_get_featurecount_stats


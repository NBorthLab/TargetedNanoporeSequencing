#!/bin/bash

# This script filters reads according to the strand and their start/end mapping position. 

# $1 scaffold
# $2 startposition of ROI (e.g. junction site)
# $3 input file: alignments.stats, output of stats_from_bam
# $4 merged.fastq file 

read1=$(expr $2 + 100) # Only reads are considered with a starting/ending position within a range of +- 100bp
read2=$(expr $2 - 100) # 100 bp to compensate for mapping uncertainty

echo "Filter reads according to their aligning position and directionality"

awk -v var="$1" -v fwd1="$read1" -v fwd2="$read2" '{if (($10=="+") && (var==$2 && $7<fwd1 && $7>fwd2)) print $1}' $3 | sort -u > junc_reads_fwd_$1\:$2.txt
# The variable "strandbias" is the endposition of the alignment
awk -v var="$1" -v rev1="$read1" -v rev2="$read2" '{strandbias= $7+$9}{if (($10=="-") && (var==$2 && strandbias<rev1 && strandbias>rev2)) print $1}' $3 | sort -u > junc_reads_rev_$1\:$2.txt

echo "Preparing fastq file of on-target reads - might take a some min"

	while read file; do
		file_name1=junc_reads_fwd_$1\:$2.txt
		grep -A 3 "$file" $4 >> ontarget_${file_name1}.fastq
	done<junc_reads_fwd_$1\:$2.txt

	while read file; do
		file_name2=junc_reads_rev_$1\:$2.txt
                grep -A 3 "$file" $4 >> ontarget_${file_name2}.fastq
        done<junc_reads_rev_$1\:$2.txt

if [[ ! -s junc_reads_fwd_$1\:$2.txt ]]
        then
        rm -rf junc_reads_fwd_$1\:$2.txt

fi
if [[ ! -s junc_reads_rev_$1\:$2.txt ]]
        then
        rm -rf junc_reads_rev_$1\:$2.txt

fi
echo "done"




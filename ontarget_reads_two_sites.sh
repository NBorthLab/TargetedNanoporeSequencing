#!/bin/bash

# ROI is here defined as the innermost (or smallest) region between gRNAs flanking a genomic region

# $1= stats_from_bam file
# $2= startpostion of innermost gRNA
# $3= endpostion of innermost gRNA
# $4= scaffold
# $5= merge.fastq file

echo "Filtering reads only mapping to ROI"

start_position1=$(expr $2 - 3000) # The "tolarance range" can be adjusted, especially when multiple gRNAs targeting the same site 
start_position2=$(expr $3 + 100)  # and are further than 3000bp apart from each other.
end_position1=$(expr $3 + 3000)   
end_position2=$(expr $2 - 100)
scaffold=$4

awk -v var5="$scaffold" -v var1="$start_position1" -v var2="$start_position2" -v var3="$end_position1" -v var4="$end_position2" 'BEGIN{OFS="\t"}$2==var5{endpos= $7+$9}{ if ( $7>var1 && $7<var2 && endpos<var3 && endpos>var4 ) print $1}' $1 > reads_${4}_${2}_${3}.txt 

echo "Preparing fastq file with on-target reads - can take a few minutes"
while read file; do 
 grep -A 3 "$file" $5 >> ontarget_reads_${4}_${2}_${3}.fastq
done<reads_${4}_${2}_${3}.txt


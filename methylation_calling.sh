#!/bin/bash

# Before this script can be used, you have to index the fastq reads: nanpolish index -d /path/to/fast5_folder/ on-target_reads.fastq

# $1= ref assembly
# $2= fastq file with ontarget reads
# $3= output name, e.g. IS_A_cd4
# $4= coordinates (Scaffold:start-end), e.g.: IS_A_cd4:1-4000

pathtoexe=/path/to/executables

# only primary alignments are allowed
$pathtoexe/minimap2 -a -x map-ont -t 30 $1 $2 | $pathtoexe/samtools view -T $1 -F 2304 -@ 30 -bS - | $pathtoexe/samtools sort -@ 30 -l 9 -o $3.sorted.no_sup_sec_ali.bam
$pathtoexe/samtools index $3.sorted.no_sup_sec_ali.bam

#call methylation
$pathtoexe/nanopolish call-methylation -t 30 -r $2 -b $3.sorted.no_sup_sec_ali.bam -g $1 -w "$4" > methylation_calls$3.tsv
$pathtoexe/calculate_methylation_frequency.py -s methylation_calls$3.tsv > methylation_frequency$3.tsv


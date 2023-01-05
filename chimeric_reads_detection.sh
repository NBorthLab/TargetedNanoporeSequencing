#!/bin/bash

# Before this script can be used, a paf file is required
# Use for minimap2 a modified ref genome as input (does include two contigs with the expression vector (EV) sequence --> cat genome.fasta EV.fasta > modified.fasta)
# E.g. minimap2 modified.fasta input.fastq > minimap2.paf

# $1= minimap2.paf
# $2= name of expression vector contig
# $3= name of expression vector contig2

# Expression vector were linearized at two different sites (contig and contig2)
echo "Filter reads mapping to expression vector"

awk -v var="$2" -v var2="$3" '$6 ~ var || $6 ~ var2 {print $1}' $1 | sort -u > reads_mapping_to_EV.txt
while read name; do
        grep $name $1 | cut -f1,5,6,8,9,12 >> plasmidreads_location_of_read.txt
done<reads_mapping_to_plasmid.txt
awk 'BEGIN{OFS="\t"; print "Readname""\t""Strand""\t""Scaffold""\t""Start_alignment""\t""End_alignment""\t""MAPQ"}{print $1,$2,$3,$4,$5,$6}' plasmidreads_location_of_read.txt > reads_mapping_to_EV_stats.tsv

rm -rf plasmidreads_location_of_read.txt





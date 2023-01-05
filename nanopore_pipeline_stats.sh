#!/bin/env bash


set -eo pipefail

trap cleanup SIGINT SIGTERM ERR EXIT

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)


usage="$(basename "${BASH_SOURCE[0]}") -r <file.fasta|file.fa|file.fna> -i <file.fastq> -f <bedfile.bed> -d <directory> [...]

Attention: A bed file with the following format must be provided: 
Chromosome target_start_coordinate target_end_coordinate name_of_target_region 
NC_048595.1     101759735       101786524       b4galt1)

This script will map nanopore reads (fastq format) to a reference genome and creates stats about the stated ROI. That includes: Chromosome, Start, End, Name, Size of ROI, Average and Median Coverage, Number of Bases mapping, Number of Reads, Mean Read Length, Mean Accuracy, Strand Bias, on-targed enrichment and the total number of aligned reads.

This script depends on: minimap2, samtools, bedtools and stats_from_bam (part of the pomoxis package). Adjust your path in the variable: path_to_executables 

Samtools includes supplementary reads (-F 256) which is necessary when dealing with chimeric reads/integration sites. If only genomic regions are from interest, excluded sup. reads (-F 2304).

Use absolute paths for arguments. E.g. -r /path/to/genome.fasta -i /path/to/fastq_directory/ -f /path/to/bedfile.bed -d /path/to/new_directory_name

Available options:

	-h		Print this help and exit
	-v		Print script debug info
	-i		location of fastq files (files are in fastq format) (required)
	-r		genome in .fasta or .fna or .fa format (required)
	-f		bedfile in .bed format (required)
	-t		number of threads (default=1)
	-d		Specify a directory name where data will be stored (e.g. nanopore_results) (required)" 



cleanup() {
trap - SIGINT SIGTERM ERR EXIT
	tmp="tmp*"
		if [[ $tmp =~ "tmp*" ]] 
		then
		rm -f tmp*
		fi
	}



rflag=false
iflag=false
direc=false
bedflag=false
threads=1


while getopts ':hvr:i:d:f:t:' option; do
  case "$option" in
    h  ) echo "$usage" >&2; exit;;
    v  ) set -x ;;
    r  ) rflag=true; REFERENCE=$OPTARG;;
    i  ) iflag=true; INPUT=$OPTARG;;
    d  ) direc=true; DIRECTORY=$OPTARG;;
    f  ) bedflag=true; BEDFILE=$OPTARG ;;
    t  ) threads=$OPTARG ;;
    \? ) echo "Invalid option: -${OPTARG}." >&2; exit 1;;
    :  ) echo "Option -$OPTARG requires an argument." >&2; exit 1;;
  esac
done
shift $(($OPTIND - 1))

if ! $rflag; then
  echo "$usage" >&2;
  echo "-r must be specified" >&2;
  exit 1;
  elif ! [[ $REFERENCE =~ fasta$ || $REFERENCE =~ fna$ || $REFERENCE =~ fa$ ]]; then
	echo "$usage" >&2;
	echo "-r wrong file fromat">&2;
	exit 1;
fi

if ! $iflag; then
  echo "$usage" >&2;
  echo "-i must be specified." >&2;
  exit 1;
fi

if ! $bedflag; then
  echo "$usage" >&2;
  echo "-f must be specified" >&2;
  exit 1;
fi

if [ "$direc" = "" ];	then
	echo "$usage" >&2;
	echo "-d must be specified" >&2;
	exit 1;
fi

if [ $threads -le 0 ]; then
	echo "$usage" >&2;
	echo "-t must be at least 1" >&2;
	exit 1;
fi


#################
##  Functions  ##
#################


size_of_ROI () {
echo "Calculating size of ROI..."

awk '{print $3 - $2}' ROI_sorted.bed > tmp_tsize.txt
paste ROI_sorted.bed tmp_tsize.txt > tmp_tsize_bed.txt
}

median_coverage () {
echo "Calculating median coverage of ROI..."

cut -f4 ROI_sorted.bed > tmp_sort_cut_bed_names.txt
$path_to_executables/samtools depth alignments/alignments.bam -b ROI_sorted.bed > tmp_coverage_median.tsv # To included 0 per-base coverage add -a
cut -f1 ROI_sorted.bed > tmp_scaffolds_ROI_sorted.bed

if (( $(cut -f1 ROI_sorted.bed | wc -l | cut -c1-2) != $(cut -f1 ROI_sorted.bed | sort -u | wc -l | cut -c1-2) )) # Check if scaffolds are unique
then

                while read file
                do
                	local variable1=$(grep "$file" ROI_sorted.bed | cut -f1)
                	local variable2=$(grep "$file" ROI_sorted.bed | cut -f2)
                	local variable3=$(grep "$file" ROI_sorted.bed | cut -f3)
                	local variable4=$(grep "$file" ROI_sorted.bed | cut -f4)

                	awk -v var1="$variable1" -v var2="$variable2" -v var3="$variable3" 'BEGIN{OFS="\t"}{ if ($1==var1 && $2>=var2 && $2<=var3) print $0}' tmp_coverage_median.tsv >> tmp_is_same_${variable4}scf.txt
                done<tmp_sort_cut_bed_names.txt
	
                for i in tmp_is_same_*
                do
			# The median cov is calculted like...
        	        sort -k3,3 -n $i | awk '{count[NR] =$3}END{if (NR % 2) {print count[(NR + 1) /2]} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0}}' > tmp_median_cov${i}.txt
                done

else
        while read names
	do
                awk -v var="$names" '$1 ~ var {print $0}' tmp_coverage_median.tsv > tmp_cov_median_${names}.txt
                sort -k3,3 -n tmp_cov_median_${names}.txt | \
		awk '{count[NR] =$3}END{if (NR % 2) {print count[(NR + 1) /2]} else {print (count[(NR / 2)] + count[(NR / 2) + 1]) / 2.0}}' > tmp_median_cov${names}.txt
        done<tmp_scaffolds_ROI_sorted.bed

fi
        cat tmp_median_cov* > tmp_m_coverage.txt
        paste tmp_tsize_bed.txt tmp_m_coverage.txt > tmp_medi_coverage.txt
}

bases_against_ROI () {
mkdir bedtools_results
echo "Calculating number of bases mapping against ROI... "

$path_to_executables/intersectBed -a alignments/alignments.bam -b ROI_sorted.bed -bed -wo > bedtools_results/overlap_bases_ROI.txt
        while read names
	do
		# If no reads are mapping to ROI, the code below will be skipped and "N.A" will be printed. 
                if ! grep $names bedtools_results/overlap_bases_ROI.txt; then echo "N.A." >> tmp_overlap${names}.unique.txt && continue; fi >> tmp_base${names}_overlap.txt
                awk '{sum_bases += $NF}END{print sum_bases/1000}' tmp_base${names}_overlap.txt > tmp_overlap${names}.unique.txt
        done<tmp_sort_cut_bed_names.txt

cat *.unique.txt > tmp_overlap.txt
paste tmp_medi_coverage.txt tmp_overlap.txt > tmp_bases.txt
}

average_coverage_ROI () {
echo "Calculating average coverage of ROI... "

paste tmp_overlap.txt tmp_tsize.txt > tmp_a_coverage.txt
# The average coverage is calculated based on the number of bases mapped to the ROI divided by the size of the ROI.
awk '{printf "%.1f\n", $1 * 1000 / $2}' tmp_a_coverage.txt > tmp_average_coverage.txt
paste tmp_bases.txt tmp_average_coverage.txt > tmp_avg_cov.txt
}

number_of_reads_ROI () {
echo "Calculating number of reads mapping against ROI... "

$path_to_executables/intersectBed -a alignments/alignments.bam -b ROI_sorted.bed -bed -c > bedtools_results/nr_reads_align_to_ROI.txt
cut -f4,16 bedtools_results/overlap_bases_ROI.txt | sort -k2,2 > tmp_meanread_sort.txt
        while read names
	do
                awk -v var="$names" '$2 ~ var {print $1}' tmp_meanread_sort.txt | sort -u | wc -l | awk '{print $1}' > tmp_number_reads${names}.align.txt
                        if [[ ! -s "tmp_number_reads${names}.align.txt" ]] # Same as before, if file is empty (e.g. no read was detected) print "N.A." and continue.
			then
                        	echo "N.A." > tmp_number_reads${names}.align.txt
                         	continue
                        fi
        done<tmp_sort_cut_bed_names.txt

cat *.align.txt > tmp_number_reads.txt
paste tmp_avg_cov.txt tmp_number_reads.txt > tmp_reads_ROI.txt
}

mean_read_length_ROI () {
echo "Calculating mean read length mapping against ROI..."

while read names
do
	awk -v var="$names" '$2 ~ var {print $1}' tmp_meanread_sort.txt | sort -u > tmp_mean_${names}.txt
        if [[ ! -s "tmp_mean_${names}.txt" ]]
	then
        	echo "N.A." > tmp_mean_read_l${names}.txt
                rm -f tmp_mean_${names}.txt
                continue
        else
                # Of each ROI the, the mean read length is calculated based on the reads mapping to these ROI.
			local var=$(cat tmp_mean_${names}.txt)
                        for i in $var; do
                                if ! grep $i alignments/alignments.stats; then echo "N.A." && continue >> tmp_mean_read_l${names}.txt; fi | \
				cut -f12 >> tmp_meanread${names}_stats.txt
                        done
        	awk '{sum+=$1}END{print sum/NR}' tmp_meanread${names}_stats.txt > tmp_mean_read_l${names}.txt
        fi
done<tmp_sort_cut_bed_names.txt

cat tmp_mean_read_l* > tmp_meanreadlength.txt
paste tmp_reads_ROI.txt tmp_meanreadlength.txt > tmp_mean_read_length.txt
}

mean_read_accuracy_ROI () {
echo "Calculating mean accuracy of reads mapping against ROI..."

while read names
do
	awk -v var="$names" '$2 ~ var {print $1}' tmp_meanread_sort.txt | sort -u > tmp_accu_${names}.txt
        	if [[ ! -s "tmp_accu_${names}.txt" ]]
		then
                	echo "N.A." > tmp_finished_meanaccu${names}.txt
                        rm -f tmp_accu_${names}.txt
                        continue
                else
                        local var=$(cat tmp_accu_${names}.txt)
                        	for i in $var
				do
                                	if ! grep $i alignments/alignments.stats; then echo "N.A." && continue >> tmp_finished_meanaccu${names}.txt; fi | \
					cut -f18 >> tmp_meanaccu${names}_stats.txt
                        	done
        		awk '{sum+=$1}END{print sum/NR}' tmp_meanaccu${names}_stats.txt > tmp_finished_meanaccu${names}.txt
                fi
done<tmp_sort_cut_bed_names.txt

cat tmp_finished_meanaccu* > tmp_meanaccuracy.txt
paste tmp_mean_read_length.txt tmp_meanaccuracy.txt > tmp_accuracy.txt
}

strand_bias_ROI () {
echo "Calculating strand bias..."

while read names
do
        awk -v var="$names" '$2 ~ var {print $1}' tmp_meanread_sort.txt | sort -u > tmp_strdk_${names}.txt
        if [[ ! -s "tmp_strdk_${names}.txt" ]]
	then
                echo "N.A." > tmp_strand${name}.txt
                continue
        else
        	local var=$(cat tmp_strdk_${names}.txt)
                	for i in $var; do
                        	if ! grep $i alignments/alignments.stats; then echo "N.A." && continue >> tmp_strand${name}.txt; fi | \
				cut -f10 >> tmp_finish_strand${names}_stats.txt # Determine how many reads are mapping to "+" or "-"
                	done
        fi
done<tmp_sort_cut_bed_names.txt

while read name
do
       local files="./tmp_finish_strand${name}_stats.txt"
                if [[ ! -s "tmp_finish_strand${name}_stats.txt" ]]
		then
                	echo "N.A." > tmp_strand${name}.txt
                        rm -f tmp_finish_strand${name}_stats.txt
                        continue
                else
                	for i in $files; do
				# Determine strand bias between 1 and -1. 1 means all reads are mapping to the sense strand, -1 all reads to anti-sense.
                        	awk '{if ( $1~/+/ ) {pos+=1} else {neg+=1}}END{c=pos-neg; d=pos+neg; print c/d}' tmp_finish_strand${name}_stats.txt > tmp_strand${name}.txt
                	done
        	fi
done<tmp_sort_cut_bed_names.txt

cat tmp_strand* > tmp_strandinfo.txt
paste tmp_accuracy.txt tmp_strandinfo.txt > tmp_target_region.txt
}

on_target_percentage_ROI () {
echo "Calculating on-target enrichment..."

totalreads=$($path_to_executables/samtools view -F 2308 -c alignments/alignments.bam) # exclude supplementary, secondary and unmapped reads
paste tmp_number_reads.txt tmp_sort_cut_bed_names.txt > tmp_ontarget_genenames.txt

while read input
do
	local ontargetreads=$(grep $input tmp_ontarget_genenames.txt |cut -f1)
        	if [ "$ontargetreads" = "0" ]
		then
                        echo "N.A." >> tmp_on_target.txt
                        continue
        	else
                	echo "scale=6 ; $ontargetreads / $totalreads *100" | bc >> tmp_on_target.txt
        	fi
done<tmp_sort_cut_bed_names.txt

paste tmp_target_region.txt tmp_on_target.txt > tmp_ontarget_percentage.txt
}


#############################
##   Skript starts here    ##
#############################

path_to_executables="/data/borth/kleitner/miniconda3/bin"

rm -f ${REFERENCE}.fai
rm -f ${REFERENCE}.mmi
rm -rf ${DIRECTORY}
mkdir ${DIRECTORY}
cd ${DIRECTORY}
cat ${INPUT}*.fastq > merge.fastq
sort -k4,4 ${BEDFILE} > ${DIRECTORY}/ROI_sorted.bed

# -F 256 to remove only secondary alignments (not supplementary alignments) - necessary to get correct stats for chimeric reads
# Use -F 2304 when only genomic reads are from interest
# stats_from_bam excludes in default secondary and supplementary alignments, add -a to include both alignment types

mkdir alignments

$path_to_executables/minimap2 -ax map-ont ${REFERENCE} merge.fastq -t ${threads} | samtools view -T ${REFERENCE} -F 256 -q 31 -@ ${threads} -bS - | samtools sort -@ ${threads} -l 9 -o alignments/alignments.bam
$path_to_executables/samtools index -@ ${threads} alignments/alignments.bam alignments/alignments.bam.bai
$path_to_executables/stats_from_bam alignments/alignments.bam -a -o alignments/alignments.stats

#Determine the size of ROI 

size_of_ROI

#Median coverage calulated with samtools

median_coverage

#Calculating total number of bases mapping against ROI

bases_against_ROI

#Calculating average coverage

average_coverage_ROI

#Counting, how many reads align to ROI (It is enough when just 1 bp of read intersects with ROI)

number_of_reads_ROI

#Mean read length targeting ROI

mean_read_length_ROI
	
#Mean accuracy of reads mapping against ROI

mean_read_accuracy_ROI

#Calcualting of strandbias. 1 or -1 implies that reads aligning solely on the plus or minus strand, respectively

strand_bias_ROI

#On-target reads Percentage

on_target_percentage_ROI

#Total number of reads 
echo "The total aligned read count is $totalreads" | tee total_number_of_reads.txt

#Create a summary table
awk 'BEGIN{print "Scaffold\tStart\tEnd\ttname\ttsize\tmedian_coverage\tkbases\tavg_coverage\tnreads\tmean_read_length\tmean_accuracy\tstrand_bias\ton-target[%]";OFS="\t"}{print $0}' tmp_ontarget_percentage.txt | tee summary_results_ROI.txt





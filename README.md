# Targeted Nanopore Sequencing

The provided scripts should help to analyze targeted nanopore sequencing data and are briefly explained below.

### Identification of junction/integration sites

##### chimeric_reads_detection.sh

Prerequisite is a modified reference genome that has the expression vector included. The vector must be linearized at two different sites and the sequence is provided as two separate contigs (If a recombinant cell line might has more than one stable integrated expression vector, do also provide those. Then the scripts must be slightly adjusted). Additionally, a paf file generated with minimap2 is necessary. As a result, a table including read names, strand, start and end aligning positions will be generated. Chimeric reads will map against the expression vector and against an endogenous position. When a chimeric read has multiple start/end position at endogenous sites, it might be an indication for a low complexity region. Reads that are only mapping against the expression vector might be short reads, a result of concatemerization events or simply the transgenic sequence between two CRISPR/Cas9 cut sites. Normally, it can be expected that multiple reads have the same start or endposition, which is an indication for a junction site. (Reads aligning to the "+" strand will have similar start positon, whereas reads aligning to "-" strand will have similar end positions). Visualizing junction/integration sites with integrative genomics viewer (IGV) helps to understand junction sites, especially more complex ones.


##### ontarget_reads_one_site.sh

Chimeric reads are grouped according to their junction sites, and a fastq file with only on-target reads will be generated. To counter mapping uncertainty, reads with a start/end alignment position within 200bp in respect to their junction site will be considered as well. The script requires a scaffold name, the junction site coordinate, an alignment.stats file (output of stats_from_bam <input.bam> -o <outputfile.stats> which is part of the ont pomoxis package) and a fastq file. 

The on-target reads were further assembled and polished. For polishing, a script of Gilpatrick et al (2020) was modified (https://github.com/timplab/Cas9Enrichment/blob/master/Brca1_assembly/Racon_medaka_polishing.sh). The junction sites were determined via blastn. 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### nanopore_pipeline_stats.sh

This script generates statistics about a beforehand specified region of interest (bed format). That includes median and average coverage, read length, on-target enrichment and more. Additional information can be found within the script. 


##### methylation_calling.sh

This script is quite self-explanatory. To identify the methylation pattern of the stable integrated expression vector, the beforehand generated assemblies of the integration sites were used as a reference assembly.

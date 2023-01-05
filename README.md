# Targeted Nanopore Sequencing

The provided scripts should help to analyze targeted nanopore sequencing data and are briefly explained below.

### Identification of junction/integration sites

##### chimeric_reads_detection.sh

Prerequisite are a modified reference genome that has the expression vector inlcuded. The vector must be lineralized at two different sites and are provided as two separate contigs (If a recombinant cell line might has more than one stable integrated expression vector, do also provide those. Then the scripts must be slightly adjusted). Additonally, a paf file generated with minimap2 is necessary. As a result a table including read names, strand, start and end aligning positions will be produced. Chimeric reads will map against the expression vector and against an endogenous position. When a chimeric read has multiple start/end position at endogenous sites, it might be an indication of a low complexity region. Reads that are only mapping against the expression vector might be short reads, a result of concatemerization events or simply the transgenic sequence between two CRISPR/Cas9 cut sites. Normally, it can be expected that mulitple reads have the same start or endposition, which is an indication for a junction site. (Reads aligning to the "+" strand will have similar start positon, whereas reads aligning to "-" strand will have similar end positions). Visualizing junction/integration sites with integrative genomics viewer (IGV) helps to understand junction sites, especially more complex ones.


##### ontarget_reads_one_site.sh

This script groups chimeric reads accoring to their junction site and generates a fastq file with only on-target reads. To counter mapping uncertainty, reads with a start/end alignment position within 200bp in respect to the junction site will be considered as well. The script requires the scaffold name, the junction site coordinate, an alignment.stats file (output of stats_from_bam <input.bam> -o <outputfile.stats> which is part of the ont pomoxis package) and a fastq file with all reads. 

Then, an assembly with the on-target reads can be generated and polished. For polishing, a script of Gilpatrick et al (2020) was modified (https://github.com/timplab/Cas9Enrichment/blob/master/Brca1_assembly/Racon_medaka_polishing.sh). The junction sites were determined via blastn. 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### nanopore_pipeline_stats.sh

This script generates statistics about a beforehand specified region of interest (bed format). That includes median and average coverage, read length, on-target enrichment and more. Additional information can be found within the script. 


##### ontarget_reads_two_sites.sh

Filters on-target reads that are enlcosed by gRNAs, and thus are only mapping to a defined region of interest. Normally, at least two gRNAs per target site are used (two upstream and two downstream in regard to a region of interest). To include all reads that are generated from all gRNAs, reads that have their start/end position 3000bp up-/downstream of the innermost gRNAs are included and count as on-target reads. Additionally, those reads must stop aligning when reaching the first gRNA after the region of interest. The script needs as input: the output of stats_from_bam, start coordinate of the innermost gRNA, endpostion of the innermost gRNA, Scaffold name and a fastq file with all reads. 


##### methylation_calling.sh

This script is quite self-explanatory. To identify the methylation pattern of the stable integrated expression vector, the beforehand generated assemblies of the integration sites were used as a reference assembly.

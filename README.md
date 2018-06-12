# Germline variant calling and genotyping with GATK best practices (2018)
Variant calling identifies SNP and Indel sites that vary from the reference genome. Genotyping determines the genotype for each individual at called variant sites. [Variant calling and genotyping from next-generation sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3593722/) provides a detailed background. This script takes a single pair of paired end fastq files from whole genome or exome sequencing data which have been previously QC'd and performs all steps necessary to produce a vcf or gvcf file containing germline SNPs and Indels.

The [GATK best practice pipeline for calling germline variants](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) is illustrated in the image below. This script `fastqToVar.pl` includes all steps up to and including "Call Variants Per-Sample" which produces a vcf file which must be filtered prior to down stream processing. Filtering has not been included in this pipeline as it is subject to user preference and needs.

<img src="https://us.v-cdn.net/5019796/uploads/editor/mz/tzm69d8e2spl.png" width="800">

## Requirements

1 Software
* [GATK](https://software.broadinstitute.org/gatk/documentation/quickstart.php)
* [BWA](https://icb.med.cornell.edu/wiki/index.php/Elementolab/BWA_tutorial)
* [Java](https://java.com/en/download/help/download_options.xml)
* [Python](https://www.python.org/)

2 User required input
* A single pair of paired end fastq files
* [Read group information](#read-group-information)

3 Known variants
* Already hardcoded into script for mouse and macaque 
* To update or use other organims, find [known variants here](ftp://ftp.ensembl.org/pub/release-92/variation/vcf/)
* If you want to use the script on a non-model organism see ["No Excuses"](https://gatkforums.broadinstitute.org/gatk/discussion/11081/base-quality-score-recalibration-bqsr)

4 Reference genome
* Already hardcoded into script for mouse and macaque
* Downloaded fasta from [Ensembl](ftp://ftp.ensembl.org/pub/release-92/fasta/)
* Indexed with: `bwa index -a fastafile.fa && samtools faidx fastafile.fa && gatk CreateSequenceDictionary -R fastafile.fa`

5 Gene intervals
* Already hardcoded into script for mouse and macaque
* Downloaded chr, start, end to gene.intervals.bed from [Ensembl biomart](https://www.ensembl.org/biomart/martview/9e094011f1f0ee298e0b004e64597103)
* Split file by chromosome with: `grep '^1\s' gene.intervals.bed > gene.intervals.1.bed # repeated for each chromosome`

To subset your dataset to quickly test change you make to the pipeline:

`zcat filename_1.fastq.gz | sed -n 1,1000p > test_1.fastq && zcat filename_2.fastq.gz | sed -n 1,1000p > test_2.fastq && gzip *.fastq`

## Single v multi sample analysis
This script `fastqToVar.pl` runs on the paired end fastq files of a single sample to output a vcf file with '-of vcf'. However, if their are multiple samples, the script can be run on each sample to output a gvcf file for each with '-of gvcf'. The multiple gvcf files can then be consolidated and used to joint call variants. This is a continuation of the [GATK best practice pipeline for calling germline variants](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) illustrated in the image above. This [tutorial](https://software.broadinstitute.org/gatk/documentation/article?id=11813) describes how to consolidate gvcfs for joint-calling. 

## Read group information
The Broad Institute provide [a detailed explanation of read groups](https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups) and why they are needed. In summary, a read group is a set of reads that were generated from a single run of a sequencing instrument, so providing read group information assists the pipeline in the identification and removal of batch effects. Read group information can be found in the fastq filename and header and also in the sample sheet or qcstats sheet if available as shown in the following example:

* fastq filename:		WTCHG_461109_50_1.fastq.gz
* fastq headers:		@K00150:286:HNGMNBBXX:5:XXXX:XXXX:XXXX 1:N:0:XXXX       

		                @instrument	    run number	flowcell ID	    lane	tile	x-pos	y-pos 	read	is filtered	    control number	index (barcode)
						@K00150:		286:		HNGMNBBXX:		5:		XXXX:	XXXX:	XXXX:	1:		N:				0: 				XXXX
						# XXXX indicates entries that differentiate reads within a read group so is not needed for read group information

* sample/qcstats sheets:

                        Index	Tag			Readgroup		    Sample name		Sample id		Library				Type			Genome			Project		Date
						50		GTCTGTCA	WTCHG_461109_50	    mpc372-2.5e		POT5490A2		106/18_MPX_10nM		SureSelectXT	mm10 	8-bp	P180007 	2018-02-08

To view the barcodes present in the fastq file: `grep '^@K00150:286' WTCHG_461109_50_2.fastq | cut -d : -f 10 | sort | uniq -c | sort -nr > barcodes.txt`

The read group information required by the script `fastqToVar.pl` can now be extracted for this example data to build the [example run command](#example-run-command):

        --READ_GROUP_NAME 		WTCHG_461109_50             # ID
        --SAMPLE_NAME 			mpc372-2.5e                	# SM
        --LIBRARY_NAME 			106/18_MPX_10nM             # LB
        --PLATFORM 				illumina                    # PL
        --PLATFORM_UNIT 		HNGMNBBXX.GTCTGTCA.5		# PU # flowcellID.[barcode|date|readgroup].lane
        --SEQUENCING_CENTER 	WTCHG                       # CN
        --RUN_DATE 				2018-02-08    				# DT

## Options

All options are required.

* -m: 'wes' or 'wgs' for whole exome or whole shotgun sequencing respectively
* -sp: currently only supports 'mouse' or 'macaque'
* -of: 'vcf' or 'gvcf' as described in [single v multi sample analysis](#single-v-multi-sample-analysis)
* -fq1: full path to the first of the paired fastq files
* -fq2: full path to the second of the paired fastq files
* -id: read group name for GATK from read group information
* -pu: platform unit for GATK from read group information
* -sm: sample name for GATK from read group information
* -pl: platform for GATK from read group information
* -cn: sequencing centre code for GATK from read group information
* -dt: run date for GATK from read group information 
* -o: full path to output folder, must end with "/"

## Example run command

`fastqToVcfOnGrid.pl /
-m wes /
-sp mouse /
-of vcf /
-fq1 WTCHG_461109_50_1.fastq.gz /
-fq2 WTCHG_461109_50_2.fastq.gz /
-id WTCHG_461109_50 /
-pu HNGMNBBXX.GTCTGTCA.5 /
-sm mpc372-2.5e -lb 106/18_MPX_10nM /
-pl illumina /
-cn WTCHG /
-dt 2018-02-08 / 
-o /path/to/output/folder/`


Once the script finishes the read group information can be printed from the BAM file with: `samtools view -H WTCHG_461109_50.bam | grep '@RG'`

        @RG	ID:WTCHG_461109_50	SM:mpc372-2.5e	LB:106/18_MPX_10nM	PL:illumina PU:HNGMNBBXX.GTCTGTCA.5 CN:WTCHG	DT:2018-02-08
About variant calling:
	- SNP calling (aka variant calling) identifies sites that vary from the reference genome
	- Genotype calling determines the genotype for each individual at called SNP site.
	- Genotype and SNP calling from next-generation sequencing data: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3593722/
	- All comands are for WGS, wheras for WES the following options should be included for BaseRecalibrator and HaplotypeCaller
		-L requires a bed file of exon start/end coordinates
		-ip 100 required to pad intervals by 100bps either end

GATK resource bundle:
	- https://software.broadinstitute.org/gatk/documentation/article.php?id=1213 
	- ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
	- https://drive.google.com/drive/folders/1aBcbV_Hlyg0wOOmZDDSBeIc0uw1r3f_w

Create index of reference genome:
	- bwa index -a bwtsw human_g1k_v37_decoy.fasta # run once to create index files

Image for gitlab:
	- https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145

About read groups
	- All reads in a read group come from the same flowcell, on the same run, from the same biological sample
	- example
		fastq filename:		WTCHG_461109_50_1.fastq.gz
		fastq headers:		@K00150:286:HNGMNBBXX:5:XXXX:XXXX:XXXX 1:N:0:XXXX       # XXXX indicates entries that differentiate between reads
		header sections:		@instrument:	run_number:	flowcell_ID:	lane:	tile:	x-pos:	y-pos 	read:	is_filtered:	control_number:	index_sequence(barcode)
								@K00150:		286:		HNGMNBBXX:		5:		XXXX:	XXXX:	XXXX:	1:		N:				0: 				XXXX
		From qcstats.html:	Index	Tag			Readgroup		Sample name		Sample id		Library				Type			Genome			Project		Date
							50		GTCTGTCA	WTCHG_461109_50	mpc372-2.5e		POT5490A2		106/18_MPX_10nM		SureSelectXT	mm10 	8-bp	P180007 	2018-02-08

		--READ_GROUP_NAME 		WTCHG_461109_50             # ID
		--PLATFORM_UNIT 		HNGMNBBXX.GTCTGTCA.5		# PU # flowcellID.barcode.lane or flowcellID.data.lane
		--SAMPLE_NAME 			mpc372-2.5e                	# SM
		--LIBRARY_NAME 			106/18_MPX_10nM             # LB
		--PLATFORM 				illumina                    # PL
		--SEQUENCING_CENTER 	WTCHG                       # CN
		--RUN_DATE 				2018-02-08    				# DT

		# The @RG line this will create:
		@RG	ID:WTCHG_461109_50	SM:mpc372-2.5e	LB:106/18_MPX_10nM	PL:illumina PU:HNGMNBBXX.GTCTGTCA.5 N:WTCHG	DT:2018-02-08

		# Useful commands:
		samtools view -H WTCHG_461109_50.bam | grep '@RG' # to print read group info from bam file
		grep '^@K00150:286' WTCHG_461109_50_2.fastq | cut -d : -f 10 | sort | uniq -c | sort -nr > indices.txt # to print index_sequences(barcodes) from unzipped fastq

Requirements:
	- known variants
		- vcf.gz and vcf.gz.tbi of known variants: ftp://ftp.ensembl.org/pub/release-92/variation/vcf/[organism]/
		- if non model organsim see "no excuses" section of https://gatkforums.broadinstitute.org/gatk/discussion/11081/base-quality-score-recalibration-bqsr
	- indexed ref genome (e.g. mouse)
		- download
			wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
		- index (remove "." if commands don't work)
			gzip -d /Mus_musculus.GRCm38.dna.toplevel.fa.gz
			bwa index -a bwtsw Mus_musculus.GRCm38.dna.toplevel.fa
			samtools faidx Mus_musculus_GRCm38.dna.toplevel.fa
			gatk CreateSequenceDictionary -R Mus_musculus.GRCm38.dna.toplevel_fa
	- gene coords bed file
			downloaded from biomart - chr, start, end

Scripts:
	- pipeline.txt
	- fastqToVcf_grid.pl

Subsetting fastq files for testing
	zcat WTCHG_461109_50_1.fastq.gz | sed -n 1,1000p > test_1.fastq
	zcat WTCHG_461109_50_2.fastq.gz | sed -n 1,1000p > test_2.fastq
	gzip *.fastq

# Single sample analysis:
https://gatkforums.broadinstitute.org/gatk/discussion/9827/single-sample-genotyping-different-workflow
https://gatkforums.broadinstitute.org/gatk/discussion/7943/single-sample-vs-multiple-samples-haplotype-caller


# command
fastqToVcfOnGrid.pl 
-m wes or wgs
-sp mouse or macaque
-of vcf or gvcf (choose vcf if a single sample, if you choose gvcf you must run GenotypeGVCF after)
-fq1 WTCHG_461109_50_1.fastq.gz 
-fq2 WTCHG_461109_50_2.fastq.gz
-id WTCHG_461109_50 
-pu HNGMNBBXX.GTCTGTCA.5 
-sm mpc372-2.5e -lb 106/18_MPX_10nM 
-pl illumina 
-cn WTCHG 
-dt 2018-02-08 
-o [outputfolder]

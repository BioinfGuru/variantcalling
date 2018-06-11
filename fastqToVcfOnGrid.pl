#!/bin/perl
# Takes paired fastq files to create a vcf file on the grid#
# Each pair of fastq files must be in their own empty folder

# example command
# /NGS/users/Kenneth/toolbox/variantcalling/fastqToVcfOnGrid.pl -m wes -sp mouse -of vcf -fq1 /NGS/users/Kenneth/jobs/michelle_simon/mm10_analysis/testdata1M/testdata1M_1.fastq.gz -fq2 /NGS/users/Kenneth/jobs/michelle_simon/mm10_analysis/testdata1M/testdata1M_2.fastq.gz -id testdata1M -pu HNGMNBBXX.GTCTGTCA.5 -sm mpc372-2.5e -lb 106/18_MPX_10nM -pl illumina -cn WTCHG -dt 2018-02-08 -o /NGS/users/Kenneth/jobs/michelle_simon/mm10_analysis/testdata1M/

use strict;
use warnings;
use Getopt::Long;

my ($m, $of, $sp, $fq1, $fq2, $id, $pu, $sm, $lb, $pl, $cn, $dt, $o, $help) = "" ;
GetOptions(
     'm=s' => \$m,
     'of=s' => \$of,     
     'sp=s' => \$sp,
 	'fq1=s' => \$fq1,
 	'fq2=s' => \$fq2,
 	'id=s' => \$id,
 	'pu=s' => \$pu,
 	'sm=s' => \$sm,
 	'lb=s' => \$lb,
 	'pl=s' => \$pl,
 	'cn=s' => \$cn,
 	'dt=s' => \$dt,
 	'o=s' => \$o,
    'help'   => \$help

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the usage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	Important: Each pair of fastq files must be in their own empty folder
     ARGUMENTS :
                    REQUIRED
                    -m\t\tmode: wgs or wes
                    -of\t\toutput format: vcf or gvcf
                    -sp\t\tspecies: mouse or macaque
                    -fq1\tfull path to first fastq file in pair e.g. \/full\/path\/to\/[read group]_1.fastq.gz
                    -fq2\tfull path to second fastq file in pair e.g. \/full\/path\/to\/[read group]_2.fastq.gz
                    -id\t\tread group (will be used to name all output files)
                    -pu\t\tplatform unit
                    -sm\t\tsample name
                    -lb\t\tlibrary name
                    -pl\t\tplatform e.g. illumina
                    -cn\t\tsequencing center, e.g. WTCHG
                    -dt\t\trun date yyyy-mm-dd
                    -o\t\tfull path to output directory
                    OPTIONAL
                    -help -> prints this help message
               \n/);
}

if ($help) { &useage ;}
if (!$m || !$of || !$sp || !$fq1 || !$fq2 || !$id || !$pu || !$sm || !$lb || !$pl || !$cn || !$dt || !$o) { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
unless (($m eq 'wgs') or ($m eq 'wes')){print "\nInvalid input for -m\nUse '-m wgs' or '-m wes'\n\n"};
unless (($sp eq 'mouse') or ($sp eq 'macaque')){print "\nInvalid input for -sp\nUse '-sp mouse' or '-sp macaque'\n\n"};
unless (($of eq 'vcf') or ($of eq 'gvcf')){print "\nInvalid input for -of\nUse '-of vcf' or '-of gvcf'\n\n"};

# create required folders
my $logs = $o."logs";
mkdir $logs;
my $temp = $o."temp";
mkdir $temp;

# get start time
my $runtimefile = $o."runtime.txt";
open (RUNTIME, ">".$runtimefile);
print (RUNTIME "Start:\t\t".`date`."\nEnd:\t\t");
close (RUNTIME);

# set gatk
my $gatk = "/usr/java/latest8/bin/java -jar /NGS/Software/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar";

# set species dependencies
my ($ref, $exons, $vcf, @chrlist);
if ($sp eq 'mouse')
     {
          $ref = '/NGS/musRefs_10/gatk/ref/ensembl92/Mus_musculus.GRCm38.dna.toplevel.fa';
          $exons = '/NGS/musRefs_10/gatk/ref/ensembl92/exons';
          $vcf = '/NGS/musRefs_10/gatk/vcf/ensembl92/mus_musculus.vcf.gz';
          @chrlist = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y);
     } elsif ($sp eq 'macaque')
          {
               $ref = '/NGS/macaqueRef/Mmul_8.0.1/ensembl/refdata/fasta/genome.fa';
               $exons = '/NGS/macaqueRef/Mmul_8.0.1/ensembl/refdata/regions';
               $vcf = '/NGS/macaqueRef/Mmul_8.0.1/ensembl/refdata/snps/macaca_mulatta.vcf';
               @chrlist = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 X Y);
          }

# create bam
system "qsub -cwd -j y -b yes -P NGS -N FqToSam.$id -o $logs -q small.q $gatk FastqToSam --FASTQ $fq1 --FASTQ2 $fq2 --OUTPUT $o/unmapped.$id.TEMP.bam --READ_GROUP_NAME $id --PLATFORM_UNIT $pu --SAMPLE_NAME $sm --LIBRARY_NAME $lb --PLATFORM $pl --SEQUENCING_CENTER $cn --RUN_DATE $dt";

# mark adapters
system "qsub -hold_jid FqToSam.$id -cwd -j y -b yes -P NGS -N MarkAdapters.$id -o $logs -q small.q $gatk MarkIlluminaAdapters --INPUT $o/unmapped.$id.TEMP.bam --OUTPUT $o/markedadapters.$id.TEMP.bam --METRICS $o/$id.markedadapters.metrics --TMP_DIR $o/temp";

# make interleaved fastq
system "qsub -hold_jid MarkAdapters.$id -cwd -j y -b yes -P NGS -N InterLeaveFq.$id -o $logs -q small.q $gatk SamToFastq --INPUT $o/markedadapters.$id.TEMP.bam --FASTQ $o/$id.interleaved.fq --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true --TMP_DIR $temp";

# map to ref
open(my $fh, '>', $temp."/command.sh");
print $fh 'bash -c "/NGS/Software/bwa/bwa mem -M -t 32 -p '.$ref.' '.$o.$id.'.interleaved.fq | '.$gatk.' MergeBamAlignment --REFERENCE_SEQUENCE '.$ref.' --UNMAPPED_BAM '.$o.'unmapped.'.$id.'.TEMP.bam --ALIGNED_BAM /dev/stdin --OUTPUT '.$o.'mapped.'.$id.'.TEMP.bam --CREATE_INDEX true --ADD_MATE_CIGAR true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --TMP_DIR '.$temp.'"';
close $fh;
system "qsub -hold_jid InterLeaveFq.$id -cwd -j y -b yes -P NGS -N Mapping.$id -o $logs -pe big 32 sh $temp/command.sh";

# mark dups
system "qsub -hold_jid Mapping.$id -cwd -j y -b yes -P NGS -N MarkDups.$id -o $logs -q small.q $gatk MarkDuplicates --INPUT $o/mapped.$id.TEMP.bam --OUTPUT $o/dedup.$id.TEMP.bam --METRICS_FILE $o/$id.dedup.metrics";

# calc bqsr
if ($m eq 'wgs')
     {
     system "qsub -hold_jid MarkDups.$id -cwd -j y -b yes -P NGS -N CalcBQSR.$id -o $logs -q small.q $gatk BaseRecalibrator --input $o/dedup.$id.TEMP.bam --reference $ref --known-sites $vcf --output $o/$id.before_bqsr.table";
     } elsif ($m eq 'wes')
          {
          system "qsub -hold_jid MarkDups.$id -cwd -j y -b yes -P NGS -N Index.$id -o $logs -q small.q samtools index $o/dedup.$id.TEMP.bam";
          system "qsub -hold_jid Index.$id -cwd -j y -b yes -P NGS -N CalcBQSR.$id -o $logs -q small.q $gatk BaseRecalibrator --input $o/dedup.$id.TEMP.bam --reference $ref --known-sites $vcf --output $o/$id.before_bqsr.table -L $exons/gene.intervals.bed -ip 100";
          }

# apply bqsr
system "qsub -hold_jid CalcBQSR.$id -cwd -j y -b yes -P NGS -N ApplyBQSR.$id -o $logs -q small.q $gatk ApplyBQSR --bqsr-recal-file $o/$id.before_bqsr.table --input $o/dedup.$id.TEMP.bam --output $o/$id.bam"; 

# get mapping stats  
system "qsub -hold_jid ApplyBQSR.$id -cwd -j y -b yes -P NGS -N FlagStats.$id -o $o/$id.flagstat.metrics -q small.q samtools flagstat $o/$id.bam";

# call variants
my @jobids;
my @vcfs;
if ($m eq 'wgs')
     {
          foreach my $chr(@chrlist)
               {
                    my $jobname = 'var.chr'.$chr.'.'.$id;
                    if ($of eq 'vcf')
                         {    
                               system "qsub -hold_jid FlagStats.$id -cwd -j y -b yes -P NGS -N $jobname -o $logs -pe big 8 $gatk HaplotypeCaller --input $o/$id.bam --output $temp/$chr.$of.gz --reference $ref --standard-min-confidence-threshold-for-calling 10.0 --dbsnp $vcf --native-pair-hmm-threads 8 -L $chr";
                         } elsif ($of eq 'gvcf')
                              {
                                   system "qsub -hold_jid FlagStats.$id -cwd -j y -b yes -P NGS -N $jobname -o $logs -pe big 8 $gatk HaplotypeCaller --input $o/$id.bam --output $temp/$chr.$of.gz --reference $ref --emit-ref-confidence GVCF --dbsnp $vcf --native-pair-hmm-threads 8 -L $chr";
                              }
                    push @jobids, $jobname;
                    push @vcfs, $temp.'/'.$chr.'.'.$of.'.gz';
               }

    } elsif ($m eq 'wes')
         {
               foreach my $chr(@chrlist)
                    {
                         system "qsub -hold_jid FlagStats.$id -cwd -j y -b yes -P NGS -N split.$chr.$id -o $logs -q small.q samtools view -b $o/$id.bam $chr -o $temp/$chr.bam";
                         system "qsub -hold_jid split.$chr.$id -cwd -j y -b yes -P NGS -N index.$chr.$id  -o $logs -q small.q samtools index $temp/$chr.bam";
                         my $jobname = 'var.chr'.$chr.'.'.$id;
                      if ($of eq 'vcf')
                         {    
                              system "qsub -hold_jid index.$chr.$id -cwd -j y -b yes -P NGS -N $jobname -o $logs -pe big 8 $gatk HaplotypeCaller --input $temp/$chr.bam --output $temp/$chr.$of.gz --reference $ref --standard-min-confidence-threshold-for-calling 10.0 --dbsnp $vcf --native-pair-hmm-threads 8 -L $exons/gene.intervals.$chr.bed -ip 100";
                         } elsif ($of eq 'gvcf')
                              {
                                   system "qsub -hold_jid index.$chr.$id -cwd -j y -b yes -P NGS -N $jobname -o $logs -pe big 8 $gatk HaplotypeCaller --input $temp/$chr.bam --output $temp/$chr.$of.gz --reference $ref --emit-ref-confidence GVCF --dbsnp $vcf --native-pair-hmm-threads 8 -L $exons/gene.intervals.$chr.bed -ip 100";
                              }                                                
                         push @jobids, $jobname;
                         push @vcfs, $temp.'/'.$chr.'.'.$of.'.gz';
                    }
         }

my $joblist = join ( ',', @jobids );
my $vcflist = join ( ' ', @vcfs );            
system "qsub -hold_jid $joblist -cwd -j y -b yes -P NGS -N concat.$id -o $logs -q small.q /NGS/Software/bcftools/bcftools_grid/installation/bin/bcftools concat -o $o/$id.$of $vcflist";
system "qsub -hold_jid concat.$id -cwd -j y -b yes -P NGS -N bgzip.$id -o $logs -q small.q bgzip $o/$id.$of";
system "qsub -hold_jid bgzip.$id -cwd -j y -b yes -P NGS -N tabix.$id -o $logs -q small.q tabix -p vcf $o/$id.$of.gz";

# tidy
system "qsub -hold_jid tabix.$id -cwd -j y -b yes -P NGS -N tidy.$id -o $logs -q small.q rm -rf $temp $o/*TEMP.ba* $o/*fq $o/*.table";

# get end time
system "qsub -hold_jid tidy.$id -cwd -j y -b yes -P NGS -N endtime.$id -o $logs -q small.q bash -c 'date >>$runtimefile'";

#!/usr/bin/perl -w

use strict;
use warnings;

#Master script for analysing whole genome paired end reads from Illumina pipleine using mosaik. The program performs QC in FASTQC, trims reads, aligns reads to human genome using Mosaik program suit and generates a coverage report.

#Programs FASTQC, filter_fastq.pl, MosaikMerge and calulate_coverage_statistics.pl kan be skipped. MosaikBuild, mosaikAlign and mosaikSort are so far required to be run. 
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
wgs_mosaik_wf.pl  -id [indir...n] -a [project ID] -s [sample ID...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata]

=head2 COMMANDS AND OPTIONS

-id/--infiledir Infile dir(s) (Mandatory: Supply whole path)

-ids/--indirscript The pipeline custom script in dir (Mandatory: Supply whole path)

-rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)

-a/--projectid The project ID (Mandatory)

-s/--sampleid The sample ID(s) (Mandatory)

-em/--email

-odf/--outdirdata The data files output directory (Supply whole path, defaults to data)

-ods/--outdirscript The script files output directory (Supply whole path, defaults to wgs_wf_scripts)

-pFQC/--fastqc Flag for running FASTQC (defaults to yes (=1))

-pFQF/--fast_fs.pl Flag for running filter_fastq.pl (defaults to yes (=1))

-pFQC2/--fastqc Flag running FASTQC after filter_fastq.pl (defaults to yes (=1))

-pMoB/--mosaikBuild Flag running MosaikBuild (defaults to yes (=1))

-pMoA/--mosaikAlign Flag running MosaikAlign (defaults to yes (=1))

-pMoS/--mosaikSort Flag running MosaikSort (defaults to yes (=1))

-pMoDup/--mosaikSortSupsnoop Flag running mosaik Dupsnoop removal of PCR duplicates (defaults to yes (=1))

-pMoM/--mosaikMerge Flag running MosaikMerge (defaults to yes (=1))

-pMoT/--mosaikText Flag running MosaikText (defaults to yes (=1))

-pCR/--calculate_coverage_statistics Flag running Calculate_coverage_statistics (defaults to yes (=1))

-pRCP/--rcovplots  Flag running rcovplots (defaults to yes (=1))

-pGZ/--gzip  Flag generating gzip sbatch (defaults to yes (=1))

-pRMMOB/--gzip  Flag generating sbatch of rm of mosaikkBuild.dat files (defaults to yes (=1))

-xcov/--xcoverage  Flag determining x-scale coverage plot in rcovplots (defaults to "30")

-crbed/--crbedfile Bed file for calculate_coverage_statistics (defaults to "CCDS-hg19.bed")

-crfasta/--crfasta Fasta file for calculate_coverage_statistics (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")

=head3 I/O

Input format ( dir/infiles.fastq )

Output format

1. Fastqc files

2. Filtered fastq reads read1 read2 and readS where one read has been filtered away

3. Mosaik.dat files

4. Mosaik_align.dat files

5. Mosaik_align_sorted.dat files

6. Mosaik_lanes_merged.dat

7. Mosaik_lanes_merged.bam

8. Calculate_coverage_results.stdout

9. Coverage plots

10. Gzipped fastq and mosaikAlign files

=head4 Dependencies

Filter_fastq.pl

1. Perl modules (Magnus Bjursell) 
#use lib '/home/magnus.bjursell/script/modules';
#use myFileIO qw(:basic);
#use myMathStat qw(max min round);

Mosaik

1. .dat files of reference
2. Jump database in script/mosaik dir

calculate_coverage_statistics.2.0.pl

1. Perl modules (Magnus Bjursell)
#use lib '/home/magnus.bjursell/script/modules';
#use myFileIO qw(:basic);
#use myMathStat qw(max min round);
#use mySeqAnalysis qw(isCanonicalChr);

R scripts
1. Average_cov_chr.R
2. Coverage_bed_pileup_genome.R
3. Coverage_hist_by_Chr.R

=cut
    
use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{wgs_mosaik_wf.pl -id [indir...n] -a [projectid] -s [sampleid] -em [e-mail] -ods [outdirscripts] -odf [outdirdata]
	       -id/--infiledir Infile dir(s), comma sep (Mandatory: Supply whole path,)
               -ids/--indirscript The pipeline custom script in dir (Mandatory: Supply whole path)
               -rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)
	       -a/--projectid The project ID  (Mandatory)
	       -s/--sampleid The sample ID(s),comma sep (Mandatory)
	       -em/--email e-mail
	       -odf/--outdirdata The data files output directory (Supply whole path, defaults to data)
	       -ods/--outdirscript The script files output directory (Supply whole path, defaults to wgs_wf_scripts)
	       -pFQC/--fastqc Flag running FASTQC (defaults to yes (=1))
	       -pFQF/--fast_fs.pl Flag running filter_fastq.pl (defaults to yes (=1))
	       -pFQC2/--fastqc Flag running FASTQC after filter_fastq.pl (defaults to yes (=1))
	       -pMoB/--mosaikBuild Flag running MosaikBuild (defaults to yes (=1))
	       -pMoA/--mosaikAlign Flag running MosaikAlign (defaults to yes (=1))
	       -pMoS/--mosaikSort Flag running MosaikSort (defaults to yes (=1))
               -pMoDup/--mosaikSortSupsnoop Flag running mosaik Dupsnoop removal of PCR duplicates (defaults to yes (=1))
	       -pMoM/--mosaikMerge Flag running MosaikMerge (defaults to yes (=1))	
	       -pMoT/--mosaikText Flag running MosaikText (defaults to yes (=1))
	       -pCR/--calculate_coverage_statistics Flag running Calculate_coverage_statistics (defaults to yes (=1))
	       -pRCP/--rcovplots  Flag running rcovplots (defaults to yes (=1))
	       -pGZ/--gzip  Flag generating gzip sbatch (defaults to yes (=1))
	       -pRMMOB/--gzip  Flag generating sbatch of rm of mosaikBuild.dat files (defaults to yes (=1))              
	       -xcov/--xcoverage  Flag determining x-scale coverage plot in rcovplots (defaults to "30")
               -crbed/--bedfile Bed file for calculate_coverage_statistics (defaults to "CCDS-hg19.bed")
               -crfasta/--crfasta Fasta file for calculate_coverage_statistics (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")
	   };
    
}

my ($aid,$em, $ids, $rd, $odf, $ods, $xcov, $crbed, $crfasta, $fnend, $filename, $filename2, $fnt, $fnt2, $help) = (0,0,0,0,0,0,30, "CCDS-hg19.bed", "Homo_sapiens.GRCh37.57.dna.concat.fa", ".sh"); #Arguments for project
my ($pFQC, $pFQF, $pFQC2, $pMoB, $pMoA, $pMoS, $pMoDup, $pMoM, $pMoT, $pCR, $pRCP, $pGZ, $pRMMOB) = (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1); #Default arguments for running programs
my (@inid,@sid); #Arrays for input dirs,sample ids
my (%infiles, %indirpath, %fqfiltInfiles, %fqfiltSInfiles, %MosBInfiles, %MosBSInfiles, %MosBAllInfiles, %lanes); 
#infiles=from platform, fgfiltInfile= from filter_fast.pl,fgfiltSInfile= Singels from filter_fast.pl, MosBInfiles for MosaikBuild, 
#MosBSInfiles for MosaikBuild singels, MosBAllInfiles for all infiles both PE and singels, Lanes for sample lanes
my ($dataHR, $dataMosB) = ({}, {}); #For reformating after filter_fastq.pl and for MosaikBuild output

GetOptions('id|infiledir:s'  => \@inid, #Comma separeted list
	   'ids|inscriptdir:s'  => \$ids, #Directory for custom scripts required by the pipeline
	   'rd|referencedir:s'  => \$rd, #directory containing references
	   'a|projectid:s'  => \$aid,
	   's|sampleid:s'  => \@sid, #Comma separeted list, one below outdirdata
	   'em|email:s'  => \$em,
	   'odf|outdirdata:s'  => \$odf, #One above sample id, must supply whole path i.e. /bubo/proj/...
	   'ods|outdirscript:s'  => \$ods,
	   'pFQC|fastqc:n' => \$pFQC,
	   'pFQF|fastqf:n' => \$pFQF,
	   'pFQC2|fastqc2:n' => \$pFQC2,
	   'pMoB|mosaikBuild:n' => \$pMoB,
	   'pMoA|mosaikAlign:n' => \$pMoA,
	   'pMoS|mosaikSort:n' => \$pMoS,
	   'pMoDup|mosaikDup:n' => \$pMoDup,
	   'pMoM|mosaikMerge:n' => \$pMoM,
	   'pMoT|mosaikText:n' => \$pMoT,
	   'pCR|cal_cov_stat:n' => \$pCR,
	   'crbed|cal_cov_bedfile:s' => \$crbed,
	   'crfasta|cal_cov_fastafile:n' => \$crfasta,
	   'pRCP|rcovplots:n' => \$pRCP,
	   'xcov|xcoverage:n' => \$xcov,
	   'pGZ|gzip:n' => \$pGZ,
	   'pRMMOB|rmmob:n' => \$pRMMOB,
	   'h|help' => \$help,
	   );

die $USAGE if( $help );

if (@inid == 0) {
   my $verbosity = 2;
 print"\n";
 pod2usage({-message => "Must supply an infile directory as comma separeted list.\n",
     -verbose => $verbosity
   });
}
if ($aid eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a project ID", "\n\n";
    die $USAGE;
}
if ( scalar(@sid) eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a sample ID as a comma separated list", "\n\n";
    die $USAGE;
}
if ( $ids eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a script dir", "\n\n";
    die $USAGE;
}
if ( $rd eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a reference dir", "\n\n";
    die $USAGE;
}
if ($odf eq 0) {
    
    print STDERR "\n";
    $odf = "/bubo/proj/$aid/private/data";
    print STDERR "Setting output data dir to: $odf", "\n\n";
}
if ($ods eq 0) {
    
    print STDERR "\n";
    $ods = "/bubo/proj/$aid/private/wgs_wf_scripts";
    print STDERR "Setting output scripts dir to: $ods", "\n\n";
}

@inid = split(/,/,join(',',@inid)); #Enables comma separated indir(s)
@sid = split(/,/,join(',',@sid)); #Enables comma separeted list of sample IDs


for (my $inputdir=0;$inputdir<scalar(@inid);$inputdir++) { #Collects inputfiles
    
    my @infiles = `cd $inid[ $inputdir ];ls *.fastq;`; #cd to input dir and collect .fastq files
    
    print STDERR "\nReads from Platform", "\n";
    print STDERR "\nSample ID", "\t", $sid[$inputdir],"\n";
    print STDERR "Inputfiles", "\n", @ { $infiles{ $sid[$inputdir] }  =[@infiles] }; #hash with sample id as key and inputfiles in dir as array 
    $indirpath{$sid[$inputdir]} = $inid[ $inputdir ];  #Catch inputdir path
}

FilterReFormat(); #Required to run filter.fastq.pl

MoSBReFormat(); #Required to format MosaikBuild outfiles correctly

#########################
###Run program part######
#########################

if ($pFQC eq 1) { #FASTQC
    
    print STDERR "\nFASTQC", "\n";
    print STDERR "Creating sbatch script FASTQC and writing script file(s) to: ", $ods,"/sampleid/fastqc/fastqc_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script FASTQC data files will be written to: ", $odf,"/sampleid/fastqc/", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Fastqc($sid[$sampleid]);
	
    }

}

if ($pFQF eq 1) { #filter_fastq.pl

    print STDERR "\nfilter_fastq.pl", "\n";

    print STDERR "Creating sbatch script Fast_fix_single and writing script file(s) to: ", $ods,"/sampleid/fastq_filter/fastq_filter_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script filter_fastq.pl data files will be written to: ", $odf,"/sampleid/fastq_filter/", "\n";


    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Fastq_filter($sid[$sampleid]);
	
    }

}

#Runs FASTQC if filter_fastq.pl were run and FASTQC were requested
if ($pFQC2 eq 1 && $pFQF eq 1) { #FASTQC after filter_fastq.pl
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	
	Fastqc_filtered($sid[$sampleid]);
	
    }
    
    print STDERR "\nFASTQC2", "\n";
    print STDERR "Creating sbatch script FASTQC2 and writing script file(s) to: ", $ods,"/sampleid/fastqc_filtered/fastqc_filtered_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script FASTQC data files will be written to: ", $odf,"/sampleid/fastqc_filtered/", "\n";
    print STDERR "Script will run after Filter_fastq sbatch have completed", "\n";
}

if ($pMoB eq 1) { #Run MosaikBuild
    
    print STDERR "\nMosaikBuild", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	
	MosaikBuild($sid[$sampleid]);
	
    }
    
    print STDERR "Creating sbatch script MosaikBuild and writing script file(s) to: ", $ods,"/sampleid/mosaik/mosaikBuild_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script MosaikBuild data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
    
    if ($pFQF eq 1) {
	
	print STDERR "Script will run after Filter_fastq sbatch have completed", "\n";
	
    }
    
}


if ($pMoA eq 1 && $pMoB eq 1) { #Run MosaikAlign but only if MosaikBuild was run previously
    
    print STDERR "\nMosaikAlign", "\n";
    
    if ($pFQF eq 0) {
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    MosaikAlignPlatform($sid[$sampleid]);
	    
	}
    }
    else {
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    MosaikAlignFilter($sid[$sampleid]);
	    
	}
    }
    print STDERR "Creating sbatch script MosaikAlign and writing script file(s) to: ", $ods,"/sampleid/mosaik/mosaikAlign_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script MosaikAlign data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
    
    if ($pFQF eq 1) {
	
	print STDERR "Script will run after Filter_fastq and MosaikBuild sbatch scripts have been completed", "\n";
	
    }
    else {
	print STDERR "Script will run after MosaikBuild sbatch scripts have been completed", "\n";
    }
    
}

if ($pMoS eq 1 && $pMoA eq 1 && $pMoB eq 1) { #Run MosaikSort but only if MosaikBuild and Align were run previously
    
    print STDERR "\nMosaikSort", "\n";

    if ($pFQF eq 0) {    
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    MosaikSortPlatform($sid[$sampleid]);
	    
	}
    }
    else {
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    MosaikSortFilter($sid[$sampleid]);
	    
	}
    }
    print STDERR "Creating sbatch script MosaikSort and writing script file(s) to: ", $ods,"/sampleid/mosaik/mosaikSort_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script MosaikSort data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
    
    if ($pFQF eq 1) {
	
	print STDERR "Script will run after Filter_fastq, MosaikBuild and MosaikAlign sbatch scripts have been completed", "\n";
	
    }
    else {
	print STDERR "Script will run after MosaikBuild and MosaikAlign sbatch scripts have been completed", "\n";
    }
}

if ($pMoM eq 1 && $pMoS eq 1 && $pMoA eq 1 && $pMoB eq 1) { #Run MosaikMerge but only if MosaikBuild, Align and Sort were run previously
    
    print STDERR "\nMosaikMerge", "\n";

    if ($pFQF eq 0) {    
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    if (scalar( @{ $MosBInfiles{$sid[$sampleid]} }) > 1) {
		
		MosaikMergePlatform($sid[$sampleid]);
		
	    }
	}
    }
    else {

	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    if (scalar( @{ $MosBAllInfiles{$sid[$sampleid]} }) > 1) {
		
		MosaikMergeFilter($sid[$sampleid]);
		
	    }
	}
    }
    print STDERR "Creating sbatch script MosaikMerge and writing script file(s) to: ", $ods,"/sampleid/mosaik/mosaikMerge_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script MosaikMerge data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge", "\n";
    
    if ($pFQF eq 1) {
	
	print STDERR "Script will run after Filter_fastq, MosaikBuild, MosaikAlign and MosaikSort sbatch scripts have been completed", "\n";
	
    }
    else {
	print STDERR "Script will run after MosaikBuild, MosaikAlign, MosaikSort sbatch scripts have been completed", "\n";
    }
}

if ($pMoT eq 1 && $pMoS eq 1 && $pMoA eq 1 && $pMoB eq 1) { #Run MosaikText but only if MosaikBuild, Align and Sort were run previously
    
    print STDERR "\nMosaikText", "\n";
    
    if ($pFQF eq 0) {    
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    MosaikTextPlatform($sid[$sampleid]);
	    
	}
    }
    else {
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    MosaikTextFilter($sid[$sampleid]);
	    
	}
    }
    print STDERR "Creating sbatch script MosaikText and writing script file(s) to: ", $ods,"/sampleid/mosaik/mosaikText_wf_sampleid.sh", "\n";
    
    if ($pFQF eq 1) {
	
	if ($pMoM eq 0) {
	    
	    print STDERR "Sbatch script MosaikText data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
	    print STDERR "Script will run after Filter_fastq, MosaikBuild, MosaikAlign and MosaikSort sbatch scripts have been completed", "\n";
	}
	else {
	    print STDERR "Sbatch script MosaikText data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge", "\n";
	    print STDERR "Script will run after Filter_fastq, MosaikBuild, MosaikAlign, MosaikSort and MosaikMerge sbatch scripts have been completed", "\n";

	}
    }
    else {
	
	if ($pMoM eq 0) {

	    print STDERR "Sbatch script MosaikText data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
	    print STDERR "Script will run after MosaikBuild, MosaikAlign and MosaikSort sbatch scripts have been completed", "\n";
	}
	else {

	    print STDERR "Sbatch script MosaikText data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge", "\n";
	    print STDERR "Script will run after MosaikBuild, MosaikAlign, MosaikSort and MosaikMerge sbatch scripts have been completed", "\n";

	}
    }
}

if ($pCR eq 1) { #Run calculate_coverage_statistics.2.0.pl
    
    print STDERR "\ncalculate_coverage_statistics.2.0.pl", "\n";
    
    if ($pFQF eq 0) {    
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    Cal_cov_statPlatform($sid[$sampleid]);
	    
	}
    }
    else {
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    Cal_cov_statFilter($sid[$sampleid]);
	    
	}
    }
    print STDERR "Creating sbatch script calculate_coverage_statistics.2.0.pl and writing script file(s) to: ", $ods,"/sampleid/mosaik/cal_cov_stat.2.0_wf_sampleid.sh", "\n";
    
    if ($pFQF eq 1) {
	
	if ($pMoM eq 0) {

	    print STDERR "Sbatch script calulate_coverage_statistics.pl data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
	    print STDERR "Script will run after Filter_fastq, MosaikBuild, MosaikAlign, MosaikSort and MosaikText sbatch scripts have been completed", "\n";
	}
	else {
	    
	    print STDERR "Sbatch script calulate_coverage_statistics.pl data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge", "\n";
	    print STDERR "Script will run after Filter_fastq, MosaikBuild, MosaikAlign, MosaikSort, MosaikMerge and MosaikText sbatch scripts have been completed", "\n";

	}
    }
    else {
	
	if ($pMoM eq 0) {

	    print STDERR "Sbatch script calulate_coverage_statistics.pl data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
	    print STDERR "Script will run after MosaikBuild, MosaikAlign, MosaikSort and MosaikText sbatch scripts have been completed", "\n";
	}
	else {
	    print STDERR "Sbatch script calulate_coverage_statistics.pl data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge", "\n";
	    print STDERR "Script will run after MosaikBuild, MosaikAlign, MosaikSort, MosaikMerge and MosaikText sbatch scripts have been completed", "\n";

	}
    }
}

if ($pRCP eq 1 && $pCR) { #Run Rcovplot scripts after calculate_coverage_statistics.p. Rscripts:  Average_cov_chr.R, Coverage_bed_pileup_genome.R,Coverage_hist_by_Chr.R  
    
    print STDERR "\nRcovplots", "\n";
    
    if ($pFQF eq 0) {    
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    RcoveragePlotsPlatform($sid[$sampleid]);
	    
	}
    }
    else {
	
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	    
	    RcoveragePlotsFilter($sid[$sampleid]);
	    
	}
    }
    print STDERR "Creating sbatch script Rcovplots and writing script file(s) to: ", $ods,"/sampleid/mosaik/rcovplots_wf_sampleid.sh", "\n";
    
    if ($pFQF eq 1) {
	
	if ($pMoM eq 0) {

	    print STDERR "Sbatch script rcovplots data files will be written to: ", $odf,"/sampleid/mosaik/coverageReport", "\n";
	    print STDERR "Script will run after Filter_fastq, MosaikBuild, MosaikAlign, MosaikSort, MosaikText, calculate_coverage_statistics sbatch scripts have been completed", "\n";
	}
	else {
	    
	    print STDERR "Sbatch script rcovplots data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge", "\n";
	    print STDERR "Script will run after Filter_fastq, MosaikBuild, MosaikAlign, MosaikSort, MosaikMerge, MosaikText and calulate_coverage_statistics sbatch scripts have been completed", "\n";

	}
    }
    else {
	
	if ($pMoM eq 0) {

	    print STDERR "Sbatch script rcovplots data files will be written to: ", $odf,"/sampleid/mosaik/", "\n";
	    print STDERR "Script will run after MosaikBuild, MosaikAlign, MosaikSort, MosaikText and calulate_coverage_statistics sbatch scripts have been completed", "\n";
	}
	else {
	    print STDERR "Sbatch script rcovplots data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge", "\n";
	    print STDERR "Script will run after MosaikBuild, MosaikAlign, MosaikSort, MosaikMerge, MosaikText and calulate_coverage_statistics sbatch scripts have been completed", "\n";

	}
    }
}

if ($pGZ eq 1) { #Run Gzip
    
    print STDERR "\nGzip for fastq files", "\n";
    print STDERR "Creating sbatch script gzip and writing script file(s) to: ", $ods,"/sampleid/gzip/gzip_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script gzip data files will be written to: ", $odf,"/sampleid/fastq/", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Gzipfastq($sid[$sampleid]);
	GzipMoA($sid[$sampleid]);
	
    }

}

if ($pRMMOB eq 1) { #Run removal of mosaikBuild files
    
    print STDERR "\nRemoval of mosaikBuild files", "\n";
    print STDERR "Creating sbatch script rmmob and writing script file(s) to: ", $ods,"/sampleid/mosaik/rmmosaikBuild_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script rmmob data files will be deleted in: ", $odf,"/sampleid/mosaik/", "\n";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Rmmobfn($sid[$sampleid]);
	
    }

}

######################
###Sub Routines#######
######################

sub Rmmobfn {
    
#Generates a sbatch script, which removes on mosaikbuild files
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the fastqc_filtered folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the fastqc_filtered script directory
    $filename = "$ods/$_[0]/mosaik/rmmosaikBuild_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (RMMOAB, ">$filename") or die "Can't write to $filename: $!\n";
    
    print RMMOAB "#! /bin/bash -l", "\n";
    print RMMOAB "#SBATCH -A ", $aid, "\n";
    print RMMOAB "#SBATCH -n 1", "\n";
    print RMMOAB "#SBATCH -C thin", "\n";
    print RMMOAB "#SBATCH -t 00:05:00", "\n";
    print RMMOAB "#SBATCH -J RMMoB_", $_[0], "\n";
    print RMMOAB "#SBATCH -e $odf/$_[0]/mosaik/info/rmmosaikBuild_$_[0].", $fnt ,".stderr.txt", "\n";
    print RMMOAB "#SBATCH -o $odf/$_[0]/mosaik/info/rmmosaikBuild_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print RMMOAB "#SBATCH --mail-type=All", "\n";
	print RMMOAB "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print RMMOAB 'echo "Running on: $(hostname)"',"\n\n";
    print RMMOAB "cd $odf/$_[0]/mosaik", "\n\n";
    print RMMOAB "#Samples", "\n\n";
    print RMMOAB 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n\n";
    
    if ( $pFQF eq 1 ) { #If run fastq_filter.pl
	
	for (my $infile=0;$infile < scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++ ) { #All mosaikAlignemnt files
	    
	    my $tempinfile = $MosBAllInfiles{$_[0]}[$infile]; #Read 1 and read 2
	    print RMMOAB "rm  ", '${inSampleDir}', "/$tempinfile", ".dat", "\n\n";
	    
	}
	print RMMOAB "wait", "\n";
    }
    else { 
	
	for (my $infile=0;$infile < scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #All mosaikAlignemnt files
	    
	    my $tempinfile = $MosBInfiles{$_[0]}[$infile]; #Read 1 and read 2
	    print RMMOAB "rm  ", '${inSampleDir}', "/$tempinfile", ".dat", "\n\n";
	    
	}
	print RMMOAB "wait", "\n";
    }
    close(RMMOAB);
    return;
}

sub GzipMoA {
    
#Generates a sbatch script and runs gzip on mosaikaligned reads
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/gzip/info;`; #Creates the fastqc_filtered folder and info data file directory
    `mkdir -p $ods/$_[0]/gzip;`; #Creates the fastqc_filtered script directory
    $filename = "$ods/$_[0]/gzip/gzipMoA_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(3*scalar( @{ $MosBInfiles{$_[0]} })); #One full lane aligned with MosaikAligned takes approx. 3 h for gzip to process, round up to nearest full hour.

    open (GZMOA, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GZMOA "#! /bin/bash -l", "\n";
    print GZMOA "#SBATCH -A ", $aid, "\n";
    print GZMOA "#SBATCH -p node -n 8", "\n";
    print GZMOA "#SBATCH -C thin", "\n";
    print GZMOA "#SBATCH -t $t:00:00", "\n";
    print GZMOA "#SBATCH -J GZMoA_", $_[0], "\n";
    print GZMOA "#SBATCH -e $odf/$_[0]/gzip/info/gzipMoA_$_[0].", $fnt ,".stderr.txt", "\n";
    print GZMOA "#SBATCH -o $odf/$_[0]/gzip/info/gzipMoA_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GZMOA "#SBATCH --mail-type=All", "\n";
	print GZMOA "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GZMOA 'echo "Running on: $(hostname)"',"\n\n";
    print GZMOA "cd $odf/$_[0]/mosaik", "\n\n";
    print GZMOA "#Samples", "\n\n";
    print GZMOA 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n\n";
    
    if ( $pFQF eq 1 ) { #If run fastq_filter.pl
    
	my $k=1;
	for (my $infile=0;$infile < scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++ ) { #All mosaikAlignemnt files
	    
	    if ($infile eq $k*8) { #Using only 8 cores
		
		print GZMOA "wait", "\n\n";
		$k=$k+1;
	    }
	    my $tempinfile = $MosBAllInfiles{$_[0]}[$infile]; #Read 1 and read 2
	    print GZMOA "gzip  ", '${inSampleDir}', "/$tempinfile", "_aligned.dat &", "\n\n";
	    
	}
	print GZMOA "wait", "\n";
    }
    else { 
	
	my $k=1;
	for (my $infile=0;$infile < scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #All mosaikAlignemnt files
	    
	    if ($infile eq $k*8) { #Using only 8 cores
		
		print GZMOA "wait", "\n\n";
		$k=$k+1;
	    }
	    my $tempinfile = $MosBInfiles{$_[0]}[$infile]; #Read 1 and read 2
	    print GZMOA "gzip  ", '${inSampleDir}', "/$tempinfile", "_aligned.dat &", "\n\n";
	    
	}
	print GZMOA "wait", "\n";
    }
    close(GZMOA);
    return;
}

sub Gzipfastq { 
    
#Generates sbatch scripts for gziping, which will not be started but can be used whenever needed.

    `mkdir -p $odf/$_[0]/gzip/info;`; #Creates the gzip folder and info data file directory
    `mkdir -p $ods/$_[0]/gzip;`; #Creates the gzip script folder
    
    $filename = "$ods/$_[0]/gzip/gzipFastq_wf_$_[0].";
    Checkfnexists($filename, $fnend);
    
    my $t = ceil(1.5*scalar( @{ $infiles{$_[0]} })); #One full lane on Hiseq takes approx. 1.5 h for gzip to process, round up to nearest full hour.
    open (GZFASTQ, ">$filename") or die "Can't write to $filename: $!\n";
    
    print GZFASTQ "#! /bin/bash -l", "\n";
    print GZFASTQ "#SBATCH -A ", $aid, "\n";
    print GZFASTQ "#SBATCH -p node -n 8", "\n";
    print GZFASTQ "#SBATCH -C thin", "\n";	
    print GZFASTQ "#SBATCH -t $t:00:00", "\n";
    print GZFASTQ "#SBATCH -J GZFQ", $_[0], "\n";
    print GZFASTQ "#SBATCH -e $odf/$_[0]/gzip/info/gzipFastq_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print GZFASTQ "#SBATCH -o $odf/$_[0]/gzip/info/gzipFastq_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print GZFASTQ "#SBATCH --mail-type=All", "\n";
	print GZFASTQ "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print GZFASTQ 'echo "Running on: $(hostname)"',"\n\n";
    print GZFASTQ "#Samples", "\n";
    print GZFASTQ "cd $indirpath{$_[0]}", "\n\n";
    print GZFASTQ 'inSampleDir="',"$indirpath{$_[0]}", '"', "\n\n";

    my $k=1;
    
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	
	if ($infile eq $k*8) { #Using only 8 cores
	    
	    print GZFASTQ "wait", "\n\n";
	    $k=$k+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";
	print GZFASTQ "gzip ", '${inSampleDir}', "/$tempinfile"," &", "\n\n";
	
    }
    print GZFASTQ "wait", "\n\n";
    return;
}

sub RcoveragePlotsFilter { 
    
#Generates sbatch scripts for R scripts:
#1. Average_cov_chr.RAverage_cov_chr.R
#2. Coverage_bed_pileup_genome.R
#3. Coverage_hist_by_Chr.R
# on files generated from calculate_coverage_statistics.2.0.pl (Filter_fastq.pl)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    
    $filename = "$ods/$_[0]/mosaik/rcovplots_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (RCovP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print RCovP "#! /bin/bash -l", "\n";
    print RCovP "#SBATCH -A ", $aid, "\n";
    print RCovP "#SBATCH -n 1 ", "\n";
    print RCovP "#SBATCH -C thin", "\n";	
    print RCovP "#SBATCH -t 00:15:00", "\n"; 
    print RCovP "#SBATCH -J RcovPlots", $_[0], "\n";
    print RCovP "#SBATCH -e $odf/$_[0]/mosaik/info/rCovplots_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print RCovP "#SBATCH -o $odf/$_[0]/mosaik/info/rCovplots_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print RCovP "#SBATCH --mail-type=All", "\n";
	print RCovP "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print RCovP 'echo "Running on: $(hostname)"',"\n\n";
    print RCovP "module load bioinfo-tools", "\n\n"; 
    print RCovP "module load R/2.12.2", "\n\n";
    print RCovP "#Samples", "\n";
    
    if ($pMoM eq 0) {

	print RCovP 'inSampleDir="',"$odf/$_[0]/mosaik/coverageReport", '"', "\n";
	print RCovP 'outSampleDir="', "$odf/$_[0]/mosaik/coverageReport", '"', "\n\n"; 
	
	for (my $infile=0;$infile<scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (filter_fastq files)
	    
	    my $tempinfile = $MosBAllInfiles{$_[0]}[$infile];
	    print RCovP "Rscript $ids/Average_cov_chr.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_per_chromosome.txt $tempinfile ", '${outSampleDir}', "\n\n";
	    print RCovP "Rscript $ids/Coverage_bed_pileup_genome.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_histogram_per_bed_and_genome.txt $tempinfile $xcov ", '${outSampleDir}', "\n\n";
	    print RCovP "Rscript $ids/Coverage_hist_by_Chr.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_histogram_per_chromosome.txt $tempinfile $xcov ", '${outSampleDir}', "\n\n";
	}
	print RCovP "wait", "\n\n";
    }
    else {

	`mkdir -p $odf/$_[0]/mosaik/mosaikMerge/coverageReport;`; #Creates the coverageReport folder 
	print RCovP 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n";
	print RCovP 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n\n";
	
	print RCovP "Rscript $ids/Average_cov_chr.R ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_S_merged.bam.coverage_per_chromosome.txt $_[0]_lanes_", @{ $lanes{$_[0]} }," ", '${outSampleDir}', "\n\n";
	print RCovP "Rscript $ids/Coverage_bed_pileup_genome.R ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_S_merged.bam.coverage_histogram_per_bed_and_genome.txt $_[0]_lanes_", @{ $lanes{$_[0]} }, " $xcov ", '${outSampleDir}', "\n\n";
	print RCovP "Rscript $ids/Coverage_hist_by_Chr.R ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_S_merged.bam.coverage_histogram_per_chromosome.txt $_[0]_lanes_", @{ $lanes{$_[0]} }, " $xcov ", '${outSampleDir}', "\n\n";
	print RCovP "wait", "\n\n";
    }
    
    return;
}

sub RcoveragePlotsPlatform { 
    
#Generates sbatch scripts for R scripts:
#1. Average_cov_chr.RAverage_cov_chr.R
#2. Coverage_bed_pileup_genome.R
#3. Coverage_hist_by_Chr.R
# on files generated from calculate_coverage_statistics.2.0.pl (Platform)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    
    $filename = "$ods/$_[0]/mosaik/rcovplots_wf_$_[0].";
    Checkfnexists($filename, $fnend);
    
    open (RCovP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print RCovP "#! /bin/bash -l", "\n";
    print RCovP "#SBATCH -A ", $aid, "\n";
    print RCovP "#SBATCH -n 1 ", "\n";
    print RCovP "#SBATCH -C thin", "\n";	
    print RCovP "#SBATCH -t 00:15:00", "\n"; 
    print RCovP "#SBATCH -J RcovPlots", $_[0], "\n";
    print RCovP "#SBATCH -e $odf/$_[0]/mosaik/info/rCovplots_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print RCovP "#SBATCH -o $odf/$_[0]/mosaik/info/rCovplots_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print RCovP "#SBATCH --mail-type=All", "\n";
	print RCovP "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print RCovP 'echo "Running on: $(hostname)"',"\n\n";
    print RCovP "module load bioinfo-tools", "\n\n"; 
    print RCovP "module load R/2.12.2", "\n\n";
    print RCovP "#Samples", "\n";
    
    if ($pMoM eq 0) {

	print RCovP 'inSampleDir="',"$odf/$_[0]/mosaik/coverageReport", '"', "\n";
	print RCovP 'outSampleDir="', "$odf/$_[0]/mosaik/coverageReport", '"', "\n\n"; 
	
	for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform files)
	    
	    my $tempinfile = $MosBInfiles{$_[0]}[$infile];
	 
	    print RCovP "Rscript $ids/Average_cov_chr.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_per_chromosome.txt $tempinfile ", '${outSampleDir}', "\n\n";
	    print RCovP "Rscript $ids/Coverage_bed_pileup_genome.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_histogram_per_bed_and_genome.txt $tempinfile $xcov ", '${outSampleDir}', "\n\n";
	    print RCovP "Rscript $ids/Coverage_hist_by_Chr.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_histogram_per_chromosome.txt $tempinfile $xcov ", '${outSampleDir}', "\n\n";
	}
	print RCovP "wait", "\n\n";
    }
    else {
	if (scalar( @{ $MosBInfiles{$_[0]} }) > 1) { #but only if there is more than one mosaikBuild file
	    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge/coverageReport;`; #Creates the coverageReport folder 
	    print RCovP 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n";
	    print RCovP 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n\n";
	    print RCovP "Rscript $ids/Average_cov_chr.R ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged.bam.coverage_per_chromosome.txt $_[0]_lanes_", @{ $lanes{$_[0]} }," ", '${outSampleDir}', "\n\n";
	    print RCovP "Rscript $ids/Coverage_bed_pileup_genome.R ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged.bam.coverage_histogram_per_bed_and_genome.txt $_[0]_lanes_", @{ $lanes{$_[0]} }, " $xcov ", '${outSampleDir}', "\n\n";
	    print RCovP "Rscript $ids/Coverage_hist_by_Chr.R ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged.bam.coverage_histogram_per_chromosome.txt $_[0]_lanes_", @{ $lanes{$_[0]} }, " $xcov ", '${outSampleDir}', "\n\n";
	    print RCovP "wait", "\n\n";
	}
	else { #Run from mosaik directory instead
	    print RCovP 'inSampleDir="',"$odf/$_[0]/mosaik/coverageReport", '"', "\n";
	    print RCovP 'outSampleDir="', "$odf/$_[0]/mosaik/coverageReport", '"', "\n\n"; 
	    
	    for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform files)
		
		my $tempinfile = $MosBInfiles{$_[0]}[$infile];
		
		print RCovP "Rscript $ids/Average_cov_chr.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_per_chromosome.txt $tempinfile ", '${outSampleDir}', "\n\n";
		print RCovP "Rscript $ids/Coverage_bed_pileup_genome.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_histogram_per_bed_and_genome.txt $tempinfile $xcov ", '${outSampleDir}', "\n\n";
		print RCovP "Rscript $ids/Coverage_hist_by_Chr.R ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam.coverage_histogram_per_chromosome.txt $tempinfile $xcov ", '${outSampleDir}', "\n\n";
	    }
	    print RCovP "wait", "\n\n";
	}
    return;
    }
}

sub Cal_cov_statPlatform { 
    
#Generates sbatch scripts and runs calculate_coverage_statistics.2.0.pl on files generated from MosaikText (Platform)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/cal_cov_stat.2.0_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(10*scalar( @{ $MosBAllInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 10 h for calculate_coverage_statistics.2.0.pl to process, round up to nearest full hour.
    
    open (CR, ">$filename") or die "Can't write to $filename: $!\n";
    
    print CR "#! /bin/bash -l", "\n";
    print CR "#SBATCH -A ", $aid, "\n";
    print CR "#SBATCH -p node -n 8 ", "\n";
    print CR "#SBATCH -C thin", "\n";
    if ($pMoM eq 0) {	
	print CR "#SBATCH -t 8:00:00", "\n";	
    }
    else{
	print CR "#SBATCH -t $t:00:00", "\n";	
    }	
    print CR "#SBATCH -J CR", $_[0], "\n";
    print CR "#SBATCH -e $odf/$_[0]/mosaik/info/cal_cov_stat.2.0_$_[0].", $fnt ,".stderr.txt", "\n";
    print CR "#SBATCH -o $odf/$_[0]/mosaik/info/cal_cov_stat.2.0_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print CR "#SBATCH --mail-type=All", "\n";
	print CR "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print CR 'echo "Running on: $(hostname)"',"\n\n";
    print CR "module load bioinfo-tools", "\n\n"; 
    print CR "module load samtools/0.1.8", "\n\n";
    print CR "#Samples", "\n";
    
    if ($pMoM eq 0) {
	
	`mkdir -p $odf/$_[0]/mosaik/coverageReport;`; #Creates the coverageReport folder

	print CR 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print CR 'outSampleDir="', "$odf/$_[0]/mosaik/coverageReport", '"', "\n\n"; 
	my $k=1;
	for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform files)
	    
	    if ($infile eq $k*8) { #Using only 8 cores
	    
		print CR "wait", "\n\n";
	    $k=$k+1;
	    }
	    my $tempinfile = $MosBInfiles{$_[0]}[$infile];
	    print CR "$ids/calculate_coverage_statistics.2.0.pl $rd/$crfasta $rd/$crbed ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam ", '${outSampleDir}', " ", '${outSampleDir}',"/$tempinfile", "_aligned_sorted_cal_cov_stat.2.0 &", "\n\n";
	}
	print CR "wait", "\n\n";
    }
    else {
	if (scalar( @{ $MosBInfiles{$_[0]} }) > 1) { #but only if there is more than one mosaikBuild file
	    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge/coverageReport;`; #Creates the coverageReport folder 
	    print CR 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge", '"', "\n";
	    print CR 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n\n";
	    print CR "$ids/calculate_coverage_statistics.2.0.pl $rd/$crfasta $rd/$crbed ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged.bam ", '${outSampleDir}' ," ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged_cal_cov_stat.2.0", "\n\n";
	    print CR "wait", "\n\n";
	}
	else {
	    `mkdir -p $odf/$_[0]/mosaik/coverageReport;`; #Creates the coverageReport folder
	    
	    print CR 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	    print CR 'outSampleDir="', "$odf/$_[0]/mosaik/coverageReport", '"', "\n\n"; 
	    
	    my $k=1;
	    for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild
		
		if ($infile eq $k*8) { #Using only 8 cores
		    
		    print CR "wait", "\n\n";
		    $k=$k+1;
		}
		
		my $tempinfile = $MosBInfiles{$_[0]}[$infile];
		print CR "$ids/calculate_coverage_statistics.2.0.pl $rd/$crfasta $rd/$crbed ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam ", '${outSampleDir}', " ", '${outSampleDir}',"/$tempinfile", "_cal_cov_stat.2.0 &", "\n\n";
		
	    }
	    print CR "wait", "\n\n";
	}
	
    }
    
    if ($pRCP eq 1) { #If run rcoverageplots.R
	
	$filename2 = "$ods/$_[0]/mosaik/rcovplots_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print CR "sbatch $filename2", "\n\n";
	print CR "wait", "\n\n";
    }
    close(CR);
    return;
}

sub Cal_cov_statFilter { 
    
#Generates sbatch scripts and runs calculate_coverage_statistics.2.0.pl on files generated from MosaikText (Filter_fastq.pl)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/cal_cov_stat.2.0_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(10*scalar( @{ $MosBAllInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 10 h for calculate_coverage_statistics.2.0.pl to process, round up to nearest full hour.
    
    open (CR, ">$filename") or die "Can't write to $filename: $!\n";
    
    print CR "#! /bin/bash -l", "\n";
    print CR "#SBATCH -A ", $aid, "\n";
    print CR "#SBATCH -p node -n 8 ", "\n";
    print CR "#SBATCH -C thin", "\n";	 
    if ($pMoM eq 0) {	
	print CR "#SBATCH -t 8:00:00", "\n";	
    }
    else{
	print CR "#SBATCH -t $t:00:00", "\n";	
    }
    print CR "#SBATCH -J CR", $_[0], "\n";
    print CR "#SBATCH -e $odf/$_[0]/mosaik/info/cal_cov_stat.2.0_$_[0].", $fnt ,".stderr.txt", "\n";
    print CR "#SBATCH -o $odf/$_[0]/mosaik/info/cal_cov_stat.2.0_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print CR "#SBATCH --mail-type=All", "\n";
	print CR "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print CR 'echo "Running on: $(hostname)"',"\n\n";
    print CR "module load bioinfo-tools", "\n\n"; 
    print CR "module load samtools/0.1.8", "\n\n";
    print CR "#Samples", "\n";
    
    if ($pMoM eq 0) {
	
	`mkdir -p $odf/$_[0]/mosaik/coverageReport;`; #Creates the coverageReport folder
	
	print CR 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print CR 'outSampleDir="', "$odf/$_[0]/mosaik/coverageReport", '"', "\n\n"; 
	
	my $k=1;
	for (my $infile=0;$infile<scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild
	    
	    if ($infile eq $k*8) { #Using only 8 cores
		
		print CR "wait", "\n\n";
		$k=$k+1;
	    }
	    
	    my $tempinfile = $MosBAllInfiles{$_[0]}[$infile];
	    print CR "$ids/calculate_coverage_statistics.2.0.pl $rd/$crfasta $rd/$crbed ", '${inSampleDir}', "/$tempinfile","_aligned_sorted.bam ", '${outSampleDir}', " ", '${outSampleDir}',"/$tempinfile", "_cal_cov_stat.2.0 &", "\n\n";
	    
	}
	print CR "wait", "\n\n";
    }
    else {
	`mkdir -p $odf/$_[0]/mosaik/mosaikMerge/coverageReport;`; #Creates the coverageReport folder 
	print CR 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge", '"', "\n";
	print CR 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n\n";
	print CR "$ids/calculate_coverage_statistics.2.0.pl $rd/$crfasta $rd/$crbed ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_S_merged.bam ", '${outSampleDir}', " ", '${outSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged_cal_cov_stat.2.0", "\n\n";
	print CR "wait", "\n\n";
    }
    
    if ($pRCP eq 1) { #If run rcoverageplots.R
	
	$filename2 = "$ods/$_[0]/mosaik/rcovplots_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print CR "sbatch $filename2", "\n\n";
	print CR "wait", "\n\n";
    }
    close(CR);
    return;
}

sub MosaikTextFilter { 
    
#Generates sbatch scripts and runs MosaikText on files generated from MosaikSort or MosaikMerge (Filter_fastq.pl)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaikMerge folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/mosaikText_wf_$_[0].";
    Checkfnexists($filename, $fnend);
    my $t = ceil(0.5*scalar( @{ $MosBAllInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 0.5 h for MosaikText to process, round up to nearest full hour.
    
    open (MosT, ">$filename") or die "Can't write to $filename: $!\n";
    
    print MosT "#! /bin/bash -l", "\n";
    print MosT "#SBATCH -A ", $aid, "\n";
    print MosT "#SBATCH -p node -n 8 ", "\n";
    print MosT "#SBATCH -C thin", "\n";	
    print MosT "#SBATCH -t $t:00:00", "\n";
    print MosT "#SBATCH -J MoT", $_[0], "\n";
    print MosT "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikText_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosT "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikText_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print MosT "#SBATCH --mail-type=All", "\n";
	print MosT "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print MosT 'echo "Running on: $(hostname)"',"\n\n";
    print MosT "module load bioinfo-tools", "\n\n"; 
    print MosT "module load mosaik-aligner/1.0.1388", "\n\n";
    print MosT "mkdir -p /scratch/mosaik_tmp", "\n";
    print MosT "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
    print MosT "#Samples", "\n";
    
    if ($pMoM eq 0) {
	
	print MosT 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print MosT 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n"; 
	
	for (my $infile=0;$infile<scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild
	    
	    my $tempinfile = $MosBAllInfiles{$_[0]}[$infile];
	    print MosT "MosaikText -in ", '${inSampleDir}', "/$tempinfile", "_aligned_sorted",'.', "dat -bam ", '${outSampleDir}', "/$tempinfile","_aligned_sorted.bam", "\n\n";
	}
	print MosT "wait", "\n\n";
    }
    else {
	
	print MosT 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge", '"', "\n";
	print MosT 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge", '"', "\n\n";
	print MosT "MosaikText -in ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_S_merged.dat -bam ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_S_merged.bam", "\n\n";
	print MosT "wait", "\n\n";
    }
    
    if ($pCR eq 1) { #If run calculate_coverage_statsistics.pl
	
	#Requires only 1 sbatch since all files will be merged
	$filename2 = "$ods/$_[0]/mosaik/cal_cov_stat.2.0_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print MosT "sbatch $filename2", "\n\n";
	print MosT "wait", "\n\n";
    }
    close(MosT);
    return;
}

sub MosaikTextPlatform { 
    
#Generates sbatch scripts and runs MosaikText on files generated from MosaikSort or MosaikMerge (Platform)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaikMerge folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/mosaikText_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(1.5*scalar( @{ $MosBInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 1.5 h for MosaikText to process, round up to nearest full hour.
    
    open (MosT, ">$filename") or die "Can't write to $filename: $!\n";
    
    print MosT "#! /bin/bash -l", "\n";
    print MosT "#SBATCH -A ", $aid, "\n";
    print MosT "#SBATCH -p node -n 8 ", "\n";
    print MosT "#SBATCH -C thin", "\n";	
    print MosT "#SBATCH -t $t:00:00", "\n";
    print MosT "#SBATCH -J MoT", $_[0], "\n";
    print MosT "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikText_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosT "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikText_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print MosT "#SBATCH --mail-type=All", "\n";
	print MosT "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print MosT 'echo "Running on: $(hostname)"',"\n\n";
    print MosT "module load bioinfo-tools", "\n\n"; 
    print MosT "module load mosaik-aligner/1.0.1388", "\n\n";
    print MosT "mkdir -p /scratch/mosaik_tmp", "\n";
    print MosT "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
    print MosT "#Samples", "\n";
    
    if ($pMoM eq 0) {
	
	print MosT 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print MosT 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n"; 
	
	for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform files)
	    
	    my $tempinfile = $MosBInfiles{$_[0]}[$infile];
	    print MosT "MosaikText -in ", '${inSampleDir}', "/$tempinfile", "_aligned_sorted",'.', "dat -bam ", '${outSampleDir}', "/$tempinfile","_aligned_sorted.bam", "\n\n";
	}
	print MosT "wait", "\n\n";
    }
    else {

	if (scalar( @{ $MosBInfiles{$_[0]} }) > 1) { #but only if there is more than one mosaikBuild file
	    print MosT 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge", '"', "\n";
	    print MosT 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge", '"', "\n\n";
	    print MosT "MosaikText -in ", '${inSampleDir}',"/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged.dat -bam ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged.bam", "\n\n";
	    print MosT "wait", "\n\n";
	}
	else {
	    print MosT 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	    print MosT 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n";
	    for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform files)
		
		my $tempinfile = $MosBInfiles{$_[0]}[$infile];
		print MosT "MosaikText -in ", '${inSampleDir}', "/$tempinfile", "_aligned_sorted",'.', "dat -bam ", '${outSampleDir}', "/$tempinfile","_aligned_sorted.bam", "\n\n";
	    }
	}
    }
    if ($pCR eq 1) { #If run calculate_coverage_statsistics.pl
	
	#Requires only 1 sbatch since all files will be merged
	$filename2 = "$ods/$_[0]/mosaik/cal_cov_stat.2.0_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print MosT "sbatch $filename2", "\n\n";
	print MosT "wait", "\n\n";
    }
    close(MosT);
    return;
}

sub MosaikMergePlatform { 
    
#Generates sbatch scripts and runs MosaikMerge on files generated from MosaikSort (Platform)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaikMerge folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/mosaikMerge_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(1.5*scalar( @{ $MosBInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 1.5 h for MosaikMerge to process, round up to nearest full hour.

    open (MosM, ">$filename") or die "Can't write to $filename: $!\n";
    
    print MosM "#! /bin/bash -l", "\n";
    print MosM "#SBATCH -A ", $aid, "\n";
    print MosM "#SBATCH -p node -n 8 ", "\n";
    print MosM "#SBATCH -C thin", "\n";	
    print MosM "#SBATCH -t $t:00:00", "\n";
    print MosM "#SBATCH -J MoM", $_[0], "\n";
    print MosM "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikMerge_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosM "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikMerge_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print MosM "#SBATCH --mail-type=All", "\n";
	print MosM "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print MosM 'echo "Running on: $(hostname)"',"\n\n";
    print MosM "module load bioinfo-tools", "\n\n"; 
    print MosM "module load mosaik-aligner/1.0.1388", "\n\n";
    print MosM "mkdir -p $odf/$_[0]/mosaik/mosaikMerge/mosaik_tmp", "\n";
    print MosM "export MOSAIK_TMP=$odf/$_[0]/mosaik/mosaikMerge/mosaik_tmp", "\n\n";
    print MosM "#Samples", "\n";
    print MosM 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
    print MosM 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge", '"', "\n\n";
    print MosM "MosaikMerge "; 	
    
    for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform files)
	
	my $tempinfile = $MosBInfiles{$_[0]}[$infile];
	print MosM "-in ", '${inSampleDir}', "/$tempinfile", "_aligned_sorted",'.', "dat ";
	
    }
    
    print MosM "-out ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_merged.dat", "\n\n";
    print MosM "wait", "\n\n";
    
    if ($pMoT eq 1) { #If run MosaikText
	
	#Requires only 1 sbatch since all files will be merged
	$filename2 = "$ods/$_[0]/mosaik/mosaikText_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print MosM "sbatch $filename2", "\n\n";
	print MosM "wait", "\n\n";
    }
    close(MosM);
    return;
}

sub MosaikMergeFilter { 
    
#Generates sbatch scripts and runs MosaikMerge on files generated from MosaikSort
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaikMerge folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/mosaikMerge_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(1.5*scalar( @{ $MosBAllInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 1.5 h for MosaikMerge to process, round up to nearest full hour.
    open (MosM, ">$filename") or die "Can't write to $filename: $!\n";
    
    print MosM "#! /bin/bash -l", "\n";
    print MosM "#SBATCH -A ", $aid, "\n";
    print MosM "#SBATCH -p node -n 8 ", "\n";
    print MosM "#SBATCH -C thin", "\n";	
    print MosM "#SBATCH -t $t:00:00", "\n";
    print MosM "#SBATCH -J MoM", $_[0], "\n";
    print MosM "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikMerge_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosM "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikMerge_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print MosM "#SBATCH --mail-type=All", "\n";
	print MosM "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print MosM 'echo "Running on: $(hostname)"',"\n\n";
    print MosM "module load bioinfo-tools", "\n\n"; 
    print MosM "module load mosaik-aligner/1.0.1388", "\n\n";
    print MosM "mkdir -p $odf/$_[0]/mosaik/mosaikMerge/mosaik_tmp", "\n";
    print MosM "export MOSAIK_TMP=$odf/$_[0]/mosaik/mosaikMerge/mosaik_tmp", "\n\n";
    print MosM "#Samples", "\n";
    print MosM 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
    print MosM 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge", '"', "\n\n";
    print MosM "MosaikMerge "; 	
    
    for (my $infile=0;$infile<scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild filter_fastq.pl
	
	my $tempinfile = $MosBAllInfiles{$_[0]}[$infile];
	print MosM "-in ", '${inSampleDir}', "/$tempinfile", "_aligned_sorted",'.', "dat ";
	
    }
    
    print MosM "-out ", '${outSampleDir}', "/$_[0]_", "lanes_", @{ $lanes{$_[0]} }, "_S_merged.dat", "\n\n";
    print MosM "wait", "\n\n";
    
    if ($pMoT eq 1) { #If run MosaikText
	
	#Requires only 1 sbatch since all files will be merged
	$filename2 = "$ods/$_[0]/mosaik/mosaikText_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print MosM "sbatch $filename2", "\n\n";
	print MosM "wait", "\n\n";
    }
    close(MosM);
    return;
}

sub MosaikSortFilter { 
    
#Generates sbatch scripts and runs MosaikSort on files generated from MosaikAlign

    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    my $k=0;
    for (my $infile=0;$infile<scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild filter_fastq.pl files
	
	$filename = "$ods/$_[0]/mosaik/mosaikSort_wf_$_[0]_$k.";
	Checkfnexists($filename, $fnend);	
	
	open (MosS, ">$filename") or die "Can't write to $filename: $!\n";
	
	print MosS "#! /bin/bash -l", "\n";
	print MosS "#SBATCH -A ", $aid, "\n";
	print MosS "#SBATCH -p node -n 8 ", "\n";
	print MosS "#SBATCH -C thin", "\n";	
	print MosS "#SBATCH -t 10:00:00", "\n";
	print MosS "#SBATCH -J MoS", "$_[0]_",$k, "\n";
	print MosS "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikSort_$_[0]_", $k,".$fnt.stderr.txt", "\n";
	print MosS "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikSort_$_[0]_", $k,".$fnt.stdout.txt", "\n";
	
	unless ($em eq 0) {
	    
	    print MosS "#SBATCH --mail-type=All", "\n";
	    print MosS "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print MosS 'echo "Running on: $(hostname)"',"\n\n";
	print MosS "module load bioinfo-tools", "\n\n"; 
	print MosS "module load mosaik-aligner/1.0.1388", "\n\n";
	print MosS "mkdir -p /scratch/mosaik_tmp", "\n";
	print MosS "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	print MosS "#Samples", "\n";
	print MosS 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print MosS 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n";

	if ($pMoDup eq 1) {  #If 1 run mosaikDupSnoop (PCR removal)
	    print MosS "#Setup shm dirs", "\n";
	    print MosS "mkdir -p /dev/shm/$_[0]/mosaik/mosaikDupSnoop", "\n";
	    print MosS 'fragdata="', "/dev/shm/$_[0]/mosaik/mosaikDupSnoop", '"', "\n\n";	
	    print MosS "#Setup dupSnoop dir", "\n";
	    print MosS "mkdir -p $odf/$_[0]/mosaik/mosaikDupSnoop/$k.$fnt", "\n";
	    print MosS 'outDupS="', "$odf/$_[0]/mosaik/mosaikDupSnoop/$k.$fnt",'"', "\n\n";
	
	    my $tempinfile = $MosBAllInfiles{$_[0]}[$infile];
	    print MosS "MosaikDupSnoop -in ", '${inSampleDir}', "/$tempinfile", "_aligned",'.', "dat -od ", '${fragdata}', "/ -afl -rmm", "\n\n";
	    print MosS "cp -f ", '${fragdata}', "/", '.', "db* ", '${outDupS}', "/", "\n\n";
	    print MosS "MosaikSort -in ", '${inSampleDir}', "/$tempinfile", "_aligned",'.', "dat -out ", '${outSampleDir}', "/$MosBAllInfiles{$_[0]}[$infile]_aligned_sorted", '.', "dat -dup " , '${fragdata}', "/ -afl -sa -rmm", "\n\n";
	}	
	else {  #No PCR duplicates removal
	    my $tempinfile = $MosBAllInfiles{$_[0]}[$infile];
	    print MosS "MosaikSort -in ", '${inSampleDir}', "/$tempinfile", "_aligned",'.', "dat -out ", '${outSampleDir}', "/$MosBAllInfiles{$_[0]}[$infile]_aligned_sorted", '.', "dat -afl -sa -rmm", "\n\n";
	}
	$k=$k+1; #Tracks nr of sbatch scripts
    }
	
    close(MosS);
    return;
}

sub MosaikSortPlatform {
#Generates sbatch scripts and runs MosaikSort on files generated from MosaikAlign
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    my $k=0;
    for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform files)
	
	$filename = "$ods/$_[0]/mosaik/mosaikSort_wf_$_[0]_$k.";
	Checkfnexists($filename, $fnend);
	
	open (MosS, ">$filename") or die "Can't write to $filename: $!\n";
	
	print MosS "#! /bin/bash -l", "\n";
	print MosS "#SBATCH -A ", $aid, "\n";
	print MosS "#SBATCH -p node -n 8 ", "\n";
	print MosS "#SBATCH -C thin", "\n";
	print MosS "#SBATCH -t 10:00:00", "\n";
	print MosS "#SBATCH -J MoS", "$_[0]_",$k, "\n";
	print MosS "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikSort_$_[0]_", $k,".$fnt.stderr.txt", "\n";
	print MosS "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikSort_$_[0]_", $k,".$fnt.stdout.txt", "\n";
	
	unless ($em eq 0) {
	    
	    print MosS "#SBATCH --mail-type=All", "\n";
	    print MosS "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print MosS 'echo "Running on: $(hostname)"',"\n\n";
	print MosS "module load bioinfo-tools", "\n\n"; 
	print MosS "module load mosaik-aligner/1.0.1388", "\n\n";
	print MosS "mkdir -p /scratch/mosaik_tmp", "\n";
	print MosS "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	print MosS "#Samples", "\n";
	print MosS 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print MosS 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n";

	if ($pMoDup eq 1) { #If 1 run mosaikDupSnoop (PCR removal) 
	    print MosS "#Setup shm dirs", "\n";
	    print MosS "mkdir -p /dev/shm/$_[0]/mosaik/mosaikDupSnoop", "\n";
	    print MosS 'fragdata="', "/dev/shm/$_[0]/mosaik/mosaikDupSnoop", '"', "\n\n";
	    print MosS "#Setup dupSnoop dir", "\n";
	    print MosS "mkdir -p $odf/$_[0]/mosaik/mosaikDupSnoop/$k.$fnt", "\n";
	    print MosS 'outDupS="', "$odf/$_[0]/mosaik/mosaikDupSnoop/$k.$fnt",'"', "\n\n"; 
		
	    my $tempinfile = $MosBInfiles{$_[0]}[$infile];
	    print MosS "MosaikDupSnoop -in ", '${inSampleDir}', "/$tempinfile", "_aligned",'.', "dat -od ", '${fragdata}', "/ -afl -rmm", "\n\n";
	    print MosS "cp -f ", '${fragdata}', "/", '.', "db* ", '${outDupS}', "/", "\n\n"; #Copy database from dev/shm
	    print MosS "MosaikSort -in ", '${inSampleDir}', "/$tempinfile", "_aligned",'.', "dat -out ", '${outSampleDir}', "/$MosBInfiles{$_[0]}[$infile]_aligned_sorted", '.', "dat -dup " , '${fragdata}', "/ -afl -sa -rmm", "\n\n";
	}
	else { #No PCR duplicates removal 
	    my $tempinfile = $MosBInfiles{$_[0]}[$infile];
	    print MosS "MosaikSort -in ", '${inSampleDir}', "/$tempinfile", "_aligned",'.', "dat -out ", '${outSampleDir}', "/$MosBInfiles{$_[0]}[$infile]_aligned_sorted", '.', "dat -afl -sa -rmm", "\n\n";
	}
	$k=$k+1; #Tracks nr of sbatch scripts
    }
    close(MosS);
    return;
}

sub MosaikAlignPlatform {

#Generates sbatch scripts for MosaikAlign on files generated from MosaikBuild
#Used when reads have been trimmed using fastq_filter.pl

    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    my $k=0;
    for (my $infile=0;$infile<scalar( @{ $MosBInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild (platform reads)
	
	$filename = "$ods/$_[0]/mosaik/mosaikAlign_wf_$_[0]_$k.";
	Checkfnexists($filename, $fnend);
	
	open (MosA, ">$filename") or die "Can't write to $filename: $!\n";
	
	print MosA "#! /bin/bash -l", "\n";
	print MosA "#SBATCH -A ", $aid, "\n";
	print MosA "#SBATCH -p node -n 8 ", "\n";
	print MosA "#SBATCH -C thin", "\n";
	print MosA "#SBATCH -t 35:00:00", "\n";
	print MosA "#SBATCH -J MoA", "$_[0]_",$k, "\n";
	print MosA "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikAlign_$_[0]_", $k, ".$fnt.stderr.txt", "\n";
	print MosA "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikAlign_$_[0]_", $k, ".$fnt.stdout.txt", "\n";
	
	unless ($em eq 0) {
	    
	    print MosA "#SBATCH --mail-type=All", "\n";
	    print MosA "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print MosA 'echo "Running on: $(hostname)"',"\n\n";
	print MosA "module load bioinfo-tools", "\n\n"; 
	print MosA "module load mosaik-aligner/1.0.1388", "\n\n";
	print MosA "mkdir -p /scratch/mosaik_tmp", "\n";
	print MosA "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	print MosA "#Samples", "\n";
	print MosA 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print MosA 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n";
	print MosA "#Reference archive", "\n";
	print MosA 'referenceArchive="', "$rd/concat.dat", '"', "\n\n";
	
	my $tempinfile = $MosBInfiles{$_[0]}[$infile];
	print MosA "MosaikAligner -in ", '${inSampleDir}', "/$tempinfile",'.', "dat -out ",'${outSampleDir}', "/$MosBInfiles{$_[0]}[$infile]_aligned", '.', "dat -ia ", '${referenceArchive}', " -hs 15 -mm 4 -mhp 100 -act 35 -bw 35 -j $rd/concat_jdb_15 -p 8", "\n\n";
	
	$k=$k+1; #Tracks nr of sbatch scripts
    }
    close(MosA);
    return;
}

sub MosaikAlignFilter { 
    
#Generates sbatch scripts for MosaikAlign on files generated from MosaikBuild
#Used when reads have been trimmed using fastq_filter.pl
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the mosaik script directory
	
    my $k=0; #tracks infile
    
    for (my $infile=0;$infile<scalar( @{ $MosBAllInfiles{$_[0]} } );$infile++) { #For all files from MosaikBuild after filt_fastq.pl
	
	$filename = "$ods/$_[0]/mosaik/mosaikAlign_wf_$_[0]_$k.";
	Checkfnexists($filename, $fnend);

	open (MosA, ">$filename") or die "Can't write to $filename: $!\n";
	
	print MosA "#! /bin/bash -l", "\n";
	print MosA "#SBATCH -A ", $aid, "\n";
	print MosA "#SBATCH -p node -n 8 ", "\n";
	print MosA "#SBATCH -C thin", "\n";
	print MosA "#SBATCH -t 35:00:00", "\n";
	print MosA "#SBATCH -J MoA", "$_[0]_",$k, "\n";
	print MosA "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikAlign_$_[0]_", $k,".$fnt.stderr.txt", "\n";
	print MosA "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikAlign_$_[0]_", $k,".$fnt.stdout.txt", "\n";
	
	unless ($em eq 0) {
	    
	    print MosA "#SBATCH --mail-type=All", "\n";
	    print MosA "#SBATCH --mail-user=$em", "\n\n";
	    
	}
	
	print MosA 'echo "Running on: $(hostname)"',"\n\n";
	print MosA "module load bioinfo-tools", "\n\n"; 
	print MosA "module load mosaik-aligner/1.0.1388", "\n\n";
	print MosA "mkdir -p /scratch/mosaik_tmp", "\n";
	print MosA "export MOSAIK_TMP=/scratch/mosaik_tmp", "\n\n";
	print MosA "#Samples", "\n";
	print MosA 'inSampleDir="',"$odf/$_[0]/mosaik", '"', "\n";
	print MosA 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n";
	print MosA "#Reference archive", "\n";
	print MosA 'referenceArchive="', "$rd/concat.dat", '"', "\n\n";  
	
	my $tempinfile = $MosBAllInfiles{$_[0]}[$infile];
	print MosA "MosaikAligner -in ", '${inSampleDir}', "/$tempinfile",'.', "dat -out ",'${outSampleDir}', "/$MosBAllInfiles{$_[0]}[$infile]_aligned", '.', "dat -ia ", '${referenceArchive}', " -hs 15 -mm 4 -mhp 100 -act 35 -bw 35 -j $rd/concat_jdb_15 -p 8", "\n\n";
	
	$k=$k+1; #Tracks nr of sbatch scripts
    }
    close(MosA);
    return;
}


sub MosaikBuild {
    
#Generates a sbatch script and runs MosaikBuild on reads from Illumina platform or fastq_filtered reads
#Starts mosaikAlign from sbatch when mosaikBuild is finished for all files and submits mosaikSort with dependency afterrok of mosaikAlign  
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/mosaikBuild_wf_$_[0].";
    Checkfnexists($filename, $fnend);
    
    my $t = ceil(0.5*scalar( @{ $fqfiltInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 0.5 h for MosaikBuild to process, round up to nearest full hour.
    
    open (MosB, ">$filename") or die "Can't write to $filename: $!\n";
    
    print MosB "#! /bin/bash -l", "\n";
    print MosB "#SBATCH -A ", $aid, "\n";
    print MosB "#SBATCH -p node -n 8 ", "\n";
    print MosB "#SBATCH -C thin", "\n";
    print MosB "#SBATCH -t $t:00:00", "\n";
    print MosB "#SBATCH -J MoB", $_[0], "\n";
    print MosB "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikBuild_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosB "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikBuild_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print MosB "#SBATCH --mail-type=All", "\n";
	print MosB "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print MosB 'echo "Running on: $(hostname)"',"\n\n";
    print MosB "module load bioinfo-tools", "\n\n"; 
    print MosB "module load mosaik-aligner/1.0.1388", "\n\n";
    print MosB 'outSampleDir="', "$odf/$_[0]/mosaik", '"', "\n\n";
    
    my $k=1;
    my $allp=0; #Required to portion out 8 files before wait and to track the MosB outfiles to correct lane
    
    if ($pFQF eq 1) { #Reads filtered
	
	print MosB 'inSampleDir="',"$odf/$_[0]/fastq_filter", '"', "\n\n";
	
	for (my $infile=0;$infile<( scalar( @{ $fqfiltInfiles{$_[0]} }) -1);$infile++) { #Files from filter_fastq.pl (Paired reads)
	    
	    if ($allp eq $k*8) { #Using only 8 cores
		
		print MosB "wait", "\n\n";
		$k=$k+1;
	    }
	    my $tempinfile = $fqfiltInfiles{$_[0]}[$infile];
	    my $tempinfile2 = $fqfiltInfiles{$_[0]}[ ($infile+1)]; #Paired read
	    $infile = $infile+1; #To correct for reading 2 files at once
	    print MosB "MosaikBuild -q ", '${inSampleDir}', "/$tempinfile -q2 ", '${inSampleDir}', "/$tempinfile2 -out ",'${outSampleDir}', "/$MosBInfiles{$_[0]}[ ($allp)]", '.', "dat -st illumina &", "\n\n";
	    $allp++; #Track nr of printed so that wait can be printed between hashes
	}
	for (my $infile=0;$infile< scalar( @{ $fqfiltSInfiles{$_[0]} });$infile++) { #Files from filter_fastq.pl (singel reads)
	    $allp++;
	    
	    if ($allp eq $k*8+1) { #Using only 8 cores
		
		print MosB "wait", "\n\n";
		$k=$k+1;
	    }
	    my $tempinfile = $fqfiltSInfiles{$_[0]}[$infile]; #Any Singels created from fastq_filter.pl
	    print MosB "MosaikBuild -q ", '${inSampleDir}', "/$tempinfile -out ",'${outSampleDir}', "/$MosBSInfiles{$_[0]}[ ($infile)]", '.', "dat -st illumina &", "\n\n";
	    
	}
	print MosB "wait", "\n\n";
	
	if ( $pMoA eq 1) {
	    
	    for (my $infile=0;$infile < ( scalar( @{ $MosBAllInfiles{$_[0]} }) );$infile++) { 
		
		$filename2 = "$ods/$_[0]/mosaik/mosaikAlign_wf_$_[0]_$infile.";
		Checkfn2exists($filename2, $fnend);
		print MosB "MosAJOB_$infile=`sbatch $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id when submitted. Enables mosaikSort to sbatch script to be initiated when mosaikAlign sbatch script is finished. 
		
	    }
	}
	if ($pMoS eq 1 ) {
	    
	    for (my $infile=0;$infile < ( scalar( @{ $MosBAllInfiles{$_[0]} }) );$infile++) {
		
		$filename2 = "$ods/$_[0]/mosaik/mosaikSort_wf_$_[0]_$infile.";
		Checkfn2exists($filename2, $fnend);
		print MosB "MosSJOB_$infile=`sbatch --dependency=afterok:",'$MosAJOB_', "$infile $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id of mosAlign
		#print MosB "wait", "\n\n";
	    }
	}
	if ($pMoM eq 1) { #If run MosaikMerge
	    
	    $filename2 = "$ods/$_[0]/mosaik/mosaikMerge_wf_$_[0].";
	    Checkfn2exists($filename2, $fnend);
	    print MosB "sbatch --dependency=afterok:";
	    for (my $infile=0;$infile < ( scalar( @{ $MosBAllInfiles{$_[0]} }) );$infile++) { #Supply job_id to wait for mosaikSort finishes all files
		
		if ($infile eq scalar( @{ $MosBAllInfiles{$_[0]} } )-1 ) { #Last job_id finish with whitespace
		    print MosB '$MosSJOB_', "$infile "
		}
		else {
		    print MosB '$MosSJOB_', "$infile:"
		}
	    }
	    print MosB $filename2, "\n\n";
	}
	if ($pMoT eq 1 && $pMoM eq 0) { #If run MosaikText and not MosaikMerge
	    
	    for (my $infile=0;$infile < ( scalar( @{ $MosBAllInfiles{$_[0]} }) );$infile++) {
		$filename2 = "$ods/$_[0]/mosaik/mosaikText_wf_$_[0].";
		Checkfn2exists($filename2, $fnend);
		print MosB "sbatch --dependency=afterok:",'$MosSJOB_', "$infile $filename2 ", "\n\n"; #Catch sbtach job_id of mosSort
		#print MosS "wait", "\n\n";
	    }
	}
    }
    else { #Reads straight from platform and no fastq_filter
	
	print MosB 'inSampleDir="', $indirpath{$_[0]}, '"', "\n\n";
	
	for (my $infile=0;$infile<( scalar( @{ $infiles{$_[0]} }) -1);$infile++) {
	    
	    if ($allp eq $k*8) { #Using only 8 cores
		
		print MosB "wait", "\n\n";
		$k=$k+1;
	    }
	    my $tempinfile = $infiles{$_[0]}[$infile];
	    chomp $tempinfile; #Removing "\n";
	    my $tempinfile2 = $infiles{$_[0]}[ ($infile+1)]; #Paired read
	    chomp $tempinfile2; #Removing "\n";
	    $infile = $infile+1; #To correctfor reading 2 files at once
	    print MosB "MosaikBuild -q ", '${inSampleDir}', "/$tempinfile -q2 ", '${inSampleDir}', "/$tempinfile2 -out ", 
	    '${outSampleDir}', "/$MosBInfiles{$_[0]}[ ($allp)]",'.',"dat -st illumina &", "\n\n";
	    $allp++; #Track nr of printed so that wait can be printed between hashes
	}
	print MosB "wait", "\n\n";
	
	if ( $pMoA eq 1) {
	    
	    for (my $infile=0;$infile < ( scalar( @{ $MosBInfiles{$_[0]} }) );$infile++) {
		$filename2 = "$ods/$_[0]/mosaik/mosaikAlign_wf_$_[0]_$infile.";
		Checkfn2exists($filename2, $fnend);
		print MosB "MosAJOB_$infile=`sbatch $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id when submitted. Enables mosaikSort to sbatch script to be initiated when mosaikAlign sbatch script is finished.
	    }
	}
	if ($pMoS eq 1) {
	    
	    for (my $infile=0;$infile < ( scalar( @{ $MosBInfiles{$_[0]} }) );$infile++) {
		$filename2 = "$ods/$_[0]/mosaik/mosaikSort_wf_$_[0]_$infile.";
		Checkfn2exists($filename2, $fnend);
		print MosB "MosSJOB_$infile=`sbatch --dependency=afterok:",'$MosAJOB_', "$infile $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id of mosAlign		
	    }
	}
	if ($pMoM eq 1) { #If run MosaikMerge 
	    
	    if (scalar( @{ $MosBInfiles{$_[0]} }) > 1) { #but only if there is more than one mosaikBuild file
		$filename2 = "$ods/$_[0]/mosaik/mosaikMerge_wf_$_[0].";
		Checkfn2exists($filename2, $fnend);
		print MosB "sbatch --dependency=afterok:";
		for (my $infile=0;$infile < ( scalar( @{ $MosBInfiles{$_[0]} }) );$infile++) { #Supply job_id to wait for mosaikSort finishes all files
		    
		    if ($infile eq scalar( @{ $MosBInfiles{$_[0]} } )-1 ) { #Last job_id finish with whitespace
			print MosB '$MosSJOB_', "$infile "
		    }
		    else {
			print MosB '$MosSJOB_', "$infile:"
		    }
		    
		}
		print MosB $filename2, "\n\n";
	    }
	    elsif($pMoT eq 1) { #Run MosaikText instead
		for (my $infile=0;$infile < ( scalar( @{ $MosBInfiles{$_[0]} }) );$infile++) {
		    $filename2 = "$ods/$_[0]/mosaik/mosaikText_wf_$_[0].";
		    Checkfn2exists($filename2, $fnend);
		    print MosB "sbatch --dependency=afterok:",'$MosSJOB_', "$infile $filename2", "\n\n"; #Catch sbtach job_id of mosSort
		}
	    }
	}
	if ($pMoT eq 1 && $pMoM eq 0) { #If run MosaikText and not MosaikMerge
	    
	    for (my $infile=0;$infile < ( scalar( @{ $MosBInfiles{$_[0]} }) );$infile++) {
		$filename2 = "$ods/$_[0]/mosaik/mosaikText_wf_$_[0].";
		Checkfn2exists($filename2, $fnend);
		print MosB "sbatch --dependency=afterok:",'$MosSJOB_', "$infile $filename2", "\n\n"; #Catch sbtach job_id of mosSort
	    }
	}
	
    }
    close(MosB);
    
    if ($pFQF eq 0) { #submit sbatch otherwise mosaikBuild will be started by filter_fastq sbatch script
	my $ret = `sbatch $filename`;
	my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
	print STDERR "Sbatch script submitted, job id: $jobID\n";
	print STDERR "To check status of job, please run \'jobinfo -j $jobID\'\n";
	print STDERR "To check status of job, please run \'squeue -j $jobID\'\n";
    }
    return;
}   


sub Fastqc_filtered {
	
#Generates a sbatch script and runs FASTQC on fastq_filetered reads
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/fastqc_filtered/info;`; #Creates the fastqc_filtered folder and info data file directory
    `mkdir -p $ods/$_[0]/fastqc_filtered;`; #Creates the fastqc_filtered script directory
    $filename = "$ods/$_[0]/fastqc_filtered/fastqc_filtered_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(0.5*scalar( @{ $fqfiltInfiles{$_[0]} })); #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.

    open (FASTQC2, ">$filename") or die "Can't write to $filename: $!\n";
    
    print FASTQC2 "#! /bin/bash -l", "\n";
    print FASTQC2 "#SBATCH -A ", $aid, "\n";
    print FASTQC2 "#SBATCH -p node -n 8 ", "\n";
    print FASTQC2 "#SBATCH -C thin", "\n";
    print FASTQC2 "#SBATCH -t $t:00:00", "\n";
    print FASTQC2 "#SBATCH -J FQC2_", $_[0], "\n";
    print FASTQC2 "#SBATCH -e $odf/$_[0]/fastqc_filtered/info/fastqc_filtered_$_[0].", $fnt ,".stderr.txt", "\n";
    print FASTQC2 "#SBATCH -o $odf/$_[0]/fastqc_filtered/info/fastqc_filtered_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print FASTQC2 "#SBATCH --mail-type=All", "\n";
	print FASTQC2 "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print FASTQC2 'echo "Running on: $(hostname)"',"\n\n";
    print FASTQC2 "cd $odf/$_[0]/fastq_filter", "\n\n";
    print FASTQC2 "module add bioinfo-tools", "\n\n";
    print FASTQC2 "module add FastQC/0.7.2", "\n\n";
    print FASTQC2 "#Samples", "\n";
    print FASTQC2 'outSampleDir="', "$odf/$_[0]/fastqc_filtered", '"', "\n\n";
    
    my $k=1;
    my $allp=0; #Required to portion out 8 files before wait
    
    for (my $infile=0;$infile<scalar( @{ $fqfiltInfiles{$_[0]} });$infile++) {
	
	if ($infile eq $k*8) { #Using only 8 cores
	    
	    print FASTQC2 "wait", "\n\n";
	    $k=$k+1;
	}
	my $tempinfile = $fqfiltInfiles{$_[0]}[$infile]; #Read 1 and read 2
	chomp $tempinfile; #Removing "\n";
	print FASTQC2 "fastqc ", $tempinfile, ' -o ${outSampleDir}';
	print FASTQC2 " &", "\n\n";
	$allp = $infile; #Track nr of printed so that wait can be printed between hashes 
	
    }
    for (my $infile=0;$infile<scalar( @{ $fqfiltSInfiles{$_[0]} });$infile++) { #Singels
	
	$allp++;
	
	if ($allp eq $k*8) { #Using only 8 cores
	    
	    print FASTQC2 "wait", "\n\n";
	    $k=$k+1;
	}
	my $tempinfile = $fqfiltSInfiles{$_[0]}[$infile]; #Any Singels
	chomp $tempinfile; #Removing "\n";
	print FASTQC2 "fastqc ", $tempinfile, ' -o ${outSampleDir}';
	print FASTQC2 " &", "\n\n";
	
    }
    print FASTQC2 "wait", "\n";
    close(FASTQC2);
    return;
}


sub Fastq_filter {
    
#Generates a sbatch script and runs Fastq_filter custom perl script. Starts FASTQC2 and MosaikBuild if enabled
#File format must be y.x.lane8_2.fastq and files must contain equal amounts of read 1 and 2. 
#$_[0] = sampleid
    
    `mkdir -p $odf/$_[0]/fastq_filter/info;`; #Creates the fastq_filter folder and info data file directory
    `mkdir -p $ods/$_[0]/fastq_filter;`; #Creates the fastq_filter script directory
   
    $filename = "$ods/$_[0]/fastq_filter/fastq_filter_wf_$_[0].";
    Checkfnexists($filename, $fnend);
  
    my $t;

    if (scalar( @{ $infiles{$_[0]} })/2 < 8 ) {# Divided by 2 since this is for both read1 and read 2.
	
	$t = 15; #One full lane on Hiseq takes approx. 10 h for Filter_fastq.pl to process, will round up to nearest full hour.
    }
    else {
	
	$t = 25; #If running more than 8 lanes at the same time
    }
    
    open (FASTQF, ">$filename") or die "Can't write to $filename: $!\n";
    
    print FASTQF "#! /bin/bash -l", "\n";
    print FASTQF "#SBATCH -A ", $aid, "\n";
    print FASTQF "#SBATCH -p node -n 8 ", "\n";
    print FASTQF "#SBATCH -C thin", "\n";    
    print FASTQF "#SBATCH -t $t:00:00", "\n";
    print FASTQF "#SBATCH -J FqFi", $_[0], "\n";
    print FASTQF "#SBATCH -e $odf/$_[0]/fastq_filter/info/fastq_filter_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print FASTQF "#SBATCH -o $odf/$_[0]/fastq_filter/info/fastq_filter_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print FASTQF "#SBATCH --mail-type=All", "\n";
	print FASTQF "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print FASTQF 'echo "Running on: $(hostname)"',"\n\n";
    print FASTQF "module load bioinfo-tools", "\n\n";
    print FASTQF "#Samples", "\n";
    print FASTQF 'outSampleDir="', "$odf/$_[0]/fastq_filter/", '"', "\n";
    print FASTQF 'infoDir="', '${outSampleDir}', 'info/"', "\n\n"; 
    
    my $k =1;

    for (my $infile=0;$infile<( scalar( @{ $infiles{$_[0]} }) -1);$infile++) {

	if ($infile eq $k*16) { #Using only 8 cores, 16 because we are supplying 2 infiles at once
	    
	    print FASTQF "wait", "\n\n";
	    $k=$k+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";
	my $tempinfile2 = $infiles{$_[0]}[ ($infile+1)]; #Paired read
	chomp $tempinfile2; #Removing "\n";
	$infile = $infile+1; #To correct for reading 2 files at once
	
	print FASTQF 'file1="', $indirpath{$_[0]}, "/$tempinfile", '"', "\n";
	print FASTQF 'file2="', $indirpath{$_[0]}, "/$tempinfile2",'"', "\n";	
	print FASTQF 'fstdout="', '${infoDir}', "fastq_filter_$_[0]_$tempinfile.stdout", '"', "\n";
	print FASTQF 'fstderr="', '${infoDir}', "fastq_filter_$_[0]_$tempinfile.stderr", '"', "\n";
	print FASTQF 'fcmdLne="', '${infoDir}', "fastq_filter_$_[0]_$tempinfile.cmdLine", '"', "\n";
	print FASTQF 'fullCmd="', $ids, "/filter_fastq.pl filter ", '${outSampleDir} ${file1} ${file2}','"', "\n";
	print FASTQF 'echo "${fullCmd} > ${fstdout} 2> ${fstderr}" > ${fcmdLne}', "\n";
	print FASTQF '${fullCmd} > ${fstdout} 2> ${fstderr} &', "\n\n";

    }

    print FASTQF "wait", "\n\n",;
    
    if ($pFQC2 eq 1) { #FASTQC after filter_fastq.pl
	
	$filename2 = "$ods/$_[0]/fastqc_filtered/fastqc_filtered_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print FASTQF "sbatch $filename2","\n\n";
	
    }
    if ($pMoB eq 1) {  #MosaikBuild after filter_fastq.pl
	
	$filename2 ="$ods/$_[0]/mosaik/mosaikBuild_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print FASTQF "sbatch $filename2","\n\n";
    }
    print FASTQF "wait", "\n";
    close(FASTQF);
    my $ret = `sbatch $filename`;
    my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
    print STDERR "Sbatch script submitted, job id: $jobID\n";
    print STDERR "To check status of job, please run \'jobinfo -j $jobID\'\n";
    print STDERR "To check status of job, please run \'squeue -j $jobID\'\n";
    return;
}

sub Fastqc {
    
#Generates a sbatch script and runs FASTQC
#$_[0] = sampleid

    `mkdir -p $odf/$_[0]/fastqc/info;`; #Creates the fastqc folder and info data file directory
    `mkdir -p $ods/$_[0]/fastqc;`; #Creates the fastqc script directory
    
    $filename = "$ods/$_[0]/fastqc/fastqc_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    my $t = ceil(0.5*scalar( @{ $infiles{$_[0]} })); #One full lane on Hiseq takes approx. 0.5 h for FASTQC to process, round up to nearest full hour.
    
    open (FASTQC, ">$filename") or die "Can't write to $filename: $!\n";
    
    print FASTQC "#! /bin/bash -l", "\n";
    print FASTQC "#SBATCH -A ", $aid, "\n";
    print FASTQC "#SBATCH -p node -n 8 ", "\n";
    print FASTQC "#SBATCH -C thin", "\n";
    print FASTQC "#SBATCH -t $t:00:00", "\n";
    print FASTQC "#SBATCH -J FQC", $_[0], "\n";
    print FASTQC "#SBATCH -e $odf/$_[0]/fastqc/info/fastqc_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print FASTQC "#SBATCH -o $odf/$_[0]/fastqc/info/fastqc_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print FASTQC "#SBATCH --mail-type=All", "\n";
	print FASTQC "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print FASTQC 'echo "Running on: $(hostname)"',"\n\n";
    print FASTQC "cd $indirpath{$_[0]}", "\n\n";
    print FASTQC "module add bioinfo-tools", "\n\n";
    print FASTQC "module add FastQC/0.7.2", "\n\n";
    print FASTQC "#Samples", "\n";
    print FASTQC 'outSampleDir="', "$odf/$_[0]/fastqc", '"', "\n\n";
    
    my $k=1;
    
    for (my $infile=0;$infile<scalar( @{ $infiles{$_[0]} });$infile++) {
	
	if ($infile eq $k*8) { #Using only 8 cores
	    
	    print FASTQC "wait", "\n\n";
	    $k=$k+1;
	}
	my $tempinfile = $infiles{$_[0]}[$infile];
	chomp $tempinfile; #Removing "\n";
	print FASTQC "fastqc ", $tempinfile, ' -o ${outSampleDir}';
	print FASTQC " &", "\n\n";

    }
    
    print FASTQC "wait", "\n";
    
    
    close(FASTQC);
    my $ret = `sbatch $filename`;
    my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
    print STDERR "Sbatch script submitted, job id: $jobID\n";
    print STDERR "To check status of job, please run \'jobinfo -j $jobID\'\n";
    print STDERR "To check status of job, please run \'squeue -j $jobID\'\n";
    return;
}

sub FilterReFormat {

#Code needed to reformat files which have not yet been created into correct format so that a sbatch script can be generated with the correct filenames. 
 
    for my $samplid ( keys %infiles ) { #For every sample id
	
	if ($pFQF eq 1) {

	    print STDERR "\nInput for FASTQC after filter_fastq.pl", "\n";
	    print STDERR "\nSample ID", "\t", $samplid,"\n";
	}
	my $k=1;
	my $itrack=0; #Needed to be able to insert an extra S.filter.fastq file in the hash{sampleid}[infile]
	for (my $i=0;$i<scalar( @ { $infiles{ $samplid } } );$i++) { #Collects inputfiles for every in dir after Fastqc_filtered and remakes format
	    
	    if ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/ ) { #Parse to filter_fastq.pl format
		
		$dataHR->{'inFile'}{'name'} = $1;
		$dataHR->{'inFile'}{'lane'} = $2;
		$dataHR->{'inFile'}{'dir0'} = $3;
	 
		$fqfiltInfiles{ $samplid }[$i]= join(".", @{$dataHR->{'inFile'}}{'name', 'lane', 'dir0'}, "filter.fastq"); #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet.
		if ($pFQF eq 1) {
		    print STDERR $fqfiltInfiles{ $samplid }[$i], "\n";
		}
	    }
	    if ($i eq $k*2-1 ) { #To print singlet reads filename ($k starts at 1 and $i at 0, hence -1)
		$fqfiltSInfiles{ $samplid }[$itrack]= join(".", @{$dataHR->{'inFile'}}{'name', 'lane'}, "S.filter.fastq"); #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet.
		if ($pFQF eq 1) {
		    print STDERR $fqfiltSInfiles{ $samplid }[$itrack], "\n"; #Store singels in separet hash
		}
		$itrack++; #Track for every lane finished
		$k=$k+1; #Multiple of *2-1
	    }
	     
	}
	
    }
    return;
}

sub MoSBReFormat {
    
#Code needed to reformat files for mosaik output, which have not yet been created into, correct format so that a sbatch script can be generated with the correct filenames. 
    
    for my $samplid ( keys %infiles ) { #For every sample id
	
	print STDERR "\nInput for MosaikAlign after MosaikBuild", "\n";
	print STDERR "\nSample ID", "\t", $samplid,"\n";
	my $k=1;
	my $itrack=0; #Needed to be able to insert an extra S.filter.fastq file in the hash{sampleid}[infile]
	for (my $i=0;$i<scalar( @ { $infiles{ $samplid } } );$i++) { #Collects inputfiles for every in dir after Fastqc_filtered and remakes format
	    if ( $infiles{$samplid}[$i] =~ /\/?([^\.\/]+\.[^\.]+)\.lane(\d+)_([12FfRr])\.fastq/ ) { #Parse to MosaikBuild output format
		
		$dataMosB->{'inFile'}{'name'}= $1;
		$dataMosB->{'inFile'}{'lane'} = $2;
		push( @ {$lanes{$samplid} }, $2);
		$MosBInfiles{ $samplid }[$itrack]= join(".", @{$dataMosB->{'inFile'}}{'name', 'lane'} ); #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet.
		$MosBSInfiles{ $samplid }[$itrack]= join(".", @{$dataMosB->{'inFile'}}{'name', 'lane'}, "S" ); #Save new format in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet.
		push (@ { $MosBAllInfiles{$samplid} }, $MosBInfiles{ $samplid }[$itrack]); #Stores all infiles
		push (@ { $MosBAllInfiles{$samplid} }, $MosBSInfiles{ $samplid }[$itrack]); #Stores all infiles
		print STDERR $MosBInfiles{ $samplid }[$itrack], ".dat", "\n"; #Note .dat is not part of hash entry
		if ( $pFQF eq 1) {

		    print STDERR $MosBSInfiles{ $samplid }[$itrack],".dat", "\n"; #Note .dat is not part of hash entry
		}
		$i++; #Skip second direction
		$itrack++; #Track for every lane finished
	    }
	    
	}
	
    }
    return;
}

sub Checkfnexists {
    
#$_[0] = complete filepath
#$_[1] = file ending

    my $fn;
    $fnt = 0; #Nr of sbatch with identical filenames
    for (my $i=0;$i<999;$i++) { #Number of possible files with the same name
	
	$fn = "$_[0]$i$_[1]"; #filename, filenr and fileending
	$fnt = $i; #Nr of sbatch with identical filenames, global variable
	if (-e $fn) { #if file exists 
	}
	else {
	    $i=999; #Exit loop
	}
	
    }
    $filename = $fn; #Transfer to global variable
    return;
}

sub Checkfn2exists {
    
#$_[0] = complete filepath
#$_[1] = file ending

    my $fn;
    for (my $i=0;$i<999;$i++) { #Number of possible files with the same name
	
	$fn = "$_[0]$i$_[1]"; #filename, filenr and fileending
	if (-e $fn) { #if file exists 
	}
	else {
	    $i=999; #Exit loop
	}
	
    }
    $filename2 = $fn; #Transfer to global variable
    return;
}

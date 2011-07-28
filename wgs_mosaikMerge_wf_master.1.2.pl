#!/usr/bin/perl -w

use strict;
use warnings;

#Master scripts for analysing whole genome paired end reads from Illumina pipleine using mosaik. The program merges mosaikSorted files and generates a coverage report.
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
wgs_mosaik_wf_master.1.0.pl  -i [infile...n] -a [project ID] -s [sample ID...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata] -oMoM [outfilemosaikMerge]

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s) (Mandatory: Supply whole path)

-ids/--indirscript The pipeline custom script in dir (Mandatory: Supply whole path)

-rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)

-a/--projectid The project ID (Mandatory)

-s/--sampleid The sample ID (Mandatory)

-em/--email

-odf/--outdirdata The data files output directory (Supply whole path, defaults to data)

-ods/--outdirscript The script files output directory (Supply whole path, defaults to wgs_wf_scripts)

-pMoM/--mosaikMerge Flag running MosaikMerge (defaults to yes (=1))

-oMoM/--outmerged Merged output file 

-pMoT/--mosaikText Flag running MosaikText (defaults to yes (=1))

-pCR/--calculate_coverage_statistics Flag running Calculate_coverage_statistics (defaults to yes (=1))
    
-pRCP/--rcovplots  Flag running rcovplots (defaults to yes (=1))
    
-xcov/--xcoverage  Flag determining x-scale coverage plot in rcovplots (defaults to "30")
    
-crbed/--crbedfile Bed file for calculate_coverage_statistics (defaults to "CCDS-hg19.bed")

-crfasta/--crfasta Fasta file for calculate_coverage_statistics (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")

=head3 I/O

Input format ( infiles_aligned_sorted_(merged).dat )

Output format:

1. Mosaik_lanes_merged.dat

2. Mosaik_lanes_merged.bam

3. Calculate_coverage_results.stdout

4. Coverage plots

=head4 Dependencies

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
	qq{wgs_mosaik_wf_master.1.0.pl -id [infile...n] -a [projectid] -s [sampleid...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata] -oMoM [outfilemosaikMerge]
	       -i/--infile Infile(s), comma sep (Mandatory: Supply whole path)
               -ids/--indirscript The pipeline custom script in dir (Mandatory: Supply whole path)
               -rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)
	       -a/--projectid The project ID (Mandatory)
	       -s/--sampleid The sample ID,comma sep (Mandatory)
	       -em/--email e-mail
	       -odf/--outdirdata The data files output directory (Supply whole path, defaults to data)
	       -ods/--outdirscript The script files output directory (Supply whole path, defaults to wgs_wf_scripts)
	       -pMoM/--mosaikMerge Flag running MosaikMerge (defaults to yes (=1))
	       -oMoM/--outmerged Merged output file
	       -pMoT/--mosaikText Flag running MosaikText (defaults to yes (=1))
	       -pCR/--calculate_coverage_statistics Flag running Calculate_coverage_statistics (defaults to yes (=1))
	       -pRCP/--rcovplots  Flag running rcovplots (defaults to yes (=1))
	       -xcov/--xcoverage  Flag determining x-scale coverage plot in rcovplots (defaults to "30")
               -crbed/--bedfile Bed file for calculate_coverage_statistics (defaults to "CCDS-hg19.bed")
               -crfasta/--crfasta Fasta file for calculate_coverage_statistics (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")
	   };
    
}

my ($aid,$em, $ids, $rd, $oMoM,$odf,$ods,$xcov, $crbed, $crfasta, $fnend, $filename, $filename2, $fnt, $fnt2, $help) = (0,0,0,0, "merged",0,0, 30, "CCDS-hg19.bed", "Homo_sapiens.GRCh37.57.dna.concat.fa", ".sh"); #Arguments for project
my ($pMoM, $pMoT, $pCR, $pRCP) = (1,1,1,1); #Default arguments for running programs
my (@infn,@sid);
my (%infiles);

GetOptions('i|infile:s'  => \@infn, #Comma separeted list
	   'ids|inscriptdir:s'  => \$ids, #Directory for custom scripts required by the pipeline
	   'rd|referencedir:s'  => \$rd, #directory containing subfolders of references
	   'a|projectid:s'  => \$aid,
	   's|sampleid:s'  => \@sid, #Comma separeted list, one below outdirdata
	   'em|email:s'  => \$em,
	   'odf|outdirdata:s'  => \$odf, #One above sample id
	   'ods|outdirscript:s'  => \$ods,
	   'pMoM|mosaikMerge:n' => \$pMoM,
	   'oMoM|outmerged:s'  => \$oMoM,
	   'pMoT|mosaikText:n' => \$pMoT,
	   'pCR|cal_cov_stat:n' => \$pCR,
	   'crbed|cal_cov_bedfile:n' => \$crbed,
	   'crfasta|cal_cov_fastafile:n' => \$crfasta,
	   'pRCP|rcovplots:n' => \$pRCP,
	   'xcov|xcoverage:n' => \$xcov,
	   'h|help' => \$help,
	   );

die $USAGE if( $help );

if (@infn == 0) {
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

@infn = split(/,/,join(',',@infn)); #Enables comma separated indir(s)
@sid = split(/,/,join(',',@sid)); #Enables comma separeted list of sample IDs


#########################
###Run program part######
#########################

if ($pMoM eq 1) {

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Mosaikmerge($sid[$sampleid]);
	
    }

    print STDERR "\nMosaikMerge", "\n";
    print STDERR "Creating sbatch script mosaikMerge and writing script file(s) to: ", $ods,"/sampleid/mosaik/mosaikMerge_mM_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script mosaikMerge data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge/", "\n";

}

if ($pMoT eq 1) { #Run mosaikText

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	MosaikText($sid[$sampleid]);
	
    }

    print STDERR "\nMosaikText", "\n";
    print STDERR "Creating sbatch script MosaikText and writing script file(s) to: ", $ods,"/sampleid/mosaik/mosaikText_mM_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script mosaikText data files will be written to: ", $odf,"/sampleid/mosaikMerge/", "\n";
print STDERR "Script will run after MosaikMerge sbatch scripts have been completed", "\n";
}

if ($pCR eq 1) { #Run  calculate_coverage_statistics.2.0.pl

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	Cal_cov_stat($sid[$sampleid]);
	
    }

    print STDERR "\ncalculate_coverage_statistics.2.0.pl", "\n";
    print STDERR "Creating sbatch script calculate_coverage_statistics.2.0.pl and writing script file(s) to: ", $ods,"/sampleid/mosaik/cal_cov_stat.2.0_mM_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script calculate_coverage_statistics.2.0.pl data files will be written to: ", $odf,"/sampleid/mosaikMerge/", "\n";
    print STDERR "Script will run after MosaikMerge and MosaikText sbatch scripts have been completed", "\n";
}

if ($pRCP eq 1 && $pCR) { #Run Rcovplot scripts after calculate_coverage_statistics.p. Rscripts:  Average_cov_chr.R, Coverage_bed_pileup_genome.R,Coverage_hist_by_Chr.R  
    
    print STDERR "\nRcovplots", "\n";    
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	
	RcoveragePlots($sid[$sampleid]);
	
    }
    
    print STDERR "Creating sbatch script Rcovplots and writing script file(s) to: ", $ods,"/sampleid/mosaik/rcovplots_mM_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script rcovplots data files will be written to: ", $odf,"/sampleid/mosaik/mosaikMerge/coverageReport", "\n";
    print STDERR "Script will run after MosaikMerge, MosaikText and calulate_coverage_statistics sbatch scripts have been completed", "\n";
    
}


######################
###Sub Routines#######
######################

sub RcoveragePlots { 
    
#Generates sbatch scripts for R scripts:
#1. Average_cov_chr.RAverage_cov_chr.R
#2. Coverage_bed_pileup_genome.R
#3. Coverage_hist_by_Chr.R
# on files generated from calculate_coverage_statistics.2.0.pl
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge/coverageReport;`; #Creates the mosaik output directory

    $filename = "$ods/$_[0]/mosaik/rcovplots_mM_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (RCovP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print RCovP "#! /bin/bash -l", "\n";
    print RCovP "#SBATCH -A ", $aid, "\n";
    print RCovP "#SBATCH -n 1 ", "\n";
    print RCovP "#SBATCH -C thin", "\n";	
    print RCovP "#SBATCH -t 00:15:00", "\n"; 
    print RCovP "#SBATCH -J RcovPlots", $_[0], "\n";
    print RCovP "#SBATCH -e $odf/$_[0]/mosaik/info/rCovplots_mM_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print RCovP "#SBATCH -o $odf/$_[0]/mosaik/info/rCovplots_mM_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print RCovP "#SBATCH --mail-type=All", "\n";
	print RCovP "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print RCovP 'echo "Running on: $(hostname)"',"\n\n";
    print RCovP "module load bioinfo-tools", "\n\n"; 
    print RCovP "module load R/2.12.2", "\n\n";
    print RCovP "#Samples", "\n";
    print RCovP 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n";
    print RCovP 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n\n";    

    print RCovP "Rscript $ids/Average_cov_chr.R ", '${inSampleDir}',"/$oMoM.bam.coverage_per_chromosome.txt $oMoM.bam ", '${outSampleDir}', "\n\n";
    print RCovP "Rscript $ids/Coverage_bed_pileup_genome.R ", '${inSampleDir}', "/$oMoM.bam.coverage_histogram_per_bed_and_genome.txt $oMoM.bam $xcov ", '${outSampleDir}', "\n\n";
    print RCovP "Rscript $ids/Coverage_hist_by_Chr.R ", '${inSampleDir}', "/$oMoM.bam.coverage_histogram_per_chromosome.txt $oMoM.bam $xcov ", '${outSampleDir}', "\n\n";
    print RCovP "wait", "\n\n";
    
    close(RCovP);
    return;
}

sub Cal_cov_stat { 
    
#Generates sbatch scripts and runs calculate_coverage_statistics.2.0.pl on files generated from MosaikText
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge/coverageReport;`; #Creates the mosaik output directory
    
    $filename = "$ods/$_[0]/mosaik/cal_cov_stat.2.0_mM_wf_$_[0].";
    Checkfnexists($filename, $fnend);
    
    my $t = ceil(8*scalar(@infn)); #One full lane on Hiseq takes approx. 8 h for calculate_coverage_statistics.2.0.pl to process, round up to nearest full hour.
    
    open (CR, ">$filename") or die "Can't write to $filename: $!\n";
    
    print CR "#! /bin/bash -l", "\n";
    print CR "#SBATCH -A ", $aid, "\n";
    print CR "#SBATCH -p node -n 8 ", "\n";
    print CR "#SBATCH -C thin", "\n";
    print CR "#SBATCH -t $t:00:00", "\n";	
    print CR "#SBATCH -J CR", $_[0], "\n";
    print CR "#SBATCH -e $odf/$_[0]/mosaik/info/cal_cov_stat.2.0_mM_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print CR "#SBATCH -o $odf/$_[0]/mosaik/info/cal_cov_stat.2.0_mM_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print CR "#SBATCH --mail-type=All", "\n";
	print CR "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print CR 'echo "Running on: $(hostname)"',"\n\n";
    print CR "module load bioinfo-tools", "\n\n"; 
    print CR "module load samtools/0.1.8", "\n\n";
    print CR "#Samples", "\n";
    print CR 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge", '"', "\n";
    print CR 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge/coverageReport", '"', "\n\n";
    print CR "$ids/calculate_coverage_statistics.2.0.pl $rd/$crfasta $rd/$crbed ", '${inSampleDir}',"/$oMoM.bam ", '${outSampleDir}' ," ", '${outSampleDir}', "/$oMoM.bam_cal_cov_stat.2.0", "\n\n";
    print CR "wait", "\n\n";
    
    
    if ($pRCP eq 1) { #If run rcoverageplots.R
	
	$filename2 = "$ods/$_[0]/mosaik/rcovplots_mM_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print CR "sbatch $filename2", "\n\n";
	print CR "wait", "\n\n";
    }
    close(CR);
    return;
}

sub MosaikText { 
    
#Generates sbatch scripts and runs MosaikText on files generated from MosaikMerge (Filter_fastq.pl)
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaikMerge folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/mosaikText_mM_wf_$_[0].";
    Checkfnexists($filename, $fnend);
    my $t = ceil(7*scalar(@infn) ); #One full lane on Hiseq takes approx. 7 h for MosaikText to process, round up to nearest full hour.
    
    open (MosT, ">$filename") or die "Can't write to $filename: $!\n";
    
    print MosT "#! /bin/bash -l", "\n";
    print MosT "#SBATCH -A ", $aid, "\n";
    print MosT "#SBATCH -p node -n 8 ", "\n";
    print MosT "#SBATCH -C thin", "\n";	
    print MosT "#SBATCH -t $t:00:00", "\n";
    print MosT "#SBATCH -J MoT", $_[0], "\n";
    print MosT "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikText_mM_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosT "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikText_mM_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
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
    print MosT 'inSampleDir="',"$odf/$_[0]/mosaik/mosaikMerge", '"', "\n";
    print MosT 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge", '"', "\n\n";

    print MosT "MosaikText -in ", '${inSampleDir}',"/$oMoM.dat -bam ", '${outSampleDir}', "/$oMoM.bam", "\n\n";
    print MosT "wait", "\n\n";
    
    
    if ($pCR eq 1) { #If run calculate_coverage_statsistics.pl
	
	#Requires only 1 sbatch since all files will be merged
	$filename2 = "$ods/$_[0]/mosaik/cal_cov_stat.2.0_mM_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print MosT "sbatch $filename2", "\n\n";
	print MosT "wait", "\n\n";
    }
    close(MosT);
    return;
}

sub Mosaikmerge {
    
#Generates sbatch scripts and runs MosaikMerge on files generated from MosaikSort
    
    `mkdir -p $odf/$_[0]/mosaik/info;`; #Creates the mosaik folder and info data file directory
    `mkdir -p $odf/$_[0]/mosaik/mosaikMerge;`; #Creates the mosaikMerge folder and info data file directory
    `mkdir -p $ods/$_[0]/mosaik;`; #Creates the mosaik script directory
    
    $filename = "$ods/$_[0]/mosaik/mosaikMerge_mM_wf_$_[0].";
    Checkfnexists($filename, $fnend);

 my $t = ceil(5*scalar(@infn) ); #One full lane on Hiseq takes approx. 3 h for MosaikMerge to process, round up to nearest full hour.    
    open (MosM, ">$filename") or die "Can't write to $filename: $!\n";

    print MosM "#! /bin/bash -l", "\n";
    print MosM "#SBATCH -A ", $aid, "\n";
    print MosM "#SBATCH -p node -n 8 ", "\n";
    print MosM "#SBATCH -C thin", "\n";	
    print MosM "#SBATCH -t $t:00:00", "\n";
    print MosM "#SBATCH -J MoM", $_[0], "\n";
    print MosM "#SBATCH -e $odf/$_[0]/mosaik/info/mosaikMerge_mM_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print MosM "#SBATCH -o $odf/$_[0]/mosaik/info/mosaikMerge_mM_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
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
    print MosM 'outSampleDir="', "$odf/$_[0]/mosaik/mosaikMerge", '"', "\n\n";
    print MosM "MosaikMerge ";
        
    for (my $infile=0;$infile<scalar(@infn);$infile++) { #For all files
	
	my $tempinfile = $infn[$infile];
	print MosM "-in $tempinfile ";
    }
   
    print MosM "-out ", '${outSampleDir}', "/$oMoM.dat ", "\n\n";
    print MosM "wait", "\n\n";

    if ($pMoT eq 1) { #If run MosaikText
	
	#Requires only 1 sbatch since all files will be merged
	$filename2 = "$ods/$_[0]/mosaik/mosaikText_mM_wf_$_[0].";
	Checkfn2exists($filename2, $fnend);
	print MosM "sbatch $filename2", "\n\n";
	print MosM "wait", "\n\n";
    }
    close(MosM);
    my $ret = `sbatch $filename`;
    my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
    print STDERR "Sbatch script submitted, job id: $jobID\n";
    print STDERR "To check status of job, please run \'jobinfo -j $jobID\'\n";
    print STDERR "To check status of job, please run \'squeue -j $jobID\'\n";
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

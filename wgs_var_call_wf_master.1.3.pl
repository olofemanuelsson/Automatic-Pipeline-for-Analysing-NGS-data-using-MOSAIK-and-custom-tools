#!/usr/bin/perl -w

use strict;
use warnings;

#Master scripts for analysing whole genome paired end reads from Illumina pipleine using samtools, and custum perl scripts. The program performs variation calls on bam files and generates output in the vcf and hsd format (haplotype).
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
wgs_var_call_wf_master.1.3.pl  -i [infile...n] -a [project ID] -s [sample ID...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata] -ids [indirscripts] -rd [referencedir]

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s) (Mandatory: Supply whole path)

-ids/--indirscript The pipeline custom script in dir (Mandatory: Supply whole path)

-rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)

-rdav/--referencesdirannov Reference(s) dir for annovar (Supply whole path, defaults to /bubo/nobackup/uppnex/annotations/annovar/humandb)

-a/--projectid The project ID (Mandatory)

-s/--sampleid The sample ID(s) (Mandatory)

-em/--email

-odf/--outdirdata The data files output directory (Supply whole path, defaults to data)

-ods/--outdirscript The script files output directory (Supply whole path, defaults to wgs_wf_scripts)

-pSTSI/--samtools Flag running samtools sort & index (defaults to yes (=1))

-pSTV/--samtools Flag running samtools view (defaults to yes (=1))

-pSTMP/--samtools Flag running samtools mpileup (defaults to yes (=1))
    
-pCBSI/--perl Flag running convert_bam_to_snps_indels_sa.pl (defaults to yes (=1))

-pVCOMP/--perl Flag running varcall_compare_vcf_MB.1.1.pl (defaults to yes (=1))

-pVCMA/--perl Flag running varcall_merge_annovar_master.1.0.pl (defaults to yes (=1))

-pANVAR/--perl Flag running annovar (defaults to yes (=1))

-pVMERGE/--perl Flag running varcall_merge_post_annovar_master.1.0.pl (defaults to yes (=1))

-STMPRef/--Flag for setting reference in samtools mpileup (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")

-annovar_dbsnp_ver/--Flag for setting the version of dbsnp in annovar (defaults to snp131)

-annovar_1000g_ver/--Flag for setting the version of 1000g in annovar (defaults to 1000g2010nov_all)

-annovar_maf_thold/--Flag for setting the maf_threshold in annovar (defaults to 0)

=head3 I/O

Input format ( infiles_aligned_sorted_(merged).dat )

Output format:

1. Mosaik_lanes_merged_sorted.bam

2. Mosaik_lanes_merged_sorted.bai

3. Mosaik_lanes_merged_sorted_chrNo.bam

4. Mosaik_lanes_merged_sorted_chrNo_var.flt.vcf

5. Mosaik_lanes_merged_sorted_chrNo_var.raw.bcf

6. Mosaik_lanes_merged_sorted_chrNo_var_MB.txt

7. Varcall_comp.vcf.annovar_master_chrNo.txt

8. Varcall_comp.vcf.annovar_master_chrNo.txt.annovar_analysis.txt

9. Annovar_master_all_subject_variants.txt

=head4 Dependencies

1. convert_bam_to_snps_indels_sa.pl

2. varcall_compare_vcf_MB.1.1.pl

3. varcall_merge_annovar_master.1.0.pl

4. varcall_merge_post_annovar_master.1.0.pl

=cut

use Pod::Usage;
use Pod::Text;
use Getopt::Long;
use POSIX;
use File::Basename;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{wgs_mosaik_wf.pl -id [infile...n] -a [projectid] -s [sampleid...n] -em [e-mail] -ods [outdirscripts] -odf [outdirdata] -ids [indirscripts] -rd [referencedir]
	       -i/--infile Infile(s), comma sep (Mandatory: Supply whole path)
               -ids/--indirscript The pipeline custom script in dir (Mandatory: Supply whole path)
               -rd/--referencesdir Reference(s) dir (Mandatory: Supply whole path)
               -rdav/--referencesdirannov Reference(s) dir for annovar (Supply whole path, defaults to /bubo/nobackup/uppnex/annotations/annovar/humandb)
	       -a/--projectid The project ID (Mandatory)
	       -s/--sampleid The sample ID,comma sep (Mandatory)
	       -em/--email e-mail
	       -odf/--outdirdata The data files output directory (Supply whole path, defaults to data)
	       -ods/--outdirscript The script files output directory (Supply whole path, defaults to wgs_wf_scripts)
	       -pSTSI/--samtools Flag running samtools sort & index (defaults to yes (=1))
	       -pSTV/--samtools Flag running samtools view (defaults to yes (=1))
	       -pSTMP/--samtools Flag running samtools mpileup (defaults to yes (=1))
	       -pCBSI/--perl Flag running convert_bam_to_snps_indels_sa.pl (defaults to yes (=1))
	       -pVCOMP/--perl Flag running varcall_compare_vcf_MB.pl (defaults to yes (=1))
               -pVCMA/--perl Flag running varcall_merge_annovar_master.1.0.pl (defaults to yes (=1))
               -pANVAR/--perl Flag running annovar (defaults to yes (=1))
               -pVMERGE/--perl Flag running varcall_merge_post_annovar_master.1.0.pl (defaults to yes (=1))
               -STMPRef/--Flag for setting reference in samtools mpileup (defaults to "Homo_sapiens.GRCh37.57.dna.concat.fa")
               -annovar_dbsnp_ver/--Flag for setting the version of dbsnp in annovar (defaults to snp131)
               -annovar_1000g_ver/--Flag for setting the version of 1000g in annovar (defaults to 1000g2010nov_all)
               -annovar_maf_thold/--Flag for setting the maf_threshold in annovar (dbsnp, 1000g etc.)
	   };
    
}

my ($aid,$em, $ids, $rd, $rdav, $odf, $ods, $fnend, $annovar_dbsnp_ver, $annovar_1000g_ver, $annovar_maf_thold, $STMPRef, $filename, $filename2, $fnt, $fnt2, $help) = (0,0,0,0, "/bubo/nobackup/uppnex/annotations/annovar/humandb", "data", "wgs_wf_scripts", ".sh", "snp131", "1000g2010nov_all", 0, "Homo_sapiens.GRCh37.57.dna.concat.fa"); #Arguments for project
my ($pSTSI, $pSTV, $pSTMP, $pCBSI,$pVCOMP, $pVCMA, $pANVAR, $pVMERGE) = (1,1,1,1,1,1,1,1); #Default arguments for running programs
my (@infn,@sid, @chr, @avcov, @jobID);
my (%infiles,%avcovfn, %dirname);

GetOptions('i|infile:s'  => \@infn, #Comma separeted list
	   'ids|inscriptdir:s'  => \$ids, #Directory for custom scripts required by the pipeline
	   'rd|referencedir:s'  => \$rd, #directory containing references
	   'rdav|referencedirannov:s'  => \$rdav, #directory containing /humandb references
	   'a|projectid:s'  => \$aid,
	   's|sampleid:s'  => \@sid, #Comma separeted list, one below outdirdata
	   'em|email:s'  => \$em,
	   'odf|outdirdata:s'  => \$odf, #One above sample id
	   'ods|outdirscript:s'  => \$ods, #One above sample id
	   'pSTSI|samtoolssort_index:n' => \$pSTSI,
	   'pSTV|samtoolsview:n' => \$pSTV,
	   'pSTMP|samtoolsmpileup:n' => \$pSTMP,
	   'pCBSI|varcallmb:n' => \$pCBSI,
	   'pVCOMP|varcallcomp:n' => \$pVCOMP, #Union between VCF and MB format
	   'pVCMA|varcallcompmerge:n' => \$pVCMA, #Merges union files per chr across all subjects
	   'pANVAR|annovar:n' => \$pANVAR, #Performs annovar filter gene, region and filter analysis
	   'pVMERGE|varcallmergeannov:n' => \$pVMERGE, #Merges annovar analysis results to one master file
	   'STMPRef|samtoolsmpilupref:s' => \$STMPRef,
	   'annovar_dbsnp_ver|annovar_dbsnp_version:n' => \$annovar_dbsnp_ver,
	   'annovar_1000g_ver|annovar__1000g_version:n' => \$annovar_1000g_ver,
	   'annovar_maf_thold|annovar_maf_threshold:n' => \$annovar_maf_thold,
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
@chr = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"); #Chr for filtering of bam file


FileReFormat(); #removes .bam ending and extracts filename

#########################
###Run program part######
#########################

if ($pSTSI eq 1) {

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	SamtoolsSortIndex($sid[$sampleid],$sampleid);
	
    }
    print STDERR "\nSamtools sort & index", "\n";
    print STDERR "Creating sbatch script samtools and writing script file(s) to: ", $ods,"/sampleid/samtools/samtools_sort_index_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script samtools sort & index data files will be written to: ", $odf,"/sampleid/samtools/", "\n";
}

if ($pSTV eq 1) {

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	SamtoolsView($sid[$sampleid]);
	
    }
    print STDERR "\nSamtools view", "\n";
    print STDERR "Creating sbatch script samtools and writing script file(s) to: ", $ods,"/sampleid/samtools/samtools_view_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script samtools view data files will be written to: ", $odf,"/sampleid/samtools", "\n";
}

if ($pSTMP eq 1 && $pSTV) { #Run samtools mpileup only if samtools view has been run

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	SamtoolsMP($sid[$sampleid]);
	
    }
    print STDERR "\nSamtools mpileup", "\n";
    print STDERR "Creating sbatch script samtools and writing script file(s) to: ", $ods,"/sampleid/samtools/samtools_mpileup_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script samtools view data files will be written to: ", $odf,"/sampleid/samtools/varcalls", "\n";
}

if ($pCBSI eq 1 ) { #Run convert_bam_to_snps_indels_sa.pl

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
    
	ConvertBamMB($sid[$sampleid]);
	
    }
    print STDERR "\nConvert_bam_to_snps_indels_sa.pl", "\n";
    print STDERR "Creating sbatch script convert_bam_to_snps_indels_sa.pl and writing script file(s) to: ", $ods,"/sampleid/convert_bam_vc/convert_bam_vc_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script samtools view data files will be written to: ", $odf,"/sampleid/convert_bam_vc", "\n";
}

if ($pVCOMP eq 1 ) { #Run varcall_compare_vcf_MB.pl
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	
	VarcomVCF_MB($sid[$sampleid]);
	
    }
    print STDERR "\nVarcall_compare_vcf_MB.pl", "\n";
    print STDERR "Creating sbatch script varcall_compare_vcf_MB.pl and writing script file(s) to: ", $ods,"/sampleid/varcall_comp/varcall_comp_vcf_MB_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script samtools view data files will be written to: ", $odf,"/sampleid/varcall_comp", "\n";
}

if ($pVCMA eq 1 ) { #Run varcall_merge_annovar_master.1.0.pl, Merges all variants across subjects per chr
    
    VCallMergeAnnovar();

    print STDERR "\nvarcall_merge_annovar_master.1.0.pl", "\n";
    print STDERR "Creating sbatch script varcall_merge_annovar_master.1.0.pl and writing script file(s) to: ", $ods,"/annovar_filter/varcall_merge_annovar_master_wf_.sh", "\n";
    print STDERR "Sbatch script varcall_merge data files will be written to: ", $odf,"/annovar_filter", "\n";
    
}

if ($pANVAR eq 1 ) { #Run annovar
    
    AnnovarFilter();
    print STDERR "\nAnnovar Analysis", "\n";
    print STDERR "Creating sbatch script annovar and writing script file(s) to: ", $ods,"/sampleid/annovar/annovar_wf_sampleid.sh", "\n";
    print STDERR "Sbatch script annovar data files will be written to: ", $odf,"/sampleid/annovar_filter", "\n";
}

if ($pVMERGE eq 1 ) { #Run varcall_merge_post_annovar_master.1.0.pl, Merges all variants across all subjects
    
    VarcallMergePostAnnovar();

    print STDERR "\nVarcall_merge_post_annovar_master.1.0.pl", "\n";
    print STDERR "Creating sbatch script varcall_merge_post_annovar_master.1.0.pl and writing script file(s) to: ", $ods,"/annovar_filter/varcall_merge_post_annovar_master_wf.sh", "\n";
    print STDERR "Sbatch script varcall_merge data files will be written to: ", $odf,"/annovar_filter", "\n";
    
}

######################
###Sub Routines#######
######################

sub VarcallMergePostAnnovar { 

#Merges all variants for all chr (and subjects) to create a master file containing all annovar information
    
    `mkdir -p $odf/annovar_filter/info;`; #Creates the annovar_filter folder and info data file directory   
    `mkdir -p $ods/annovar_filter;`; #Creates the annovar_filter script folder
    
    $filename = "$ods/annovar_filter/varcall_merge_post_annovar_master_wf.";
    Checkfnexists($filename, $fnend);
    
    open (VMERGE, ">$filename") or die "Can't write to $filename: $!\n";
    
    print VMERGE "#! /bin/bash -l", "\n";
    print VMERGE "#SBATCH -A ", $aid, "\n";
    print VMERGE "#SBATCH -p node -n 1", "\n";
    print VMERGE "#SBATCH -C thin", "\n";	
    print VMERGE "#SBATCH -t 00:30:00", "\n";
    print VMERGE "#SBATCH -J VMERGE", "\n";
    print VMERGE "#SBATCH -e $odf/annovar_filter/info/varcall_merge_post_annovar_master_wf.", $fnt ,".stderr.txt", "\n";
    print VMERGE "#SBATCH -o $odf/annovar_filter/info/varcall_merge_post_annovar_master_wf.", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print VMERGE "#SBATCH --mail-type=All", "\n";
	print VMERGE "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print VMERGE 'echo "Running on: $(hostname)"',"\n\n";
    print VMERGE "#Samples", "\n";
    print VMERGE 'inSampleDir="',"$odf/annovar_filter", '"', "\n"; # All variants for all subjects have been merged here by chr
    print VMERGE 'inSamplePrefix="',"/varcall_comp.vcf.annovar_master_", '"', "\n";
    print VMERGE 'outSampleDir="', "$odf/annovar_filter", '"', "\n\n";
    
    print VMERGE "perl $ids/varcall_merge_post_annovar_master.1.0.pl -i ";
    
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr
	
	
	if ($chr eq scalar(@chr)-1) {
	    
	    print VMERGE '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ";
	}
	else {
	    print VMERGE '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt,";	
	    
	}
    }
    print VMERGE "-nos ", scalar(@sid), " -o ", '${outSampleDir}', "/annovar_master_all_subject_variants.txt", "\n\n";
#scalar(@sid) required for correct number of columns
    close(VMERGE);   
    return;
}

sub AnnovarFilter { 

#Filter SNVs by gene, region and databases
#Works on the varcall_comp.vcf.annovar_master_$chr.txt files where all variants per subject is present. Prints new files for each analysis

    `mkdir -p $odf/annovar_filter/info;`; #Creates the annovar folder and info data file directory   
    `mkdir -p $ods/annovar_filter;`; #Creates the annovar script folder
   
    $filename = "$ods/annovar_filter/annovar_master_filter_wf.";
    Checkfnexists($filename, $fnend);
    open (ANVARF, ">$filename") or die "Can't write to $filename: $!\n";
    
    print ANVARF "#! /bin/bash -l", "\n";
    print ANVARF "#SBATCH -A ", $aid, "\n";
    print ANVARF "#SBATCH -p node -n 8", "\n";
    print ANVARF "#SBATCH -C thin", "\n";	
    print ANVARF "#SBATCH -t 7:00:00", "\n";
    print ANVARF "#SBATCH -J ANVARMASTF", "\n";
    print ANVARF "#SBATCH -e $odf/annovar_filter/info/annovar_master_filter_wf.", $fnt ,".stderr.txt", "\n";
    print ANVARF "#SBATCH -o $odf/annovar_filter/info/annovar_master_filter_wf.", $fnt ,".stdout.txt", "\n";    
    unless ($em eq 0) {
	
	print ANVARF "#SBATCH --mail-type=All", "\n";
	print ANVARF "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print ANVARF 'echo "Running on: $(hostname)"',"\n\n";
    print ANVARF "module load bioinfo-tools", "\n\n"; 
    print ANVARF "module add annovar/2011.05.06", "\n\n";
    print ANVARF 'inRefDir="', $rdav, '"', "\n\n";
    print ANVARF "#Samples", "\n";
    print ANVARF 'inSampleDir="',"$odf/annovar_filter", '"', "\n"; # All variants for all subjects have been merged here by chr
    print ANVARF 'outSampleDir="', "$odf/annovar_filter", '"', "\n\n"; #All subjects   
    print ANVARF 'inSamplePrefix="',"/varcall_comp.vcf.annovar_master_", '"', "\n";    
    print ANVARF 'outSamplePrefix="',"/varcall_comp.vcf.annovar_master_", '"', "\n";   
    
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	print ANVARF "annovar/annotate_variation.pl -geneanno -buildver hg19 -dbtype refgene ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inRefDir} &', "\n\n";
	print ANVARF "annovar/annotate_variation.pl -regionanno -dbtype mce46way -buildver hg19 ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inRefDir} &', "\n\n"; #SNVs in conserved regions	
	print ANVARF "annovar/annotate_variation.pl -regionanno -dbtype segdup -buildver hg19 ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inRefDir} &', "\n\n"; #Finding SNVs in segmental duplications   
	print ANVARF "annovar/annotate_variation.pl -filter -buildver hg19 -dbtype $annovar_1000g_ver --maf_threshold $annovar_maf_thold ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inRefDir} &', "\n\n"; # SNVs in 1000g data 
	print ANVARF "annovar/annotate_variation.pl -filter -buildver hg19 -dbtype $annovar_dbsnp_ver --maf_threshold $annovar_maf_thold ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inRefDir} &', "\n\n"; # SNVs in dbsnp131 data  
	print ANVARF "annovar/annotate_variation.pl -filter -buildver hg19 -dbtype generic -genericdbfile hg19_cg46.txt --maf_threshold $annovar_maf_thold ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inRefDir} &', "\n\n"; #SNVs in complete genomics data
	print ANVARF "annovar/annotate_variation.pl -filter -buildver hg19 -dbtype avsift ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inSampleEnd} ${inRefDir} &', "\n\n"; #Finding benign SNVS
	print ANVARF "annovar/annotate_variation.pl -filter -buildver hg19 -dbtype ljb_pp2 ", '${inSampleDir}', '${inSamplePrefix}', "$chr[$chr].txt ", '${inSampleEnd} ${inRefDir} &', "\n\n"; #Removing benign SNVS using polyphen 2
	print ANVARF "wait", "\n\n";    
    }
        
    if ($pVMERGE eq 1) {  #varcall_merge_post_annovar_master.1.0.pl after annovar analysis
	$filename2 ="$ods/annovar_filter/varcall_merge_post_annovar_master_wf.";
	Checkfn2exists($filename2, $fnend);
	print ANVARF "sbatch $filename2","\n\n"; 
    }  
    close(ANVARF);
    return;
}
sub VCallMergeAnnovar { 
    
#Merges all variants for all subjects per chr 
    
    `mkdir -p $odf/annovar_filter/info;`; #Creates the annovar folder and info data file directory   
    `mkdir -p $ods/annovar_filter;`; #Creates the annovar script folder
    
    $filename = "$ods/annovar_filter/varcall_merge_annovar_master_wf.";
    Checkfnexists($filename, $fnend);
    
    open (VCMA, ">$filename") or die "Can't write to $filename: $!\n";
    
    print VCMA "#! /bin/bash -l", "\n";
    print VCMA "#SBATCH -A ", $aid, "\n";
    print VCMA "#SBATCH -p node -n 8", "\n";
    print VCMA "#SBATCH -C thin", "\n";	
    print VCMA "#SBATCH -t 00:05:00", "\n";
    print VCMA "#SBATCH -J VARMANNMAS", "\n";
    print VCMA "#SBATCH -e $odf/annovar_filter/info/varcall_merge_annovar_master_wf.", $fnt ,".stderr.txt", "\n";
    print VCMA "#SBATCH -o $odf/annovar_filter/info/varcall_merge_annovar_master_wf.", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print VCMA "#SBATCH --mail-type=All", "\n";
	print VCMA "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print VCMA 'echo "Running on: $(hostname)"',"\n\n";
    print VCMA "#Samples", "\n";
    print VCMA 'outSampleDir="', "$odf/annovar_filter", '"', "\n\n"; #All subjects  
    print VCMA 'outSamplePrefix="',"/varcall_comp.vcf.annovar_master", '"', "\n";   
    
    my $k=1;
   
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $k*8) { #Using only 8 cores
	    
	    print VCMA "wait", "\n\n";
	    $k=$k+1;
	}
	print VCMA "perl $ids/varcall_merge_annovar_master.1.0.pl -i "; #Calls script and infile(s) flag
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sample ids
	    
	    my $tempinfile = $avcovfn{$sid[$sampleid]};
	    
	    if ($sampleid eq scalar(@sid)-1) { #Ensure that the last print occurs with blankspace in the end
		
		print VCMA "$odf/$sid[$sampleid]/varcall_comp/$tempinfile", "_", "$chr[$chr]", "_varcall_comp.vcf.annovar.txt ";
	    }
	    else {
		print VCMA "$odf/$sid[$sampleid]/varcall_comp/$tempinfile", "_", "$chr[$chr]", "_varcall_comp.vcf.annovar.txt,";	
		
	    }
	}
	print VCMA "-s "; #Sampleid
	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {
	    
	    my $tempinfile = $avcovfn{$sid[$sampleid]};
	    
	    if ($sampleid eq scalar(@sid)-1) {
		
		print VCMA "$sid[$sampleid] ";
	    }
	    else {
		print VCMA "$sid[$sampleid],";	
		
	    }
	}
	
	print VCMA '-o ${outSampleDir}', '${outSamplePrefix}',"_", "$chr[$chr].txt &", "\n\n"; #Merges all variants from all subject per chromosome
	
    }
    print VCMA "wait", "\n\n";
    if ($pANVAR eq 1) {  #Run annovar_filter after varcall_merge_annovar_master.1.0.pl
	$filename2 ="$ods/annovar_filter/annovar_master_filter_wf.";
	Checkfn2exists($filename2, $fnend);
	print VCMA "sbatch $filename2","\n\n"; 
    }
    close(VCMA);
    return;
}

sub VarcomVCF_MB { 

#Compare and merge MB-format with VCF format to VCF format (greedy list)

    `mkdir -p $odf/$_[0]/varcall_comp/info;`; #Creates the varcall_comp folder and info data file directory   
    `mkdir -p $ods/$_[0]/varcall_comp;`; #Creates the varcall_comp script folder
   

    $filename = "$ods/$_[0]/varcall_comp/varcall_comp_vcf_MB_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (VCOMP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print VCOMP "#! /bin/bash -l", "\n";
    print VCOMP "#SBATCH -A ", $aid, "\n";
    print VCOMP "#SBATCH -p node -n 8", "\n";
    print VCOMP "#SBATCH -C thin", "\n";	
    print VCOMP "#SBATCH -t 00:10:00", "\n"; #10:00:00", "\n";
    print VCOMP "#SBATCH -J VCOMP_", $_[0], "\n";
    print VCOMP "#SBATCH -e $odf/$_[0]/varcall_comp/info/varcall_comp_vcf_MB_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print VCOMP "#SBATCH -o $odf/$_[0]/varcall_comp/info/varcall_comp_vcf_MB_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print VCOMP "#SBATCH --mail-type=All", "\n";
	print VCOMP "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print VCOMP 'echo "Running on: $(hostname)"',"\n\n";
    print VCOMP "#Samples", "\n";
    print VCOMP 'inSampleDir="',"$odf/$_[0]/samtools/var_call", '"', "\n";
    print VCOMP 'inSampleDir_MB="',"$odf/$_[0]/convert_bam_vc", '"', "\n";
    print VCOMP 'outSampleDir="', "$odf/$_[0]/varcall_comp", '"', "\n\n";   
    
    my $k=1;
    my $tempinfile = $avcovfn{$_[0]};
    
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $k*8) { #Using only 8 cores
	    
	    print VCOMP "wait", "\n\n";
	    $k=$k+1;
	}
	print VCOMP "perl $ids/varcall_compare_vcf_MB.1.1.pl ", '${inSampleDir}', "/$tempinfile","_sorted_", "$chr[$chr]", "_var.flt.vcf ", '${inSampleDir_MB}', "/$tempinfile","_sorted_", "$chr[$chr]", "_var_MB.txt -o ",'${outSampleDir}', "/$tempinfile", "_", "$chr[$chr]", "_varcall_comp.vcf &", "\n\n";
	
    }
    print VCOMP "wait", "\n\n";      
    close(VCOMP);
    return;
}

sub ConvertBamMB { 

    `mkdir -p $odf/$_[0]/convert_bam_vc/info;`; #Creates the convert_bam_vc folder and info data file directory    
    `mkdir -p $ods/$_[0]/convert_bam_vc;`; #Creates the convert_bam_vc script folder
   

    $filename = "$ods/$_[0]/convert_bam_vc/convert_bam_vc_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (CBMB, ">$filename") or die "Can't write to $filename: $!\n";
    
    print CBMB "#! /bin/bash -l", "\n";
    print CBMB "#SBATCH -A ", $aid, "\n";
    print CBMB "#SBATCH -p node -n 8", "\n";
    print CBMB "#SBATCH -C thin", "\n";	
    print CBMB "#SBATCH -t 10:00:00", "\n";
    print CBMB "#SBATCH -J CBMB_", $_[0], "\n";
    print CBMB "#SBATCH -e $odf/$_[0]/convert_bam_vc/info/convert_bam_vc_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print CBMB "#SBATCH -o $odf/$_[0]/convert_bam_vc/info/convert_bam_vc_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print CBMB "#SBATCH --mail-type=All", "\n";
	print CBMB "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print CBMB 'echo "Running on: $(hostname)"',"\n\n";
    print CBMB "module load bioinfo-tools", "\n\n"; 
    print CBMB "module load samtools/0.1.12-10", "\n\n"; #Requires newer version for bcftools
    print CBMB "#Reference Archive", "\n";
    print CBMB 'referenceArchive="', "$rd/$STMPRef", '"', "\n\n"; 
    print CBMB "#Samples", "\n";
    print CBMB 'inSampleDir="',"$odf/$_[0]/samtools", '"', "\n";
    print CBMB 'outSampleDir="', "$odf/$_[0]/convert_bam_vc", '"', "\n\n";   

    my $tempinfile = $avcovfn{$_[0]};
    my $k=1;
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $k*8) { #Using only 8 cores
	    
	    print CBMB "wait", "\n\n";
	    $k=$k+1;
	}
	print CBMB "$ids/convert_bam_to_snps_indels_sa.pl ", '${referenceArchive} ', '${inSampleDir}', "/$tempinfile","_sorted_", "$chr[$chr].bam > ", '${outSampleDir}', "/$tempinfile","_sorted_", "$chr[$chr]", "_var_MB.txt &", "\n\n";
	
    }
    $k=1; #Resetting for new infile
    print CBMB "wait", "\n\n";
      
    close(CBMB);   
    return;
}

sub SamtoolsMP { 
    
#Mpileup
#$_[0] = sampleid

    `mkdir -p $odf/$_[0]/samtools/info;`; #Creates the samtools folder and info data file directory
    `mkdir -p $odf/$_[0]/samtools/var_call;`; #Creates the  samtools var_call folder    
    `mkdir -p $ods/$_[0]/samtools;`; #Creates the samtools script folder
   

    $filename = "$ods/$_[0]/samtools/samtools_mpileup_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (STMP, ">$filename") or die "Can't write to $filename: $!\n";
    
    print STMP "#! /bin/bash -l", "\n";
    print STMP "#SBATCH -A ", $aid, "\n";
    print STMP "#SBATCH -p node -n 8", "\n";
    print STMP "#SBATCH -C thin", "\n";	
    print STMP "#SBATCH -t 10:00:00", "\n"; #10:00:00", "\n";
    print STMP "#SBATCH -J STMP_", $_[0], "\n";
    print STMP "#SBATCH -e $odf/$_[0]/samtools/info/samtools_mpileup_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print STMP "#SBATCH -o $odf/$_[0]/samtools/info/samtools_mpileup_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print STMP "#SBATCH --mail-type=All", "\n";
	print STMP "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print STMP 'echo "Running on: $(hostname)"',"\n\n";
    print STMP "module load bioinfo-tools", "\n\n"; 
    print STMP "module load samtools/0.1.12-10", "\n\n"; #Requires newer version for bcftools
    print STMP "#Reference Archive", "\n";
    print STMP 'referenceArchive="', "$rd/$STMPRef", '"', "\n\n"; 
    print STMP "#Samples", "\n";
    print STMP 'inSampleDir="',"$odf/$_[0]/samtools", '"', "\n";
    print STMP 'outSampleDir="', "$odf/$_[0]/samtools/var_call", '"', "\n\n";   
    
    my $tempinfile = $avcovfn{$_[0]};
    my $k=1;
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $k*8) { #Using only 8 cores
	    
	    print STMP "wait", "\n\n";
	    $k=$k+1;
	}
	print STMP "samtools mpileup -ugf ", '${referenceArchive} ', '${inSampleDir}', "/$tempinfile","_sorted_", "$chr[$chr].bam | bcftools view -bvcg - > ", '${outSampleDir}', "/$tempinfile","_sorted_", "$chr[$chr]", "_var.raw.bcf &", "\n\n";
	
    }
    
    $k=1; #Resetting for new infile
    print STMP "wait", "\n\n";
    
    print STMP "#Samples", "\n";
    print STMP 'inSampleDir="',"$odf/$_[0]/samtools/var_call", '"', "\n";
    print STMP 'outSampleDir="', "$odf/$_[0]/samtools/var_call", '"', "\n\n";
    
    $tempinfile = $avcovfn{$_[0]};
    @avcov = `cut -f5 $odf/$_[0]/mosaik/mosaikMerge/coverageReport/$avcovfn{$_[0]}.bam.coverage_per_chromosome.txt`; #Collects average coverage of each chromosome for bcftools -D, should be 2x avg. cov.
    shift @avcov; #Removes header from array 

    foreach my $value (@avcov) { $value = $value * 2; } #Multiply all entries by 2

    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr	    
	
	if ($chr eq $k*8) { #Using only 8 cores
	    
	    print STMP "wait", "\n\n";
	    $k=$k+1;
	}
	chomp($avcov[$chr]);
	print STMP "bcftools view ", '${inSampleDir}', "/$tempinfile","_sorted_", "$chr[$chr]", "_var.raw.bcf | vcfutils.pl varFilter -D$avcov[$chr] > ", '${outSampleDir}', "/$tempinfile","_sorted_", "$chr[$chr]", "_var.flt.vcf &", "\n\n";
	
    }
    $k=1; #Resetting for new infile
    print STMP "wait", "\n\n";
    close(STMP);   
    return;
}

sub SamtoolsView { 
    
#View

    `mkdir -p $odf/$_[0]/samtools/info;`; #Creates the samtools folder and info data file directory
    `mkdir -p $ods/$_[0]/samtools;`; #Creates the samtools folder and info data file directory

    $filename = "$ods/$_[0]/samtools/samtools_view_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (STV, ">$filename") or die "Can't write to $filename: $!\n";
    
    print STV "#! /bin/bash -l", "\n";
    print STV "#SBATCH -A ", $aid, "\n";
    print STV "#SBATCH -p node -n 8", "\n";
    print STV "#SBATCH -C thin", "\n";	
    print STV "#SBATCH -t 10:00:00", "\n"; 
    print STV "#SBATCH -J STV_", $_[0], "\n";
    print STV "#SBATCH -e $odf/$_[0]/samtools/info/samtools_view_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print STV "#SBATCH -o $odf/$_[0]/samtools/info/samtools_view_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print STV "#SBATCH --mail-type=All", "\n";
	print STV "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print STV 'echo "Running on: $(hostname)"',"\n\n";
    print STV "module load bioinfo-tools", "\n\n"; 
    print STV "module load samtools/0.1.12-10", "\n\n";
    print STV "#Samples", "\n";
    
    my $k=1;
    my $tempinfile = $avcovfn{$_[0]};
    
    for (my $chr=0;$chr<scalar(@chr);$chr++) { #For all chr
	
	if ($chr eq $k*8) { #Using only 8 cores
	    
	    print STV "wait", "\n\n";
	    $k=$k+1;
	}
	print STV "samtools view -b -o $odf/$_[0]/samtools/$tempinfile", "_sorted_", "$chr[$chr].bam $odf/$_[0]/samtools/$tempinfile", "_sorted.bam $chr[$chr] &", "\n\n";
    }
    print STV "wait", "\n\n";
    close(STV);   
    return;
}

sub SamtoolsSortIndex { 
#Sort, index
#$_[0]= $sampleid
#$_[1]= $sampleidcounter

    `mkdir -p $odf/$_[0]/samtools/info;`; #Creates the samtools folder and info data file directory
    `mkdir -p $ods/$_[0]/samtools;`; #Creates the samtools folder and info data file directory

    $filename = "$ods/$_[0]/samtools/samtools_sort_index_wf_$_[0].";
    Checkfnexists($filename, $fnend);

    open (STSI, ">$filename") or die "Can't write to $filename: $!\n";
    
    print STSI "#! /bin/bash -l", "\n";
    print STSI "#SBATCH -A ", $aid, "\n";
    print STSI "#SBATCH -p node -n 8", "\n";
    print STSI "#SBATCH -C thin", "\n";	
    print STSI "#SBATCH -t 20:00:00", "\n";
    print STSI "#SBATCH -J STSI_", $_[0], "\n";
    print STSI "#SBATCH -e $odf/$_[0]/samtools/info/samtools_sort_index_wf_$_[0].", $fnt ,".stderr.txt", "\n";
    print STSI "#SBATCH -o $odf/$_[0]/samtools/info/samtools_sort_index_wf_$_[0].", $fnt ,".stdout.txt", "\n";
    
    unless ($em eq 0) {
	
	print STSI "#SBATCH --mail-type=All", "\n";
	print STSI "#SBATCH --mail-user=$em", "\n\n";
	
    }
    
    print STSI 'echo "Running on: $(hostname)"',"\n\n";
    print STSI "module load bioinfo-tools", "\n\n"; 
    print STSI "module load samtools/0.1.12-10", "\n\n";
    print STSI "#Samples", "\n";
    
    my $tempinfile = $infiles{$_[0]}; #Whole path
    my $tempoutfile = $avcovfn{$_[0]}; #Just filename
    print STSI "samtools sort $tempinfile.bam $odf/$_[0]/samtools/$tempoutfile","_sorted", "\n\n";
    print STSI "samtools index $odf/$_[0]/samtools/$tempoutfile","_sorted.bam", "\n\n";   
    print STSI "wait", "\n\n";

    if ($_[1] eq scalar(@sid)-1 ) {#For last sampleID

	for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sampleIDs
   
	    if ($pSTV eq 1) {  #Samtools View after samtools sort & index
	
		$filename2 ="$ods/$sid[$sampleid]/samtools/samtools_view_wf_$sid[$sampleid].";
		Checkfn2exists($filename2, $fnend);
		if ($sampleid eq scalar(@sid)-1 ) { #Go ahead and submit since samToolsSort & Index sbatch script is finished.
		    print STSI "STVJOB=`sbatch $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id when submitted. Enables SamToolsView sbatch script to be initiated when samToolsSort & Index sbatch script is finished.
		}
		else { #Wait for the initial samToolsSort & Index to be finished for the other sampleids
		    print STSI "STVJOB=`sbatch --dependency=afterok:$jobID[$sampleid]"," $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id when submitted. Enables SamToolsView sbatch script to be initiated when samToolsSort & Index sbatch script is finished.
		}
	    }
	    if ($pSTMP eq 1) {  #Samtools mpileup after samtools sort, index and view
		
		$filename2 ="$ods/$sid[$sampleid]/samtools/samtools_mpileup_wf_$sid[$sampleid].";
		Checkfn2exists($filename2, $fnend);
		print STSI "STMPJOB=`sbatch --dependency=afterok:",'$STVJOB'," $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id of Samtools mpileup
	    }
	    if ($pCBSI eq 1) {  #convert_bam_to_snps_indels_sa.pl after samtools sort, index and view
		$filename2 ="$ods/$sid[$sampleid]/convert_bam_vc/convert_bam_vc_wf_$sid[$sampleid].";
		Checkfn2exists($filename2, $fnend);
		print STSI "CBSIJOB=`sbatch --dependency=afterok:",'$STVJOB', " $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id of convert_bam_to_snps_indels_sa.pl
	    }
	    if ($pVCOMP eq 1 ) {  #varcall_compare_vcf_MB.pl after samtools sort, index, view, mpileup and convert_bam_to_snps_indels_sa.pl
		
		$filename2 ="$ods/$sid[$sampleid]/varcall_comp/varcall_comp_vcf_MB_wf_$sid[$sampleid].";
		Checkfn2exists($filename2, $fnend);
		print STSI "VCOMPJOB_$sampleid=`sbatch --dependency=afterok:",'$STMPJOB', ":", '$CBSIJOB', " $filename2 | egrep -o -e", '"\b[0-9]+$"`', "\n\n"; #Catch sbatch job_id varcall_comp_vcf_MB per sampleid
	    }
	}
	if ($pVCMA eq 1) {  #varcall_merge_annovar_master.1.0.pl after varcall_compare_vcf_MB.1.1.pl for all subjects
	    $filename2 ="$ods/annovar_filter/varcall_merge_annovar_master_wf.";
	    Checkfn2exists($filename2, $fnend);
	    print STSI "sbatch --dependency=afterok:";
	    
	    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) { #For all sampleIDs
		
		if ($sampleid eq scalar(@sid)-1 ) {
		    print STSI '$VCOMPJOB_',"$sampleid ";
		}
		else{
		    print STSI '$VCOMPJOB_',"$sampleid:";
		} 
	    }
	    print STSI $filename2,"\n\n";
	}
    }
    close(STSI);
    my $ret = `sbatch $filename`;
    my ($jobID) = ($ret =~ /Submitted batch job (\d+)/);
    push(@jobID, $jobID);
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


sub FileReFormat {
    
#Code needed to reformat files name removing .bam and keeping stub. 

    for (my $infile=0;$infile<scalar(@infn);$infile++) {  

	if ( basename($infn[$infile]) =~ /(\S+)\.bam$/ ) { #Parse to .bam format just filename
	   $avcovfn{$sid[$infile]} = $1; #Stores infile per sample id
	    
	}
	if ( $infn[$infile] =~ /(\S+)\.bam$/ ) { #Parse to .bam format whole path
	    
	    $infiles{$sid[$infile]} = $1;
	    $dirname{$sid[$infile]} = dirname($infn[$infile]); #Stores infile dir   
	}	
    }
    return;
}

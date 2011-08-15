#Adds depth-info for non-var sites of every sample and writes new varcall "master" file with added non-var depth info.
#
#NB: prefix of bamfile needs to exactly match samplename and have _lanes directly appended to the samplename prefix.
#
#USAGE: perl get.siteinfo.pl -var /proj/b2011075/melanoma_exomseq/mosaik_pipe/outdata/run2/annovar_filter/annovar_master_all_subject_variants.txt -bam /proj/b2011075/melanoma_exomseq/mosaik_pipe/outdata/217_1/samtools/217_1_lanes_45_S_merged_sorted.bam -ods /proj/b2011075/melanoma_exomseq/mosaik_pipe/outscripts/run2/annovar_filter -a b2010035 -em daniel.edsgard@scilifelab.se

use strict;
use File::Basename;
use File::Spec::Functions;
use Getopt::Long;

my ($varfile, @sample_bamfiles, $aid, $em, $outscriptdir);
GetOptions('var|varfile:s', => \$varfile,
	   'bam|infile:s'  => \@sample_bamfiles, #Comma separeted list
	   'a|projectid:s'  => \$aid,
	   'em|email:s'  => \$em,
	   'ods|outdirscript:s'  => \$outscriptdir,
	  );

#  my $varfile = shift @_;
#  my $sample_bamfiles_arrayref = shift @_;
#  my $outscriptdir = shift @_;
#  my $aid = shift @_;
#  my $em = shift @_;

#  my @sample_bamfiles = @{$sample_bamfiles_arrayref};

my $outdir = dirname($varfile);

my $sbatch_file = catfile($outscriptdir, 'addnonvarinfo.sbatch');
open(ADDNONVAR, '>', $sbatch_file) or die("Couldn't open $sbatch_file: $!\n");

print ADDNONVAR "#! /bin/bash -l", "\n";
print ADDNONVAR "#SBATCH -A ", $aid, "\n";
print ADDNONVAR "#SBATCH -p core -n 1", "\n";
print ADDNONVAR "#SBATCH -t 2:00:00", "\n"; #10:00:00", "\n";
print ADDNONVAR "#SBATCH -J ADDNONVAR", "\n";
print ADDNONVAR "#SBATCH -e $outdir/addnonvar.stderr.txt", "\n";
print ADDNONVAR "#SBATCH -o $outdir/addnonvar.stdout.txt", "\n";
    
unless ($em eq 0) {	
  print ADDNONVAR "#SBATCH --mail-type=All", "\n";
  print ADDNONVAR "#SBATCH --mail-user=$em", "\n\n";	
}
    
print ADDNONVAR 'echo "Running on: $(hostname)"',"\n\n";
print ADDNONVAR 'module load bioinfo-tools', "\n"; 
print ADDNONVAR 'module load samtools/0.1.12-10', "\n";
print ADDNONVAR 'module load R/2.13.0' . "\n\n";
  

#For every sample print all the sites where it didnt have any variation but any of the other samples had.
my $cmd = 'Rscript get.sites.R ' . $varfile . ' ' . $outdir;
print ADDNONVAR $cmd . "\n";


#Call mpileup for every sample: mpileup -l jsample.sitefile.
#samtools mpileup -l $sitesfile $infile >file.tmp.
#my @nonvarsites_files = <$outdir/*.nonvarsites.txt>;  
foreach my $jsample_bamfile (@sample_bamfiles) {

  #get sample id
  #prefix of bamfile needs to exactly match samplename and have _lanes directly appended to the samplename prefix.
  $jsample_bamfile =~ m/.*\/([^\/]*)_lanes[^\/]*\.bam$/;
  my $jsample = $1;
    
  my $jsample_sitesfile = catfile($outdir, $jsample . '.nonvarsites.txt');
  my $jsample_nonvars_depth_file = catfile($outdir, $jsample . '.nonvarsites.info');
  my $cmd = 'samtools mpileup -l ' . $jsample_sitesfile . ' ' . $jsample_bamfile . ' >' . $jsample_nonvars_depth_file;
  print ADDNONVAR $cmd . "\n";
}


#Add depth to varcall file for samples with no called variant.
my $varfile_addeddepth = $varfile . '.addeddepth';
my $cmd = 'Rscript depthadd2varcallfile.R ' . $varfile . ' ' . $varfile_addeddepth;
print ADDNONVAR $cmd . "\n";

close(ADDNONVAR);

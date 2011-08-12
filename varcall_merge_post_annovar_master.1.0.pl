#!/usr/bin/perl - w

use strict;
use warnings;

#Compiles variation calls for annovar analysis files
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS

varcall_merge_post_annovar_master.1.0.pl -i infile1..n -o outfile.txt

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s)

-nos/--numberofsubjects Tracks the nuber of subjects to make sure that all info stays in proper column

-o/--outfile The output file (defaults to annovar_master.txt)

=head3 I/O

Input format (VCF/Custom )

Output format (tab separate list)

=cut

use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
	qq{varcall_merge_post_annovar_master.1.0.pl -i [queryfile.txt...n] -nos No_of_subjects -o outfile.txt
	       -i/--infile Infile(s), comma sep
               -nos/--numberofsubjects Tracks the nuber of subjects to make sure that all info stays in proper column
	       -o/--outfile The output file (defaults to annovar_master.txt)
	   };
    
}

my ($of, $nos, $arraypos, $help) = ("annovar_master.txt", 0, 0);
my (@infn,@sid, @chr);
@chr = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"); #Chr for ordering master_output file

GetOptions('i|infile:s'  => \@infn, #Comma separeted list
	   'nos|noofsubjects:s'  => \$nos,
	   'o|outfile:s'  => \$of,
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
if ($nos eq 0) {
 print STDERR "\n", "Should specify how many subjects that are included or you migt overwrite prescence of subjects.\n";
 print STDERR "\n", "Setting number of subjects to 2, change by using flag -nos.\n";
 $nos = 2;
}

@infn = split(/,/,join(',',@infn)); #Enables comma separated indir(s)

$arraypos = 7+$nos; #Chr,start,stop,Ref,Variant, Location and Gene/Na is always the same. Makes sure that ana anlysis always ends up in the same coulmn independent if a variant is present or not in db etc. Changes after each read of new file.

my (%allVariants, %allVariants_chr, %allVariants_chr_unique, %allVariants_chr_sorted);
my (@allVariants, @allVariants_unique, @allVariants_sorted);

#Read annotations
for (my $inputfiles=0;$inputfiles<scalar(@infn);$inputfiles++) {
    
     ReadVarCAnnoVMaster($infn[$inputfiles]);
     ReadAnnoVarFunc($infn[$inputfiles]);
     ReadAnnoExoVarFunc($infn[$inputfiles], $arraypos);
     ReadAnnoConservedR($infn[$inputfiles], $arraypos);
     ReadAnnoSegDup($infn[$inputfiles], $arraypos);
     ReadAnno1000g($infn[$inputfiles], $arraypos);
     ReadAnnoDbsnp($infn[$inputfiles], $arraypos);
     ReadAnnocg46($infn[$inputfiles], $arraypos);
     ReadAnnoAvsift($infn[$inputfiles], $arraypos);
     ReadAnnoPP2($infn[$inputfiles], $arraypos);
     unless ( $inputfiles eq scalar(@infn)-1 ) {
	 $arraypos = 7+$nos; #Reset for new infile
     }
}

SortAllVariants();
WriteMasterAnnoV($of);

#add site depth info among samples where the site was not called as a variant
#AddNonVarSiteInfo($of, \@sample_indexedbamfiles, $ods, $aid, $em);

sub AddNonVarSiteInfo{
#NB: Prefix of bamfile needs to exactly match samplename (as printed in the master annotation txt file) and the samplename prefix directly appended by '_lanes'.

#SEE: get.siteinfo.pl:
#  perl get.siteinfo.pl -var /proj/b2011075/melanoma_exomseq/mosaik_pipe/outdata/run2/annovar_filter/annovar_master_all_subject_variants.txt -bam /proj/b2011075/melanoma_exomseq/mosaik_pipe/outdata/217_1/samtools/217_1_lanes_45_S_merged_sorted.bam -ods /proj/b2011075/melanoma_exomseq/mosaik_pipe/outscripts/run2/annovar_filter -a b2010035 -em daniel.edsgard@scilifelab.se

  return;
}

sub SortAllVariants {
#Creates an array of all position which are unique and in sorted ascending order
    
    for my $chr (keys %allVariants)  { #For all chr
	
	for my $pos (keys %{ $allVariants{$chr} } )  { #For all pos
	    
	    for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants
		#print $chr,"\t", $pos,"\t";
		#print $variant,"\t";
		push ( @{$allVariants_chr{$chr} },$pos );
		#push(@allVariants, $pos); #All non-overlapping positions per chr but two variants at the same chr,pos yields two entries in @allVariants.
		#for (my $i=0;$i<scalar( @{$allVariants{$chr}{$pos}{$variant} });$i++) {
		#  print $allVariants{$chr}{$pos}{$variant}[$i], "\t";
		#}
		#print "\n";
	    }
	}
	my %seen = (); @{$allVariants_chr_unique{$chr} } = grep { ! $seen{$_} ++ } @{$allVariants_chr{$chr} }; #Unique entries only 
	@{$allVariants_chr_sorted{$chr} } = sort { $a <=> $b } @{ $allVariants_chr_unique{$chr} }; #Sorts keys to be able to print sorted table later 
    }  
    print STDERR "Sorted all non overlapping entries\n";
	
}

sub ReadAnnoPP2 {
#Reads annovar.hg19_ljb_pp2_dropped, higher scores means more deleterious predictions
#$_[0] = filename
#$_[1] = $arraypos
    
    open(ACR, "<$_[0].hg19_ljb_pp2_dropped") or die "Can't open $_[0].hg19_ljb_pp2_dropped:$!, \n";    
    
    while (<ACR>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    if ($allVariants{$temp[2]}{$temp[3]}{$temp[6]}) { #If present in ljb_pp2 genomes db add info
		$allVariants{$temp[2]}{$temp[3]}{$temp[6]}[$_[1]] = $temp[0] .=";$temp[1]";
		#push(@{ $allVariants{$temp[2]}{$temp[3]}{$temp[6]} }, $temp[0] .=";$temp[1]"); # Hash{chr}{pos}{variant}, all variants, adds AVSIFT (SIFT_score>0.05 means benign)
	    }
	}
    } 	
    $arraypos++;
    close(ACR);
    print STDERR "Finished Reading Infile $_[0].hg19_ljb_pp2_dropped","\n";
    return;
}

sub ReadAnnoAvsift {
#Reads annovar.hg19_avsift_dropped
#$_[0] = filename
#$_[1] = $arraypos
    
    open(ACR, "<$_[0].hg19_avsift_dropped") or die "Can't open $_[0].hg19_avsift_dropped:$!, \n";    
    
    while (<ACR>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    if ($allVariants{$temp[2]}{$temp[3]}{$temp[6]}) { #If present in avsift db add info
		$allVariants{$temp[2]}{$temp[3]}{$temp[6]}[$_[1]] = $temp[0] .=";$temp[1]";
		#push(@{ $allVariants{$temp[2]}{$temp[3]}{$temp[6]} }, $temp[0] .=";$temp[1]"); # Hash{chr}{pos}{variant}, all variants, adds AVSIFT (SIFT_score>0.05 means benign)
	    }
	    
	}
    } 	
    $arraypos++;
    close(ACR);
    print STDERR "Finished Reading Infile $_[0].hg19_avsift_dropped","\n";
    return;
}

sub ReadAnnocg46 {
#Reads annovar.hg19_generic_dropped (cg46)

  my $db = 'cg46';

  my $inputfileprefix = shift @_;
  my $local_arraypos = shift @_;
  my $inputfile = $inputfileprefix . '.hg19_generic_dropped';
    
    open(ACR, '<', $inputfile) or die "Can't open $inputfile:$!, \n";    
    
    while (defined(my $line = <ACR>)) {
        chomp $line;
	
	if ($line =~ m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if ($line =~ m/(\S+)/) {
	    
	    my @cols = split("\t", $line);	    #Loads variant calls
	    my ($maf, $chr, $start, $obsallele) = ($cols[1], $cols[2], $cols[3], $cols[6]);

	    if ($allVariants{$chr}{$start}{$obsallele}) { #If present in complete genomics 46 genomes db add info
		$allVariants{$chr}{$start}{$obsallele}[$local_arraypos] = $db . ';' . $maf;
		#push(@{ $allVariants{$cols[2]}{$start}{$obsallele} }, "cg46"); # Hash{chr}{pos}{variant}, all variants, adds COMPLETE GENOMICS 46 GENOMES presence
	    }
	}
    } 	
    $arraypos++;
    close(ACR);
    print STDERR "Finished Reading Infile $inputfile","\n";
    return;
}

sub ReadAnnoDbsnp {
#Reads annovar.hg19_snp131_dropped
    
  my $inputfileprefix = shift @_;
  my $local_arraypos = shift @_;
  my $inputfile = $inputfileprefix . '.hg19_snp131_dropped';

  open(ACR, '<', $inputfile) or die "Can't open $inputfile:$!, \n";    
    
    while (defined(my $line = <ACR>)) {
        chomp $line;
	
	if ($line =~ m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if ($line =~ m/(\S+)/) {
	    
	    my @cols = split("\t", $line);	    #Loads variant calls
	    my ($dbver, $rs_id, $chr, $start, $obsallele) = ($cols[0], $cols[1], $cols[2], $cols[3], $cols[6]);

	    if ($allVariants{$chr}{$start}{$obsallele}) { #If present in dbsnp add info
		$allVariants{$chr}{$start}{$obsallele}[$local_arraypos] = $dbver;
		$allVariants{$chr}{$start}{$obsallele}[$local_arraypos + 1] = $rs_id;
		#push(@{ $allVariants{$temp[2]}{$temp[3]}{$temp[6]} }, $temp[0], $temp[1] ); # Hash{chr}{pos}{variant}, all variants, adds DBSNP presence and ID
	    }
	}
    } 	
    $arraypos = $arraypos + 2;
    close(ACR);
    print STDERR "Finished Reading Infile $inputfile","\n";
    return;
}

sub ReadAnno1000g {
#Reads annovar.hg19_ALL.sites.2010_11_dropped
    
  my $inputfileprefix = shift @_;
  my $local_arraypos = shift @_;
  my $inputfile = $inputfileprefix . '.hg19_ALL.sites.2010_11_dropped';

    open(ACR, '<', $inputfile) or die "Can't open $inputfile:$!, \n";    
    
    while (defined(my $line = <ACR>)) {
        chomp $line;
	
	if ($line =~ m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if ($line =~ m/(\S+)/) {
	    
	    my @cols = split("\t", $line);	    #Loads variant calls
	    my ($dbver, $maf, $chr, $start, $obsallele) = ($cols[0], $cols[1], $cols[2], $cols[3], $cols[6]);

	    if ($allVariants{$chr}{$start}{$obsallele}) { #If present in 1000g db add info
		$allVariants{$chr}{$start}{$obsallele}[$local_arraypos] = $dbver . ';' . $maf;
		#push(@{ $allVariants{$temp[2]}{$temp[3]}{$temp[6]} }, $temp[0]); # Hash{chr}{pos}{variant}, all variants, adds 1000G presence
	    }
	}
    } 	
    $arraypos++;
    close(ACR);
    print STDERR "Finished Reading Infile $inputfile","\n";
    return;
}

sub ReadAnnoSegDup {
#Reads annovar.hg19_genomicSuperDups
#$_[0] = filename
#$_[1] = $arraypos
    
    open(ACR, "<$_[0].hg19_genomicSuperDups") or die "Can't open $_[0].hg19_genomicSuperDups:$!, \n";    
    
    while (<ACR>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    if ($allVariants{$temp[2]}{$temp[3]}{$temp[6]}) { #If present in segmental duplications add info
		$allVariants{$temp[2]}{$temp[3]}{$temp[6]}[$_[1]] = $temp[0].=";$temp[1]";
#		push(@{ $allVariants{$temp[2]}{$temp[3]}{$temp[6]} }, $temp[0].=";$temp[1]" ); # Hash{chr}{pos}{variant}, all variants, adds SEGDUP;SCORE=\d+;NAME=chr\d+
	    }
	}
    } 	
    $arraypos++;
    close(ACR);
    print STDERR "Finished Reading Infile $_[0].hg19_genomicSuperDups","\n";
    return;
}

sub ReadAnnoConservedR {
#Reads annovar..hg19_phastConsElements46way file
#$_[0] = filename
#$_[1] = $arraypos
    open(ACR, "<$_[0].hg19_phastConsElements46way") or die "Can't open $_[0].hg19_phastConsElements46way:$!, \n";    
    
    while (<ACR>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    if ($allVariants{$temp[2]}{$temp[3]}{$temp[6]}) { #If present in conserved region add info
		$allVariants{$temp[2]}{$temp[3]}{$temp[6]}[$_[1]] = $temp[0].=";$temp[1]";
	#	push(@{ $allVariants{$temp[2]}{$temp[3]}{$temp[6]} }, $temp[0].=";$temp[1]" ); # Hash{chr}{pos}{variant}, all variants, adds CONSERVED(db);SCORE=\d+;NAME=lod=\d+
	    }
	}
    } 	
    $arraypos++;
    close(ACR);
    print STDERR "Finished Reading Infile $_[0].hg19_phastConsElements46way","\n";
    return;
}

sub ReadAnnoExoVarFunc {
#Reads annovar.exonic_variant_function file
#$_[0] = filename
#$_[1] = $arraypos
    
    open(AVF, "<$_[0].exonic_variant_function") or die "Can't open $_[0].exonic_variant_function:$!, \n";    
    
    while (<AVF>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    if ( $allVariants{$temp[3]}{$temp[4]}{$temp[7]}) {
		$allVariants{$temp[3]}{$temp[4]}{$temp[7]}[$_[1]] = $temp[1];
		$allVariants{$temp[3]}{$temp[4]}{$temp[7]}[$_[1]+1] = $temp[2];
		#push(@{ $allVariants{$temp[3]}{$temp[4]}{$temp[7]} }, $temp[1], $temp[2] ); # Hash{chr}{pos}{variant}, all variants, adds SYNONYMOUS and AMINO ACID CHANGE
	    }
	} 		
    }
    $arraypos = $arraypos+2;
    close(AVF);
    print STDERR "Finished Reading Infile $_[0].exonic_variant_function","\n";
    return;
}

sub ReadAnnoVarFunc {
#Reads annovar.variant_function file
#$_[0] = filename
    
    open(AVF, "<$_[0].variant_function") or die "Can't open $_[0].variant_function:$!, \n";    
    
    while (<AVF>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    push(@{ $allVariants{$temp[2]}{$temp[3]}{$temp[6]} }, $temp[0], $temp[1] ); # Hash{chr}{pos}{variant}, all variants, adds LOCATION and GENE
	    
	}
    } 	
    close(AVF);
    print STDERR "Finished Reading Infile $_[0].variant_function","\n";
    return;
}

sub ReadVarCAnnoVMaster {
#Reads varcall_comp.vcf.annovar_master_$chr.txt file
#$_[0] = filename
    
    open(VCAVM, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<VCAVM>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    $allVariants{$temp[0]}{$temp[1]}{$temp[4]} = [@temp]; # Hash{chr}{pos}{variant}, all variants non overlapping and array [chr->unknown] All info starting from chr	    
	}
    } 	
    close(VCAVM);
    print STDERR "Finished Reading Infile $_[0]","\n";
    return;
}

sub WriteMasterAnnoV {
    
#Prints tab separated file
#Serves as master file output from annovar analysis merged over all subjects and chr
#All variants non-overlapping
    
    my $filename = shift;
    open (ANVAR, ">$filename") or die "Can't write to $filename: $!\n"; 
    for (my $chrpos=0; $chrpos<scalar( @chr ); $chrpos++) { #For all chr positions
	
	my $chr = $chr[$chrpos];
	if ($allVariants{$chr}) {
	    
	    for (my $i=0;$i<scalar( @{$allVariants_chr_sorted{$chr} } );$i++)  { #For all pos per chr	
		
		my $pos = $allVariants_chr_sorted{$chr}[$i]; #pos keys to hash from sorted arrray
		
		for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants
		    
		    print ANVAR "chr$chr","\t", $pos,"\t", $allVariants{$chr}{$pos}{$variant}[2], "\t", $allVariants{$chr}{$pos}{$variant}[3], "\t", $variant,"\t";   #$allVariants{$chr}{$pos}{$variant}[2] = stop position, #$allVariants{$chr}{$pos}{$variant}[3] = reference nucleotide
		    #for (my $k=5;$k<scalar( @{$allVariants{$chr}{$pos}{$variant} } );$k++)  { #For all entries
		    for (my $k=5; $k < scalar( $arraypos ); $k++)  { #For all entries indepent if true or not
			
			if ($allVariants{$chr}{$pos}{$variant}[$k]){
			    print ANVAR $allVariants{$chr}{$pos}{$variant}[$k], "\t";
			}
			else {  print ANVAR "-","\t";
			}
		    }
		    print ANVAR "\n";
		}
	    }
	}
    }
    close (ANVAR);
    print STDERR "Finished Writing Master file","\n";
    return;
}

#printing all in sorted order
#for my $chr (keys %allVariants)  { #For all chr
#    for (my $i=0;$i<scalar( @{$allVariants_chr_sorted{$chr} } );$i++)  { #For all pos
#	for my $variant (keys % { $allVariants{$chr}{$allVariants_chr_sorted{$chr}[$i]} })  { #For all variants
#	    print $chr,"\t", $allVariants_chr_sorted{$chr}[$i],"\t";
#	    print $variant,"\t";
#	    for (my $k=0;$k<scalar( @{$allVariants{$chr}{$allVariants_chr_sorted{$chr}[$i]}{$variant} });$k++) {
#		print $allVariants{$chr}{$allVariants_chr_sorted{$chr}[$i]}{$variant}[$k], "\t";
#	    }
#	    print "\n";
#	}
 #   }
#} 
#    for my $chr (keys %allVariants)  { #For all chr
#	for (my $i=0;$i<scalar( @allVariants_sorted );$i++)  { #For all pos
#	    for my $variant (keys % { $allVariants{$chr}{$allVariants_sorted[$i]} })  { #For all variants
#		print $chr,"\t", $allVariants_sorted[$i],"\t";
#		print $variant,"\t";
#		for (my $k=0;$k<scalar( @{$allVariants{$chr}{$allVariants_sorted[$i]}{$variant} });$k++) {
#		    print $allVariants{$chr}{$allVariants_sorted[$i]}{$variant}[$k], "\t";
#		}
#		print "\n";
#	    }
#	}
 #  } 

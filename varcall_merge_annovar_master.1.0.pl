#!/usr/bin/perl - w

use strict;
use warnings;

#Compares and merges variation calls for annovar.txt files per chr
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS

varcall_merge_annovar_master.1.0.pl -i infile1..n -s sampleid1..n -o outfile.txt

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s)

-s/--sampleid The sample ID(s)

-o/--outfile The output file (defaults to varcall_merge_annovar_master.chr.txt)

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
	qq{varcall_merge_annovar_master.1.0.pl -i [queryfile.txt...n] -s sampleid..n -o outfile.txt
	       -i/--infile Infile(s), comma sep
	       -s/--sampleid The sample ID(s), comma sep
	       -o/--outfile The output file (defaults to varcall_merge_annovar_master.chr.txt)
	   };
    
}

my ($of, $help) = ("varcall_merge_annovar_master");
my (@infn,@sid);

GetOptions('i|infile:s'  => \@infn, #Comma separeted list
	   's|sampleid:s'  => \@sid, #Comma separeted list
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
if ( scalar(@sid) eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a sample ID as a comma separated list", "\n\n";
    die $USAGE;
}

@infn = split(/,/,join(',',@infn)); #Enables comma separated indir(s)
@sid = split(/,/,join(',',@sid)); #Enables comma separeted list of sample IDs

my (%allVariants, %subjectVariants);
my (@allVariants, @allVariants_unique, @allVariants_sorted);

    
for (my $inputfiles=0;$inputfiles<scalar(@infn);$inputfiles++) {
    
    ReadVarCVCF($infn[$inputfiles], $inputfiles);    #Reads all variants for all infiles
}

SortAllVariants(); #Creates an sorted array of positions with unique entries
WriteMasterAnnoV($of);   #Writes all variants per chr to disc

sub SortAllVariants {
#Creates an array of all position which are unique and in sorted ascending order
    
    for my $chr (keys %allVariants)  { #For all chr

	for my $pos (keys %{ $allVariants{$chr} } )  { #For all pos
	    	    
	    for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants

		push(@allVariants, $pos); #All non-overlapping positions per chr but two variants at the same chr,pos yields two entries in @allVariants.
		print "\n";
	    }
	}
    }    
    my %seen = (); @allVariants_unique = grep { ! $seen{$_} ++ } @allVariants; #Unique entries only
    @allVariants_sorted = sort { $a <=> $b } @allVariants_unique; #Sorts keys to be able to print sorted table later
    print STDERR "Sorted all non overlapping entries\n";
   
}

sub ReadVarCVCF {
#Reads varcall_comp.vcf.annovar.txt file
#$_[0] = filename
#$_[1] = Nr of inputfile
    
    open(AVF, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<AVF>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    $subjectVariants{$_[1]}{$temp[0]}{$temp[1]}{$temp[4]} = [@temp]; # Hash{inputfile}{chr}{pos}{variant} of array[chr->unknown] All info starting from chr
	    $allVariants{$temp[0]}{$temp[1]}{$temp[4]} = [@temp]; # Hash{chr}{pos}{variant}, all variants non overlapping and [DP=X...FQ=Y/PV4=Z] 
	    
	}
    } 	
    close(AVF);
    print STDERR "Finished Reading Infile $_[0]","\n";
    return;
}

sub WriteMasterAnnoV {
    
#Prints tab separated file with chr, start, stop, ref, subjects_found DP..PV4
#Serves as master file input to annovar analysis merged over all subjects per chr
#All variants non-overlapping
    
    my $filename = shift;
    my $postracker=0;

    for my $chr (keys %allVariants)  { #For all chr
	open (ANVAR, ">$filename") or die "Can't write to $filename: $!\n"; 
	for (my $i=0;$i<scalar( @allVariants_sorted );$i++)  { #For all pos
	    
	    my $pos = $allVariants_sorted[$i]; #pos keys to hash from sorted arrray
	  
	    for my $variant (keys % { $allVariants{$chr}{$pos} })  { #For all variants
		
		print ANVAR $chr,"\t", $pos,"\t", $allVariants{$chr}{$pos}{$variant}[2], "\t", $allVariants{$chr}{$pos}{$variant}[3], "\t", $variant,"\t";		#$allVariants{$chr}{$pos}{$variant}[2] = stop position, $allVariants{$chr}{$pos}{$variant}[3] = reference nucleotide

		for my $infiles (keys %subjectVariants) { #Loop over all infiles
		    
		    if ($subjectVariants{$infiles}{$chr}{$pos}{$variant}) { #If present in infile
			print ANVAR $sid[$infiles];
			if ( $subjectVariants{$infiles}{$chr}{$pos}{$variant}[7] =~m/AF1=(\d+.\d+|\d+)/ ) { #Allel freq
			    if ($1 > 0.85) { #Homozygous
				print ANVAR "Hom;AF1=$1;";
			    }
			    else {
				print ANVAR "Het;AF1=$1;";
			    }
			}
			if ( $subjectVariants{$infiles}{$chr}{$pos}{$variant}[7] =~m/(DP=\d+)/ ) { #Read Depth
			    print ANVAR "$1";
			}
			print ANVAR "\t";
		    }
		    else { print ANVAR "-","\t"; #To keep the same number of tabs independent of present or not 
		    }
		}
		print ANVAR "\n";
	    }
	}
	close (ANVAR);
	return;
    }
}

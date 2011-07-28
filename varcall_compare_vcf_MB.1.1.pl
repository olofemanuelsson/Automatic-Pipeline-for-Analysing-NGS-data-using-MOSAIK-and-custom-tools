#!/usr/bin/perl - w

use strict;
use warnings;

#Compiles and compares variation calls for vcf-format files and M.Bjur format per subject
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS

varcall_comp.pl infile1.vcf infile2_MB outfile

=head2 COMMANDS AND OPTIONS

-o/--outfile The output file (defaults to STDIN)

=head3 I/O

Input format (1st. VCF/Custom 2nd. (MBjur) )

Output format (txt)

=cut

use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
qq{varcall_compare_vcf_MB.1.1.pl < [infile1.vcf] [infile2.MB]  -o outfile.vcf
-o/--outfile The output file (defaults to varcall_comp_out.vcf)
};

}
my ($of, $help) = ("varcall_comp_out.vcf");

GetOptions('o|outfile:s'  => \$of,
'h|help' => \$help,
);

die $USAGE if( $help );

my ($infileVCF, $infilevarCallMB) = @ARGV;

if (@ARGV == 0) {
    my $verbosity = 2;
  print"\n";
  pod2usage({-message => "Must supply an infile.\n",
	     -verbose => $verbosity
	    });
}

my (%chrMTVCF, %chrMTVCFNoRefM, %allVCF);
my ($chr, $matchcount, $matchtrack, $nomatchcount) = 0;
my (@chrMTVCF, @chrMTVCF_sorted, @VCFheader);

if (@ARGV != 1) {

    
    ReadVarCallMB ($infilevarCallMB);
    
    ReadVCF($infileVCF);
   
    sortkeys();
    
    WriteVarCalls($of);
    
    Writehaplotype($of);

    WriteAnnovar($of);

    #WritePolyPhen2($of);
} 

else {
    
    ReadVCF($infileVCF);
    
}


sub ReadVarCallMB {
    
    open(VCFMB, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<VCFMB>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }		
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$');	    #'Loads variant calls
	    shift (@temp); #Removes empty string (MT)
	    if ($temp[2] =~ m/\*|\+|\-/) { #Catch indels. Presently indel called by VCFMB will be ignored. 
		#print $temp[2], "\n";
	    }
	    else {
		$temp[2] =~ /(\w)\=>(\w)\W(\d+)/; #Parses X=>Y and read depth		 
		my $snv = $2;
		my $DP = "DP=";
		$DP .=$3;
		$DP .=";";
		
		if ($2 eq "R") {
		    
		if ($1 eq "A") {
		    $snv = "G";
		}
		else {
		    $snv = "A";
		}	  
		}	
		if ($2 eq "Y") {
		    
		    if ($1 eq "C") {
			$snv = "T";
		    }
		    else {
			$snv = "C";
		    }	  
		}
		if ($2 eq "S") {
		    
		    if ($1 eq "G") {
			$snv = "C";
		    }
		    else {
			$snv = "G";
		    }
		}	
		if ($2 eq "W") {
		    
		    if ($1 eq "A") {
			$snv = "T";
		    }
		    else {
			$snv = "A";
		    }	  
		}
		if ($2 eq "K") {
		    
		    if ($1 eq "G") {
			$snv = "T";
		    }
		    else {
			$snv = "G";
		    }
		}	
		if ($2 eq "M") {
		    
		if ($1 eq "A") {
		    $snv = "C";
		}
		else {
		    $snv = "A";
		}
		}			
		my @temp2 = ($temp[1], ".", $1, $snv, "-", ".", $DP);
		if ($2 =~/[RYSWKM]/) {
		    #push(@temp2, "het")
		    $temp2[6] .= "AF1=0.5";
		}
		#else {
		#	$temp2[6] .= "AF1=1";
		#   }
		$chrMTVCF{$temp[1]} = [@temp2]; #0 based
	    }
	} 
    }
    close(VCFMB);
    print STDERR "Finished Reading infile (VCFMB)","\n";
    return;    	
    
    #$temp[2] =~ /(\w)\=>(\w)\W(\d+)\;([ACGT]):(\d+),([ACGT]):(\d+)/;
    #print STDERR $4, "\t", $5, "\t", $6, "\t", $7, "\n"; 
}


sub ReadVCF {
    
    open(VCF, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<VCF>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	
	if (/##/) {		# Loads VCF header into array (required for IGV)
	    
	    push( @VCFheader, $_);
	    next;
    	}
	if (/#/) {		# Loads Table names into array
	    
	    push( @VCFheader, $_);
	    next;
	    
    	}	
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads variant calls
	    $chr = shift (@temp); #Removes chr string
	    $temp[3] =~ /(\w+)(,)(\w+)/;
	    if ($2) {
		$temp[3] = $1;
		$chrMTVCFNoRefM{$temp[0]} = [@temp];
		$temp[3] = $3;
	    }
	    $chrMTVCF{$temp[0]} = [@temp];	    
	}
    } 	
    close(VCF);
    print STDERR "Finished Reading infile (VCF)","\n";
    return;
}

sub sortkeys {
    
    for my $keys (keys %chrMTVCF)  {
	
	push(@chrMTVCF, $keys);	#Stores keys for sorting later
    }
    
    @chrMTVCF_sorted = sort { $a <=> $b } @chrMTVCF; #Sorts keys to be able to print sorted table later
    print STDERR "Sorted VCF Keys", "\n";   
}

sub Writehaplotype { #For haplogrep format
    
    my $filename = shift;
    print STDERR "NOTE: SampeID in hsd format cannot contain backslash", "\n";
    open (HAPCALL, ">$filename.hsd") or die "Can't write to $filename.hsd: $!\n";
    
    print HAPCALL "SampleId", "\t", "Range", "\t",	"Haplogroup", "\t",	"Polymorphisms", "\n";
    print HAPCALL $of, "\t", $chrMTVCF_sorted[0],"-",$chrMTVCF_sorted[ scalar(@chrMTVCF_sorted)-1 ], "\t", "\t";
    
    for (my $i=0;$i<scalar(@chrMTVCF_sorted);$i++) { #For all sorted keys
	
	my $keys = $chrMTVCF_sorted[$i]; #keys to hash from sorted arrray
	
	if ($chrMTVCF{$keys}[6] =~ m/INDEL/ || $chrMTVCF{$keys}[3] =~ m/,/ ) {
	    
	}
	else {
	    print HAPCALL $chrMTVCF{$keys}[0];
	    print HAPCALL $chrMTVCF{$keys}[3],"\t";
	}
    }
}

sub WriteVarCalls {
    
    
    my $filename = shift;
    
    open (VARCALL, ">$filename") or die "Can't write to $filename: $!\n";
    
    for (my $i=0;$i<scalar(@VCFheader);$i++) {
	
	print VARCALL $VCFheader[$i],"\n";
	
    }
    
    for (my $i=0;$i<scalar(@chrMTVCF_sorted);$i++) { #For all sorted keys
	
	my $keys = $chrMTVCF_sorted[$i]; #keys to hash from sorted arrray
	
	print VARCALL $chr, "\t";	#Adds chr name again
	
	for (my $i=0;$i<scalar( @{$chrMTVCF{$keys} });$i++) {
	    
	    print VARCALL $chrMTVCF{$keys}[$i],"\t";
	}
	print VARCALL "\n";
    }
    close (VARCALL);
    return;
}

sub WriteAnnovar {
    
    
    my $filename = shift;
    
    open (ANNOVAR, ">$filename.annovar.txt") or die "Can't write to $filename: $!\n";
    
    for (my $i=0;$i<scalar(@chrMTVCF_sorted);$i++) { #For all sorted keys
	
	my $keys = $chrMTVCF_sorted[$i]; #keys to hash from sorted arrray
	
	print ANNOVAR $chr, "\t";	#Adds chr name again
	
	if($chrMTVCFNoRefM{$keys}) { #For A,T in vcf, currently is not accepted
	    
	    for (my $k=0;$k<scalar( @{$chrMTVCFNoRefM{$keys} });$k++) {
		print ANNOVAR $chrMTVCFNoRefM{$keys}[$k],"\t";
		
		if ($k eq 0) { #print again and jump ID
		    
		    print ANNOVAR ($chrMTVCFNoRefM{$keys}[$k]+length($chrMTVCFNoRefM{$keys}[2])-1),"\t";
		    $k++;
		}
		
	    }
	    print ANNOVAR "\n";
	    print ANNOVAR $chr, "\t";	#Adds chr name again
	}
	for (my $i=0;$i<scalar( @{$chrMTVCF{$keys} });$i++) {
	    
	    print ANNOVAR $chrMTVCF{$keys}[$i],"\t";
	    
	    if ($i eq 0) { #print again and jump ID
		
		print ANNOVAR ($chrMTVCF{$keys}[$i]+length($chrMTVCF{$keys}[2])-1),"\t";
		$i++;
	    }
	}
	print ANNOVAR "\n";
    }
    close (ANNOVAR);
    return;
}

sub WritePolyPhen2 {
    
    
    my $filename = shift;
    
    open (PP, ">$filename.pp.txt") or die "Can't write to $filename: $!\n";
    
    for (my $i=0;$i<scalar(@chrMTVCF_sorted);$i++) { #For all sorted keys
	
	my $keys = $chrMTVCF_sorted[$i]; #keys to hash from sorted arrray
	
	if ($chrMTVCF{$keys}[6] =~ m/INDEL/ || $chrMTVCF{$keys}[3] =~ m/,/ ) {
	    
	}
	else {
	    print PP "chr$chr:", $keys, " ";	#Adds chr name again
	    print PP $chrMTVCF{$keys}[2], "/", $chrMTVCF{$keys}[3];
	    print PP "\n";	    
	}
    }
    close (PP);
    return;
}


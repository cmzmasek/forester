#!/usr/bin/perl -W

# extractSWISS-PROT.pl
# --------------------
#
# Copyright (C) 2001 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Created: 09/25/01
# Author: Christian M. Zmasek
# zmasek@genetics.wustl.edu
# http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Last modified 05/23/02

# Purpose. Extracts ID, DE, and species from a "sprot.dat" file.
#          The output is used by "rio.pl".
#          If a species list (format: SWISS-PROT-code=full name) is supplied:
#          only sequences from species found in this list are written to
#          outfile (recommended).
#
# Usage.   extractSWISS-PROT.pl <infile> <outfile> [species list] 

# Remark.  Need to re-run this if species in species tree or species list
#          are added/changed or if a new version of Pfam is used!!


use strict;


my $VERSION       = "1.001";
my $infile        = "";
my $outfile       = "";
my $speciesfile   = "";
my $return_line   = "";
my $read          = 0;
my $ac            = "";
my $de            = "";
my $os            = "";
my %Species_names = (); # SWISS-PROT name -> "".
my $i             = 0;

if ( @ARGV != 2 && @ARGV != 3 ) {
    &errorInCommandLine();
}

$infile  = $ARGV[ 0 ];
$outfile = $ARGV[ 1 ];

if ( @ARGV == 3 ) {
    $speciesfile = $ARGV[ 2 ];
    unless ( ( -s $speciesfile ) && ( -f $speciesfile ) && ( -T $speciesfile ) ) {
        die "\n$0: <<$speciesfile>> does not exist, is empty, or is not a plain textfile.\n\n";
    }
    &readSpeciesNamesFile( $speciesfile );
}

if ( -e $outfile ) {
    die "\n$0: <<$outfile>> already exists.\n\n";
}
unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n$0: <<$infile>> does not exist, is empty, or is not a plain textfile.\n\n";
}

open( IN, "$infile" ) || die "\n$0: Cannot open file <<$infile>>: $!\n";
open( OUT, ">$outfile" ) || die "\n$0: Cannot create file <<$outfile>>: $!\n";

print OUT "# extractTrembl.pl version: $VERSION\n"; 
print OUT "# trembl.dat file: $infile\n";
print OUT "# output file    : $outfile\n";
print OUT "# species file   : $speciesfile\n";
print OUT "# date           : ".`date`."\n\n";

$read = 0; 

while ( $return_line = <IN> ) {
    if ( $return_line =~ /^ID\s+(\S+)/ ) {
        $ac = $1;
        $read = 1;
        if ( $ac =~ /[A-Z0-9]+_([A-Z0-9]+)/ ) {
            $os = $1;
        }
        else {
            die "\n$0: Unexpected format: $ac.\n\n";
        }
        if ( $speciesfile ne "" ) {
            unless ( exists( $Species_names{ $os } ) ) {
                $read = 0;
                $ac = "";
                $de = "";
                $os = "";
                next;
            }
        }
    }
    elsif ( $return_line =~ /^DE\s+(.+)/ && $read == 1 ) {
        if ( $de ne "" ) {
            $de .= " ".$1;
        }
        else {
            $de = $1;
        }
    }
    elsif ( $return_line =~ /^\/\// && $read == 1 ) {
        $read = 0;
        print OUT "$ac;$de;$os\n";
        $ac = "";
        $de = "";
        $os = "";
        $i++;
    }
}

close( IN ); 

print OUT "\n # $i entries.\n";

close( OUT );

exit( 0 );



# Reads in species file.
# Format: SWISS-PROT=full name (e.g. "BACSU=Bacillus subtilis")
# Lines beginning with "#" are ignored.
# One argument: species file-name
# Last modified: 04/24/01
sub readSpeciesNamesFile {
    my $infile = $_[ 0 ];
    my $return_line = "";
    my $sp          = "";
    my $full        = "";

    unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
        die "\n$0: readSpeciesNamesFile: <<$infile>> does not exist, is empty, or is not a plain textfile.\n";
    }

    open( IN_RSNF, "$infile" ) || die "\n$0: Cannot open file <<$infile>>: $!\n";
    while ( $return_line = <IN_RSNF> ) {
        if ( $return_line !~ /^\s*#/ && $return_line =~ /(\S+)=(.+)/ ) {
            $sp   = $1;
            $full = $2;
            $full =~ s/^\s+//;
            $full =~ s/\s+$//;
            $Species_names{ $sp } = "";
        }
    }
    close( IN_RSNF );

    return;
}



sub errorInCommandLine {
    print "\n";
    print " extractSWISS-PROT.pl  $VERSION\n";
    print " --------------------\n";
    print "\n";
    print " Christian Zmasek (zmasek\@genetics.wustl.edu)\n";
    print "\n";
    print " Purpose. Extracts ID, DE, and species from a \"sprot.dat\" file.\n";
    print "          The resulting output is used by \"rio.pl\".\n";
    print "          If a species list (format: SWISS-PROT-code=full name) is supplied:\n";
    print "          only sequences from species found in this list are written to\n";
    print "          outfile (recommended).\n";
    print "\n";
    print " Usage.   extractSWISS-PROT.pl <infile> <outfile> [species list]\n";
    print "\n";
    print " Remark.  Need to re-run this if species in species tree or species list\n";
    print "          are added/changed or if a new version of Pfam is used!!\n";
    print "\n\n";
    exit( -1 );
}

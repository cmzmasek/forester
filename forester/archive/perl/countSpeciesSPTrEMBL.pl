#!/usr/bin/perl -W

# countSpeciesSPTrEMBL.pl
# -----------------------
#
# Copyright (C) 2003 Christian M. Zmasek
# All rights reserved
#
# Created: 02/27/03
# Last modified: 02/27/03
# Author: Christian M. Zmasek
# zmasek@genetics.wustl.edu
# http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Last modified 05/23/02

# Purpose. Counts species in SWISS-PROT and TrEMBL.
#
# Usage.   countSpeciesSPTrEMBL.pl <path/to/trembl.dat> <path/to/sprot.dat> <outfile> 
#


use strict;


my $VERSION       = "1.000";
my $infile_sp     = "";
my $infile_tr     = "";
my $outfile       = "";

my $return_line   = "";
my $read          = 0;
my $os            = "";
my %species_count = (); # full name -> count.


if ( @ARGV != 3 ) {
    &errorInCommandLine();
}

$infile_tr = $ARGV[ 0 ];
$infile_sp = $ARGV[ 1 ];
$outfile   = $ARGV[ 2 ];



if ( -e $outfile ) {
    die "\n$0: <<$outfile>> already exists.\n\n";
}
unless ( ( -s $infile_tr ) && ( -f $infile_tr ) && ( -T $infile_tr ) ) {
    die "\n$0: <$infile_tr>> does not exist, is empty, or is not a plain textfile.\n\n";
}
unless ( ( -s $infile_sp ) && ( -f $infile_sp ) && ( -T $infile_sp ) ) {
    die "\n$0: <<$infile_sp>> does not exist, is empty, or is not a plain textfile.\n\n";
}

open( IN_TR, "$infile_tr" ) || die "\n$0: Cannot open file <<$infile_tr>>: $!\n";
open( IN_SP, "$infile_sp" ) || die "\n$0: Cannot open file <<$infile_sp>>: $!\n";
open( OUT, ">$outfile" ) || die "\n$0: Cannot create file <<$outfile>>: $!\n";


$read = 0; 

while ( $return_line = <IN_TR> ) {
    if ( $return_line =~ /^AC\s+(\S+);/ ) {
        $read = 1;
    }
    elsif ( $return_line =~ /^OS\s+(.+)\.\s*$/ && $read == 1 ) {
        $os = $1;
        $os =~ s/\(.+\)//g;
        $os =~ s/^\s+//;
        $os =~ s/\s+$//;
        $os =~ s/\.$//;
        if ( exists( $species_count{ $os } ) ) {
            $species_count{ $os } = $species_count{ $os } + 1;   
        }
        else {
            $species_count{ $os } = 1;
        }
        print "$os\n";
    }
    elsif ( $return_line =~ /^\/\// && $read == 1 ) {
        $read = 0;
        $os = "";
    }
}

close( IN_TR ); 

$read = 0;
$os = "";
$return_line =  "";

while ( $return_line = <IN_SP> ) {
    if ( $return_line =~ /^ID\s+(\S+)/ ) {
        $read = 1;
    }
    elsif ( $return_line =~ /^OS\s+(.+)\s*$/ && $read == 1 ) {
        $os = $1;
        $os =~ s/\(.+//g;
        $os =~ s/^\s+//;
        $os =~ s/\s+$//;
        $os =~ s/\.$//;
        $read = 0;
        if ( exists( $species_count{ $os } ) ) {
            $species_count{ $os } = $species_count{ $os } + 1;   
        }
        else {
            $species_count{ $os } = 1;
        }
        print "$os\n";
    }
    elsif ( $return_line =~ /^\/\// && $read == 1 ) {
        $read = 0;
        $os = "";
    }
}

close( IN_SP );


foreach my $species ( sort { $species_count{ $b } <=> $species_count{ $a } } keys %species_count ) {
    print OUT "$species: $species_count{$species}\n";
}


print "\n\nDone!\n\n";

close( OUT );

exit( 0 );






sub errorInCommandLine {
    print "\n";
    print " countSpeciesSPTrEMBL.pl $VERSION\n";
    print " -----------------------\n";
    print "\n";
    print " Christian Zmasek (zmasek\@genetics.wustl.edu)\n";
    print "\n";
    print " Purpose. Counts species in SWISS-PROT and TrEMBL.\n";
    print "\n";
    print " Usage.   countSpeciesSPTrEMBL.pl <path/to/trembl.dat> <path/to/sprot.dat> <outfile>\n";
    print "\n";
    exit( -1 );
}

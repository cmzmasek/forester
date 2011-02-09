#!/usr/bin/perl -W

# extractSpecies.pl
# ----------------
#
# Copyright (C) 2003 Christian M. Zmasek
# All rights reserved
#
# Created: 09/03/03
# Author: Christian M. Zmasek
# zmasek@genetics.wustl.edu
# http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Last modified 03/12/04 (Added gg)

# Purpose. Adds species information to a file describing a phylogenetic
#          tree in the following format (by way of example):
#          "((ceINX_CE33055:0.02883,cbINX_CB09748:0.02934):0.36899[&&NHX:B=100],..."
#          ce stands for "CAEEL". The hash  %SPECIES needs to be set accordingly.
#  


use strict;


my %SPECIES = ( 
                "dm" => "DROME",
                "ag" => "ANOGA",
                "ce" => "CAEEL",
                "cb" => "CAEBR",
                "ci" => "CIOIN",
                "fr" => "FUGRU",
                "gg" => "CHICK",
                "rn"  => "RAT",
                "mm"  => "MOUSE",
                "hs"  => "HUMAN"
               ); 
                

my $infile        = "";
my $outfile       = "";
my $intree        = "";
my $return_line   = "";

if ( @ARGV != 1 && @ARGV != 2 ) {
    &errorInCommandLine();
}

$infile = $ARGV[ 0 ];

if ( @ARGV == 1 ) {
    $outfile = $infile;
    $outfile =~ s/\.nhx$//;
    $outfile .= "_species.nhx";
}

if ( @ARGV == 2 ) {
    $outfile = $ARGV[ 1 ];
}




if ( -e $outfile ) {
    die "\n$0: <<$outfile>> already exists.\n\n";
}
unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n$0: <<$infile>> does not exist, is empty, or is not a plain textfile.\n\n";
}

open( IN, "$infile" ) || die "\n$0: Cannot open file <<$infile>>: $!\n";
open( OUT, ">$outfile" ) || die "\n$0: Cannot create file <<$outfile>>: $!\n";

while ( $return_line = <IN> ) {
    $return_line =~ s/\s+//g;
    $return_line =~ s/\+/_/g;

    $intree .= $return_line;

}

close( IN ); 

while ( ( my $short, my $long ) = each ( %SPECIES ) ) {
    
    while ( $intree =~ /[(),]($short[^\[]+?)[(),]/ ) {
        
        my $name_and_length = $1;
        
        print "$name_and_length   ->   $name_and_length\[\&\&NHX:S=$long\]\n";
        
        $intree =~ s/$name_and_length/$name_and_length\[&&NHX:S=$long\]/;
    
    }

}

print OUT $intree;

close( OUT );

print "\n\nDone!\n\n";

exit( 0 );



sub errorInCommandLine {
    print "\n";
    print "extractSpecies.pl infile [outfile]";
    print "\n\n";
    exit( -1 );
}

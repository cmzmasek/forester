#!/usr/bin/perl -W

# $Id: gs_aa_extract.pl,v 1.2 2008/03/09 00:11:50 cmzmasek Exp $

# This extracts the AA sequences from GENSCAN output files
# Copyright (C) 2008-2009 Christian M. Zmasek
# All rights reserved
# Created 2007-07-28 in Winterthur, Switzerland by CMZ

# Usage: gs_aa_extract.pl <genscan-output infile> <outfile>

use strict;

if ( scalar( @ARGV ) != 2 ) {
    print "\ngs_aa_extract.pl <genscan-output infile> <outfile>\n\n";
    exit( -1 );
}

my $infile  = $ARGV[ 0 ];
my $outfile = $ARGV[ 1 ];

if ( -e $outfile) {
    die "\n$0: \"$outfile\" already exists.\n\n";
}
unless( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n$0: cannot read from \"$infile\".\n\n";
}

open( IN, "$infile" ) || die "\n$0: Cannot open file \"$infile\": $!\n";
open( OUT, ">$outfile" ) || die "\n$0: Cannot create file \"$outfile\": $!\n";

my $line = "";
my $desc = "";

while ( $line = <IN> ) {
    if ( $line =~ /^>/ ) {
        $desc = $line;
    }
    elsif ( $line =~ /^[A-Z]+$/ ) {
        if ( length( $desc ) > 0 ) {
            print OUT $desc;
            $desc = "";
        }
        print OUT $line;
    }
}

close( OUT );

print( "\nOK\n" );

exit( 0 );


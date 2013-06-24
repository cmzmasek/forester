#!/usr/bin/perl -W

# $Id: pf_cutoff_extract.pl,v 1.4 2009/11/11 02:28:19 cmzmasek Exp $

# This extracts GA, TC, or NC score cutoff values from
# Pfam HMM files (GA1, TC1, NC1)
# Copyright (C) 2008-2009 Christian M. Zmasek
# All rights reserved
# Created 2007-08-01 in Winterthur, Switzerland by CMZ

# Usage: pf_cutoff_extract.pl <Pfam HMM file> <GA|TC|NC> <outfile>

use strict;

if ( scalar( @ARGV ) != 3 ) {
    print "\npf_cutoff_extract.pl <Pfam HMM file> <GA|TC|NC> <outfile>\n\n";
    exit( -1 );
}

my $infile      = $ARGV[ 0 ];
my $cutoff_type = uc( $ARGV[ 1 ] );
my $outfile     = $ARGV[ 2 ];

my $GA = "GA";
my $TC = "TC";
my $NC = "NC"; 

if ( -e $outfile ) {
    die "\n$0: \"$outfile\" already exists.\n\n";
}
unless( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n$0: cannot read from \"$infile\".\n\n";
}
unless( $cutoff_type eq $GA || $cutoff_type eq $TC || $cutoff_type eq $NC ) {
    die "\n$0: illegal value \"$cutoff_type\" for cutoff type.\n\n";
}

open( IN, "$infile" ) || die "\n$0: Cannot open file \"$infile\": $!\n";
open( OUT, ">$outfile" ) || die "\n$0: Cannot create file \"$outfile\": $!\n";

my $line        = "";
my $name        = "";
my $line_number = 0;
my $n           = 0;

while ( $line = <IN> ) {
    $line_number++;
    if ( $line =~ /^NAME\s+(.+)/ ) {
        if ( length( $name ) > 0 ) {
            die "\n$0: Unexpected line $line at line $line_number: $!\n";
        }
        $name = $1;
    }
    elsif ( $line =~ /^$cutoff_type\s+(\S+)\s+[^;]+/ ) {
        if ( length( $name ) < 1 ) {
            die "\n$0: Unexpected line $line at line $line_number: $!\n";
        }
        $n++;
        print OUT "$name $1\n";
        $name = "";
    }
    elsif ( $line =~ /\/\// ) {
       $name = ""; 
    }
}

close( OUT ) || die "\n$0: Cannot close file \"$outfile\": $!\n";;

print( "\nExtracted $n $cutoff_type" . "1 values to \"$outfile\"\n" );
print( "\nOK\n" );

exit( 0 );


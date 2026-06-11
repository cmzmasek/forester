#!/usr/bin/perl -W

# forester -- software libraries and applications
# for evolutionary biology and genomics.
# Copyright (C) 2026 Christian M. Zmasek
# All rights reserved
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact: czmasek at jcvi dot org

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


#!/usr/bin/perl -W

# rio_slave.pl
# ------------
#
# Copyright (C) 2002 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Created: 01/18/02
# Author: Christian M. Zmasek
# zmasek@genetics.wustl.edu
# http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Last modified: 02/20/02


# Arguments:

# 0: first block in multiple alignment to process
# 1: last block in multiple alignment to process
# 2: name of resampled alignment, inc. query
# 3: matrix number
# 4: name of query
# 5: PWD file
# 6: seed for random number generator for neighbor
# 7: node number
# 8: temp dir


use strict;

use FindBin;
use lib $FindBin::Bin;
use rio_module;

if ( @ARGV != 9 ) {
    &dieWithUnexpectedError( "argument count is off" );
}

my $start       = $ARGV[ 0 ];
my $end         = $ARGV[ 1 ];
my $align       = $ARGV[ 2 ];
my $matrix_n    = $ARGV[ 3 ];
my $name        = $ARGV[ 4 ];
my $pwd_file    = $ARGV[ 5 ];
my $seed        = $ARGV[ 6 ];
my $number      = $ARGV[ 7 ];
my $temp_dir    = $ARGV[ 8 ];

my $b           = 0;
my $outfile     = "";
my $mytemp_dir  = $temp_dir."/dir_".$number;

mkdir(  $mytemp_dir, 0700 )
|| &dieWithUnexpectedError( "Could not create \"$mytemp_dir\"" );

unless ( ( -e $mytemp_dir ) && ( -d $mytemp_dir ) ) {
    &dieWithUnexpectedError( "\"$mytemp_dir\" does not exist, or is not a directory" );
}
 
       
&executePuzzleDQObootstrapped( $align, $matrix_n );
       
system( "mv", $align.".dist", $mytemp_dir."/DISTs_TO_QUERY" )
&& &dieWithUnexpectedError( "could not mv" );

unlink( $align );
  
sleep( 2 );

&dividePWDfile( $pwd_file,
                $mytemp_dir."/DIVIDED",
                $start,
                $end );

&addDistsToQueryToPWDfile( $mytemp_dir."/DIVIDED",
                           $mytemp_dir."/DISTs_TO_QUERY",
                           $mytemp_dir."/PWD_INC_QUERY",
                           $name );
                            
unlink( $mytemp_dir."/DIVIDED" );                           

$b = $end - $start + 1;

chdir ( $mytemp_dir ) 
|| &dieWithUnexpectedError( "Could not chdir to \"$mytemp_dir\"" );

&executeNeighbor( $mytemp_dir."/PWD_INC_QUERY",
                  $b,
                  1, # randomize input order
                  $seed,
                  1 ); # lower-triangular data matrix


unlink( "outfile", $mytemp_dir."/PWD_INC_QUERY", $mytemp_dir."/DISTs_TO_QUERY" );

system( "mv", "outtree", "../MAKETREEOUT".$MULTIPLE_TREES_FILE_SUFFIX.$number )
&& &dieWithUnexpectedError( "could not mv" );   

sleep( 1 );

chdir( ".." )
|| &dieWithUnexpectedError( "Could not chdir to \"..\"" );

rmdir( $mytemp_dir ) || &dieWithUnexpectedError( "Could not delete \"$mytemp_dir\"" );

$outfile = "FINISHED_$number";

open( OUT, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
close( OUT );

exit( 0 );




sub dividePWDfile {
    my $pwd_file        = $_[ 0 ];
    my $outfile         = $_[ 1 ];
    my $start           = $_[ 2 ]; # e.g. 0
    my $end             = $_[ 3 ]; # e.g. 9
    
    my $c               = 0;
    my $write           = 0;
    my $return_line     = "";
    
    &testForTextFilePresence( $pwd_file );
    
    open( IN_PWD, "$pwd_file" ) || &dieWithUnexpectedError( "Cannot open file \"$pwd_file\"" );
    open( OUT_PWD, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
    
    while ( $return_line = <IN_PWD> ) {
        if ( $return_line =~ /^\s*(\d+)\s*$/ ) {
            if ( $c >= $start && $c <= $end ) {
                $write = 1;
            }
            elsif ( $c > $end ) {
                last;
            }
            $c++;
        }
        if ( $write == 1 ) {
            print OUT_PWD $return_line;
        }    
    }
    
    close( IN_PWD );
    close( OUT_PWD );
    
    return;
 
} ## dividePWDfile








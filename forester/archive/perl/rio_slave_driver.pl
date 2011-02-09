#!/usr/bin/perl -W

# rio_slave_driver.pl
# -------------------
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


# 0: block size
# 1: number of blocks which have a size of block size + 1
# 2: name of resampled alignment, inc. query
# 3: matrix number
# 4: name of query
# 5: PWD file
# 6: temp dir
# 7: seed for random number generator for neighbor
# 8...: list of node names




use strict;

use FindBin;
use lib $FindBin::Bin;
use rio_module;

if ( @ARGV < 9 ) {
    &dieWithUnexpectedError( "argumnet count off" ); 
}


my $block_size    = shift( @ARGV );
my $larger_blocks = shift( @ARGV );
my $align         = shift( @ARGV );
my $matrix_n      = shift( @ARGV );
my $name          = shift( @ARGV );
my $pwd_file      = shift( @ARGV );
my $temp_dir      = shift( @ARGV );
my $seed          = shift( @ARGV );
my @nodelist      = @ARGV;
my $start         = 0;   
my $end           = 0;
my $x             = 0;
my $node          = "";


$start = 0;

if ( $larger_blocks > 0 ) {
    $end = $block_size;
}
else {
    $end = $block_size - 1;
}

for ( $x = 0; $x < scalar( @nodelist ); $x++ ) {
    my $child_pid;
    $node = $nodelist[ $x ];
	    
    if ( !defined( $child_pid = fork() ) ) {
        &dieWithUnexpectedError( "cannot fork" );
    }
    elsif ( $child_pid ) {
        # I'm the parent, forking off $nodelist number of children
    }
    else {	
        exec( "ssh",
              $node,
              "/usr/bin/perl",
              $RIO_SLAVE,
              $start,
              $end,
              $align.$x,
              $matrix_n,
              $name,
              $pwd_file,
              $seed,
              $x,
              $temp_dir )
              || &dieWithUnexpectedError( "could not \"exec ssh $node /usr/bin/perl $RIO_SLAVE\"" );

    }
    $larger_blocks--;   
    if ( $larger_blocks > 0 ) {
        $start += ( $block_size + 1 );
        $end   += ( $block_size + 1 );
    }
    elsif ( $larger_blocks == 0 ) {
        $start += ( $block_size + 1 );
        $end   += $block_size;
    }
    else {
        $start += $block_size;
        $end   += $block_size;
    }
}
      
exit( 0 );

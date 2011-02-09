#!/usr/bin/perl -w

# bootstrap_cz.pl
# ---------------
# Copyright (C) 1999-2003 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek 
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Created: 05/17/01
#
# Last modified 08/26/03
#
# Purpose:
# Bootstrap resamples an alignment in PHYLIP sequential format <bootstraps>
# times.
# Amino acid sequences must only be represented by uppercase letters (A-Z)
# and '-'.
# In mode 0 it saves the positions which it used to create the
# bootstrapped alignment into <positions outfile>.
# Mode 1 allows to recreate exactly the same boostrapped alignment
# by reading in a <positions infile>.
# Sequence names are normalized to $LENGTH_OF_NAME characters.
# The output alignment is in PHYLIP's sequential or interleaved format.
# (These two are the same in this case, since all the seqs will be one
# line in length (no returns in seq).)
# 
# Usage:
# bootstrap_cz.pl <mode (0 or 1)> <bootstraps> <alignment infile> 
# <alignment outfile> <positions out- (mode 0) or in-file (mode 1)>
# [random number seed (mode 0 only)]
# 

use strict;
use FindBin;
use lib $FindBin::Bin;

use rio_module;

my $VERSION        = "2.001";

my $modus          = -1;  # 0 to create pos. file, 1 to use premade pos. file
my $bootstraps     = -1;
my $infile         = "";
my $outalign_file  = "";
my $positions_file = "";
my $seed           = -1;


$modus          = $ARGV[ 0 ];
$bootstraps     = $ARGV[ 1 ];
$infile         = $ARGV[ 2 ];
$outalign_file  = $ARGV[ 3 ];
$positions_file = $ARGV[ 4 ];
$seed           = $ARGV[ 5 ];

if ( @ARGV != 5 && @ARGV != 6 ) {
    &printUsage();
    exit( -1 );
}

if ( $modus != 0 && $modus != 1 ) {
    &printUsage();
    exit( -1 );
}

if ( $modus == 0 && @ARGV != 6 ) {
    &printUsage();
    exit( -1 );
}

if ( $modus == 1 && @ARGV != 5 ) {
    &printUsage();
    exit( -1 );
}

if ( $bootstraps < 1 ) {
    &printUsage();
    exit( -1 );
}

if ( $seed && $seed < 0 ) {
    &printUsage();
    exit( -1 );
}


unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n\nbootstrap_cz.pl: \"$infile\" does not exist, is empty, or is not a plain textfile.\n\n";
}
if ( -e $outalign_file ) {
    die "\n\nbootstrap_cz.pl: \"$outalign_file\" already exists.\n\n";
}

if ( $modus == 0 ) {
    if ( -e $positions_file ) {
        die "\n\nbootstrap_cz.pl: \"$positions_file\" already exists.\n\n";
    }
}
else {
    unless ( ( -s $positions_file ) && ( -f $positions_file ) && ( -T $positions_file ) ) {
        die "\n\nbootstrap_cz.pl: \"$positions_file\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
}

if ( $modus == 0 ) {
    &bootstrap( $modus, $bootstraps, $infile, $outalign_file, $positions_file, $seed );
}
else {
    &bootstrap( $modus, $bootstraps, $infile, $outalign_file, $positions_file );
}


exit( 0 );



# Methods
# -------


# Five/six arguemnts:
# 1. Mode: 0 to create pos. file, 1 to use premade pos. file
# 2. bootstraps
# 3. Alignment infile name
# 4. Outfile name
# 5. file name for positions file (created if mode is 0, read if mode is 1)
# [6. If modus is 0: seed for random number generator]
#
# This method is very similar to method "pfam2phylip" "in makeTree.pl".
#
# Last modified: 05/17/01
#
sub bootstrap { 

    my $modus           = $_[ 0 ];
    my $bootstraps      = $_[ 1 ];
    my $infile          = $_[ 2 ];
    my $outalign_file   = $_[ 3 ];
    my $positions_file  = $_[ 4 ];
  
    my @seq_name        = ();
    my @seq_array       = ();
    my @random_numbers  = ();
    my $return_line     = "";
    my $seq             = "";
    my $x               = 0;
    my $y               = 0;
    my $seq_no          = 0;
    my $original_length = 0;
    my $max_x           = 0;
    my $n               = 0;
    my $i               = 0;
    my $random          = -1;
    my $length          = 0;
    my $number_of_seqs  = 0;
    my $number_of_colm  = 0;


    # Checks the arguments
    # --------------------
 
    if ( $modus == 0 ) {
        if ( !$_[ 5 ] ) {
            die "\n\n$0: bootstrap: Failed to give a seed for random number generator.\n\n";
        }
        srand( $_[ 5 ] );
    } 
    elsif( $modus == 1 ) {
        if ( $_[ 5 ] ) {
            die "\n\n$0: bootstrap: Must not give a seed for random number generator.\n\n";
        }
        unless ( ( -s $positions_file ) && ( -f $positions_file ) && ( -T $positions_file ) ) {
            die "\n\n$0: bootstrap: <<$positions_file>> does not exist, is empty, or is not a plain textfile.\n\n";
        }
    }
    else {
        die "\n\n$0: bootstrap: modus must be either 0 or 1.\n\n";
    }

    unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
        die "\n\n$0: bootstrap: <<$infile>> does not exist, is empty, or is not a plain textfile.\n\n";
    }



    # Reads in the alignment
    # ----------------------

    open( IN, "$infile" ) || die "\n$0: bootstrap: Cannot open file <<$infile>>: $!";
    while ( $return_line = <IN> ) {
    
        if ( $return_line =~ /^\s*(\d+)\s+(\d+)/ ) {
            $number_of_seqs = $1;
            $number_of_colm = $2;
        }
        elsif ( $return_line =~ /^(\S+)\s+(\S+)/ ) {
            $seq_name[ $seq_no ] = substr( $1, 0, $LENGTH_OF_NAME );
            $seq = $2;
            if ( $original_length == 0 ) {
                $original_length = length( $seq );
            }
            elsif ( $original_length != length( $seq ) ) {
                die "\n\n$0: Sequences do not have the same length.\n\n";
            }
            for ( $x = 0; $x < $original_length; $x++ ) {
                $seq_array[ $x ][ $seq_no ] = substr( $seq, $x, 1 );
            }
            $seq_no++;
        }
    }
    close( IN );

    if ( ( $number_of_seqs != $seq_no ) 
    || ( $number_of_colm != $original_length ) ) {
        die "\n\n$0: Number of sequences or number of columns are inconsisten with the values given in the alignment.\n\n";
    }

    # Adusts the length of the names to $LENGTH_OF_NAME
    # -------------------------------------------------
  
    for ( $y = 0; $y < $seq_no; $y++ ) {
        $length = length( $seq_name[ $y ] );
        for ( $i = 0; $i <= ( $LENGTH_OF_NAME - $length - 1 ); $i++ ) {
	        $seq_name[ $y ] .= " ";
        }
    }


   
    # Bootstraps $bootstraps times and writes the outputfiles
    # -------------------------------------------------------
    
    open( OUT, ">$outalign_file" ) || die "\n\n$0: bootstrap: Cannot create file <<$outalign_file>>: $!";
    if ( $modus == 0 ) {
        open( OUT_P, ">$positions_file" ) || die "\n\n$0: bootstrap: Cannot create file <<$positions_file>>: $!";
    }
    else {
        open( IN_P, "$positions_file" ) || die "\n\n$0: bootstrap: Cannot open file <<$positions_file>>: $!";
    }

    for ( $n = 0; $n < $bootstraps; $n++ ) {
        
        if ( $modus == 0 ) {
            for ( $x = 0; $x < $original_length; $x++ ) {
                $random = int( rand( $original_length ) );
                print OUT_P "$random ";
                $random_numbers[ $x ] = $random;
            }
            print OUT_P "\n";
        }
        else {
            $return_line = <IN_P>;
            if ( !$return_line || $return_line !~ /\d/ ) {
                die "\n\n$0: bootstrap: <<$positions_file>> seems too short or otherwise unsuitable.\n\n";
            }
            $return_line =~ s/^\s+//;
            $return_line =~ s/\s+$//;
            @random_numbers = split( /\s+/, $return_line );
            if ( scalar( @random_numbers ) != $original_length ) {
               die "\n\n$0: bootstrap: <<$positions_file>> seems not to correspond to <<$infile>>.\n\n";
            }
        }       

        print OUT " $seq_no  $original_length\n";

        for ( $y = 0; $y < $seq_no; $y++ ) {
            print OUT "$seq_name[ $y ]";
            
            for ( $x = 0; $x < $original_length; $x++ ) {
                $random = $random_numbers[ $x ];
                if ( !$seq_array[ $random ][ $y ] || $seq_array[ $random ][ $y ] !~ /[A-Z]|-/ ) {
                    die "\n\n$0: Sequence must be represented by uppercase letters A-Z and \"-\" only.\n\n";
                }    
                print OUT $seq_array[ $random ][ $y ];
            }
            print OUT "\n";
        }
    }
    
    close( OUT );

    if ( $modus == 0 ) {
        print OUT_P "\n";
        close( OUT_P );  
    }
    else {
        close( IN_P );
    }

    return;

} ## bootstrap



sub printUsage {
    print "\n";
    print " bootstrap_cz.pl  $VERSION\n";
    print " ---------------\n";
    print "\n";
    print " Christian Zmasek (zmasek\@genetics.wustl.edu)\n";
    print "\n";
    print " Purpose:\n";
    print " Bootstrap resamples an alignment in PHYLIP sequential format\n";
    print " <bootstraps> times.\n";
    print " In mode 0 it saves the positions which it used to create the\n";
    print " bootstrapped alignment into <positions outfile>.\n";
    print " Mode 1 allows to recreate exactly the same boostrapped alignment\n";
    print " by reading in a <positions infile>.\n";
    print " Sequence names are normalized to $LENGTH_OF_NAME characters.\n";
    print " The output alignment is in PHYLIP's sequential or interleaved format.\n";
    print " (These two are the same in this case, since all the seqs will be one\n";
    print " line in length (no returns in seq).)\n";
    print "\n";
    print " Usage:\n";
    print " bootstrap_cz.pl <mode (0 or 1)> <bootstraps> <alignment infile>\n";
    print " <alignment outfile> <positions out (mode 0) or infile (mode 1)>\n";
    print " [random number seed (mode 0 only)]\n";
    print "\n";
} ## printUsage


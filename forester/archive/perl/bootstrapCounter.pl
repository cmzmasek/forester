#!/usr/bin/perl -w
#
# bootstrapCounter.pl
# -------------------
#
# Copyright (C) 2001 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek 
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Created: 04/04/01
#
# Last modified 08/16/01
#
#
# Objective. Determines the distribution of top orthology bootstrap values
#            of a Xrio.pl output file.
#
# Usage. "bootstrapCounter.pl <infile = Xrio.pl-output> <outfile>"
#
# Important. The result of this is meaningful ONLY if the thresholds 
#            for output of the RIO analysis are set to zero (L=0 R=0).
#
# Format for infile: 
# ... 
#
# # ############################################################################
# # Annotation: B0511.6 CE17345   helicase (ST.LOUIS) TR:O61815 protein_id:AAC17654.1
# # HMM       : ABC_tran
# # score     : -59.6
# # E-value   : 1.1
# # Query has not been aligned (score lower than gathering cutoff).
# # ############################################################################
# 
# 
# # ############################################################################
# # Annotation: B0511.7 CE17346    (ST.LOUIS) TR:O61817 protein_id:AAC17655.1
# # HMM       : FHA
# # score     : 71.6
# # E-value   : 1.7e-17
# RIO - Resampled Inference of Orthologs
# Version: 1.000
# ------------------------------------------------------------------------------
# Alignment file: /tmp/Xriopl9846081980/Full-FHA
# Alignment     : FHA domain
# HMM           : FHA
# Query file    : /tmp/Xriopl9846081980/__queryfile__
# ==============================================================================
# 
# Query         : CE17346.FHA_CAEEL/45-114
# 
# Number (in %) of observed orthologies (o) and super orthologies (s) to query
# in bootstrapped trees, evolutionary distance to query:
#  
# Sequence              Description                                                                  # o[%] s[%]  distance
# --------              -----------                                                                  ---- ----  --------
# YC67_MYCTU/308-372    -                                                                              20   14  1.577840
# FRAH_ANASP/204-277    FRAH PROTEIN.                                                                  17   16  1.532670
# ABA2_NICPL/557-633    ZEAXANTHIN EPOXIDASE PRECURSOR (EC 1.14.-.-).                                  14   11  1.885700
# ABA2_LYCES/563-639    ZEAXANTHIN EPOXIDASE PRECURSOR (EC 1.14.-.-).                                  14   11  2.140000
# 
# 
# 
# Distance values (based on ML branch length values on consensus tree)
# --------------------------------------------------------------------
# Given the thresholds for distance calculations:
# No sequence is considered orthologous to query.
#
# ... 
                                                        


use strict;

my $VERSION            = 0.200;

my $infile             = "";       
my $outfile            = "";
my $return_line        = "";
my $results            = 0;
my $o_bootstraps       = 0;
my $s_bootstraps       = 0;
my @o_bootstraps_array = ();
my @s_bootstraps_array = ();
my $total              = 0;
my $i                  = 0;


if ( @ARGV != 2 ) {
    &errorInCommandLine();
    exit ( -1 ); 
}

$infile  = $ARGV[ 0 ];
$outfile = $ARGV[ 1 ];

if ( -e $outfile ) {
    die "\n$0: <<$outfile>> already exists.\n";
}
unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n$0: <<$infile>> does not exist, is empty, or is not a plain textfile.\n";
}


open( IN, "$infile" ) || die "\n$0: Cannot open file <<$infile>>: $!\n";

$results = 0;
for ( $i = 0; $i <= 100; ++$i ) {
    $s_bootstraps_array[ $i ] = $o_bootstraps_array[ $i ] = 0;
}

while ( $return_line = <IN> ) {

    if ( $return_line =~ /^\s*--------\s+/ ) {
        $results = 1;
    }
    elsif ( $return_line =~ /^\s*Distance\s+values\s+/i ) {
        $results = 0;
    }
    elsif ( $results == 1 && $return_line =~ /^\s*!NO\s+ORTHOLOGS/ ) {
        $o_bootstraps_array[ 0 ]++;
        $s_bootstraps_array[ 0 ]++;
        $total++;
        $results = 0;
    }
    elsif ( $results == 1 && $return_line =~ /(\S+)\s+(\S+)\s+\S+\s*$/ ) {
        $o_bootstraps = $1;
        $s_bootstraps = $2;
        $results = 0;
        if ( $o_bootstraps > 100 || $s_bootstraps > 100 
        || $o_bootstraps < 0 ) {
            print "o bootstraps: $o_bootstraps\n";
            print "s bootstraps: $s_bootstraps\n";
            die "\n\n$0: Error: Boostrap value(s) out of range.\n\n";
        }
        
        $total++;
        $o_bootstraps_array[ $o_bootstraps ]++;
        $s_bootstraps_array[ $s_bootstraps ]++;
        
    }
}

close( IN );


open( OUT, ">$outfile" ) || die "\n$0: Cannot create file \"$outfile\": $!\n";

print OUT "bootstrapCounter.pl version: $VERSION\n\n";
print OUT "Distribution of top bootstrap values\n\n";
print OUT "Input file : $infile\n";
print OUT "Output file: $outfile\n";
print OUT "Date       : ".`date`."\n";
print OUT "Total: $total\n\n";
print OUT "top-orthology-bootstraps vs. count:\n";
for ( $i = 0; $i < @o_bootstraps_array; ++$i ) {
    print OUT "$i $o_bootstraps_array[ $i ]\n";
}
print OUT "\n\ntop-super-orthology-bootstraps vs. count:\n";
for ( $i = 0; $i < @s_bootstraps_array; ++$i ) {
    print OUT "$i $s_bootstraps_array[ $i ]\n";
}
close( OUT );

print( "\nDone.\n\n" );

exit( 0 );



sub errorInCommandLine {
    print "\n";
    print " bootstrapCounter.pl version: $VERSION\n";
    print " Usage: \"bootstrapCounter.pl <infile = Xrio.pl-output> <outfile>\"\n";
    print " Important: The result of this is meaningful ONLY if the thresholds\n"; 
    print " for output of the RIO analysis are set to zero (L=0 R=0).\n";
    print "\n";
    exit( -1 );
}



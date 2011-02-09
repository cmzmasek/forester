#!/usr/bin/perl -w
#
# bootstrapSelector.pl
# --------------------
#
# Copyright (C) 2001 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek 
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Created: 04/06/01
#
# Last modified 09/24/01
#
#
# Objective. Selection of RIO analysis results with top ortholgy
#            bootstrap values greater or less than a threshold.
#
# Usage: "bootstrapSelector.pl <threshold options> <infile = Xrio.pl-output> <outfile>"
# Options: "l" for "less or equal" ("grater or equal" is default)
#          "c" for "all hits must meet threshold in case of
#              multiple copies of the same domain in the query"
#              (default: "at least one")
# Example: "bootstrapSelector.pl 95lc OUTFILE_At_1 At_1_out" 
#
# Important. The result of this is meaningful ONLY if the thresholds 
#            for output of the RIO analysis are set to zero (L=0 R=0).
#
#
# Format for infile: 
#
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

my $VERSION              = 1.000;
my $threshold            = 0;
my $infile               = "";       
my $outfile              = "";
my $summary_outfile      = "";
my $return_line          = "";
my $identifier           = "";
my $top1                 = "";
my $analysis_performed   = 0;
my $reading              = 0;
my $i                    = 0;
my @lines                = ();
my $larger               = 1;
my $complete             = 0;
my $total                = 0;

if ( @ARGV != 3 ) {
    &errorInCommandLine();
    exit ( -1 ); 
}

$threshold  = $ARGV[ 0 ];
$infile     = $ARGV[ 1 ];
$outfile    = $ARGV[ 2 ];
$summary_outfile = $outfile.".short";

if ( -e $outfile ) {
    die "\n$0: <<$outfile>> already exists.\n";
}
unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n$0: <<$infile>> does not exist, is empty, or is not a plain textfile.\n";
}


if ( $threshold =~ /l/ ) {
    $larger = 0;
    $threshold =~ s/l//;
}
if ( $threshold =~ /c/ ) {
    $complete = 1;
    $threshold =~ s/c//;
}

open( IN, "$infile" ) || die "\n$0: Cannot open file <<$infile>>: $!\n";

open( OUT, ">$outfile" ) || die "\n$0: Cannot create file \"$outfile\": $!\n";
open( OUT_SUMMARY, ">$summary_outfile" ) || die "\n$0: Cannot create file \"$summary_outfile\": $!\n";

print OUT "bootstrapSelector.pl version: $VERSION\n\n";
print OUT "Selection of RIO analysis results with top ortholgy\n";
print OUT "bootstrap values greater or less than a threshold.\n\n";
if ( $larger == 1 ) {
    print OUT "Threshold  : Grater than or equal to $threshold\n";
}
else {
    print OUT "Threshold  : Less than or equal to $threshold\n";
}
print OUT "In case of multiple copies of the same domain in the query:\n";
if ( $complete == 1 ) {
    print OUT "All hits must meet threshold.\n";
}
else {
    print OUT "At least one hit must meet threshold.\n";
}
print OUT "Input file       : $infile\n";
print OUT "Output file      : $outfile\n";
print OUT "Output file short: $summary_outfile\n";
print OUT "Date             : ".`date`."\n\n\n";

while ( $return_line = <IN> ) {
    
    if ( $return_line =~ /^\s*# Annotation:\s*(.+)/ ) {
        $identifier = $1;
        $identifier = substr( $identifier, 0, 60); 
        $analysis_performed = 0;
        $reading = 1;
        $i = 0;
        @lines = ();
    }

    if ( $reading == 1 && $return_line =~ /^\s*RIO/ ) {
        $analysis_performed = 1;
    }
    
    if ( $reading == 1 
    && $return_line =~ /^\s*# ####################################/ ) {
        if ( $analysis_performed == 1 ) {
            &analyze();
        }
        $reading = 0;
    }

    if ( $reading == 1 ) {
        $lines[ $i++ ] = $return_line;
    }
}

close( IN );

print OUT "\n\nTotal: $total\n";

close( OUT );
close( OUT_SUMMARY );

print "\nTotal: $total\n";
print "Done.\n\n";

exit( 0 );


sub analyze {
    my $j            = 0;
    my $results      = 0;
    my $o_bootstraps = 0;
    $top1            = "";
 
    for ( $j = 0; $j < $i; $j++ ) {

        if ( $lines[ $j ] =~ /^\s*--------\s+/ ) {
            $results = 1;
        }
        elsif ( $lines[ $j ] =~ /^\s*Distance\s+values\s+/i ) {
            $results = 0;
        }
        elsif ( $results == 1
        && ( $lines[ $j ] =~ /\S+\s+\S+\s+\S+\s*$/ 
        || $lines[ $j ] =~ /^\s*!NO\s+ORTHOLOGS/ ) ) {

            if ( $lines[ $j ] =~ /^\s*!NO\s+ORTHOLOGS/ ) {
                $o_bootstraps = 0;
            }
            else {
                $lines[ $j ] =~ /(\S+)\s+\S+\s+\S+\s*$/;
                $o_bootstraps = $1;
                if ( $top1 eq "" ) {
                    $top1 = $lines[ $j ];
                    $top1 =~ s/\n//;
                    $top1 =~ s/\s{2,}/ /g;
                }
            }

            $results = 0;

            if ( $o_bootstraps > 100 || $o_bootstraps < 0 ) {
                print "o bootstraps: $o_bootstraps\n";
                die "\n\n$0: Error: Boostrap value(s) out of range.\n\n";
            }
            
            if ( $larger == 1 ) {
                if ( $complete != 1 && $o_bootstraps >= $threshold ) {
                    &writeout();
                    $total++;
                    return;
                }
                elsif ( $complete == 1 && $o_bootstraps < $threshold ) {
                    return;
                }
            }
            else {
                if ( $complete != 1 && $o_bootstraps <= $threshold ) {
                    &writeout();
                    $total++;
                    return;
                }
                elsif ( $complete == 1 && $o_bootstraps > $threshold ) {
                    return;
                }
            }
        }
    }
    if ( $complete == 1 ) {
        &writeout();
        $total++;
    }
    return;
}



sub writeout {
    my $j = 0;
    print OUT "# ############################################################################\n";
    for ( $j = 0; $j < $i; ++$j ) {
        print OUT "$lines[ $j ]";
    }
    print OUT "# ############################################################################\n\n\n";
    print OUT_SUMMARY "$identifier [top 1: $top1]\n\n";
}



sub errorInCommandLine {
    print "\n";
    print " bootstrapCounter.pl version: $VERSION\n";     
    print " Usage: \"bootstrapSelector.pl <threshold options> <infile = Xrio.pl-output> <outfile>\"\n";
    print " Options: \"l\" for \"less or equal\" (\"grater or equal\" is default)\n";
    print "          \"c\" for \"all hits must meet threshold in case of\n";
    print "          multiple copies of the same domain in the query\"\n";
    print "          (default: \"at least one\")\n";
    print " Example:\n";
    print " \"bootstrapSelector.pl 95lc OUTFILE_At_1 At_1_out\"\n\n";
    print " Important: The result of this is meaningful ONLY if the thresholds\n"; 
    print " for output of the RIO analysis are set to zero (L=0 R=0).\n\n";
    exit( -1 );
}



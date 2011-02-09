#!/usr/bin/perl -W

# mt.pl
# -----
#
# Copyright (C) 2003 Christian M. Zmasek
# All rights reserved
#
# Author: Christian M. Zmasek
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Version: 1.000
# Created on: 09/05/03
# Last modified: 09/05/03
#
#
#
# Calculates trees based on all alignments/files in a given directory using
# makeTree.pl.
#
#

use strict;
use FindBin;
use lib $FindBin::Bin;
use rio_module2;


my $PREPROCESSING_COMMAND = "";
my $PERFORM_PREPROCESSING = 0;

my $POSTPROCESSING_COMMAND = "/nfs/dm3/homedir1/czmasek/RIO1.24/perl/extractSpecies.pl";
my $PERFORM_POSTPROCESSING = 1;


my $MY_TEMP_DIR           = $TEMP_DIR_DEFAULT;   # $TEMP_DIR_DEFAULT is inherited 
                                                 # from rio_module.pm 




my $options               = ""; # Options for makeTree.pl, see makeTree.pl.
                               

my $suffix                = "";
my $use_suffixes          = 0;
my $input_dir             = "";
my $output_dir            = "";

my $i                     = 0;
my $filename              = "";
my @filenames             = ();






# Analyzes the options:
# ---------------------

unless ( @ARGV == 3 || @ARGV == 4 ) {
    &printUsage();
}

$options    = $ARGV[ 0 ]; 
$input_dir  = $ARGV[ 1 ];
$output_dir = $ARGV[ 2 ];

if ( @ARGV == 3 ) {
    $use_suffixes = 0;
}
elsif ( @ARGV == 4 ) {
    $use_suffixes = 1;
    $suffix = $ARGV[ 3 ];
}


$input_dir   = &addSlashAtEndIfNotPresent( $input_dir ); 
$output_dir  = &addSlashAtEndIfNotPresent( $output_dir );
$MY_TEMP_DIR = &addSlashAtEndIfNotPresent( $MY_TEMP_DIR );




# This adds a "-" before the options for makeTree:
# ------------------------------------------------
unless ( $options =~ /^-/ ) {
    $options = "-".$options;
}





# This creates the temp file:
# --------------------------

my $time = time;
my $ii   = 0;

my $temp_file = $MY_TEMP_DIR."mt".$time.$ii;

while ( -e $temp_file ) {
    $ii++;
    $temp_file = $MY_TEMP_DIR."mt".$time.$ii;
}



opendir( DIR, $input_dir ) || error( "Cannot open directory \"$input_dir\": $!" );

$i = 0;

while( defined( $filename = readdir( DIR ) ) ) {
    if ( $filename =~ /^\.\.?$/ ) {
        next;
    }
    if ( $use_suffixes == 1 && $filename !~ /$suffix$/ ) {
        next;
    }
   
    $filenames[ $i ] = $filename;
    $i++;
}

close( DIR );

$i = 0;

FOREACH: foreach $filename ( @filenames ) {

    # If the corresponding tree seems to already exists, do next one.
    if ( -e "$output_dir$filename.nhx" ) {
        next FOREACH;
    }

    print "\n\n\n\n";
    print "MT.PL\n";
    print "working on: $filename\n";
    
    print "[tree calculation $i]\n";
    print "=====================================================================\n\n\n";


    unlink( "$output_dir$filename.aln",
            "$output_dir$filename.log",
            "$output_dir$filename.nbd"  );

    print( "MT.PL: executing:\n" );
    
    my $inputfile = $input_dir.$filename;
    
    my $outputfilename = "";
    
    if ( $use_suffixes == 1 ) {
        $outputfilename = $output_dir . $filename;
        $outputfilename =~ s/$suffix$//;
        $outputfilename =~ s/\.$//;
        $outputfilename .= ".nhx";
    }
    else {
        $outputfilename = $output_dir . $filename . ".nhx";
    }
    
    
 
    if ( $PERFORM_PREPROCESSING == 1 ) {
         my $pre_command = "$PREPROCESSING_COMMAND";
    
         print( "$pre_command\n" );
         system( $pre_command ) && &error( "Could not execute \"$pre_command\"" );
    }
 
    $MAKETREE = "/nfs/dm3/homedir1/czmasek/RIO1.24/perl/makeTree2.pl"; # <<<<<<<<<<<<<<<<<<<<<<<-------------------~~~~~~~~~~~~~~~~~~~~~~~
    
    my $command = "$MAKETREE $options $inputfile $outputfilename";
    
    print( "$command\n" );
    system( $command ) && &error( "Could not execute \"$command\"" );
   
   
   
    if ( $PERFORM_POSTPROCESSING == 1 ) {
         my $post_command = "$POSTPROCESSING_COMMAND $outputfilename";
    
         print( "$post_command\n" );
         system( $post_command ) && &error( "Could not execute \"$post_command\"" );
    }
 
 
   
    $i++;

}



print( "\n\n\nMT.PL: Done!\n" );

exit( 0 );






sub error{

    my $text = $_[ 0 ];

    print( "\nxt.pl: ERROR:\n" );
    print( "$text\n\n" );

    exit( -1 );

}




sub printUsage {
    print "\n";
    print " mt.pl\n";
    print " _____\n";
    print " \n";
    print " Copyright (C) 2003 Christian M. Zmasek\n";
    print " All rights reserved\n";
    print "\n";
    print " Author: Christian M. Zmasek\n";
    print " zmasek\@genetics.wustl.edu\n";
    print " http://www.genetics.wustl.edu/eddy/forester/\n";
    print "\n";
    print "\n";
    print " Purpose\n";
    print " -------\n";
    print "\n";
    print " Tree construction using makeTree.pl on all alignments/files\n";
    print " in a given directory.\n"; 
    print "\n";
    print "\n";
    print " Usage\n";
    print " -----\n";
    print "\n"; 
    print "   mt.pl <options for makeTree.pl> <input directory: aligments> <output\n";
    print "         directory> [suffix for alignments to be used in input directory]\n";
    print "\n";
    print "   If a suffix is given, it will be removed for the output files.\n";
    print "\n";
    print "\n";
    print " Example\n";
    print " -------\n";
    print "\n";
    print "   \"mt.pl NS21UTRB100DX alignments/ trees/ .aln\"\n";
    print "\n";
    print "\n";
    print "\n";
    exit( -1 );

}

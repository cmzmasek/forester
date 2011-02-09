#!/usr/bin/perl -W

# pfam2pwd.pl
# -----------
# Copyright (C) 1999-2002 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek 
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Created: 05/17/01
#
# Last modified 02/20/03
#
#
#  See RIO_INSTALL on how to use this program.
#  ------------------------------------------
#

use strict;

use FindBin;
use lib $FindBin::Bin;
use rio_module;

my $VERSION = "3.002";



# =============================================================================
# =============================================================================
#
# THESE VARIABLES NEED TO BE SET BY THE USER
# ------------------------------------------
#


# Pfam alignments to calculate pairwise distances from:
# -----------------------------------------------------
my $MY_PFAM_FULL_DIRECTORY = "/path/to/Pfam/Full/"; # must end with "/"



# This file lists all the alignments for which to calculate pairwise distances
# from. If left empty, ALL the alignments in $MY_PFAM_FULL_DIRECTORY
# will be used:
# ----------------------------------------------------------------------------
my $ALGNS_TO_USE_LIST_FILE = "";



# This is _VERY IMPORTANT_. It determines the species whose sequences 
# are being used (sequences from species not listed in $MY_SPECIES_NAMES_FILE
# are ignored). Normally, one would use the same list as RIO uses
# ($SPECIES_NAMES_FILE in "rio_module.pm") -- currently "tree_of_life_bin_1-6.nhx".
#
# For certain large families (such  as protein kinases), one might use a
# species file which contains less species in order to be able to finish
# the calculations in reasonable time.
# For example, to exclude most mammals, use:
# my $MY_SPECIES_NAMES_FILE = $PATH_TO_FORESTER."data/species/tree_of_life_bin_1-6_species_list_NO_RAT_MONKEYS_APES_SHEEP_GOAT_HAMSTER"
# (to only use sequences from SWISS-PROT add this line:
# $TREMBL_ACDEOS_FILE = $PATH_TO_FORESTER."data/NO_TREMBL";)
# ----------------------------------------------------------------------------
my $MY_SPECIES_NAMES_FILE  = $SPECIES_NAMES_FILE;



# This is were the output goes (must end with "/")
# ------------------------------------------------
my $MY_RIO_PWD_DIRECTORY   = "/path/to/pfam2pwd_out/pwd/";
my $MY_RIO_BSP_DIRECTORY   = "/path/to/pfam2pwd_out/bsp/";
my $MY_RIO_NBD_DIRECTORY   = "/path/to/pfam2pwd_out/nbd/";
my $MY_RIO_ALN_DIRECTORY   = "/path/to/pfam2pwd_out/aln/";
my $MY_RIO_HMM_DIRECTORY   = "/path/to/pfam2pwd_out/hmm/";



# A directory to create temporary files in:
# -----------------------------------------
my $MY_TEMP_DIR            = "/tmp/"; # must end with "/"



# Alignments in which the number of sequences after pruning (determined
# by "$MY_SPECIES_NAMES_FILE") is lower than this, are ignored
# (no calculation of pwds):
# ------------------------------------------------------------------
my $MIN_SEQS  = 5;



# Alignments in which the number of sequences after pruning (determined
# by "$MY_SPECIES_NAMES_FILE") is greater than this, are ignored
# (no calculation of pwds):
# ------------------------------------------------------------------
my $MAX_SEQS  = 700;



# Seed for the random number generator for bootstrapping (must be 4n+1):
# ---------------------------------------------------------------------
my $MY_SEED   = 85; 



# This is used to choose the model to be used for the (ML)
# distance calculation:
# IMPORTANT: "$MY_MATRIX_FOR_PWD" in "rio_module.pm" needs to
# have the same value, when the pwds calculated are going to
# be used for RIO!
# 0 = JTT
# 2 = BLOSUM 62
# 3 = mtREV24
# 5 = VT
# 6 = WAG
# PAM otherwise
# --------------------------------------------------------
my $MY_MATRIX = 2; 



#
# End of variables which need to be set by the user.
#
# =============================================================================
# =============================================================================








my $too_small             = 0;
my $too_large             = 0;
my $i                     = 0;
my $seqs                  = 0;
my $filename              = "";
my $tmp_dir               = "";
my $current_dir           = "";
my $return_line           = "";
my @filenames             = ();
my @too_small_names       = ();
my @too_large_names       = ();
my %Species_names_hash    = ();
my %AC_OS                 = (); # AC -> species name
my %AC_DE                 = (); # AC -> description
my %ALGNS_TO_USE          = (); # name of alignment -> ""
my $use_algns_to_use_list = 0;
my $LOGFILE               = "00_pfam2pwd_LOGFILE";
   $HMMBUILD              = $HMMBUILD." --amino";


&createTempdir();


&startLogfile();


opendir( DIR, $MY_PFAM_FULL_DIRECTORY ) || die "\n\n$0: Cannot open directory $MY_PFAM_FULL_DIRECTORY: $!\n\n";
$i = 0;
while( defined( $filename = readdir( DIR ) ) ) {
    if ( $filename =~ /^\.\.?$/ ) {
        next;
    }
    $filenames[ $i ] = $filename;
    $i++;
}
close( DIR );


&readSpeciesNamesFile( $MY_SPECIES_NAMES_FILE );

&readTrEMBL_ACDEOS_FILE();

if ( defined( $ALGNS_TO_USE_LIST_FILE ) && $ALGNS_TO_USE_LIST_FILE =~ /\w/ ) {
    $use_algns_to_use_list = 1;
    &readListFile();
}


$current_dir = `pwd`;
$current_dir =~ s/\s//;
chdir ( $tmp_dir ) 
|| die "\n\n$0: Unexpected error: Could not chdir to <<$tmp_dir>>: $!";

$i = 0;

FOREACH_ALIGN: foreach $filename ( @filenames ) {

    # If the corresponding pwd, positions, and aln files seem to already exists, do next one.
    if ( ( -e $MY_RIO_PWD_DIRECTORY.$filename.$SUFFIX_PWD ) 
    &&   ( -e $MY_RIO_BSP_DIRECTORY.$filename.$SUFFIX_BOOT_STRP_POS ) 
    &&   ( -e $MY_RIO_NBD_DIRECTORY.$filename.$SUFFIX_PWD_NOT_BOOTS )
    &&   ( -e $MY_RIO_ALN_DIRECTORY.$filename.$ALIGN_FILE_SUFFIX )
    &&   ( -e $MY_RIO_HMM_DIRECTORY.$filename.$SUFFIX_HMM ) ) {
        next FOREACH_ALIGN;
    }

    if ( $use_algns_to_use_list == 1 && !exists( $ALGNS_TO_USE{ $filename } ) ) {
        next FOREACH_ALIGN;
    }


    $seqs = &removeSeqsFromPfamAlign( $MY_PFAM_FULL_DIRECTORY.$filename,
                                      "REM_SEQ_OUTFILE",
                                      1 );
    if ( $seqs < $MIN_SEQS ) { 
        unlink( "REM_SEQ_OUTFILE" );
        $too_small_names[ $too_small++ ] = $filename;
        next FOREACH_ALIGN;
    }
    elsif ( $seqs > $MAX_SEQS ) {
        unlink( "REM_SEQ_OUTFILE" );
        $too_large_names [ $too_large++ ] = $filename;
        next FOREACH_ALIGN;
    }


    print "\n\n\n";
    print " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    print " $i: $filename ($seqs seqs)\n";
    print " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    print "\n";

    # If one of the two file exists from a previous (interrupted) run.
    unlink( $MY_RIO_PWD_DIRECTORY.$filename.$SUFFIX_PWD );
    unlink( $MY_RIO_BSP_DIRECTORY.$filename.$SUFFIX_BOOT_STRP_POS );
    unlink( $MY_RIO_NBD_DIRECTORY.$filename.$SUFFIX_PWD_NOT_BOOTS );
    unlink( $MY_RIO_ALN_DIRECTORY.$filename.$ALIGN_FILE_SUFFIX );
    unlink( $MY_RIO_HMM_DIRECTORY.$filename.$SUFFIX_HMM );
    
    
    &executeHmmbuild( "REM_SEQ_OUTFILE",
                      $MY_RIO_ALN_DIRECTORY.$filename.$ALIGN_FILE_SUFFIX,
                      "hmm" );

    if ( unlink( "hmm" ) != 1 ) {
        die "\n\n$0: Unexpected error: Could not delete <<hmm>>: $!";
    }

    if ( unlink( "REM_SEQ_OUTFILE" ) != 1 ) {
        die "\n\n$0: Unexpected error: Could not delete <<REM_SEQ_OUTFILE>>: $!";
    }
    
    executeHmmbuildHand( $MY_RIO_ALN_DIRECTORY.$filename.$ALIGN_FILE_SUFFIX,
                         $MY_RIO_HMM_DIRECTORY.$filename.$SUFFIX_HMM );

    system( $HMMCALIBRATE, $MY_RIO_HMM_DIRECTORY.$filename.$SUFFIX_HMM )
    && die "\n\n$0: Could not execute \"$HMMCALIBRATE $MY_RIO_HMM_DIRECTORY.$filename.$SUFFIX_HMM\": $!";

    &pfam2phylipMatchOnly( $MY_RIO_ALN_DIRECTORY.$filename.$ALIGN_FILE_SUFFIX, "infile" );
    
    &executePuzzle( "infile", $MY_MATRIX );
   
    system( "mv", "infile.dist", $MY_RIO_NBD_DIRECTORY.$filename.$SUFFIX_PWD_NOT_BOOTS ) 
    && die "\n\n$0: Unexpected error: $!";

    &executeBootstrap( "infile",
                       $BOOTSTRAPS,
                       "BOOTSTRAPPED_ALGN",
                       $MY_RIO_BSP_DIRECTORY.$filename.$SUFFIX_BOOT_STRP_POS,
                       $MY_SEED );
    
    if ( unlink( "infile" ) != 1 ) {
        die "\n\n$0: Unexpected error: Could not delete <<infile>>: $!";
    }

    
    &executePuzzleBootstrapped( "BOOTSTRAPPED_ALGN", $MY_MATRIX );

    ##if ( unlink( "outfile" ) != 1 ) {
    ##    die "\n\n$0: Unexpected error: Could not delete <<outfile>>: $!";
    ##}


    system( "mv", "BOOTSTRAPPED_ALGN".".dist", $MY_RIO_PWD_DIRECTORY.$filename.$SUFFIX_PWD )
    && die "\n\n$0: Unexpected error: $!\n\n";

    if ( unlink( "BOOTSTRAPPED_ALGN" ) != 1 ) {
        die "\n\n$0: Unexpected error: Could not delete <<BOOTSTRAPPED_ALGN>>: $!";
    }
    
    $i++;

} ## End of FOREACH_ALIGN loop.


chdir( $current_dir ) 
|| die "\n\n$0: Unexpected error: Could not chdir to <<$current_dir>>: $!";

rmdir( $tmp_dir );

&finishLogfile();

print "\n\n\n";
print( "pfam2pwd.pl: Done.\n" );
print( "Successfully calculated $i pairwise distance files.\n" );
print( "Too large alignments (>$MAX_SEQS): $too_large\n" );
print( "Too small alignments (<$MIN_SEQS): $too_small\n" );
print( "See the logfile \"$MY_RIO_PWD_DIRECTORY".$LOGFILE."\"\n" );
print "\n\n\n";

exit( 0 );






# Methods
# -------



# Three arguments:
# 1. Stockholm alignment
# 2. Outalignment
# 3. Outhmm
# Returns the options used.
# Last modified: 06/26/01
sub executeHmmbuild {

    my $full         = $_[ 0 ];
    my $outalignment = $_[ 1 ];
    my $outhmm       = $_[ 2 ];
    my $options      = "";

    unless ( ( -s $full ) && ( -f $full ) && ( -T $full ) ) {
        die "\n\n$0: \"$full\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
    
    $options = getHmmbuildOptionsFromPfam( $full );

    $options =~ s/-f//;
    $options =~ s/-g//;
    $options =~ s/-s//;
    $options =~ s/-F//;
    $options =~ s/-A//;
    $options =~ s/-o\s+\S+//;
    $options =~ s/(\s|^)[^-]\S+/ /g;

    if ( $options =~ /--prior/ ) {
        my $basename = basename( $full );
        $basename .= ".PRIOR";
        $options =~ s/--prior/--prior $PRIOR_FILE_DIR$basename/;
    }

    # Remove for versions of HMMER lower than 2.2.
    if ( $options =~ /--informat\s+\S+/ ) {
        $options =~ s/--informat\s+\S+/-/; 
    }
    
    system( "$HMMBUILD $options -o $outalignment $outhmm $full" )
    && die "\n\n$0: Could not execute \"$HMMBUILD $options -o $outalignment $outhmm $full\".\n\n";
 
    return $options;

} ## executeHmmbuild.


# Two arguments:
# 1. Stockholm alignment
# 2. Outhmm
# Returns the options used.
# Last modified: 06/26/01
sub executeHmmbuildHand {

    my $full         = $_[ 0 ];
    my $outhmm       = $_[ 1 ];
    my $options      = "";

    unless ( ( -s $full ) && ( -f $full ) && ( -T $full ) ) {
        die "\n\n$0: \"$full\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
    
    $options = getHmmbuildOptionsFromPfam( $full );

    $options =~ s/-f//;
    $options =~ s/-g//;
    $options =~ s/-s//;
    $options =~ s/-F//;
    $options =~ s/-A//;
    $options =~ s/-o\s+\S+//;
    $options =~ s/(\s|^)[^-]\S+/ /g;

    if ( $options =~ /--prior/ ) {
        my $basename = basename( $full );
        $basename .= ".PRIOR";
        $options =~ s/--prior/--prior $PRIOR_FILE_DIR$basename/;
    }

    # Remove for versions of HMMER lower than 2.2.
    if ( $options =~ /--informat\s+\S+/ ) {
        $options =~ s/--informat\s+\S+/-/; 
    }
    
    system( "$HMMBUILD --hand $options $outhmm $full" )
    && die "\n\n$0: Could not execute \"$HMMBUILD -- hand $options $outhmm $full\".\n\n";
 
    return $options;

} ## executeHmmbuildHand.



# One argument:
# Pfam align name.
# Last modified: 02/26/01
sub getHmmbuildOptionsFromPfam {

    my $infile      = $_[ 0 ];
    my $return_line = "";
    my $result      = "";

    unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
        die "\n\n$0: \"$infile\" does not exist, is empty, or is not a plain textfile.\n\n";
    }

    open( GHO, $infile ) || die "\n\n$0: Unexpected error: Cannot open file <<$infile>>: $!";
    while ( $return_line = <GHO> ) {
        if ( $return_line =~ /^\s*#.*hmmbuild\s+(.+)\s*$/ ) {
            $result = $1;
            close( GHO );
            return $result;
        }
    }
    close( GHO );
    return $result;

} ## getHmmbuildOptionsFromPfam



# Similar to the method with the same name in "rio.pl".
# Removes sequences from a Pfam flat file.
# Adds species to TrEMBL seqs.
# It can remove all sequences not from species listed in a species names file.
# It can remove all sequences which do not have a SWISS-PROT name (XXXX_XXXXX)
# Three arguments:
# 1. Pfam flat file name
# 2. outfile name
# 3. 1 to remove TrEMBL seqs with "(FRAGMENT)" in their DE line.
# Returns the number of sequences in the resulting alignment.
# If a query name is given, it returns -1 if query is not found in alignment,
# -10 if the name is not unique.
# Last modified: 05/24/02
sub removeSeqsFromPfamAlign {
    my $infile                   = $_[ 0 ];
    my $outfile                  = $_[ 1 ];
    my $remove_frags             = $_[ 2 ];
    my $return_line              = "";
    my $saw_sequence_line        = 0;
    my $number_of_seqs           = 0;
    my $DE                       = "";
    my $OS                       = "";
    my $AC                       = "";
    my $i                        = 0;
    my $length                   = 0;
    my $seq_name                 = "";
    my $seq                      = "";
   
    
    open( OUT_RNSP, ">$outfile" ) || die "\n\n$0: Unexpected error: Cannot create file \"$outfile\": $!";
    open( IN_RNSP, "$infile" ) || die "\n\n$0: Unexpected error: Cannot open file <<$infile>>: $!";
    while ( $return_line = <IN_RNSP> ) {

        if ( $saw_sequence_line == 1
        && !&containsPfamNamedSequence( $return_line )
        && !&isPfamCommentLine( $return_line ) ) {
            # This is just for counting purposes.
            $saw_sequence_line = 2;
        }
        if ( &isPfamSequenceLine( $return_line ) ) { 
            if ( $saw_sequence_line == 0 ) {
                $saw_sequence_line = 1;
            }
            $return_line =~ /^\s*(\S+)\s+(\S+)/;
            $seq_name = $1;
            $seq      = $2;
            if ( !&startsWithSWISS_PROTname( $return_line ) ) {
                $seq_name =~ /^(\S+)\//;
                $AC = $1;
                unless( exists( $AC_OS{ $AC } ) ) {
                    #ACs not present in "ACDEOS" file.
                    next;
                }
                $OS = $AC_OS{ $AC };
                if ( !$OS || $OS eq "" ) {
                    die "\n\n$0: Unexpected error: species for \"$AC\" not found.\n\n";
                }   
                unless( exists( $Species_names_hash{ $OS } ) ) {
                    next;
                }
                if ( $remove_frags == 1 ) {
                    $DE = $AC_DE{ $AC };
                    if ( $DE && $DE =~ /\(FRAGMENT\)/ ) {
                        next;
                    }
                }
                $seq_name =~ s/\//_$OS\//;
            }
            else {
                if ( $return_line =~ /_([A-Z0-9]{1,5})\// ) {
                    unless( exists( $Species_names_hash{ $1 } ) ) {
                        next;
                    }
                }
                # remove everything whose species cannot be determined.
                else {
                    next;
                }
            }
            $length = length( $seq_name );
            for ( $i = 0; $i <= ( $LENGTH_OF_NAME - $length - 1 ); $i++ ) {
	            $seq_name .= " ";
            }
            $return_line = $seq_name.$seq."\n";
        }

        if ( !&isPfamCommentLine( $return_line ) ) {
            print OUT_RNSP $return_line;
        }

        if ( $saw_sequence_line == 1 ) {
            $number_of_seqs++;
        }
    } ## while ( $return_line = <IN_RNSP> )
    close( IN_RNSP );
    close( OUT_RNSP );
    
    return $number_of_seqs;

} ## removeSeqsFromPfamAlign







# Reads in (SWISS-PROT) species names from a file.
# Names must be separated by newlines.
# Lines beginning with "#" are ignored.
# A possible "=" and everything after is ignored.
# One argument: species-names-file name
# Last modified: 04/24/01
sub readSpeciesNamesFile {
    my $infile      = $_[ 0 ];
    my $return_line = "";
    my $species     = "";

    unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
        die "\n\n$0: Error: \"$infile\" does not exist, is empty, or is not a plain textfile.\n\n";
    }

    open( IN_RSNF, "$infile" ) || die "\n\n$0: Unexpected error: Cannot open file <<$infile>>: $!\n\n";
    while ( $return_line = <IN_RSNF> ) {
        if ( $return_line !~ /^\s*#/ && $return_line =~ /(\S+)/ ) {
            $species = $1;
            $species =~ s/=.+//;
            $Species_names_hash{ $species } = "";
        }
    }
    close( IN_RSNF );

    return;
} ## readSpeciesNamesFile



# Last modified: 02/21/03
sub readTrEMBL_ACDEOS_FILE {

    my $return_line = ""; 

    unless ( ( -s $TREMBL_ACDEOS_FILE  ) && ( -f $TREMBL_ACDEOS_FILE  ) && ( -T $TREMBL_ACDEOS_FILE ) ) {
        die "\n\n$0: Error: \"$TREMBL_ACDEOS_FILE\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
    # Fill up (huge) hashs.
    open( HH, "$TREMBL_ACDEOS_FILE" ) || die "\n\n$0: Unexpected error: Cannot open file <<$TREMBL_ACDEOS_FILE>>: $!\n\n";
    while ( $return_line = <HH> ) {
       
        if ( $return_line =~ /(\S+);([^;]*);(\S+)/ ) {
             $AC_OS{ $1 } = $3;
             $AC_DE{ $1 } = $2;
        }
    }
    close( HH ); 
} ## readTrEMBL_ACDEOS_FILE


# Last modified: 02/21/03
sub readListFile {

    my $return_line = "";
    
    unless ( ( -s $ALGNS_TO_USE_LIST_FILE ) && ( -f $ALGNS_TO_USE_LIST_FILE ) && ( -T $ALGNS_TO_USE_LIST_FILE ) ) {
        die "\n\n$0: Error: \"$ALGNS_TO_USE_LIST_FILE\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
    # Fill up hash.
    open( LF, "$ALGNS_TO_USE_LIST_FILE" ) || die "\n\n$0: Unexpected error: Cannot open file <<$ALGNS_TO_USE_LIST_FILE>>: $!\n\n";
    while ( $return_line = <LF> ) {
        if ( $return_line =~ /^\s*(\S+)\s*$/ ) {  # just a list
            $ALGNS_TO_USE{ $1 } = "";
        }
        elsif ( $return_line =~ /^\s*\S+\s+\S+\s+(\S+)/ ) { # "changes" list from Pfam
            $ALGNS_TO_USE{ $1 } = "";
        }
       
    }
    close( LF ); 

} ## readListFile



# Five arguments:
# 1. Name of inputfile
# 2. Bootstraps
# 2. Name of output alignment file
# 3. Name of output positions file
# 4. Seed for random number generator
#
# Last modified: 06/23/01
sub executeBootstrap {
    my $infile     = $_[ 0 ];
    my $bootstraps = $_[ 1 ];
    my $outalign   = $_[ 2 ];
    my $positions  = $_[ 3 ];
    my $seed       = $_[ 4 ];

    system( "$BOOTSTRAP_CZ_PL 0 $bootstraps $infile $outalign $positions $seed" )
    && die "\n\n$0: executeBootstrap:\nCould not execute \"$BOOTSTRAP_CZ_PL 0 $bootstraps $infile $outalign $positions $seed\".\n\n";

} ## executeBootstrap




# Last modified: 05/22/02
sub createTempdir {

    my $ii = 0;
    my $time = time;

    $tmp_dir = $MY_TEMP_DIR.$time.$ii;
   
    while ( -e $tmp_dir ) {
        $ii++;
        $tmp_dir = $MY_TEMP_DIR.$time.$ii;
    }

    mkdir(  $tmp_dir, 0777 )
    || die "\n\n$0: Unexpected error: Could not create <<$tmp_dir>>: $!\n\n";

    unless ( ( -e $tmp_dir ) && ( -d $tmp_dir ) ) {
        die "\n\n$0: Unexpected error: failed to create <<$tmp_dir>>.\n\n";
    }

} ## createTempdir



# Last modified: 05/17/01
sub startLogfile {
    if ( -e $MY_RIO_PWD_DIRECTORY.$LOGFILE ) {
        print "\npfam2pwd.pl:\n";
        print "logfile $MY_RIO_PWD_DIRECTORY"."$LOGFILE already exists\n";
        print "rename it or place it in another directory\n";
        exit( -1 ); 
    } 

    open( L, ">$MY_RIO_PWD_DIRECTORY".$LOGFILE )
    || die "\n\n$0: startLogfile: Cannot create logfile: $!\n\n";
    print L "Min seqs           : $MIN_SEQS\n";
    print L "Max seqs           : $MAX_SEQS\n";
    print L "Seed               : $MY_SEED\n";
    print L "TrEMBL ACDEOS file : $TREMBL_ACDEOS_FILE\n";
    print L "Species names file : $MY_SPECIES_NAMES_FILE\n";
    print L "Pfam directory     : $MY_PFAM_FULL_DIRECTORY\n";
    print L "PWD outputdirectory: $MY_RIO_PWD_DIRECTORY\n";
    print L "BSP outputdirectory: $MY_RIO_BSP_DIRECTORY\n";
    print L "NBD outputdirectory: $MY_RIO_NBD_DIRECTORY\n";
    print L "ALN outputdirectory: $MY_RIO_ALN_DIRECTORY\n";
    print L "HMM outputdirectory: $MY_RIO_HMM_DIRECTORY\n";
    print L "Start date         : ".`date`;
    if ( $MY_MATRIX == 0 ) {
        print L "Matrix             : JTT\n";
    }
    elsif ( $MY_MATRIX == 2 ) {
        print L "Matrix             : BLOSUM 62\n";
    }
    elsif ( $MY_MATRIX == 3 ) {
        print L "Matrix             : mtREV24\n";
    }
    elsif ( $MY_MATRIX == 5 ) {
        print L "Matrix             : VT\n";
    }
    elsif ( $MY_MATRIX == 6 ) {
        print L "Matrix             : WAG\n";
    }
    elsif ( $MY_MATRIX == 7 ) {
        print L "Matrix             : auto\n";
    } 
    else {
        print L "Matrix             : PAM\n";
    } 
} ## startLogfile



# Last modified: 05/17/01
sub finishLogfile {
    my $j = 0;
    print L "\n\n";
    print L "Successfully calculated $i pairwise distance files.\n";
    print L "Too large alignments (>$MAX_SEQS): $too_large\n";
    print L "Too small alignments (<$MIN_SEQS): $too_small\n";
    print L "Finish date        : ".`date`."\n\n";
    
    print L "List of the $too_large alignments which were ignored because they\n";
    print L "contained too many sequences (>$MAX_SEQS) after pruning:\n";
    for ( $j = 0; $j < $too_large; ++$j ) { 
        print L "$too_large_names[ $j ]\n";
    } 
    print L "\n\n";
    print L "List of the $too_small alignments which were ignored because they\n";
    print L "contained not enough sequences (<$MIN_SEQS) after pruning:\n";
    for ( $j = 0; $j < $too_small; ++$j ) { 
        print L "$too_small_names[ $j ]\n";
    } 
    print L "\n";
    close( L );
} ## finishLogfile





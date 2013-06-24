# $Id: forester.pm,v 1.26 2010/12/13 19:00:22 cmzmasek Exp $
#
# FORESTER -- software libraries and applications
# for evolutionary biology research and applications.
#
# Copyright (C) 2007-2009 Christian M. Zmasek
# Copyright (C) 2007-2009 Burnham Institute for Medical Research
# All rights reserved
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
#
# Contact: phylosoft @ gmail . com
#     WWW: www.phylosoft.org/forester
#
#
#


package forester;
use strict;
require Exporter;

our $VERSION = 1.000;

our @ISA    = qw( Exporter );

our @EXPORT = qw( executeConsense
                  executePhyloPl
                  executePuzzleDQO
                  executePuzzleDQObootstrapped
                  pfam2phylipMatchOnly
                  startsWithSWISS_PROTname
                  isPfamSequenceLine
                  isPfamCommentLine
                  containsPfamNamedSequence
                  isRFline
                  executeProtpars
                  setModelForPuzzle
                  setRateHeterogeneityOptionForPuzzle
                  setParameterEstimatesOptionForPuzzle
                  executePuzzleBootstrapped
                  executePuzzle
                  executeFastme
                  executeNeighbor
                  executeFitch
                  executeBionj
                  executeWeighbor
                  executePhyml
                  executeHmmfetch
                  addDistsToQueryToPWDfile
                  testForTextFilePresence
                  exitWithWarning
                  dieWithUnexpectedError
                  addSlashAtEndIfNotPresent
                  $LENGTH_OF_NAME
                  $MIN_NUMBER_OF_AA
                  $TREMBL_ACDEOS_FILE
                  $SWISSPROT_ACDEOS_FILE
                  $SPECIES_NAMES_FILE
                  $SPECIES_TREE_FILE_DEFAULT
                  $MULTIPLE_TREES_FILE_SUFFIX
                  $LOG_FILE_SUFFIX
                  $ALIGN_FILE_SUFFIX
                  $TREE_FILE_SUFFIX
                  $ADDITION_FOR_RIO_ANNOT_TREE
                  $SUFFIX_PWD
                  $SUFFIX_BOOT_STRP_POS
                  $MULTIPLE_PWD_FILE_SUFFIX
                  $SUFFIX_PWD_NOT_BOOTS
                  $SUFFIX_HMM
                  $MATRIX_FOR_PWD 
                  $RIO_PWD_DIRECTORY
                  $RIO_BSP_DIRECTORY
                  $RIO_NBD_DIRECTORY
                  $RIO_ALN_DIRECTORY
                  $RIO_HMM_DIRECTORY
                  $PFAM_FULL_DIRECTORY 
                  $PFAM_SEED_DIRECTORY 
                  $PRIOR_FILE_DIR
                  $PFAM_HMM_DB
                  $FORESTER_JAR
                  $SEQBOOT
                  $NEIGHBOR
                  $PROTPARS
                  $CONSENSE
                  $PROML
                  $PHYLIP_VERSION
                  $PUZZLE
                  $PUZZLE_VERSION
                  $FASTME
                  $FASTME_VERSION
                  $BIONJ
                  $BIONJ_VERSION
                  $WEIGHBOR
                  $WEIGHBOR_VERSION
                  $RAXML
                  $RAXML_VERSION
                  $PHYML
                  $PHYML_VERSION
                  $HMMALIGN
                  $HMMSEARCH
                  $HMMBUILD
                  $HMMFETCH
                  $SFE
                  $HMMCALIBRATE
                  $P7EXTRACT 
                  $MULTIFETCH 
                  $BOOTSTRAP_CZ
                  $BOOTSTRAP_CZ_PL
                  $SUPPORT_TRANSFER
                  $SUPPORT_STATISTICS
                  $NEWICK_TO_PHYLOXML
                  $PHYLO_PL
                  $RIO_PL
                  $DORIO
                  $PUZZLE_DQO
                  $BOOTSTRAPS
                  $PATH_TO_FORESTER
                  $JAVA
                  $NODE_LIST
                  $RIO_SLAVE_DRIVER
                  $RIO_SLAVE
                  $TEMP_DIR_DEFAULT
                  $EXPASY_SPROT_SEARCH_DE
                  $EXPASY_SPROT_SEARCH_AC    
 );




# =============================================================================
# =============================================================================
#
# THESE VARIABLES ARE ENVIRONMENT DEPENDENT, AND NEED TO BE SET ACCORDINGLY
# BY THE USER
# -------------------------------------------------------------------------
#

# For using just "phylo_pl.pl", only the following variables need to be set 
# $JAVA
# $FORESTER_JAR
# $TEMP_DIR_DEFAULT
# $SEQBOOT
# $CONSENSE 
# $PUZZLE
# $FASTME
# $NEIGHBOR
# $FITCH
# $BIONJ
# $WEIGHBOR
# $PHYML
# $PROTPARS

# Software directory:
# ---------------------

our $SOFTWARE_DIR              = "/home/czmasek/SOFTWARE/";


# Java virtual machine:
# ---------------------
our $JAVA                      = $SOFTWARE_DIR."JAVA/jdk1.6.0_03/bin/java";


# Where all the temporary files can be created:
# ---------------------------------------------
our $TEMP_DIR_DEFAULT          = "/tmp/";


# Programs from Joe Felsenstein's PHYLIP package:
# -----------------------------------------------
our $SEQBOOT                   = $SOFTWARE_DIR."PHYLIP/phylip-3.68/src/seqboot";
our $NEIGHBOR                  = $SOFTWARE_DIR."PHYLIP/phylip-3.68/src/neighbor";
our $PROTPARS                  = $SOFTWARE_DIR."PHYLIP/phylip-3.68/src/protpars";
our $PROML                     = $SOFTWARE_DIR."PHYLIP/phylip-3.68/src/proml";
our $FITCH                     = $SOFTWARE_DIR."PHYLIP/phylip-3.68/src/fitch";
our $CONSENSE                  = $SOFTWARE_DIR."PHYLIP/phylip-3.68/src/consense";
our $PHYLIP_VERSION            = "3.68";

# TREE-PUZZLE:
# ------------
our $PUZZLE                    = $SOFTWARE_DIR."TREE_PUZZLE/tree-puzzle-5.2/src/puzzle";
our $PUZZLE_VERSION            = "5.2";

# FASTME:
# -----------------------------------------------------
our $FASTME                    = $SOFTWARE_DIR."FASTME/fastme2.0/fastme";
our $FASTME_VERSION            = "2.0";

# BIONJ:
# -----------------------------------------------------
our $BIONJ                    = $SOFTWARE_DIR."BIONJ/bionj";
our $BIONJ_VERSION            = "[1997]";

# WEIGHBOR:
# -----------------------------------------------------
our $WEIGHBOR                 = $SOFTWARE_DIR."WEIGHBOR/Weighbor/weighbor";
our $WEIGHBOR_VERSION         = "1.2.1";

# PHYML:
# -----------------------------------------------------
our $PHYML                    = $SOFTWARE_DIR."PHYML/phyml_v2.4.4/exe/phyml_linux";
our $PHYML_VERSION            = "2.4.4";

# RAXML:
# -----------------------------------------------------
our $RAXML                    = $SOFTWARE_DIR."RAXML/RAxML-7.0.4/raxmlHPC";
our $RAXML_VERSION            = "7.0.4";


# forester.jar. This jar file is currently available at: http://www.phylosoft.org 
# -------------------------------------------------------------------------------

our $FORESTER_JAR             = $SOFTWARE_DIR."FORESTER/DEV/forester/forester/java/forester.jar";



# End of variables which need to be set by the user for using "phylo_pl.pl".














# Tool from forester.jar to transfer support values:
# -------------------------------------------------
our $SUPPORT_TRANSFER          = $JAVA." -cp $FORESTER_JAR org.forester.application.support_transfer";



# Tool from forester.jar for simple statistics for support values:
# ----------------------------------------------------------------
our $SUPPORT_STATISTICS          = $JAVA." -cp $FORESTER_JAR org.forester.application.support_statistics";


# Tool from forester.jar to transfer nh to phyloXML:
# -------------------------------------------------
our $NEWICK_TO_PHYLOXML          = $JAVA." -cp $FORESTER_JAR org.forester.application.phyloxml_converter";



# FORESTER itself (currently not needed for "phylo_pl.pl"):
# ---------------------------------------------------------
our $PATH_TO_FORESTER          = ""; 


# Pfam data (not needed for phylo_pl.pl):
# --------------------------------------
our $PFAM_FULL_DIRECTORY       = "/path/to/Pfam/Full/"; 
our $PFAM_SEED_DIRECTORY       = "/path/to/Pfam/Seed/";
our $PFAM_HMM_DB               = "/path/to/Pfam/Pfam_ls"; # Need to run "hmmindex" on this
                                                          # to produce .ssi file.
                                                          # Then, for example
                                                          # "setenv HMMERDB /home/rio/pfam-6.6/"


$PATH_TO_FORESTER = &addSlashAtEndIfNotPresent( $PATH_TO_FORESTER );


# Description lines and species from SWISS-PROT and TrEMBL (not needed for phylo_pl.pl):
# -------------------------------------------------------------------------------------
our $TREMBL_ACDEOS_FILE        = $PATH_TO_FORESTER."data/trembl22_ACDEOS_1-6";
                                 
our $SWISSPROT_ACDEOS_FILE     = $PATH_TO_FORESTER."data/sp40_ACDEOS_1-6";



# Names of species which can be analyzed and analyzed 
# against (must also be in tree $SPECIES_TREE_FILE_DEFAULT).
# By using a list with less species, RIO analyses become faster
# but lose phylogenetic resolution. 
# For many purposes, list "tree_of_life_bin_1-6_species_list"
# in "data/species/" might be sufficient:
# (not needed for phylo_pl.pl)
# --------------------------------------------------------------
our $SPECIES_NAMES_FILE        = $PATH_TO_FORESTER."data/species/tree_of_life_bin_1-6_species_list";



# A default species tree in NHX format.
# For many purposes, tree "tree_of_life_bin_1-6.nhx"
# in "data/species/" might be fine:
# (not needed for phylo_pl.pl)
# --------------------------------------------------
our $SPECIES_TREE_FILE_DEFAULT = $PATH_TO_FORESTER."data/species/tree_of_life_bin_1-6.nhx";



# Data for using precalculated distances:
# (not needed for phylo_pl.pl)
# ---------------------------------------
our $MATRIX_FOR_PWD            = 2;  # The matrix which has been used for the pwd in $RIO_PWD_DIRECTORY.
                                     # 0=JTT, 1=PAM, 2=BLOSUM 62, 3=mtREV24, 5=VT, 6=WAG.
           
our $RIO_PWD_DIRECTORY         = $PATH_TO_FORESTER."example_data/";  # all must end with "/"
our $RIO_BSP_DIRECTORY         = $PATH_TO_FORESTER."example_data/";
our $RIO_NBD_DIRECTORY         = $PATH_TO_FORESTER."example_data/";
our $RIO_ALN_DIRECTORY         = $PATH_TO_FORESTER."example_data/";
our $RIO_HMM_DIRECTORY         = $PATH_TO_FORESTER."example_data/";



#
# End of variables which need to be set by the user.
#
# =============================================================================
# =============================================================================





$TEMP_DIR_DEFAULT    = &addSlashAtEndIfNotPresent( $TEMP_DIR_DEFAULT ); 
$PFAM_FULL_DIRECTORY = &addSlashAtEndIfNotPresent( $PFAM_FULL_DIRECTORY );
$PFAM_SEED_DIRECTORY = &addSlashAtEndIfNotPresent( $PFAM_SEED_DIRECTORY );



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These variables should normally not be changed:
#

our $PRIOR_FILE_DIR            = $PATH_TO_FORESTER."data/priors_for_hmmbuild/"; 
                                 # Directory containing dirichlet prior
                                 # files needed for certain aligments     
                                 # by hmmbuild (e.g. Collagen).





# TREE-PUZZLE:
our $PUZZLE_DQO                = $PATH_TO_FORESTER."puzzle_dqo/src/puzzle";

# HMMER:
our $HMMALIGN                  = $PATH_TO_FORESTER."hmmer/binaries/hmmalign";
our $HMMSEARCH                 = $PATH_TO_FORESTER."hmmer/binaries/hmmsearch";
our $HMMBUILD                  = $PATH_TO_FORESTER."hmmer/binaries/hmmbuild";
our $HMMFETCH                  = $PATH_TO_FORESTER."hmmer/binaries/hmmfetch";
our $SFE                       = $PATH_TO_FORESTER."hmmer/binaries/sfetch";
our $HMMCALIBRATE              = $PATH_TO_FORESTER."hmmer/binaries/hmmcalibrate";

our $P7EXTRACT                 = $PATH_TO_FORESTER."perl/p7extract.pl";
our $MULTIFETCH                = $PATH_TO_FORESTER."perl/multifetch.pl";


# RIO/FORESTER:
our $BOOTSTRAP_CZ              = $PATH_TO_FORESTER."C/bootstrap_cz";
our $BOOTSTRAP_CZ_PL           = $PATH_TO_FORESTER."perl/bootstrap_cz.pl";
#our $SUPPORT_TRANSFER         = $JAVA." -cp $PATH_TO_FORESTER"."java forester.tools.transfersBranchLenghts";
#our $SUPPORT_TRANSFER         = $JAVA." -cp /home/czmasek/SOFTWARE/FORESTER/forester3/forester.jar org.forester.tools.SupportTransfer";

our $PHYLO_PL                  = $PATH_TO_FORESTER."perl/phylo_pl.pl";
our $RIO_PL                    = $PATH_TO_FORESTER."perl/rio.pl";
our $DORIO                     = $JAVA." -cp $PATH_TO_FORESTER"."java forester.tools.DoRIO";
# parallel RIO:
our $RIO_SLAVE_DRIVER          = $PATH_TO_FORESTER."perl/rio_slave_driver.pl";
our $RIO_SLAVE                 = $PATH_TO_FORESTER."perl/rio_slave.pl";
our $NODE_LIST                 = $PATH_TO_FORESTER."data/node_list.dat";

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


our $BOOTSTRAPS         = 100;
our $MIN_NUMBER_OF_AA   = 20;  # After removal of gaps, if less, gaps are not removed.
our $LENGTH_OF_NAME     = 10;




our $MULTIPLE_TREES_FILE_SUFFIX  = ".mlt";
our $LOG_FILE_SUFFIX             = ".log";
our $ALIGN_FILE_SUFFIX           = ".aln";
our $TREE_FILE_SUFFIX            = ".nhx";
our $ADDITION_FOR_RIO_ANNOT_TREE = ".rio";
our $SUFFIX_PWD                  = ".pwd";
our $MULTIPLE_PWD_FILE_SUFFIX    = ".mpwd";
our $SUFFIX_BOOT_STRP_POS        = ".bsp";
our $SUFFIX_PWD_NOT_BOOTS        = ".nbd";
our $SUFFIX_HMM                  = ".hmm";

our $EXPASY_SPROT_SEARCH_DE      = "http://www.expasy.org/cgi-bin/sprot-search-de?";
our $EXPASY_SPROT_SEARCH_AC      = "http://www.expasy.org/cgi-bin/sprot-search-ac?";



# One argument: input multiple trees file
# Last modified: 07/05/01
sub executeConsense {
    my $in  = $_[ 0 ];

    &testForTextFilePresence( $in );
   
    system( "$CONSENSE >/dev/null 2>&1 << !
$in
Y
!" ) 
    && &dieWithUnexpectedError( "Could not execute \"$CONSENSE $in\"" ); 
    
    return;
}



# Four arguments:
# 1. options ("-" is not necessary)
# 2. alignment or pwd file
# 3. outfile
# 4. temp dir
# Last modified: 07/05/01
sub executePhyloPl {

    my $opts         = $_[ 0 ]; 
    my $B            = $_[ 1 ];
    my $C            = $_[ 2 ];
    my $D            = $_[ 3 ];
    
    &testForTextFilePresence( $B );

    $opts = "-".$opts;

    system( "$PHYLO_PL $opts $B $C $D" )
    && &dieWithUnexpectedError( "Could not execute \"$PHYLO_PL $opts $B $C $D\"" ); 
    
} ## executePhyloPl




# Two arguments:
# 1. Name of inputfile
# 2. matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
#    5 = VT; 6 = WAG; 7 = auto; PAM otherwise 
sub executePuzzleDQO {
    my $in            = $_[ 0 ];
    my $matrix_option = $_[ 1 ];
    my $mat           = "";
    
    &testForTextFilePresence( $in );

    $mat = setModelForPuzzle( $matrix_option );

    system( "$PUZZLE_DQO $in >/dev/null 2>&1 << !$mat
y
!" )
    && &dieWithUnexpectedError( "Could not execute \"$PUZZLE_DQO\"" );
    
    return;

} ## executePuzzleDQO




# Two arguments:
# 1. Name of inputfile
# 2. matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
#    5 = VT; 6 = WAG; 7 = auto; PAM otherwise
# Last modified: 01/28/02
sub executePuzzleDQObootstrapped {
    my $in            = $_[ 0 ];
    my $matrix_option = $_[ 1 ];
    

    my $l             = 0;
    my $slen          = 0;
    my $counter       = 0; 
    my $mat           = "";
    my $a             = "";
    my @a             = ();
    
    &testForTextFilePresence( $in );
   
    open( GRP, "<$in" ) || &dieWithUnexpectedError( "Cannot open file \"$in\"" );
    while( <GRP> ) { 
        if ( $_ =~ /^\s*\d+\s+\d+\s*$/ ) { 
            $counter++; 
        } 
    }
    close( GRP ); 

    $l   = `cat $in | wc -l`;
    $slen   = $l / $counter;

    system( "split -$slen $in $in.splt." )
    && &dieWithUnexpectedError( "Could not execute \"split -$slen $in $in.splt.\"" );
 
    @a = <$in.splt.*>;

    $mat = setModelForPuzzle( $matrix_option );

    foreach $a ( @a ) {
                
        system( "$PUZZLE_DQO $a >/dev/null 2>&1 << !$mat
y 
!" )
        && &dieWithUnexpectedError( "Could not execute \"$PUZZLE_DQO $a\"" );

        system( "cat $a.dist >> $in.dist" )
        && &dieWithUnexpectedError( "Could not execute \"cat outdist >> $in.dist\"" );
  
        unlink( $a, $a.".dist" );
    }
    
    return;

} ## executePuzzleDQObootstrapped



# Transfers a Pfam (SELEX) alignment to a 
# PHYLIP sequential style alignment.
# It only writes "match columns" as indicated by the
# "# RF" line ('x' means match).
#
# Three arguments:
# 1. infile name
# 2. outfile name
# 3. 1 to NOT ensure that match states contain only 'A'-'Z' or '-'
#
# Returns the number of match states (=length of output alignment),
#         the length of the input alignment,
#         the number of seqs in the input alignment  
#
# Last modified: 07/07/01
#
sub pfam2phylipMatchOnly { 

    my $infile          = $_[ 0 ];
    my $outfile         = $_[ 1 ];
    my $ne              = $_[ 2 ];
    my @seq_name        = ();
    my @seq_array       = ();
    my $return_line     = "";
    my $seq             = "";
    my $x               = 0;
    my $y               = 0;
    my $i               = 0;
    my $x_offset        = 0;
    my $max_x           = 0;
    my $rf_y            = 0;
    my $number_colum    = 0;
    my $not_ensure      = 0;
    my $saw_rf_line     = 0;
    
    if ( $ne && $ne == 1 ) {
        $not_ensure = 1;
    }

    &testForTextFilePresence( $infile );

    open( INPP, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );

    # This reads in the first block. It reads in the seq names.
    while ( 1 ) {
        if ( &isPfamSequenceLine( $return_line ) ) {
            $return_line =~ /^(\S+)\s+(\S+)/;
            $seq_name[ $y ] = substr( $1, 0, $LENGTH_OF_NAME );
            $seq = $2;
            for ( $x = 0; $x < length( $seq ); $x++ ) {
                $seq_array[ $x ][ $y ] = substr( $seq, $x, 1 );
            }            
            $y++;
        }
        elsif ( &isRFline( $return_line ) ) {
            $saw_rf_line = 1;
            $return_line =~ /\s+(\S+)\s*$/;
            $seq = $1;
            $x_offset = length( $seq );
            $rf_y = $y;
            for ( $x = 0; $x < $x_offset; $x++ ) {
                $seq_array[ $x ][ $rf_y ] = substr( $seq, $x, 1 );
            }  
            last;
        }

        $return_line = <INPP>;

        if ( !$return_line ) {
             &dieWithUnexpectedError( "Alignment not in expected format (no RF line)" );
        }
    }

    if ( $saw_rf_line != 1 ) {
         &dieWithUnexpectedError( "Alignment not in expected format (no RF line)" );
    }    

    $y = 0;
    $max_x = 0;

    # This reads all blocks after the 1st one.
    while ( $return_line = <INPP> ) {
        if ( &isPfamSequenceLine( $return_line ) ) {
            $return_line =~ /^\S+\s+(\S+)/;
            $seq = $1;
            for ( $x = 0; $x < length( $seq ); $x++ ) {
                $seq_array[ $x + $x_offset ][ $y % $rf_y ] = substr( $seq, $x, 1 );
            }           
            $y++;            
        } 
        elsif ( &isRFline( $return_line ) ) {
            if ( $y != $rf_y ) {
                &dieWithUnexpectedError( "Alignment not in expected format" );
            }

            $return_line =~ /\s+(\S+)\s*$/;
            $seq = $1;
            $max_x = length( $seq );
           
            for ( $x = 0; $x < length( $seq ); $x++ ) {
                $seq_array[ $x + $x_offset ][ $rf_y ] = substr( $seq, $x, 1 );
            }
  
            $y = 0;
            $x_offset = $x_offset + $max_x;
            $max_x = 0;
        }
    }
    
    close( INPP );

    # Counts the match states, and hence the number of aa in the alignment:
    for ( $x = 0; $x < $x_offset; $x++ ) {
        if ( !$seq_array[ $x ][ $rf_y ] ) {
            &dieWithUnexpectedError( "Alignment not in expected format" );
        }
        if ( $seq_array[ $x ][ $rf_y ] eq 'x' ) {
            $number_colum++;
        }
    }

    # Writes the file:

    open( OUTPP, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
    print OUTPP "$rf_y $number_colum\n";
    for ( $y = 0; $y < $rf_y; $y++ ) {
        print OUTPP "$seq_name[ $y ]";
        for ( $i = 0; $i < ( $LENGTH_OF_NAME - length( $seq_name[ $y ] ) ); $i++ ) {
            print OUTPP " ";
        }
        for ( $x = 0; $x < $x_offset; $x++ ) {
            if ( $seq_array[ $x ][ $rf_y ] eq 'x' ) {
                if ( !$seq_array[ $x ][ $y ] ) {
                    &dieWithUnexpectedError( "Alignment not in expected format" );
                }
                if ( $not_ensure != 1 && $seq_array[ $x ][ $y ] !~ /[A-Z]|-/ ) {
                    &dieWithUnexpectedError( "Alignment not in expected format (match states must only contain 'A'-'Z' or '-')" );
                }
                print OUTPP "$seq_array[ $x ][ $y ]";
            }    
        }
        print OUTPP "\n";
    }  
    close( OUTPP );

    return $number_colum, $x_offset, $rf_y;

} ## pfam2phylipMatchOnly



# Returns whether the argument (a String) 
# starts with a SWISS-PROT name (SEQN_SPECI).
# Last modified: 06/21/01
sub startsWithSWISS_PROTname {
    return ( $_[ 0 ] =~ /^[A-Z0-9]{1,4}_[A-Z0-9]{1,5}/ );
}



# Returns whether the argument starts with XXX.. XXXXX.. and the first
# character is not a "#". 
# Last modified: 06/21/01
sub isPfamSequenceLine {
    return( !&isPfamCommentLine( $_[ 0 ] ) 
    && &containsPfamNamedSequence( $_[ 0 ] ) );
}



# Returns whether the argument does start with a "#".
# Last modified: 06/21/01
sub isPfamCommentLine {
    return ( $_[ 0 ] =~ /^#/ );
}



# Returns whether the argument starts with XXX XXXXX. 
# Last modified: 06/21/01
sub containsPfamNamedSequence {
    return ( $_[ 0 ] =~ /^\S+\s+\S+/ );
}


# Returns whether the argument starts with XXX XXXXX. 
# Last modified: 06/21/01
sub isRFline {
    return ( $_[ 0 ] =~ /^#.*RF/ );
}



# Three arguments:
# 1. pairwise distance file
# 2. number of bootstraps
# 3. initial tree: BME, GME or NJ
# Last modified: 2008/12/31
sub executeFastme {
    my $inpwd    = $_[ 0 ];
    my $bs       = $_[ 1 ];
    my $init_opt = $_[ 2 ];
      
    &testForTextFilePresence( $inpwd );
    my $command = "";
    if ( $bs > 0 ) {
        $command = "$FASTME -b $init_opt -i $inpwd -n $bs -s b";
    }
    else {
        $command = "$FASTME -b $init_opt -i $inpwd -s b";
    }    
    print $command;
    
    system( $command );
    

} ## executeFastme


# Four arguments:
# 1. pairwise distance file
# 2. number of bootstraps
# 3. seed for random number generator
# 4. lower-triangular data matrix? 1: yes; no, otherwise
sub executeNeighbor {
    my $inpwd  = $_[ 0 ];
    my $bs     = $_[ 1 ];
    my $s      = $_[ 2 ];
    my $l      = $_[ 3 ];
    my $multi  = "";
    my $lower  = "";
    
    &testForTextFilePresence( $inpwd );
   
    if (  $bs >= 2 ) {
        $multi = "
M
$bs
$s";
    }
    if ( $l == 1 ) {
        $lower = "
L"; 
    }

    system( "$NEIGHBOR >/dev/null 2>&1 << !
$inpwd$multi$lower
2
3
Y
!" )
    && &dieWithUnexpectedError( "Could not execute \"$NEIGHBOR $inpwd$multi$lower\"" );
   
} ## executeNeighbor


# Seven arguments:
# 1. pairwise distance file
# 2. number of bootstraps
# 3. seed for random number generator
# 4. number of jumbles for input order 
# 5. lower-triangular data matrix? 1: yes; no, otherwise
# 6. FM for Fitch-Margoliash, ME for ME# 6.
# 7. 1 to use globale rearragements
sub executeFitch {
    my $inpwd            = $_[ 0 ];
    my $bs               = $_[ 1 ];
    my $s                = $_[ 2 ];
    my $j                = $_[ 3 ];
    my $l                = $_[ 4 ];
    my $m                = $_[ 5 ];
    my $use_global_rearr = $_[ 6 ];
    my $jumble = "";
    my $multi  = "";
    my $lower  = "";
    my $method = "";
   
    my $global = "";
    if ( $use_global_rearr == 1 ) { 
        $global = "
G"; 
    }
    
    &testForTextFilePresence( $inpwd );

    if ( $m eq "FM" ) {
        $method = "";
    }
    elsif ( $m eq "ME" ) {
        $method = "
D"; 
    }
    else {
        &dieWithUnexpectedError( "method for FITCH must be either FM or ME" );    
    }

    if ( $j >= 1 ) {
        $jumble = "
J
$s
$j"; 
    }
   
    if (  $bs >= 2 ) {
        $multi = "
M
$bs
$s";
    }
    if ( $l == 1 ) {
        $lower = "
L"; 
    }

    # jumble must be set BEFORE multi!
    system( "$FITCH 2>&1 << !
$inpwd$method$global$jumble$multi$lower
3
Y
!" )
    && &dieWithUnexpectedError( "Could not execute \"$FITCH $inpwd$method$global$jumble$multi$lower\"" );
    # 3: Do NOT print out tree
  
} ## executeFitch



# Two arguments:
# 1. pairwise distance file
# 2. outfile
sub executeBionj {
    my $inpwd    = $_[ 0 ];
    my $out      = $_[ 1 ];
      
    &testForTextFilePresence( $inpwd );
    my $command = "$BIONJ $inpwd $out";
      
    system( $command )
    && &dieWithUnexpectedError( $command );
    
} 

# Four arguments:
# 1. (effective) sequence length
# 2. (effective) number of bases
# 3. pairwise distance file
# 4. outfile
sub executeWeighbor {
    my $L = $_[ 0 ];
    my $b = $_[ 1 ];
    my $i = $_[ 2 ];
    my $o = $_[ 3 ];
      
    &testForTextFilePresence( $i );
    my $command = "$WEIGHBOR -L $L -b $b -i $i -o $o";
      
    system( $command )
    && &dieWithUnexpectedError( $command );
    
} 

# Six arguments:
# 1. DNA or Amino-Acids sequence filename (PHYLIP format)
# 2. number of data sets to analyse (ex:3)
# 3. Model: JTT | MtREV | Dayhoff | WAG | VT | DCMut | Blosum62 (Amino-Acids)
# 4. number of relative substitution rate categories (ex:4), positive integer
# 5. starting tree filename (Newick format), your tree filename | BIONJ for a distance-based tree
# 6. 1 to estimate proportion of invariable sites, otherwise, fixed proportion "0.0" is used 
# PHYML produces several results files :
# <sequence file name>_phyml_lk.txt : likelihood value(s)
# <sequence file name>_phyml_tree.txt : inferred tree(s)
# <sequence file name>_phyml_stat.txt : detailed execution stats 
sub executePhyml {
    my $sequences = $_[ 0 ]; 
    my $data_sets = $_[ 1 ]; 
    my $model     = $_[ 2 ]; 
    my $nb_categ  = $_[ 3 ];  
    my $tree      = $_[ 4 ];
    my $estimate_invar_sites = $_[ 5 ];
    
    if ( $data_sets < 1 ) {
        $data_sets = 1
    }
    
    my $invar          = "0.0";   # proportion of invariable sites,
                                  # a fixed value (ex:0.0) | e to get the maximum likelihood estimate
    if ( $estimate_invar_sites == 1 ) {
        $invar = "e";
    }
    
    my $data_type      = "1";     # 0 = DNA | 1 = Amino-Acids
    my $format         = "i";     # i = interleaved sequence format | s = sequential
    my $bootstrap_sets = "0";     # number of bootstrap data sets to generate (ex:2)
                                  # only works with one data set to analyse
  
    my $alpha          = "e";     # gamma distribution parameter,
                                  # a fixed value (ex:1.0) | e to get the maximum likelihood estimate
   
    my $opt_topology   = "y";     # optimise tree topology ? y | n
    my $opt_lengths    = "y";     # optimise branch lengths and rate parameters ? y | n
      
    if ( $data_sets > 1 ) {
        # No need to calc branch lengths for bootstrapped analysis
        $opt_lengths = "n";
    } 
      
    &testForTextFilePresence( $sequences );
    my $command = "$PHYML $sequences $data_type $format $data_sets $bootstrap_sets $model $invar $nb_categ $alpha $tree $opt_topology $opt_lengths";
      
    print( "\n$command\n");  
      
    system( $command )
    && &dieWithUnexpectedError( $command );
    
} 




# Four arguments:
# 1. name of alignment file (in correct format!)
# 2. number of bootstraps
# 3. jumbles: 0: do not jumble; >=1 number of jumbles
# 4. seed for random number generator
sub executeProtpars {
    my $align  = $_[ 0 ];
    my $bs     = $_[ 1 ];
    my $rand   = $_[ 2 ];
    my $s      = $_[ 3 ];
    my $jumble = "";
    my $multi  = "";
    
   
    &testForTextFilePresence( $align );

    if ( $bs > 1 && $rand < 1 ) {
        $rand = 1;
    }

    if ( $rand >= 1 ) {
        $jumble = "
J
$s
$rand"; 
    }
   
    if (  $bs > 1 ) {
        $multi = "
M
D
$bs";
    }
   
    system( "$PROTPARS  2>&1 << !
$align$jumble$multi
Y
!" )
    && &dieWithUnexpectedError( "Could not execute \"$PROTPARS $align$jumble$multi\"" );
    # 3: Do NOT print out tree
    
      
    return;

} ## executeProtpars



# "Model of substitution" order for DQO TREE-PUZZLE 5.0:
# Auto
# m -> Dayhoff (Dayhoff et al. 1978)
# m -> JTT (Jones et al. 1992)
# m -> mtREV24 (Adachi-Hasegawa 1996)
# m -> BLOSUM62 (Henikoff-Henikoff 92)
# m -> VT (Mueller-Vingron 2000) 
# m -> WAG (Whelan-Goldman 2000)
# m -> Auto
# One argument:
# matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
# 5 = VT; 6 = WAG; 7 = auto; PAM otherwise
# Last modified: 07/07/01
sub setModelForPuzzle {
    my $matrix_option = $_[ 0 ];
    my $matr          = "";

    if ( $matrix_option == 0 ) { # JTT
        $matr = "
m
m";
    }
    elsif ( $matrix_option == 2 ) { # BLOSUM 62
        $matr = "
m
m
m
m";   
    }
    elsif ( $matrix_option == 3 ) { # mtREV24
        $matr = "
m
m
m";
    }
    elsif ( $matrix_option == 5 ) { # VT 
        $matr = "
m
m
m
m
m";
    }
    elsif ( $matrix_option == 6 ) { # WAG
        $matr = "
m
m
m
m
m
m";
    }
    elsif ( $matrix_option == 7 ) { # auto
        $matr = "";
    }          
    else { # PAM
        $matr = "
m"       
    }   

    return $matr;

} ## setModelForPuzzle

# One argument:
# Model of rate heterogeneity:
#    1 for "8 Gamma distributed rates"
#    2 for "Two rates (1 invariable + 1 variable)"
#    3 for "Mixed (1 invariable + 8 Gamma rates)"
#    otherwise: Uniform rate
# Last modified: 09/08/03 
sub setRateHeterogeneityOptionForPuzzle {
    my $rate_heterogeneity_option  = $_[ 0 ];
    my $opt                        = "";

    if ( $rate_heterogeneity_option == 1 ) {
        $opt = "
w";
    }
    elsif ( $rate_heterogeneity_option == 2 ) {
        $opt = "
w
w";   
    }
    elsif ( $rate_heterogeneity_option == 3 ) {
        $opt = "
w
w
w";
    }
    else {
        $opt = "";       
    }   

    return $opt;
} ## setRateHeterogeneityOptionForPuzzle 


# One argument:
# Parameter estimates: 1 for "Exact (slow)"; "Approximate (faster)" otherwise
# Last modified: 09/08/03 
sub setParameterEstimatesOptionForPuzzle {
    my $parameter_estimates_option  = $_[ 0 ];
    my $opt                         = "";

    if ( $parameter_estimates_option == 1 ) {
        $opt = "
e";
    }
    else {
        $opt = "";       
    }   

    return $opt;
} ## setParameterEstimatesOptionForPuzzle 



# three/four/five arguments:
# 1. Name of inputfile
# 2. matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
#    5 = VT; 6 = WAG; 7 = auto; PAM otherwise
# 3. Number of sequences in alignment
# 4. Parameter estimates: 1 for "Exact (slow)"; "Approximate (faster)" otherwise
# 5. Model of rate heterogeneity:
#    1 for "8 Gamma distributed rates"
#    2 for "Two rates (1 invariable + 1 variable)"
#    3 for "Mixed (1 invariable + 8 Gamma rates)"
#    otherwise: Uniform rate
sub executePuzzleBootstrapped {
    my $in                         = $_[ 0 ];
    my $matrix_option              = $_[ 1 ];
    my $number_of_seqs             = $_[ 2 ];
    my $parameter_estimates_option = $_[ 3 ];
    my $rate_heterogeneity_option  = $_[ 4 ];

    my $l             = 0;
    my $slen          = 0;
    my $counter       = 0; 
    my $mat           = "";
    my $est           = "";
    my $rate          = "";
    my $a             = "";
    my @a             = ();
    
    &testForTextFilePresence( $in );
   
    open( GRP, "<$in" ) || die "\n\n$0: Unexpected error: Cannot open file <<$in>>: $!";
    while( <GRP> ) { 
        if ( $_ =~ /^\s*\d+\s+\d+\s*$/ ) { 
            $counter++; 
        } 
    }
    close( GRP ); 

    $l   = `cat $in | wc -l`;
    $slen   = $l / $counter;

    system( "split --suffix-length=4 -$slen $in $in.splt." )
    && die "\n\n$0: executePuzzleDQObootstrapped: Could not execute \"split --suffix-length=4 -$slen $in $in.splt.\": $!";
    
    @a = <$in.splt.*>;
   
    $mat = setModelForPuzzle( $matrix_option );
     if ( $parameter_estimates_option ) {
        $est = &setParameterEstimatesOptionForPuzzle( $parameter_estimates_option );
    }
    if ( $rate_heterogeneity_option ) {
        $rate = &setRateHeterogeneityOptionForPuzzle( $rate_heterogeneity_option );
    }
    
    my $k="";
    if (  $number_of_seqs <= 257 ) {
        $k = "k";
    }

    foreach $a ( @a ) {
        print "-".$a."\n";        
        system( "$PUZZLE $a << !
$k
k
k$mat$est$rate
y 
!" )
        && die "$0: Could not execute \"$PUZZLE $a\"";

        system( "cat $a.dist >> $in.dist" )
        && die "$0: Could not execute \"cat outdist >> $in.dist\"";
  
        unlink( $a, $a.".dist", $a.".puzzle" );
    }
    
    return;

} ## executePuzzleBootstrapped





# three/four/five arguments:
# 1. Name of inputfile
# 2. Matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
#    5 = VT; 6 = WAG; 7 = auto; PAM otherwise
# 3. Number of sequences in alignment
# 4. Parameter estimates: 1 for "Exact (slow)"; "Approximate (faster)" otherwise
# 5. Model of rate heterogeneity:
#    1 for "8 Gamma distributed rates"
#    2 for "Two rates (1 invariable + 1 variable)"
#    3 for "Mixed (1 invariable + 8 Gamma rates)"
#    otherwise: Uniform rate
sub executePuzzle {
    my $in                         = $_[ 0 ];
    my $matrix_option              = $_[ 1 ];
    my $number_of_seqs             = $_[ 2 ];
    my $parameter_estimates_option = $_[ 3 ];
    my $rate_heterogeneity_option  = $_[ 4 ];
    my $mat                        = "";
    my $est                        = "";
    my $rate                       = "";
    
    &testForTextFilePresence( $in );

    $mat = &setModelForPuzzle( $matrix_option );
    if ( $parameter_estimates_option ) {
        $est = &setParameterEstimatesOptionForPuzzle( $parameter_estimates_option );
    }
    if ( $rate_heterogeneity_option ) {
        $rate = &setRateHeterogeneityOptionForPuzzle( $rate_heterogeneity_option );
    }
    
    my $k="";
    if (  $number_of_seqs <= 257 ) {
        $k = "k";
    }


    system( "$PUZZLE $in << !
$k
k
k$mat$est$rate
y
!" )
    && die "$0: Could not execute \"$PUZZLE\"";
    
    return;

} ## executePuzzle




# Preparation of the pwd file
sub addDistsToQueryToPWDfile {
    my $pwd_file          = $_[ 0 ];
    my $disttoquery_file  = $_[ 1 ];
    my $outfile           = $_[ 2 ];
    my $name_of_query     = $_[ 3 ];
    my $name_of_query_    = ""; 
    my $return_line_pwd   = "";
    my $return_line_dq    = "";
    my $num_of_sqs        = 0;
    my $block             = 0;
    my $name_from_pwd     = "X";
    my $name_from_dq      = "Y";
    my @dists_to_query    = ();
    my $i                 = 0;
    
    &testForTextFilePresence( $pwd_file );
    &testForTextFilePresence( $disttoquery_file );
    
    $name_of_query_ = $name_of_query;
    for ( my $j = 0; $j <= ( $LENGTH_OF_NAME - length( $name_of_query ) - 1 ); ++$j ) {
        $name_of_query_ .= " ";
    }

    open( OUT_AD, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
    open( IN_PWD, "$pwd_file" ) || &dieWithUnexpectedError( "Cannot open file \"$pwd_file\"" );
    open( IN_DQ, "$disttoquery_file" ) || &dieWithUnexpectedError( "Cannot open file \"$disttoquery_file\"" );
    
    W: while ( $return_line_pwd = <IN_PWD> ) {
        

        if ( $return_line_pwd =~ /^\s*(\d+)\s*$/ ) {
            $num_of_sqs = $1;
            $num_of_sqs++;
            if ( $block > 0 ) {
                print OUT_AD "$name_of_query_  ";
                for ( my $j = 0; $j < $i; ++$j ) {
                    print OUT_AD "$dists_to_query[ $j ]  ";
                }
                print OUT_AD "0.0\n";
            }
            print OUT_AD "  $num_of_sqs\n";
            $block++;
            @dists_to_query = ();
            $i = 0;
        }

        if ( $block == 1 
        && $return_line_pwd =~ /^\s*(\S+)\s+\S+/ ) {
            $name_from_pwd = $1;
            
            if ( !defined( $return_line_dq = <IN_DQ> ) ) {
                &dieWithUnexpectedError( "\"$disttoquery_file\" seems too short" );
            }
            
            if ( $return_line_dq !~ /\S/ ) {
                if ( !defined( $return_line_dq = <IN_DQ> ) ) {
                    &dieWithUnexpectedError( "\"$disttoquery_file\" seems too short" );
                }
            }
            $return_line_dq =~ /^\s*(\S+)\s+(\S+)/;
            $name_from_dq = $1;
            $dists_to_query[ $i++ ] = $2;
           
            
            if ( $name_from_pwd ne $name_from_dq ) {
                &dieWithUnexpectedError( "Order of sequence names in \"$pwd_file\" and \"$disttoquery_file\" is not the same" );
            }
            print OUT_AD $return_line_pwd;
          
        } 
        elsif ( $block > 1 
        && $return_line_pwd =~ /^\s*(\S+)\s+\S+/ ) {
            $name_from_pwd = $1;
            if ( !defined( $return_line_dq = <IN_DQ> ) ) {
                &dieWithUnexpectedError( "\"$disttoquery_file\" seems too short" );
            }
            if ( $return_line_dq !~ /\S/ ) {
                if ( !defined( $return_line_dq = <IN_DQ>) ) {
                    &dieWithUnexpectedError( "\"$disttoquery_file\" seems too short" );
                }
            }
            $return_line_dq =~ /^\s*\S+\s+(\S+)/;
            $dists_to_query[ $i++ ] = $1;
            print OUT_AD $return_line_pwd;
        }
    }
    print OUT_AD "$name_of_query_  ";
    for ( my $j = 0; $j < $i; ++$j ) {
        print OUT_AD "$dists_to_query[ $j ]  ";
    }
    print OUT_AD "0.0\n";

    close( OUT_AD );
    close( IN_PWD );
    close( IN_DQ );
    return $block;
    
} ## addDistsToQueryToPWDfile




# Three arguments:
# 1. HMMER model db
# 2. name of HMM
# 3. outputfile name
# Last modified: 02/27/01
sub executeHmmfetch {

    my $db      = $_[ 0 ];
    my $name    = $_[ 1 ];
    my $outfile = $_[ 2 ];
    
    system( "$HMMFETCH $db $name > $outfile" )
    && &dieWithUnexpectedError( "Could not execute \"$HMMFETCH $db $name > $outfile\"" );
    return;

} ## executeHmmfetch



# Checks wether a file is present, not empty and a plain textfile.
# One argument: name of file.
# Last modified: 07/07/01
sub testForTextFilePresence {
    my $file = $_[ 0 ];
    unless ( ( -s $file ) && ( -f $file ) && ( -T $file ) ) {
        dieWithUnexpectedError( "File \"$file\" does not exist, is empty, or is not a plain textfile" );
    }
} ## testForTextFilePresence


# Last modified: 02/21/03
sub addSlashAtEndIfNotPresent {
    my $filename = $_[ 0 ];
    $filename =~ s/\s+//g;
    unless ( $filename =~ /\/$/ ) {
       $filename = $filename."/";
    }
    return $filename;
} ## addSlashAtEndIfNotPresent



# Last modified: 02/15/02
sub exitWithWarning {
    
    my $text = $_[ 0 ];
    if ( defined( $_[ 1 ] ) && $_[ 1 ] == 1 ) {
        print( "<H4 class=\"error\">user error</H4>\n" );
        print( "<P>\n" );
        print( "<B>$text</B>\n" );
        print( "</P>\n" );
        print( "<P> &nbsp</P>\n" );
    }
    else {
        print( "\n\n$text\n\n" );
    }
    
    exit( 0 );

} ## exit_with_warning



# Last modified: 02/15/02
sub dieWithUnexpectedError {

    my $text = $_[ 0 ];

    die( "\n\n$0:\nUnexpected error (should not have happened):\n$text\n$!\n\n" );

} ## dieWithUnexpectedError



1;

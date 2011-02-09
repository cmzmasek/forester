# Copyright (C) 2002-2003 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek
# zmasek@genetics.wustl.edu
# http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Last modified 03/13/03


package rio_module2;
use strict;
require Exporter;

our $VERSION = 3.20;

our @ISA    = qw( Exporter );

our @EXPORT = qw( executeConsense
                  executeMakeTree
                  executePuzzleDQO
                  executePuzzleDQObootstrapped
                  pfam2phylipMatchOnly
                  startsWithSWISS_PROTname
                  isPfamSequenceLine
                  isPfamCommentLine
                  containsPfamNamedSequence
                  isRFline
                  executeNeighbor
                  executeProtpars
                  setModelForPuzzle
                  setRateHeterogeneityOptionForPuzzle
                  setParameterEstimatesOptionForPuzzle
                  executePuzzleBootstrapped
                  executePuzzle
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
                  $SEQBOOT
                  $NEIGHBOR
                  $PROTPARS
                  $CONSENSE
                  $PUZZLE 
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
                  $TRANSFERSBRANCHLENGHTS
                  $MAKETREE
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



# RIO itself:
# -----------
our $PATH_TO_FORESTER          = "/nfs/dm3/homedir1/czmasek/RIO1.24/"; 


# Java virtual machine:
# ---------------------
our $JAVA                      = "/usr/local/java/jdk/bin/java";



# Where all the temporary files can be created:
# ---------------------------------------------
our $TEMP_DIR_DEFAULT          = "/tmp/";

      

# Pfam data:
# ----------
our $PFAM_FULL_DIRECTORY       = "/path/to/Pfam/Full/"; 
our $PFAM_SEED_DIRECTORY       = "/path/to/Pfam/Seed/";
our $PFAM_HMM_DB               = "/path/to/Pfam/Pfam_ls"; # Need to run "hmmindex" on this
                                                          # to produce .ssi file.
                                                          # Then, for example
                                                          # "setenv HMMERDB /home/rio/pfam-6.6/"


$PATH_TO_FORESTER = &addSlashAtEndIfNotPresent( $PATH_TO_FORESTER );


# Description lines and species from SWISS-PROT and TrEMBL:
# ---------------------------------------------------------
our $TREMBL_ACDEOS_FILE        = $PATH_TO_FORESTER."data/trembl22_ACDEOS_1-6";
                                 
our $SWISSPROT_ACDEOS_FILE     = $PATH_TO_FORESTER."data/sp40_ACDEOS_1-6";



# Names of species which can be analyzed and analyzed 
# against (must also be in tree $SPECIES_TREE_FILE_DEFAULT).
# By using a list with less species, RIO analyses become faster
# but lose phylogenetic resolution. 
# For many purposes, list "tree_of_life_bin_1-6_species_list"
# in "data/species/" might be sufficient:
# --------------------------------------------------------------
our $SPECIES_NAMES_FILE        = $PATH_TO_FORESTER."data/species/tree_of_life_bin_1-6_species_list";



# A default species tree in NHX format.
# For many purposes, tree "tree_of_life_bin_1-6.nhx"
# in "data/species/" might be fine:
# --------------------------------------------------
our $SPECIES_TREE_FILE_DEFAULT = $PATH_TO_FORESTER."data/species/tree_of_life_bin_1-6.nhx";



# Data for using precalculated distances:
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


# PHYLIP:
our $SEQBOOT                   = $PATH_TO_FORESTER."phylip_mod/exe/seqboot";
our $NEIGHBOR                  = $PATH_TO_FORESTER."phylip_mod/exe/neighbor";
our $PROTPARS                  = $PATH_TO_FORESTER."phylip_mod/exe/protpars";
our $CONSENSE                  = $PATH_TO_FORESTER."phylip_mod/exe/consense";

# TREE-PUZZLE:
our $PUZZLE                    = $PATH_TO_FORESTER."puzzle_mod/src/puzzle";
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
our $TRANSFERSBRANCHLENGHTS    = $JAVA." -cp $PATH_TO_FORESTER"."java forester.tools.transfersBranchLenghts";
our $MAKETREE                  = $PATH_TO_FORESTER."perl/makeTree.pl";
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
our $LENGTH_OF_NAME     = 26;




our $MULTIPLE_TREES_FILE_SUFFIX  = ".mlt";
our $LOG_FILE_SUFFIX             = ".log";
our $ALIGN_FILE_SUFFIX           = ".aln";
our $TREE_FILE_SUFFIX            = ".nhx";
our $ADDITION_FOR_RIO_ANNOT_TREE = ".rio";
our $SUFFIX_PWD                  = ".pwd";
our $SUFFIX_BOOT_STRP_POS        = ".bsp";
our $SUFFIX_PWD_NOT_BOOTS        = ".nbd";
our $SUFFIX_HMM                  = ".hmm";

our $EXPASY_SPROT_SEARCH_DE      = "http://www.expasy.org/cgi-bin/sprot-search-de?";
our $EXPASY_SPROT_SEARCH_AC      = "http://www.expasy.org/cgi-bin/sprot-search-ac?";



# One argument: input multiple trees file
# Last modified: 07/05/01
sub executeConsense {
    my $in  = $_[ 0 ];

    &testForTextFilePresence( "$in" );
   
    system( "$CONSENSE >/dev/null 2>&1 << !
$in
Y
!" ) 
    && &dieWithUnexpectedError( "Could not execute \"$CONSENSE \"" ); 
    
    return;
}



# Four arguments:
# 1. options ("-" is not necessary)
# 2. alignment or pwd file
# 3. outfile
# 4. temp dir
# Last modified: 07/05/01
sub executeMakeTree {

    my $opts         = $_[ 0 ]; 
    my $B            = $_[ 1 ];
    my $C            = $_[ 2 ];
    my $D            = $_[ 3 ];
    
    &testForTextFilePresence( $B );

    $opts = "-".$opts;

    system( "$MAKETREE $opts $B $C $D" )
    && &dieWithUnexpectedError( "Could not execute \"$MAKETREE $opts $B $C $D\"" ); 
    
} ## executeMakeTree




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




# Five arguments:
# 1. pairwise distance file
# 2. number of bootstraps
# 3. randomize_input_order: 0: do not randomize input order; >=1 jumble
# 4. seed for random number generator
# 5. lower-triangular data matrix? 1: yes; no, otherwise
# Last modified: 06/08/01
sub executeNeighbor {
    my $inpwd  = $_[ 0 ];
    my $bs     = $_[ 1 ];
    my $rand   = $_[ 2 ];
    my $s      = $_[ 3 ];
    my $l      = $_[ 4 ];
    my $jumble = "";
    my $multi  = "";
    my $lower  =  "";
    
    
    &testForTextFilePresence( $inpwd );

    if ( $rand >= 1 ) {
        $jumble = "
J
$s"; 
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


    system( "$NEIGHBOR >/dev/null 2>&1 << !
$inpwd$jumble$multi$lower
2
3
Y
!" )
    && &dieWithUnexpectedError( "Could not execute \"$NEIGHBOR $inpwd$jumble$multi$lower\"" );
    # 3: Do NOT print out tree
    
      
    return;

} ## executeNeighbor



# Four arguments:
# 1. name of alignment file (in correct format!)
# 2. number of bootstraps
# 3. jumbles: 0: do not jumble; >=1 number of jumbles
# 4. seed for random number generator
# Last modified: 03/13/04
sub executeProtpars {
    my $alin   = $_[ 0 ];
    my $bs     = $_[ 1 ];
    my $rand   = $_[ 2 ];
    my $s      = $_[ 3 ];
    my $jumble = "";
    my $multi  = "";
    
   
    &testForTextFilePresence( $alin );

    if ( $bs >= 2 && $rand < 1 ) {
         $rand = 1;
    }

    if ( $rand >= 1 ) {
        $jumble = "
J
$s
$rand"; 
    }
   
    if (  $bs >= 2 ) {
        $multi = "
M
D
$bs";
    }
   


    system( "$PROTPARS  2>&1 << !
$alin$jumble$multi
I
3
Y
!" )
    && &dieWithUnexpectedError( "Could not execute \"$PROTPARS $alin$jumble$multi\"" );
    # 3: Do NOT print out tree
    # I: Interleaved
      
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



# Two/three/four arguments:
# 1. Name of inputfile
# 2. matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
#    5 = VT; 6 = WAG; 7 = auto; PAM otherwise
# 3. Parameter estimates: 1 for "Exact (slow)"; "Approximate (faster)" otherwise
# 4. Model of rate heterogeneity:
#    1 for "8 Gamma distributed rates"
#    2 for "Two rates (1 invariable + 1 variable)"
#    3 for "Mixed (1 invariable + 8 Gamma rates)"
#    otherwise: Uniform rate
# Last modified: 09/08/03 (added 3rd and 4th parameter)
sub executePuzzleBootstrapped {
    my $in                         = $_[ 0 ];
    my $matrix_option              = $_[ 1 ];
    my $parameter_estimates_option = $_[ 2 ];
    my $rate_heterogeneity_option  = $_[ 3 ];

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

	system( "split -$slen $in $in.splt." )
    && die "\n\n$0: executePuzzleDQObootstrapped: Could not execute \"split -$slen $in $in.splt.\": $!";
    
    @a = <$in.splt.*>;
   
    $mat = setModelForPuzzle( $matrix_option );
     if ( $parameter_estimates_option ) {
        $est = &setParameterEstimatesOptionForPuzzle( $parameter_estimates_option );
    }
    if ( $rate_heterogeneity_option ) {
        $rate = &setRateHeterogeneityOptionForPuzzle( $rate_heterogeneity_option );
    }

    foreach $a ( @a ) {
        print "-".$a."\n";        
        system( "$PUZZLE $a << !
k
k$mat$est$rate
y 
!" )
        && die "$0: Could not execute \"$PUZZLE $a\"";

        system( "cat $a.dist >> $in.dist" )
        && die "$0: Could not execute \"cat outdist >> $in.dist\"";
  
        unlink( $a, $a.".dist", $a.".tree" );
    }
    
    return;

} ## executePuzzleBootstrapped





# Two/three/four arguments:
# 1. Name of inputfile
# 2. Matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
#    5 = VT; 6 = WAG; 7 = auto; PAM otherwise
# 3. Parameter estimates: 1 for "Exact (slow)"; "Approximate (faster)" otherwise
# 4. Model of rate heterogeneity:
#    1 for "8 Gamma distributed rates"
#    2 for "Two rates (1 invariable + 1 variable)"
#    3 for "Mixed (1 invariable + 8 Gamma rates)"
#    otherwise: Uniform rate
# Last modified: 09/08/03 (added 3rd and 4th parameter)
sub executePuzzle {
    my $in                         = $_[ 0 ];
    my $matrix_option              = $_[ 1 ];
    my $parameter_estimates_option = $_[ 2 ];
    my $rate_heterogeneity_option  = $_[ 3 ];
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
    

    system( "$PUZZLE $in << !
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

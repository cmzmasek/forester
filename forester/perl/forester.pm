# FORESTER -- software libraries and applications
# for evolutionary biology research and applications.
#
# Copyright (C) 2020 Christian M. Zmasek
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



package forester;
use strict;
require Exporter;

our $VERSION = 1.000;

our @ISA    = qw( Exporter );

our @EXPORT = qw( executeConsense
                  executePhyloPl
                  executeProtpars
                  setModelForPuzzle
                  setRateHeterogeneityOptionForPuzzle
                  setParameterEstimatesOptionForPuzzle
                  executePuzzleBootstrapped
                  executePuzzle
                  executeFastme
                  executeNeighbor
                  executeFitch
                  testForTextFilePresence
                  exitWithWarning
                  dieWithUnexpectedError
                  addSlashAtEndIfNotPresent
                  $LENGTH_OF_NAME
                  $MIN_NUMBER_OF_AA
                  $MULTIPLE_TREES_FILE_SUFFIX
                  $LOG_FILE_SUFFIX
                  $ALIGN_FILE_SUFFIX
                  $TREE_FILE_SUFFIX
                  $SUFFIX_PWD
                  $MULTIPLE_PWD_FILE_SUFFIX
                  $SUFFIX_PWD_NOT_BOOTS
                  $MATRIX_FOR_PWD 
                  $PRIOR_FILE_DIR
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
                  $RAXMLNG
                  $RAXMLNG_VERSION
                  $SFE
                  $SUPPORT_TRANSFER
                  $SUPPORT_STATISTICS
                  $NEWICK_TO_PHYLOXML
                  $PHYLO_PL
                  $BOOTSTRAPS
                  $PATH_TO_FORESTER
                  $JAVA
                  $TEMP_DIR_DEFAULT
                  
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
# $PROTPARS
# $RAXMLNG

# Software directory:
# ---------------------

our $SOFTWARE_DIR              = "/Users/czmasek/SOFT/";


# Java virtual machine:
# ---------------------
our $JAVA                      = "/usr/bin/java";


# Where all the temporary files can be created:
# ---------------------------------------------
our $TEMP_DIR_DEFAULT          = "/tmp/";


# Programs from Joe Felsenstein's PHYLIP package:
# -----------------------------------------------
our $SEQBOOT                   = $SOFTWARE_DIR."phylip-3.697/src/seqboot";
our $NEIGHBOR                  = $SOFTWARE_DIR."phylip-3.697/src/neighbor";
our $PROTPARS                  = $SOFTWARE_DIR."phylip-3.697/src/protpars";
our $PROML                     = $SOFTWARE_DIR."phylip-3.697/src/proml";
our $FITCH                     = $SOFTWARE_DIR."phylip-3.697/src/fitch";
our $CONSENSE                  = $SOFTWARE_DIR."phylip-3.697/src/consense";
our $PHYLIP_VERSION            = "3.697";

# TREE-PUZZLE:
# ------------
our $PUZZLE                    = $SOFTWARE_DIR."tree-puzzle-5.3.rc16-macosx/src/puzzle";
our $PUZZLE_VERSION            = "5.3.rc16";

# FASTME:
# -----------------------------------------------------
our $FASTME                    = $SOFTWARE_DIR."fastme-2.1.6.4/src/fastme";
our $FASTME_VERSION            = "2.1.6.4";

# RAXML:
# -----------------------------------------------------
our $RAXMLNG_VERSION           = "1.2.0";
our $RAXMLNG                   = "/Users/czmasek//SOFT/RAXML/raxml-ng";


# forester.jar. 
# --------------------------------------------------------------------------------------------------------------------

our $FORESTER_JAR              = "/Users/czmasek/IdeaProjects/forester/forester/java/forester.jar";



# End of variables which need to be set by the user for using "phylo_pl.pl".





# Tool from forester.jar to transfer support values:
# -------------------------------------------------
our $SUPPORT_TRANSFER            = $JAVA." -cp $FORESTER_JAR org.forester.application.support_transfer";



# Tool from forester.jar for simple statistics for support values:
# ----------------------------------------------------------------
our $SUPPORT_STATISTICS          = $JAVA." -cp $FORESTER_JAR org.forester.application.support_statistics";


# Tool from forester.jar to transfer nh to phyloXML:
# -------------------------------------------------
our $NEWICK_TO_PHYLOXML          = $JAVA." -cp $FORESTER_JAR org.forester.application.phyloxml_converter";



# FORESTER itself (currently not needed for "phylo_pl.pl"):
# ---------------------------------------------------------
our $PATH_TO_FORESTER          = ""; 



$PATH_TO_FORESTER = &addSlashAtEndIfNotPresent( $PATH_TO_FORESTER );




#
# End of variables which need to be set by the user.
#
# =============================================================================
# =============================================================================




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These variables should normally not be changed:



our $BOOTSTRAPS         = 100;
our $MIN_NUMBER_OF_AA   = 20;  # After removal of gaps, if less, gaps are not removed.
our $LENGTH_OF_NAME     = 10;
my  $FASTME_T_OPTION    = "-T 6";



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



# Three arguments:
# 1. pairwise distance file
# 2. number of bootstraps
# 3. output name
# Last modified: 2020/05/19
sub executeFastme {
    my $inpwd  = $_[ 0 ];
    my $bs     = $_[ 1 ];
    my $output = $_[ 2 ];  
    &testForTextFilePresence( $inpwd );
    #########
    my $inpwd_reformated = $inpwd."_";
    open(FH, '<', $inpwd) or die $!;
    open(W, '>', $inpwd_reformated ) or die $!;
    my $prev = '';
    while(<FH>) {
        my $line = $_;
        chomp($line);
        if ( $line =~ /^\s+\d+\s*$/ ) {
            if ( length( $prev) > 0 ) {
                print W $prev;
                print W "\n";
                $prev = '';
            }
            print W $line;
            print W ( "\n" );
       }
       elsif ( $line =~ /^\S+/ ) {
           if ( length( $prev) > 0 ) {
               print W $prev;
               print W "\n";
           }
           $prev = $line;
       }
       elsif ( $line =~ /^\s+(.+)/ ) {
           $prev = $prev."  ".$1;
       }
    }
    if ( length( $prev) > 0 ) {
        print W $prev;
        print W "\n";
    }
    close(FH);
    close(W);
    ########
    
    my $command = "";
    if ( $bs > 1 ) {
        $command = "$FASTME -n -s -i $inpwd_reformated -D $bs -o $output -v 2";
    }
    else {
        $command = "$FASTME -n -s -i $inpwd_reformated -o $output -v 2";
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


# "Model of substitution" order for TREE-PUZZLE 5.2:
# For amino acids:
# Auto
# m -> Dayhoff (Dayhoff et al. 1978)
# m -> JTT (Jones et al. 1992)
# m -> mtREV24 (Adachi-Hasegawa 1996)
# m -> BLOSUM62 (Henikoff-Henikoff 92)
# m -> VT (Mueller-Vingron 2000) 
# m -> WAG (Whelan-Goldman 2000)
# m -> Auto
#
# For nucleotides:
# HKY (Hasegawa et al. 1985)
# m -> TN (Tamura-Nei 1993)
# m -> GTR (e.g. Lanave et al. 1980)
# m -> SH (Schoeniger-von Haeseler 1994)
# m -> HKY (Hasegawa et al. 1985)
#
# One argument:
# matrix option:
# 0 = JTT
# 2 = BLOSUM 62
# 3 = mtREV24
# 5 = VT
# 6, 13 = WAG
# 7 = auto
# 9 = HKY [na]
# 10 = TN [na]
# 11 = GTR [na]
# 12 = SH [na]
#     PAM otherwise
# Last modified: 17/04/26
sub setModelForPuzzle {
    my $matrix_option = $_[ 0 ];
    my $matr          = "";

    if ( $matrix_option == 0 || $matrix_option == 11 ) { # JTT or GTR
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
    elsif ( $matrix_option == 3 || $matrix_option == 12) { # mtREV24 or SH
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
    elsif ( $matrix_option == 6 || $matrix_option == 13 ) { # WAG
        $matr = "
m
m
m
m
m
m";
    }
    elsif ( $matrix_option == 7 || $matrix_option == 9 ) { # auto or HKY
        $matr = "";
    } 
    elsif ( $matrix_option == 10 ) { # TN
        $matr = "
m"       
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

    system( "split -a 4 -$slen $in $in.splt." )
    && die "\n\n$0: executePuzzleDQObootstrapped: Could not execute \"split -a 4 -$slen $in $in.splt.\": $!";
    
    @a = <$in.splt.*>;
   
    $mat = setModelForPuzzle( $matrix_option );
     if ( $parameter_estimates_option ) {
        $est = &setParameterEstimatesOptionForPuzzle( $parameter_estimates_option );
    }
    if ( $rate_heterogeneity_option ) {
        $rate = &setRateHeterogeneityOptionForPuzzle( $rate_heterogeneity_option );
    }
    
    #my $k="";
    #if (  $number_of_seqs <= 257 ) {
    #    $k = "k";
    #}

    foreach $a ( @a ) {
        print "-".$a."\n";        
        system( "$PUZZLE $a << !
k
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
    
    #my $k="";
    #if (  $number_of_seqs <= 257 ) {
    #    $k = "k";
    #}


    system( "$PUZZLE $in << !
k
k
k$mat$est$rate
y
!" )
    && die "$0: Could not execute \"$PUZZLE\"";
    
    return;

} ## executePuzzle




# Last modified: 02/21/03
sub addSlashAtEndIfNotPresent {
    my $filename = $_[ 0 ];
    $filename =~ s/\s+//g;
    unless ( $filename =~ /\/$/ ) {
       $filename = $filename."/";
    }
    return $filename;
} ## addSlashAtEndIfNotPresent


# Checks whether a file is present, not empty and a plain textfile.
# One argument: name of file.
# Last modified: 07/07/01
sub testForTextFilePresence {
    my $file = $_[ 0 ];
    unless ( ( -s $file ) && ( -f $file ) && ( -T $file ) ) {
        dieWithUnexpectedError( "File \"$file\" does not exist, is empty, or is not a plain textfile" );
    }
} ## testForTextFilePresence



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



# Last modified: 2020/01/23
sub dieWithUnexpectedError {

    my $text = $_[ 0 ];

    die( "\n\n$0:\nUnexpected forester error (should not have happened):\n$text\n$!\n\n" );

} ## dieWithUnexpectedError



1;

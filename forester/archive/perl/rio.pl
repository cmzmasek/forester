#!/usr/bin/perl -W

# rio.pl
# ------
#
# Copyright (C) 2000-2002 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Created: 11/25/00
# Author: Christian M. Zmasek
# zmasek@genetics.wustl.edu
# http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Last modified 09/06/03
#

#
#
# Available at:  http://www.genetics.wustl.edu/eddy/forester/
# RIO webserver: http://www.rio.wustl.edu/
#
# Reference:
# Zmasek C.M. and Eddy S.R. (2002)
# RIO: Analyzing proteomes by automated phylogenomics using
# resampled inference of orthologs.
# BMC Bioinformatics 3:14
# http://www.biomedcentral.com/1471-2105/3/14/
#
# It is highly recommended that you read this paper before
# installing and/or using RIO. (Included in the RIO 
# distribution as PDF: "RIO.pdf".)
# 
#
# Before rio.pl can be used, some variables in rio_module.pm need to be set, 
# as described in RIO_INSTALL.
#
# Usage: rio.pl <Mode: 1, 2, 3, or 4> <tagged arguments, single letter arguments>
# -----
#
#
# Examples:
# --------
# % RIO1.1/perl/rio.pl 1 A=aconitase Q=RIO1.1/LEU2_HAEIN N=QUERY_HAEIN O=out1 p I C E
# 
# % RIO1.1/perl/rio.pl 2 A=aconitase N=LEU2_LACLA/5-449 O=out2 p I C E
#
# % RIO1.1/perl/rio.pl 3 A=/path/to/my/pfam/Full/aconitase H=aconitase Q=RIO1.1/LEU2_HAEIN N=QUERY_HAEIN O=out3 p I C E
#
# % RIO1.1/perl/rio.pl 4 A=/path/to/my/pfam/Full/aconitase N=LEU2_LACLA/5-449 O=out4 p I C E
#
# % RIO1.1/perl/rio.pl 3 A=/path/to/my/pfam/Full/aconitase b=/path/to/my/pfam/Seed/aconitase Q=RIO1.1/LEU2_HAEIN N=QUERY_HAEIN O=out5 p I C E
#
#
# Modes:
# ------
#
# 1: RIO analysis based on precalculated pairwise distances
#    alignment does not contain query sequence
#
# 2: RIO analysis based on precalculated pairwise distances
#    alignment does contain query sequence
#
# 3: RIO analysis based on Pfam alignments,
#    alignment does not contain query sequence
#
# 4: RIO analysis based on Pfam alignments,
#    alignment does contain query sequence
#
#
#
# Tagged arguments:
# -----------------
#
# No "G=", "H=", "F=", "T=", "a=", "b=", "s", "f" in modes 1 and 2.
#
#
# A=<String> Pfam alignment name (mandatory). This specifies the alignment
#            against which the RIO analysis is to be performed.
#            In modes 1 and 2: Pfam model (alignment) name
#                              (e.g. "A=aconitase").
#            In modes 3 and 4: Pfam alignment path/name
#                              (e.g. "A=/path/to/your/pfam/Full/aconitase").
#
# Q=<String> Path/name of file containing the query sequence
#            (in FASTA format or raw sequence) (mandatory in modes 1 and 3).
#
# N=<String> Query name (mandatory). This must include the SWISS-PROT code
#            for the species of the query after a "_" (e.g. "N=QUERY_HAEIN").
#            If the query sequence is already in the alignment (modes 2 and 4)
#            the complete name needs to be specified -- including "/xxx-xxx".
#
# O=<String> Output file path/name (mandatory).
#
# T=<char>   Model for pairwaise distance calculation:
#            J=JTT, B=BLOSUM 62, M=mtREV24, V=VT, W=WAG, P=PAM.
#            BLOSUM 62 is default.
#            (Not in modes 1 and 2; these modes use $MATRIX_FOR_PWD instead.)
#
#            In modes 1 and 3, a HMM is needed to align the query sequence to
#            the alignment and either one of the following options must be
#            employed:
# H=<String> HMM name: This uses hmmfetch to retrieve a HMM from
#            $PFAM_HMM_DB.
# F=<String> HMM file: This directly reads the HMM from a file.
#
# S=<String> Species tree file path/name (in NHX format) (optional).
#            If not specified, $SPECIES_TREE_FILE_DEFAULT is used.
#
# G=<String> Species names file (optional). Only sequences associated with
#            species found in this file are used.
#            In the species names file, individual species names must be
#            separated by newlines and lines starting with "#" are ignored.    
#            While only sequences associated with species found in the species
#            tree ("S=") are used for the actual RIO analysis, this allows to 
#            remove sequences prior to tree calculation (which is the most
#            time consuming step).  
#
# P=<int>    Sort priority (default is 12):
#            0 : Ortholog  
#            1 : Ortholog,         Super ortholog
#            2 : Super ortholog,   Ortholog
#            3 : Ortholog,         Distance
#            4 : Distance,         Ortholog
#            5 : Ortholog,         Super ortholog,   Distance  
#            6 : Ortholog,         Distance,         Super ortholog
#            7 : Super ortholog,   Ortholog,         Distance
#            8 : Super ortholog,   Distance,         Ortholog
#            9 : Distance,         Ortholog,         Super ortholog
#           10 : Distance,         Super ortholog,   Ortholog
#           11 : Ortholog,         Subtree neighbor, Distance 
#           12 : Ortholog,         Subtree neighbor, Super ortholog,   Distance (default)
#           13 : Ortholog,         Super ortholog,   Subtree neighbor, Distance
#           14 : Subtree neighbor, Ortholog,         Super ortholog,   Distance
#           15 : Subtree neighbor, Distance,         Ortholog,         Super ortholog   
#           16 : Ortholog,         Distance,         Subtree neighbor, Super ortholog 
#           17 : Ortholog,         Subtree neighbor, Distance,         Super ortholog 
#
# a=<int>    Bootstraps for tree construction (not in modes 1 and 2).
#            Default is 100.   
#
# L=<int>    Threshold for orthologies for output. Default is 0.
# v=<int>    Threshold for ultra-paralogies for output. Default is 50.
#
# U=<int>    Threshold for orthologies for distance calculation. Default is 60.
#
# X=<int>    In case of more than one putative orthologs:
#            number of sd the distance query - LCA has to differ
#            from the mean to generate a warning. Default is 2.
#
# Y=<int>    In case of no putative orthologs:
#            number of sd the distance query - root has to differ
#            from mean to generate a warning. Default is 2.
#
# Z=<double> In case of one putative ortholog:
#            threshold for factor between the two distances to their
#            LCA (larger/smaller) to generate a warning. Default is 2.
#
# B=<int>    Threshold for subtree-neighborings. Default is 0.
#
# b=<String> Build HMM from seed alignment with "hmmbuild -s" (optional).
#            This is to prevent from finding multiple domains per sequence
#            (i.e. prevents "cutting" the query sequence). Give path/name to
#            Seed with this.
#
# j=<String> Name for temporary directory (optional).
#
# y=<int>    Seed for random number generator. Default is 41.
#
# I          Create and save a rooted, with duplication vs speciation,
#            and orthology information annotated gene tree.
#            If precalculated distances are used (modes 1 and 2): this gene
#            tree is a NJ tree calculated based on the non-bootstrap resampled
#            (original) pairwise distances.
#            If precalculated distances are not used (modes 3 and 4): this gene
#            is a consenus tree with ML branch length values and is also
#            annotated with bootstrap values for each node.
#
#            Options for output:
# p          Output ultra-paralogs.
# D          Description from SWISS-PROT and TrEMBL.
# C          Complete description from SWISS-PROT and TrEMBL.
# E          118 character output instead of 78 character output.
#
# K          Keep intermediate files (they will go into the same directory 
#            as the output file, their names are the same as of the output
#            file, with various suffixes added).
#
# s          Ignore non SWISS-PROT sequences (i.e. sequences from TrEMBL)
#            in the Pfam alignment.
#
# f          Try to ignore TrEMBL "fragments" (sequences with "fragment" in
#            their description).
#
# +          Parallel, use machines listed in file $NODE_LIST.
#
# x          RIO used as web server -- HTML output.
#
#
#
#
# History
# -------
# 09/06/03: Removal of minor bug. Only create consenus tree with ML branch length
#           values if "I" option used (in modes 3 or 4) -- the problem/bug was that
#           this tree was always created whether "I" was used or not. 
#


use strict;

use FindBin;
use lib $FindBin::Bin;
use Net::Ping;
use rio_module;

use File::Basename;


my $VERSION                           = "5.010";

my $E_VALUE_THRESHOLD                 = 0.01; # For HMMSEARCH.
my $SORT_DEFAULT                      = 12;
my $THRESHOLD_ORTHOLOGS_DEFAULT       = 0;
my $THRESHOLD_SN_DEFAULT              = 0;
my $THRESHOLD_ORTHOLOGS_DEFAULT_DC    = 60;
my $T_ULTRA_PARALOGS_DEFAULT          = 50;
my $WARN_NO_ORTHOS_DEFAULT            = 2;
my $WARN_MORE_THAN_ONE_ORTHO_DEFAULT  = 2;
my $WARN_ONE_ORTHO_DEFAULT            = 2;
my $MIN_NUMBER_OF_SEQS_IN_ALN         = 4;
my $BOOSTRAPS_FOR_MAKETREE_DEFAULT    = 100;
my $SEED_FOR_MAKETREE_DEFAULT         = 41;
my $MATRIX_DEFAULT                    = 2; # 2=BLOSUM62

my $DO_RIO_TEMP_OUTFILE               = "DoRIO_OUTFILE";
my $TEMP_HMM_FILE                     = "HMMFILE";

my $DEFAULT_OPTIONS_FOR_MAKETREE      = "XR";


# I/O files, names:
my $alignment                      = "";
my $hmm_file                       = "";
my $hmm_name                       = "";
my $seqX_file                      = "";
my $species_tree_file              = "";
my $outfile                        = "";
my $outfile_annot_nhx_tree         = "";
my $query_name                     = "";
my $multiple_trees_file            = "";
my $distance_matrix_file           = "";
my $maketree_out_tree_file         = "";
my $seed_aln_for_hmmbuild          = "";
my $temp_dir                       = "";
my $bsp_file                       = "";
my $pwd_file                       = "";
my $nbd_file                       = "";
my $output_dir                     = "";
my $species_names_file             = " "; # Must be " ".
my $options_for_makeTree           = "";


# multiple choice options:
my $mode                           = 0;
my $sort                           = $SORT_DEFAULT;
my $matrix_n                       = $MATRIX_DEFAULT; # 0=JTT 1=PAM 2=BLOSUM62 3=mtREV24 5=VT 6=WAG




# yes/no options:
my $description                    = 0;
my $complete_description           = 0;
my $long_output                    = 0;
my $keep                           = 0;
my $non_sp                         = 1; # 0 to remove non SP seqs.
my $safe_nhx                       = 0;
my $no_frags                       = 0;
my $output_ultraparalogs           = 0;
my $parallel                       = 0;
my $output_HTML                    = 0;


# numerical options:
my $warn_no_orthos                 = $WARN_NO_ORTHOS_DEFAULT;
my $warn_more_than_one_ortho       = $WARN_MORE_THAN_ONE_ORTHO_DEFAULT;
my $warn_one_ortho                 = $WARN_ONE_ORTHO_DEFAULT;
my $boostraps_for_makeTree         = $BOOSTRAPS_FOR_MAKETREE_DEFAULT;
my $seed_for_makeTree              = $SEED_FOR_MAKETREE_DEFAULT;
my $t_orthologs                    = $THRESHOLD_ORTHOLOGS_DEFAULT;
my $t_sn                           = $THRESHOLD_SN_DEFAULT;
my $t_orthologs_dc                 = $THRESHOLD_ORTHOLOGS_DEFAULT_DC;
my $t_ultra_paralogs               = $T_ULTRA_PARALOGS_DEFAULT;


# internal variables:
my $print_header_for_orthologies = 0;
my $print_header_for_s_paralogs  = 0;
my $length_of_alignment          = 0;
my $length_of_orig_alignment     = 0;
my $time                         = 0;
my $ii                           = 0;
my $j                            = 0;
my $jj                           = 0;
my $number_of_seqs_in_aln        = 0;
my $f                            = 0;
my $saw_distance_values          = 0;
my $saw_ultra_paralogs           = 0;
my $bootstraps                   = 0;
my $ext_nodes_in_trees_analyzed  = 0;
my $time_total                   = 0;
my $time_tree_calc               = 0;
my $time_tree_calcT              = 0;
my $time_rio                     = 0;
my $time_rioT                    = 0;
my $time_dqopuzzle               = 0;
my $time_dqopuzzleT              = 0;
my $time_addingdists             = 0;
my $time_addingdistsT            = 0;
my $processors                   = 0;
my $block_size                   = 0;
my $larger_blocks                = 0;
my $printed_ultra_paralogs       = 0;
 
my $dorio_outfile                = "";
my $options_for_DoRIO            = "";
my $ortho_name                   = "";
my $orthos                       = 0;
my $s_orthos                     = 0;
my $subtree_neighbors            = 0;
my $dist                         = 0;
my $s_para_name                  = "";
my $s_paras                      = 0;   
my $sort_priority                = "";
my $return_line                  = "";
my $matrix                       = "";
my $command_line                 = "";
my $command_line_for_hmmbuild    = "";
my $current_dir                  = "";
my @complete_names               = ();
my @temp_array                   = ();
my %Species_names_hash           = ();
my %AC_DE                        = (); # AC => DE from "ACDEOS" TrEMBL file.
my %SP_AC_DE                     = (); # ID => DE from "ACIDOS" SWISS-PROT file.
my %names_in_pwd_file            = ();
my @nodelist = ();

my $start_date                   = `date`;




# This analyzes the options:
# --------------------------

$time_total = time;

if ( @ARGV < 4 ) {
    &printHelp();
}

$command_line = "$0 ";
for ( $j = 0; $j < @ARGV; ++$j ) {
    $command_line .= "$ARGV[ $j ] ";
}

&analyzeCommandLine( @ARGV );

if ( $species_tree_file eq "" ) {
    $species_tree_file = $SPECIES_TREE_FILE_DEFAULT;
}

&CheckArguments;

$options_for_makeTree = $DEFAULT_OPTIONS_FOR_MAKETREE;
$options_for_makeTree .= "S".$seed_for_makeTree;


if ( $mode == 1 || $mode == 2 ) {
            
    if ( $mode == 1 ) {
        $hmm_file = $RIO_HMM_DIRECTORY.$alignment.$SUFFIX_HMM;
        $bsp_file = $RIO_BSP_DIRECTORY.$alignment.$SUFFIX_BOOT_STRP_POS;
        &userErrorCheckForTextFileExistence( $hmm_file );
        &userErrorCheckForTextFileExistence( $bsp_file );
    }

    $pwd_file  = $RIO_PWD_DIRECTORY.$alignment.$SUFFIX_PWD;
    $nbd_file  = $RIO_NBD_DIRECTORY.$alignment.$SUFFIX_PWD_NOT_BOOTS;
    $alignment = $RIO_ALN_DIRECTORY.$alignment.$ALIGN_FILE_SUFFIX;
    &userErrorCheckForTextFileExistence( $pwd_file );
    &userErrorCheckForTextFileExistence( $nbd_file );
    &userErrorCheckForTextFileExistence( $alignment );
    $no_frags                = 0;
    $non_sp                  = 1;
       
    $options_for_makeTree .= "F";
}
elsif ( $mode == 3 || $mode == 4 ) {
    if ( $safe_nhx == 1 ) {
        $options_for_makeTree .= "U";
    }
    else {
        $options_for_makeTree .= "#";
    }  
    $options_for_makeTree .= "D"; # To calc. and keep pairwise distances.
    $options_for_makeTree .= "B".$boostraps_for_makeTree;
      
}

if ( $output_HTML == 1 ) { 
    $| = 1;
    $complete_description = 1;
    $long_output          = 1;

}

if ( $mode == 1 || $mode == 3 || $mode == 4 ) {

    if ( $mode == 1 ) {
        $matrix_n = $MATRIX_FOR_PWD;
    }
    
    if ( $matrix_n == 0 ) {
        $options_for_makeTree .= "J";
        $matrix = "JTT (Jones et al. 1992)";
    }
    elsif ( $matrix_n == 1 ) {  # PAM is makeTree's default.
        $matrix = "PAM (Dayhoff et al. 1978)";
    }    
    elsif ( $matrix_n == 2 ) {
        $options_for_makeTree .= "L";
        $matrix = "BLOSUM 62 (Henikoff-Henikoff 92)";
    }
    elsif ( $matrix_n == 3 ) {
        $options_for_makeTree .= "M";
        $matrix = "mtREV24 (Adachi-Hasegawa 1996)";
    }
    elsif ( $matrix_n == 5 ) {
        $options_for_makeTree .= "T";
        $matrix = "VT (Mueller-Vingron 2000)";
    }
    elsif ( $matrix_n == 6 ) {
        $options_for_makeTree .= "W";
        $matrix = "WAG (Whelan-Goldman 2000)";
    }
    else {
        &dieWithUnexpectedError( "Failed sanity check" );
    }
}


# This creates the temp directory:
# --------------------------------

$ii = 0;

$time = time;

if ( $temp_dir eq "" ) { 
    $temp_dir = $TEMP_DIR_DEFAULT.$time.$ii;
}
else {
    $temp_dir = $temp_dir.$ii;
}

while ( -e $temp_dir ) {
    $ii++;
    $temp_dir =  $TEMP_DIR_DEFAULT.$time.$ii;
}

mkdir(  $temp_dir, 0700 ) || &dieWithUnexpectedError( "Could not create \"$temp_dir\"" );

unless ( ( -e $temp_dir ) && ( -d $temp_dir ) ) {
    &dieWithUnexpectedError( "\"$temp_dir\" does not exist, or is not a directory" );
}



# The analysis starts here:
# -------------------------

$dorio_outfile = $temp_dir."/".$DO_RIO_TEMP_OUTFILE;

$output_dir  = dirname( $outfile );

unless ( ( -e $output_dir ) && ( -d $output_dir ) ) {
    &userError( "Outfile directory (\"$output_dir\") does not exist,\n or is not a directory." );
}

if ( $mode == 1 || $mode == 3 ) {
    $query_name = substr( $query_name, 0, $LENGTH_OF_NAME - 10 );
}





if ( $mode == 1 || $mode == 3 ) {

    # Prepares the query file:
    # ------------------------
    $query_name = &seqFile2CleanedUpFastaFile( $seqX_file,
                                               "$temp_dir/QUERY_SEQ",
                                               $query_name );
    if ( $query_name eq "" ) {
        &userError( "Query file \"$seqX_file\") does not appear to contain a valid name\n and/or \"-N\" option has not been used." );
    }

    if ( $mode == 3 ) {
        # Prepares the HMM:
        # -----------------
        if ( $hmm_file eq "" ) {
            $hmm_file = $temp_dir."/".$TEMP_HMM_FILE;
            if ( $hmm_name ne "" ) {
                &executeHmmfetch( $PFAM_HMM_DB, $hmm_name, $hmm_file );
            }
            elsif ( $seed_aln_for_hmmbuild ne "" ) {
                $command_line_for_hmmbuild = &executeHmmbuild( $seed_aln_for_hmmbuild, $hmm_file );
            }
        }
       
    }
}




# This might remove non SWISS PROT seqs, TreMBL fragments,
# and seqs from species not in $species_names_file from the alignment:
# --------------------------------------------------------------------
if ( $mode == 3 || $mode == 4 ) { 
    #if ( $do_not_removeSeqsFromPfamAlign != 1 ) {

    if ( $mode == 3 ) {
        &removeSeqsFromPfamAlign( $alignment,
                                  $temp_dir."/ALIGN2",
                                  " ",
                                  $species_names_file,
                                  $non_sp,
                                  $no_frags );
    }
    else {
        &removeSeqsFromPfamAlign( $alignment,
                                  $temp_dir."/ALIGN2",
                                  $query_name,
                                  $species_names_file,
                                  $non_sp,
                                  $no_frags );
    }
    
}



# If necessary, this aligns the query to the pfam alignment
# using hmmsearch, p7extract.pl, multifetch.pl, and hmmalign
# from the HMMER package:
# ----------------------------------------------------------
if ( $mode == 1 || $mode == 3 ) {
    if ( $mode == 1 ) {
       
        $f = &alignWithHmmalign( $alignment,
                                 $temp_dir."/QUERY_SEQ",
                                 $hmm_file,
                                 $temp_dir."/HMMALIGNOUT",
                                 1 ); # --mapali

        
    }
    else {
      
        $f = &alignWithHmmalign( $temp_dir."/ALIGN2",
                                 $temp_dir."/QUERY_SEQ",
                                 $hmm_file,
                                 $temp_dir."/HMMALIGNOUT",
                                 0 ); # --withali
        
    }
    if ( $f != 1 ) { 
        if ( $alignment =~ /.+\/(.+)/ ) {
            $alignment = $1;
        }
        if ( $alignment =~ /(.+)\..+/ ) {
            $alignment = $1;
        }
        &cleanUpTempDir();
        if ( $output_HTML == 1 ) {
            &exitWithWarning( "query sequence does not contain sufficient similarity to the \"$alignment\" domain", 1 );
        }
        else {
            &exitWithWarning( "Query sequence does not contain sufficient similarity to the \"$alignment\" domain" );
        }
    }


    # In case query contains more than one of the same domain:

    @complete_names = &getCompleteName( $temp_dir."/HMMALIGNOUT", $query_name );

    if ( @complete_names < 1 ) {
        &dieWithUnexpectedError( "Could not find \"$query_name in $temp_dir"."/HMMALIGNOUT\"" );
    }
}
elsif ( $mode == 2 || $mode == 4 ) {
    # Here, this is just for checking:
    if ( $mode == 2 ) {
        @complete_names = &getCompleteName( $alignment, $query_name );
    } 
    elsif ( $mode == 4 ) {
        @complete_names = &getCompleteName( $temp_dir."/ALIGN2", $query_name );
    }
    if ( @complete_names < 1 ) {
        &dieWithUnexpectedError( "Could not find \"$query_name in $temp_dir"."/HMMALIGNOUT\"" );
    }
    @complete_names = ();
    $complete_names[ 0 ] = $query_name;
}

if ( $parallel == 1 ) {
    &readInNodesList();
    &pingNodes();
    $processors = scalar( @nodelist );
    if ( $processors < 2 ) {
        $parallel = 0;
    }
    if ( $processors > $BOOTSTRAPS  ) {
        $processors = $BOOTSTRAPS;
    }
    else {
       $block_size = int $BOOTSTRAPS / $processors;
       $larger_blocks = $BOOTSTRAPS - ( $block_size * $processors ); # number of blocks which have a size of
                                                                     # block_size + 1
    
   }
}


# This opens the output file:
# ---------------------------
if ( $output_HTML != 1 ) {
    open( OUT, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
}

# This starts printing to the output file:
# ----------------------------------------
&printHeader();



# This loop goes through the different domains of the query
# which aligned to the alignment (in modes 2 and 4 this can
# obviously be only one):
# -----------------------------------------------------------
for ( $jj = 0; $jj < @complete_names; ++$jj ) {
    
    if ( $mode == 1 ) {
        # Moves the query to the last line(s) of the alignment.
        # Removes other querie domains $complete_names[i]-- for which i != $jj
        # --------------------------------------------------------------------
       
        &moveToLast( $complete_names[ $jj ],       
                     $temp_dir."/HMMALIGNOUT",
                     $temp_dir."/MOVETOLASTOUT",
                     \@complete_names );

    }
    
    if ( $mode == 1 || $mode == 3 ) {
        if ( $mode == 1 ) {
            @temp_array = &pfam2phylipMatchOnly( $temp_dir."/MOVETOLASTOUT",
                                                 $temp_dir."/ALIGN2_PHYLIP_MO",
                                                 0 );
        }
        else {
            @temp_array = &pfam2phylipMatchOnly( $temp_dir."/HMMALIGNOUT",
                                                 $temp_dir."/ALIGN2",
                                                 1 );
        }
        $length_of_alignment      = $temp_array[ 0 ];
        $length_of_orig_alignment = $temp_array[ 1 ];
        $number_of_seqs_in_aln    = $temp_array[ 2 ];
    }
    elsif ( $mode == 2 || $mode == 4 ) {

        $query_name = $complete_names[ 0 ];
    
        if ( $mode == 4 ) {
            if ( !&startsWithSWISS_PROTname( $query_name ) ) {
                # Query is not a SWISS-PROT sequence.
                $query_name = &getCompleteNameForTrEMBLquerySeq( $temp_dir."/ALIGN2",
                                                                 $query_name );
            }

            $number_of_seqs_in_aln = &countSeqsInPfamAlign( $temp_dir."/ALIGN2" );
        }    
        else { 
            if ( !&startsWithSWISS_PROTname( $query_name ) ) {
                # Query is not a SWISS-PROT sequence.
                $query_name = &getCompleteNameForTrEMBLquerySeq( $alignment,
                                                                 $query_name );
            }
            $number_of_seqs_in_aln = &countSeqsInPfamAlign( $alignment );
        }


     
    } 
    
    if ( $number_of_seqs_in_aln < $MIN_NUMBER_OF_SEQS_IN_ALN ) {
        &cleanUpTempDir();
        if ( $output_HTML == 1 ) {
            &exitWithWarning( "Removal of sequences resulted in an alignment with less than $MIN_NUMBER_OF_SEQS_IN_ALN sequences ($number_of_seqs_in_aln)", 1 );
        }
        else {
            &exitWithWarning( "Removal of sequences resulted in an alignment with less than $MIN_NUMBER_OF_SEQS_IN_ALN sequences ($number_of_seqs_in_aln)" );
        }
    }


    if ( $mode == 1 ) {

        unlink( $temp_dir."/ALIGN2_BOOTSTRAPPED" );
        
        if ( $parallel == 1 ) {
                &executeBootstrap_cz( $BOOTSTRAPS,
                                      $bsp_file,
                                      $temp_dir."/ALIGN2_PHYLIP_MO",
                                      $temp_dir."/ALIGN2_BOOTSTRAPPED",
                                      $processors );
              
        }
        else {

            &executeBootstrap_cz( $BOOTSTRAPS,
                                  $bsp_file,
                                  $temp_dir."/ALIGN2_PHYLIP_MO",
                                  $temp_dir."/ALIGN2_BOOTSTRAPPED" );
          
        }

       
        $current_dir = `pwd`;
        $current_dir =~ s/\s//;

        chdir ( $temp_dir ) || &dieWithUnexpectedError( "Could not chdir to \"$temp_dir\"" );


        if ( $parallel == 1 ) {

            my $number       = 0;
            my $all_finished = 0;

            system( $RIO_SLAVE_DRIVER,
                    $block_size,
                    $larger_blocks,
                    $temp_dir."/ALIGN2_BOOTSTRAPPED",
                    $matrix_n,
                    $complete_names[ $jj ],
                    $pwd_file,
                    $temp_dir,
                    $seed_for_makeTree,
                    @nodelist ) 
            && &dieWithUnexpectedError( "Could not execute \"$RIO_SLAVE_DRIVER\"" );
           
            while ( $all_finished != 1 ) {
                for ( $number = 0; $number < $processors; $number++ ) {
                    unless ( -e "FINISHED_$number" ) {
	    	        $number = -1;
                    }
                }
                $all_finished = 1;
            }

            sleep( 1 );

            system( "mv",
                    "MAKETREEOUT".$MULTIPLE_TREES_FILE_SUFFIX."0",
                    "MAKETREEOUT".$MULTIPLE_TREES_FILE_SUFFIX )
            && &dieWithUnexpectedError( "$!" );

            for ( $number = 1; $number < $processors; $number++ ) {
                system( "cat MAKETREEOUT$MULTIPLE_TREES_FILE_SUFFIX$number >> MAKETREEOUT$MULTIPLE_TREES_FILE_SUFFIX" )
                && &dieWithUnexpectedError( "$!" );
                if ( unlink( "MAKETREEOUT$MULTIPLE_TREES_FILE_SUFFIX$number" ) != 1 ) {
                    &dieWithUnexpectedError( "Could not delete \"MAKETREEOUT$MULTIPLE_TREES_FILE_SUFFIX$number" );
                }
            }

            # Sanity check: Counts ";" in "MAKETREEOUT$MULTIPLE_TREES_FILE_SUFFIX".
            if ( `grep -c ';' MAKETREEOUT$MULTIPLE_TREES_FILE_SUFFIX` != $BOOTSTRAPS ) {
                &dieWithUnexpectedError( "\"MAKETREEOUT$MULTIPLE_TREES_FILE_SUFFIX\" does not contain $BOOTSTRAPS \";\"" );
            } 

            for ( $number = 0; $number < $processors; $number++ ) {
                if ( unlink( "FINISHED_$number" ) != 1 ) {
                    &dieWithUnexpectedError( "Could not delete \"FINISHED_$number\"" );
                }
            }

            &executeConsense( "MAKETREEOUT".$MULTIPLE_TREES_FILE_SUFFIX );
            unlink( "outfile", "intree" );

            system( "mv", "outtree", "MAKETREEOUT.nhx" )
            && &dieWithUnexpectedError( "$!" );


        } 
        else {
            $time_dqopuzzle = time; #time
            &executePuzzleDQObootstrapped( "ALIGN2_BOOTSTRAPPED", $matrix_n );
            $time_dqopuzzle = time - $time_dqopuzzle; #time
            $time_dqopuzzleT += $time_dqopuzzle; #time

            system( "mv", "ALIGN2_BOOTSTRAPPED.dist", "DISTs_TO_QUERY" )
            && &dieWithUnexpectedError( "$!" );
        }
       
       
        &executePuzzleDQO( "ALIGN2_PHYLIP_MO", $matrix_n );
        
        unlink( "ALIGN2_PHYLIP_MO" );

        system( "mv", "ALIGN2_PHYLIP_MO.dist", "DIST_TO_QUERY" )
        && &dieWithUnexpectedError( "$!" );
       
        if ( $parallel != 1 ) {
            $time_addingdists = time;
            &addDistsToQueryToPWDfile( $pwd_file,
                                       $temp_dir."/DISTs_TO_QUERY",
                                       $temp_dir."/PWD_INC_QUERY",
                                       $complete_names[ $jj ] );
       

            $time_addingdists = time - $time_addingdists;   
            $time_addingdistsT += $time_addingdists;
        }
        &addDistsToQueryToPWDfile( $nbd_file,
                                   $temp_dir."/DIST_TO_QUERY",
                                   $temp_dir."/NBD_INC_QUERY",
                                   $complete_names[ $jj ] );
        
    }

    if ( $mode == 2 ) {  
        $current_dir = `pwd`;
        $current_dir =~ s/\s//;
        chdir ( $temp_dir ) 
        || &dieWithUnexpectedError( "Could not chdir to \"$temp_dir\"" );

    }


    if ( $parallel != 1 ) { 
        unlink( $temp_dir."/MAKETREEOUT".$TREE_FILE_SUFFIX );
    }

    $time_tree_calc = time;

    # This calculates the trees
    # -------------------------

    if ( $mode == 1 || $mode == 2 ) {

        if ( $mode == 1 ) { 

            &executeNeighbor( $temp_dir."/NBD_INC_QUERY",
                              0,
                              0,
                              0,
                              1 );

            unlink( "outfile" );
            system( "mv", "outtree", "NBD_NJ_TREE" )
            && &dieWithUnexpectedError( "$!" );
            if ( $parallel != 1 ) {   
                &executeMakeTree( $options_for_makeTree,
                                  $temp_dir."/PWD_INC_QUERY",
                                  $temp_dir."/MAKETREEOUT".$TREE_FILE_SUFFIX,
                                  $temp_dir."/maketree_tempdir" );
            }

        }
        else {
            &executeNeighbor( $nbd_file,
                              0,
                              0,
                              0,
                              1 );

            unlink( "outfile" );
            system( "mv", "outtree", "NBD_NJ_TREE" )
            && &dieWithUnexpectedError( "$!" ); 

            &executeMakeTree( $options_for_makeTree,
                              $pwd_file,
                              $temp_dir."/MAKETREEOUT".$TREE_FILE_SUFFIX,
                              $temp_dir."/maketree_tempdir" );
         
        }

        chdir( $current_dir ) 
        || &dieWithUnexpectedError( "Could not chdir to \"$current_dir\"" );   


    }
    elsif ( $mode == 3 || $mode == 4 ) {
        &executeMakeTree( $options_for_makeTree,
                          $temp_dir."/ALIGN2",
                          $temp_dir."/MAKETREEOUT".$TREE_FILE_SUFFIX,
                          $temp_dir."/maketree_tempdir" );

        unlink( $temp_dir."/MAKETREEOUT".$ALIGN_FILE_SUFFIX );
    }
  

    $time_tree_calc = time - $time_tree_calc;
    $time_tree_calcT += $time_tree_calc;

    if  ( $keep == 1 ) {
        
        system( "cp", $temp_dir."/MAKETREEOUT".$TREE_FILE_SUFFIX,            $outfile.$TREE_FILE_SUFFIX );
        system( "cp", $temp_dir."/MAKETREEOUT".$LOG_FILE_SUFFIX,             $outfile.$LOG_FILE_SUFFIX );
        system( "cp", $temp_dir."/MAKETREEOUT".$MULTIPLE_TREES_FILE_SUFFIX,  $outfile.$MULTIPLE_TREES_FILE_SUFFIX );
        if ( $mode == 1 || $mode == 2 ) {
            system( "cp", $temp_dir."/NBD_NJ_TREE", $outfile."-NJ".$TREE_FILE_SUFFIX );
        }
        
    }

    unlink( $temp_dir."/ALIGN2" );

    $multiple_trees_file    = $temp_dir."/MAKETREEOUT".$MULTIPLE_TREES_FILE_SUFFIX;
    $maketree_out_tree_file = $temp_dir."/MAKETREEOUT".$TREE_FILE_SUFFIX;
    $distance_matrix_file   = $temp_dir."/MAKETREEOUT".$SUFFIX_PWD_NOT_BOOTS;

    
    if ( $mode == 1 || $mode == 3 ) {
        $query_name = $complete_names[ $jj ];
    }
 
    $options_for_DoRIO            = ""; 

    # This will result in saving of the annotated consenus tree:
    # ----------------------------------------------------------
    if ( $safe_nhx == 1 ) {
        my $number = $jj + 1;
        if ( @complete_names > 1 ) {
            $outfile_annot_nhx_tree = $outfile.$ADDITION_FOR_RIO_ANNOT_TREE."-".$number.$TREE_FILE_SUFFIX;
        }
        else {
            $outfile_annot_nhx_tree = $outfile.$ADDITION_FOR_RIO_ANNOT_TREE.$TREE_FILE_SUFFIX;
        }
    }



    if ( $sort > 2 ) {
        if ( $mode == 3 || $mode == 4 ) {
            $options_for_DoRIO .= " D=".$distance_matrix_file;
        }
        elsif ( $mode == 1 ) {
            $options_for_DoRIO .= " d=".$temp_dir."/DIST_TO_QUERY";
        }
        elsif ( $mode == 2 ) {
            $options_for_DoRIO .= " D=".$nbd_file;
        }
    }
    $options_for_DoRIO .= " M=".$multiple_trees_file;
    $options_for_DoRIO .= " 'N=".$query_name."'";
    $options_for_DoRIO .= " S=".$species_tree_file;
    $options_for_DoRIO .= " O=".$dorio_outfile;
    $options_for_DoRIO .= " P=".$sort;
    $options_for_DoRIO .= " L=".$t_orthologs;
    $options_for_DoRIO .= " B=".$t_sn;
    $options_for_DoRIO .= " U=".$t_orthologs_dc;
    $options_for_DoRIO .= " X=".$warn_more_than_one_ortho;
    $options_for_DoRIO .= " Y=".$warn_no_orthos;
    $options_for_DoRIO .= " Z=".$warn_one_ortho;
    
    if ( $mode == 1 || $mode == 2 ) {
        $options_for_DoRIO .= " T=".$temp_dir."/NBD_NJ_TREE";
        $options_for_DoRIO .= " t=".$maketree_out_tree_file;
    }
    elsif ( $mode == 3 || $mode == 4 ) {
        if ( $safe_nhx == 1 ) { # Added 09/04/03.
            $options_for_DoRIO .= " T=".$maketree_out_tree_file;
        }
    }    

    if ( $safe_nhx == 1 ) {
        $options_for_DoRIO .= " I";
    }
    if ( $output_ultraparalogs == 1 ) {
        $options_for_DoRIO .= " p";
        $options_for_DoRIO .= " v=".$t_ultra_paralogs;
    }

    $time_rio  = time;

    &executeDoRIO( $options_for_DoRIO );

    $time_rio  = time - $time_rio;
    $time_rioT += $time_rio;

    unless ( ( -s $dorio_outfile ) && ( -f $dorio_outfile ) && ( -T $dorio_outfile ) ) {
        close( OUT );
        unlink( $outfile );
        &dieWithUnexpectedError( "failure during execution of RIO (no output generated)" );
    }

    if ( $safe_nhx == 1 ) {
        system( "mv",
                 $temp_dir."/".$DO_RIO_TEMP_OUTFILE.$ADDITION_FOR_RIO_ANNOT_TREE.$TREE_FILE_SUFFIX,
                 $outfile_annot_nhx_tree ) 
        && &dieWithUnexpectedError( "$!" );
    }


    open( IN, "$dorio_outfile" )
    || &dieWithUnexpectedError( "Cannot open file \"$dorio_outfile\"" );

    $saw_distance_values          = 0;
    $saw_ultra_paralogs           = 0;
    $printed_ultra_paralogs       = 0;
    $print_header_for_orthologies = 1;
    $print_header_for_s_paralogs  = 1;
     
   


    # This generates the report
    # -------------------------

    W: while ( $return_line = <IN> ) {

        if ( $return_line =~ /distance values:/i ) {
            $saw_distance_values = 1;
            &printTitleForDistanceValues();
        }
        elsif ( $return_line =~ /ultra paralogs/i ) {
            $saw_ultra_paralogs = 1;
        }
        elsif ( $return_line =~ /^mean bootstrap/i ) {
            &printMeanBootstraps();
        }
        elsif ( $return_line =~ /sort priority\s*:\s*(.+)/i ) {
            $sort_priority = $1;
        }
        elsif ( $return_line =~ /ext nodes\s*:\s*(.+)/i ) {
            $ext_nodes_in_trees_analyzed = $1 - 1; # One seq is query.
        }
        elsif ( $return_line =~ /bootstraps\s*:\s*(\S+)/i ) { 
            if ( $jj == @complete_names - 1 ) {
                $bootstraps = $1;
                if ( $output_HTML == 1 ) { 
                   $| = 1;
                }
                &printOptions();
                last W;
            }
        }
        elsif ( $saw_distance_values != 1 
        && $saw_ultra_paralogs != 1
        && $return_line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(\S*)/ ) {
            $ortho_name        = $1;
            $orthos            = $2;
            $subtree_neighbors = $3;
            $s_orthos          = $4;
            $dist              = $5;

            if ( $print_header_for_orthologies == 1 ) {
                &printHeaderForOrthologies();
                $print_header_for_orthologies = 0;
            }
            &printOrthologies();
        }
        elsif ( $saw_distance_values != 1
        && $saw_ultra_paralogs != 1
        && $return_line =~ /^\s*-\s*$/  ) {
            $ortho_name = "-";
            $orthos     = 0;
            $s_orthos   = 0;
            $dist       = 0;
            if ( $print_header_for_orthologies == 1 ) {
                &printHeaderForOrthologies();
                $print_header_for_orthologies = 0;
            }
            &printOrthologies();
        }
        elsif ( $output_ultraparalogs == 1
        && $saw_ultra_paralogs == 1
        && $return_line =~ /(\S+)\s+(\S+)\s+(\S+)/ ) {
            $s_para_name = $1;
            $s_paras     = $2; 
            $dist        = $3;
            if ( $print_header_for_s_paralogs == 1 ) {
                &printHeaderForSparalogs();
                $print_header_for_s_paralogs = 0;
            }
            &printUltraParlogs();
            $printed_ultra_paralogs = 1;
        }
        elsif ( $output_ultraparalogs == 1
        && $saw_ultra_paralogs == 1
        && $return_line =~ /^\s*-\s*$/ ) {
            &printNoUltraParalogs();
        }
        elsif ( $return_line =~ /Bootstraps/ ) {
            $saw_distance_values = 0;
        }
        elsif ( $saw_distance_values == 1 && $saw_ultra_paralogs != 1 ) {
            &printDistanceValues();
        }
       
    }
    close( IN );
    
} # End of for loop going through possible 
  # multiple matches to the same alignment/model.

if ( $output_HTML != 1 ) {
    close( OUT );
}

&cleanUpTempDir();

if ( $output_HTML != 1 ) {
    print( "\n\nrio.pl successfully terminated.\nOutput written to: $outfile\n\n" );
}

exit( 0 );









# ===========================================================
#                                                     Methods
# -----------------------------------------------------------




# ----------------------------------------------------------- 
#                                       Parallization related
# -----------------------------------------------------------



# Last modified: 02/02/02
sub readInNodesList {

    &testForTextFilePresence( $NODE_LIST );

    open( NIN, "$NODE_LIST" ) || &dieWithUnexpectedError( "Cannot open file \"$NODE_LIST\"" );
    
    while ( <NIN> ) {
        if ( $_ =~ /(\S+)/ ) {
            push( @nodelist, $1 );
        }
    }
    close( NIN );
    return;
}



# Last modified: 02/02/02
sub pingNodes {
    my @temp_node_list = ();
    my $p = Net::Ping->new( "tcp", 2 ); # or "udp"
    my $n = "";

    foreach $n ( @nodelist ) {
        if ( defined( $p->ping( $n ) ) ) {
            push( @temp_node_list, $n );
        }
    }
    @nodelist = ();
    @nodelist = @temp_node_list;
    return;

}




# ----------------------------------------------------------- 
#                                              Output related
# -----------------------------------------------------------


# Last modified: 03/07/01
sub printHeader {
    
    if ( $output_HTML != 1 ) {
        print OUT "RIO - Resampled Inference of Orthologs\n";
        print OUT "Version: $VERSION\n";
        print OUT "------------------------------------------------------------------------------\n";

        print OUT "Pfam alignment file                     : $alignment\n";
        if ( $mode == 3 ) {
            print OUT "Pfam alignment description              : ".&getDescriptionFromPfam( $alignment )."\n";
        }
        if ( $mode == 1 || $mode == 2 ) {
            print OUT "Bootstrapped pairwise distances file    : $pwd_file\n";
            print OUT "Not bootstrapped pairwise distances file: $nbd_file\n";
            print OUT "Bootstrap positions file                : $bsp_file\n";
        }
        if ( $mode == 1 || $mode == 3 ) { 
            if ( $seed_aln_for_hmmbuild ne "" ) {
                print OUT "HMM                                     : built based on $seed_aln_for_hmmbuild\n";
            }
            elsif ( $hmm_name ne "" ) {
                print OUT "HMM                                     : $hmm_name\n";
            }
            else {
                print OUT "HMM                                     : $hmm_file\n";
            }
            print OUT "Query file                              : $seqX_file\n";
        }
        print OUT "==============================================================================\n\n";
    }
    
} ## printHeader




# Last modified: 03/07/01
sub printHeaderForOrthologies {
    
    if ( $output_HTML != 1 ) {
        if ( $jj > 0 ) {
            print OUT "\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
        }

        print OUT "Query         : $query_name\n\n";
    
        if ( @complete_names > 1 ) {
            my $size = @complete_names;
            my $number = $jj + 1;
            print OUT "More than one region of the query were aligned to the profile HMM.\n";
            print OUT "This is for domain #$number out of $size.\n\n";
        }
    
        print OUT "Number (in %) of observed orthologies (o), \"subtree-neighborings\" (n),\n"; 
        print OUT "and super-orthologies (s) to query in bootstrapped trees, evolutionary\n";
        print OUT "distance to query (as average number of aa replacements per residue):\n\n";
        if ( $long_output != 1 ) {
            print OUT "Sequence                Description                        o[%] n[%] s[%]  distance\n";
            print OUT "--------                -----------                        ---- ---- ----  --------\n";
        }
        else {
            print OUT "Sequence                Description                                                                o[%] n[%] s[%]  distance\n";
            print OUT "--------                -----------                                                                ---- ---- ----  --------\n";
        }
    }
    else {
        if ( $jj > 0 ) {
            print "</TABLE>\n";
            print "<P> &nbsp </P>\n";
            print "<HR NOSHADE COLOR=\"#CCCCCC\">\n";
            print "<P> &nbsp </P>\n";
        }
    
        if ( @complete_names > 1 ) {
            my $size = @complete_names;
            my $number = $jj + 1;
            print "<P>More than one region of the query were aligned to the profile HMM. \n";
            print "This is for domain #$number out of $size.</P>\n";
        }
        print "<P>Query         : $query_name</P>\n";
        print "<H4 class = \"title\">Orthologies, subtree-neighborings, super-orthologies</H4>\n";

        print "<P>Number (in %) of observed orthologies (o), \"subtree-neighborings\" (n), \n"; 
        print "and super-orthologies (s) to query in bootstrapped trees, evolutionary \n";
        print "distance to query (as average number of aa replacements per residue):</P>\n";
        if ( $ortho_name ne "-" ) {
        print "<TABLE BORDER=\"0\" CELLPADDING=\"1\">\n";
        
        print "<TR VALIGN=\"TOP\"> <TD NOWRAP> <B>Sequence</B> </TD><TD NOWRAP> <B>Description</B> </TD><TD NOWRAP ALIGN=\"RIGHT\"> <B>o[%]</B> </TD><TD NOWRAP ALIGN=\"RIGHT\"> <B>n[%]</B> </TD><TD NOWRAP ALIGN=\"RIGHT\"> <B>s[%]</B> </TD><TD NOWRAP> &nbsp <B>distance</B> </TD> </TR>\n";
        }
    }

} ## printHeaderForOrthologies



# Last modified: 10/15/01
sub printHeaderForSparalogs {

    if ( $output_HTML != 1 ) {
        print OUT "\nUltra-paralogs\n";
        print OUT "--------------\n";
        print OUT "Number (in %) of observed ultra-paralogies (up) to query\n";
        print OUT "in bootstrapped trees, evolutionary distance to query (as average number\n";
        print OUT "of aa replacements per residue):\n\n";
        if ( $long_output != 1 ) {
            print OUT "Sequence                Description                       up[%]  distance\n";
            print OUT "--------                -----------                       -----  --------\n";
        }
        else {
            print OUT "Sequence                Description                                                               up[%]  distance\n";
            print OUT "--------                -----------                                                               -----  --------\n";
        }
    }
    else {
        print "<H4 class = \"title\">Ultra-paralogs</H4>\n";
        print "<P>Number (in %) of observed ultra-paralogies (up) to query \n";
        print "in bootstrapped trees, evolutionary distance to query (as average number \n";
        print "of aa replacements per residue):</P>\n";
        print "<TABLE BORDER=\"0\" CELLPADDING=\"1\">\n";
        print "<TR VALIGN=\"TOP\"> <TD NOWRAP> <B>Sequence</B> </TD><TD NOWRAP> <B>Description</B> </TD><TD NOWRAP ALIGN=\"RIGHT\"> <B>up[%]</B> </TD><TD NOWRAP> &nbsp <B>distance</B> </TD> </TR>\n";
       
   }
   
} ## printHeaderForSparalogs



# Last modified: 03/07/01
sub printOrthologies {
    my @cut   = ();
    my $i     = 0;
    my $descp = "";
    $orthos   = &roundToInt( $orthos );
    $s_orthos = &roundToInt( $s_orthos );

    if ( $sort > 10 ) {
        $subtree_neighbors = &roundToInt( $subtree_neighbors );
    }
    
    if ( ( $description == 1 || $complete_description == 1 )
    && $ortho_name ne "-" ) {
        
        if ( $non_sp != 1 ) {
            if ( &startsWithSWISS_PROTname( $ortho_name ) ) {
                $descp = &getDescriptionFromSWISSPROT_ACDEOSfile( $SWISSPROT_ACDEOS_FILE, $ortho_name );
            }
            else {
                $descp = "-";
            }
        }
        else {
            if ( &startsWithSWISS_PROTname( $ortho_name ) ) {
                $descp = &getDescriptionFromSWISSPROT_ACDEOSfile( $SWISSPROT_ACDEOS_FILE, $ortho_name );
            }
            else {
                $descp = &getDescriptionFromTrEMBL_ACDEOSfile( $TREMBL_ACDEOS_FILE, $ortho_name );
            }
        }

        if ( $output_HTML != 1 ) {
            if ( $long_output == 1 ) {
                @cut = &cutDescription( $descp, 73 );
            }
            else {
                @cut = &cutDescription( $descp, 33 );
            }
            $descp = $cut[ 0 ];
        }
    }
    if ( $descp eq "" ) {
        $descp = "-";
    }

    if ( $output_HTML != 1 ) {

	if ( $ortho_name eq "-" ) {
	    print OUT "\nNO ORTHOLOGS in alignment with the current thresholds for output\n";
	}
	elsif ( $dist ne "-" ) {
	    if ( $long_output == 1 ) {
		print OUT sprintf "%-24.24s%-74.74s%5s%5s%5s%10.6f", $ortho_name,$descp,$orthos,$subtree_neighbors,$s_orthos,$dist;
	    }
	    else {
		print OUT sprintf "%-24.24s%-34.34s%5s%5s%5s%10.6f", $ortho_name,$descp,$orthos,$subtree_neighbors,$s_orthos,$dist;
	    }
	}
	else {
	    if ( $long_output == 1 ) {
		print OUT sprintf "%-24.24s%-74.74s%5s%5s%5s%10.10s", $ortho_name,$descp,$orthos,$subtree_neighbors,$s_orthos,$dist;
	    }
	    else {
		print OUT sprintf "%-24.24s%-34.34s%5s%5s%5s%10.10s", $ortho_name,$descp,$orthos,$subtree_neighbors,$s_orthos,$dist;
	    }
	}
	if ( $complete_description == 1 ) {
	    for ( $i = 1; $i < @cut; ++$i ) {
		print OUT "\n";
		if ( $long_output == 1 ) {
		    print OUT sprintf "                        %-74.74s", $cut[ $i ];
		}
		else {
		    print OUT sprintf "                        %-34.34s", $cut[ $i ];
		}
	    }
        } 
	print OUT "\n";
    }
    else {
        if ( $ortho_name eq "-" ) {
            print "<H4 class = \"warnings\">NO ORTHOLOGS in alignment with the current thresholds for output</H4>\n";
        }
        else {
            $ortho_name = &replaceNameWithLinkToExpasy( $ortho_name );
            print "<TR VALIGN=\"TOP\"> <TD NOWRAP> $ortho_name </TD><TD> $descp </TD><TD NOWRAP ALIGN=\"RIGHT\"> $orthos </TD><TD NOWRAP ALIGN=\"RIGHT\"> $subtree_neighbors </TD><TD NOWRAP ALIGN=\"RIGHT\"> $s_orthos </TD><TD NOWRAP> &nbsp $dist </TD> </TR>\n";
        }
    }

} ## printOrthologies



sub replaceNameWithLinkToExpasy {
    my $name = $_[ 0 ];

    if ( $name =~ /(.+)_(.+)\/(.+)/ ) {
        my $desc    = $1;
        my $spec    = $2;
        my $numbers = $3;
        if ( length( $desc ) <= 4 ) {
            $name = "<A HREF=\"".$EXPASY_SPROT_SEARCH_DE.$desc."_".$spec."\" TARGET=\"_blank\">".$desc."_".$spec."</A>\/".$numbers;
        }
        else {
            $name = "<A HREF=\"".$EXPASY_SPROT_SEARCH_AC.$desc."\" TARGET=\"_blank\">".$desc."</A>_".$spec."\/".$numbers;
        }
    }
    
    return $name;
   
} ## replaceNameWithLinkToExpasy




# Last modified: 10/15/01
sub printUltraParlogs {
    my @cut   = ();
    my $i     = 0;
    my $descp = "";
    $s_paras  = &roundToInt( $s_paras );

    if ( ( $description == 1 || $complete_description == 1 )
    && $s_para_name ne "-" ) {
        
        if ( $non_sp != 1 ) {
            if ( &startsWithSWISS_PROTname( $s_para_name ) ) {
                $descp = &getDescriptionFromSWISSPROT_ACDEOSfile( $SWISSPROT_ACDEOS_FILE, $s_para_name );
            }
            else {
                $descp = "-";
            }
        }
        else {
            if ( &startsWithSWISS_PROTname( $s_para_name ) ) {
                $descp = &getDescriptionFromSWISSPROT_ACDEOSfile( $SWISSPROT_ACDEOS_FILE, $s_para_name );
            }
            else {
                $descp = &getDescriptionFromTrEMBL_ACDEOSfile( $TREMBL_ACDEOS_FILE, $s_para_name );
            }
        }
        
        if ( $output_HTML != 1 ) {
            if ( $long_output == 1 ) {
                @cut = &cutDescription( $descp, 73 );
            }
            else {
                @cut = &cutDescription( $descp, 33 );
            }
            $descp = $cut[ 0 ];
        }
    }
    if ( $descp eq "" ) {
        $descp = "-";
    }

    if ( $output_HTML != 1 ) {

        if ( $dist ne "-" ) {
            if ( $long_output == 1 ) {
                print OUT sprintf "%-24.24s%-74.74s%5s%10.6f", $s_para_name,$descp,$s_paras,$dist;
            }
            else {
                print OUT sprintf "%-24.24s%-34.34s%5s%10.6f", $s_para_name,$descp,$s_paras,$dist;
            }
        }
        else {
            if ( $long_output == 1 ) {
                print OUT sprintf "%-24.24s%-74.74s%5s%10.10s", $s_para_name,$descp,$s_paras,$dist;
            }
            else {
                print OUT sprintf "%-24.24s%-34.34s%5s%10.10s", $s_para_name,$descp,$s_paras,$dist;
            }
        }
        if ( $complete_description == 1 ) {
            for ( $i = 1; $i < @cut; ++$i ) {
                print OUT "\n";
                if ( $long_output == 1 ) {
                    print OUT sprintf "                        %-74.74s", $cut[ $i ];
                }
                else {
                    print OUT sprintf "                        %-34.34s", $cut[ $i ];
                }
            }
        } 
        print OUT "\n";

    }
    else {
        $s_para_name = &replaceNameWithLinkToExpasy( $s_para_name );
        print "<TR VALIGN=\"TOP\"> <TD NOWRAP> $s_para_name </TD><TD> $descp </TD><TD NOWRAP ALIGN=\"RIGHT\"> $s_paras </TD><TD NOWRAP> &nbsp $dist </TD> </TR>\n";
    }
   
} ## printUltraParlogs



sub printNoUltraParalogs {
    if ( $output_HTML != 1 ) { 
        print OUT "\nUltra-paralogs\n";
        print OUT "--------------\n";
        print OUT "\nNO ULTRA-PARALOGS in alignment with the current threshold of $t_ultra_paralogs%\n";
    }
    else {
        print "<H4 class = \"title\">Ultra-paralogs</H4>\n";
        print "<H4 class = \"warnings\">NO ULTRA-PARALOGS in alignment with the current threshold of $t_ultra_paralogs%</H4>\n";
    }
} ## printNoUltraParalogs



# Called by method "printOrthologies".
# Last modified: 02/27/01
sub cutDescription {
    my $line = $_[ 0 ];
    my $size = $_[ 1 ];
    my @cut  = ();
    my $i    = 0;
   
    while ( ( length( $line ) ) > $size ) {
        $cut[ $i++ ] = substr( $line, 0, $size );
        $line = substr( $line, $size );
    }
    $cut[ $i++ ] = $line;
    return @cut;
} ## cutDescription




# Last modified: 02/27/01
sub printTitleForDistanceValues {
    if ( $output_HTML != 1 ) {
        if ( $mode == 1 || $mode == 2 ) {
            print OUT "\n\nDistance values (based on NJ tree of original alignment)\n";
            print OUT "--------------------------------------------------------\n";
        }
        elsif ( $mode == 3 || $mode == 4 ) {
            print OUT "\n\nDistance values (based on ML branch length values on consensus tree)\n";
            print OUT "--------------------------------------------------------------------\n";
        }
    }
    else {
        print "<H4 class = \"title\">Distance values (based on NJ tree of original alignment)</H4>\n";
    }
    
} ## printTitleForDistanceValues




# Last modified: 02/27/01
sub printDistanceValues {
    if ( $output_HTML != 1 ) {
        print OUT "$return_line";
    }
    else {
        chomp( $return_line );
        if ( $return_line =~ /WARNING/ ) {
            $return_line =~ s/\+\/-/ &plusmn /;
            $return_line =~ s/\*/ &times /;
            print "<H4 class = \"warnings\">$return_line</H4>\n";
        }
        elsif ( $return_line =~ /lca\s+is/i ) {
            print "<P class = \"nomargins\">$return_line</P>\n";
        }
        elsif ( $return_line =~ /orthologous/i ) {
            print "<P class = \"nomargins\">$return_line</P>\n";
        }
        elsif ( $return_line =~ /distance\s+of\s+query/i ) {
            print "<TABLE BORDER=\"0\" CELLPADDING=\"1\">\n";
        }
        if ( $return_line =~ /(.+)=(.+)/ ) {
            print "<TR VALIGN=\"TOP\"><TD>$1</TD><TD> = $2</TD></TR>\n";
        }
        if ( $return_line =~ /sum\s+/i || $return_line =~ /distance\s+of\s+ortholog\s+to\s+LCA/i ) {
            print "</TABLE>\n";
        }
    }
} ## printDistanceValues




# Last modified: 02/27/01
sub printMeanBootstraps {
    if ( $output_HTML != 1 ) {
        print OUT "\n\n$return_line";
    }
    else {
        chomp( $return_line );
        $return_line =~ s/\+\/-/ &plusmn /;
        print "</TABLE>\n";
        print "<P>$return_line</P>\n";
    }
} ## printMeanBootstraps




# Last modified: 02/12/02
sub printOptions {

    if ( $output_HTML != 1 ) {
        print OUT "\n\n\n==============================================================================\n";
	    if ( $number_of_seqs_in_aln >= $MIN_NUMBER_OF_SEQS_IN_ALN ) {
	        print OUT "RIO options\n";
	        print OUT "-----------\n";
	        print OUT "Mode                                                      : ";
	        if ( $mode == 1 ) { 
                print OUT "precalc. pwd files with alignment not containing query (1)\n";
	        }
	        elsif ( $mode == 2 ) { 
                print OUT "precalc. pwd files with alignment containing query (2)\n";
	        }
	        elsif ( $mode == 3 ) { 
                print OUT "alignment not containing query (3)\n";
	        }
	        elsif ( $mode == 4 ) { 
		    print OUT "alignment containing query (4)\n";
	        }
	        print OUT "Bootstraps                                                : $bootstraps\n";
	        print OUT "Species tree                                              : $species_tree_file\n";
	        if ( $safe_nhx == 1 ) {
		    if ( $mode == 3 || $mode == 4 ) { 
		        if ( @complete_names > 1 ) {
			    $outfile_annot_nhx_tree =~ s/-\d+\.nhx/-X.nhx/;
			    print OUT "Saved annotated consensus trees (ML branch lengths)       : $outfile_annot_nhx_tree\n";
		        }
		        else {
			    print OUT "Saved annotated consensus tree (ML branch lengths)        : $outfile_annot_nhx_tree\n";
		        }
		    }
		    elsif ( $mode == 1 || $mode == 2 ) { 
		        if ( @complete_names > 1 ) {
			    $outfile_annot_nhx_tree =~ s/-\d+\.nhx/-X.nhx/;
			    print OUT "Saved annotated NJ trees (based on original alignment)    : $outfile_annot_nhx_tree\n";
		        }
		        else {
			    print OUT "Saved annotated NJ tree (based on original alignment)     : $outfile_annot_nhx_tree\n";
		        }
		    }
	        }
	        print OUT "Threshold for output for orthologies          (L=)        : $t_orthologs\n";
            print OUT "Threshold for output for subtree-neighborings (B=)        : $t_sn\n";
	        print OUT "Threshold for distance calc for orthologies   (U=)        : $t_orthologs_dc\n";

	        print OUT "When to generate warnings:\n";
	        print OUT "More than one ortholog:  diff. in standard deviations (X=): $warn_more_than_one_ortho\n";
	        print OUT "No  orthologs         :  diff. in standard deviations (Y=): $warn_no_orthos\n";
	        print OUT "One ortholog          :  factor                       (Z=): $warn_one_ortho\n";
	        if ( $output_ultraparalogs == 1 ) {
                print OUT "Output ultra-paralogs (p)\n";
                print OUT "Threshold for ultra-paralogies (v=)                       : $t_ultra_paralogs\n";
	        }
	        print OUT "Sort priority: $sort_priority\n";
	    }

	    print OUT "\nOptions for the calculation of the phylgenetic trees\n";
	    print OUT "----------------------------------------------------\n";
	    if ( $mode == 1 ) {
	        print OUT "Model for pairwise distance calculations           : $matrix\n";
	    }
	    elsif ( $mode == 3 || $mode == 4 ) {
	        print OUT "Model for pairwise dist and ML branch length calc. : $matrix\n";
	    }
	    if ( $mode == 1 || $mode == 3 || $mode == 4 ) {
	        print OUT "Columns in alignment used for tree calc            : $length_of_alignment\n";
	        print OUT "Columns in original alignment                      : $length_of_orig_alignment\n";
	    }
	    print OUT "Sequences in alignment used for trees (incl query) : $number_of_seqs_in_aln\n";

	    if ( $mode == 3 || $mode == 4 ) { 
	        print OUT "Removed non-SWISS-PROT sequences                   : ";
	        if ( $non_sp == 1 ) {
                print OUT "no\n";
	        }
	        else {
                print OUT "yes\n";
	        }
	        if ( $non_sp == 1 ) {
                print OUT "Removed \"TrEMBL fragments\"                         : ";
		        if ( $no_frags == 1 ) {
		            print OUT "yes\n";
		        }
		        else {
		            print OUT "no\n";
		        }
	        }
	    }
	    if ( $mode == 1 || $mode == 2 ) {
	        print OUT "Prgrm to calc. branch lengths for distance values  : PHYLIP NEIGHBOR (NJ)\n";
	    }
	    elsif ( $mode == 3 || $mode == 4 ) {
	        print OUT "Prgrm to calc branch lengths for distance values   : TREE-PUZZLE\n";
	    }
	    if ( $seed_aln_for_hmmbuild ne "" ) {
	        print OUT "HMM was built with hmmbuild using options          : $command_line_for_hmmbuild\n";
	    }
	    if ( ( $mode == 3 || $mode == 4 ) && $species_names_file =~ /\S/ ) {
	        print OUT "File listing species used for tree calculation (G=): $species_names_file\n";
	    }
	    print OUT "Seed for random number generator                   : $seed_for_makeTree\n";
	    print OUT "Options for makeTree                               : $options_for_makeTree\n";    

	    $time_total = time - $time_total;

	    print OUT "\nTime and date\n";
	    print OUT "-------------\n";
	    if ( $mode == 1 ) {
	        print OUT "Time requirement dqo puzzle          : $time_dqopuzzleT s\n";
	    }

	    print OUT "Time requirement for tree calculation: $time_tree_calcT s\n";
	    print OUT "Time requirement for SDI and RIO     : $time_rioT s\n";
	    print OUT "Total time requirement               : $time_total s\n";
	    print OUT "Date started                         : $start_date";
	    print OUT ( "Date finished                        : ".`date` );

	    print OUT "\nCommand line\n";
	    print OUT "------------\n";
	    print OUT "$command_line\n";
	    if ( $parallel == 1 ) {
	        print OUT "\nProcessors used: @nodelist\n";
	    }
    }
    else {
        if ( $printed_ultra_paralogs == 1 ) {
            print "</TABLE>\n";
        }
        if ( $species_tree_file =~ /.+\/(.+)/ ) {
            $species_tree_file = $1;
        }
        print "<H4 class = \"title\">Options</H4>\n";
        print "<TABLE BORDER=\"0\" CELLPADDING=\"1\">\n";
        print "<TR VALIGN=\"TOP\"><TD> Bootstraps: </TD><TD> $bootstraps </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Species tree: </TD><TD> $species_tree_file </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Threshold for output for orthologies: </TD><TD> $t_orthologs </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Threshold for output for subtree-neighborings: </TD><TD> $t_sn </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Threshold for distance calc for orthologies: </TD><TD> $t_orthologs_dc </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> When to generate warnings </TD><TD> </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> More than one ortholog [diff in standard deviations]: </TD><TD> $warn_more_than_one_ortho </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> No  orthologs [diff in standard deviations]: </TD><TD> $warn_no_orthos </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> One ortholog [factor]: </TD><TD> $warn_one_ortho </TD></TR>\n";
        if ( $output_ultraparalogs == 1 ) {
            print "<TR VALIGN=\"TOP\"><TD> Output ultra-paralogs </TD><TD> </TD></TR>\n";
            print "<TR VALIGN=\"TOP\"><TD> Threshold for ultra-paralogies: </TD><TD> $t_ultra_paralogs </TD></TR>\n";
        }
        print "<TR VALIGN=\"TOP\"><TD> Sort priority: </TD><TD> $sort_priority </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Model for pairwise distance calculations:</TD><TD> $matrix </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Columns in alignment used for tree calc: </TD><TD> $length_of_alignment </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Columns in original alignment: </TD><TD> $length_of_orig_alignment </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Sequences in alignment used for trees (incl query): </TD><TD> $number_of_seqs_in_aln </TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Seed for random number generator: </TD><TD> $seed_for_makeTree </TD></TR>\n";    
        print "</TABLE>\n";
    
        $time_total = time - $time_total;

        print "<P> &nbsp </P>\n";
        print "<TABLE BORDER=\"0\" CELLPADDING=\"1\">\n";
        print "<TR VALIGN=\"TOP\"><TD> Time requirement: </TD><TD> $time_total s</TD></TR>\n";
        print "<TR VALIGN=\"TOP\"><TD> Date started: </TD><TD> $start_date </TD></TR>\n";
        print ( "<TR VALIGN=\"TOP\"><TD> Date finished: </TD><TD> ".`date`." </TD></TR>\n" ); 
        if ( $parallel == 1 ) {
            print "<TR VALIGN=\"TOP\"><TD> Number of processors used: </TD><TD> ".scalar( @nodelist )." </TD></TR>\n";
        }
        print "</TABLE>\n";
    }

} ## printOptions










# -----------------------------------------------------------
#                                 Execution of other programs
# -----------------------------------------------------------





# Two arguments:
# 1. seed
# 2. outfile
# Returns the options used.
# Last modified: 05/11/01
sub executeHmmbuild {

    my $seed    = $_[ 0 ];
    my $outfile = $_[ 1 ];
    my $options = "";

    &testForTextFilePresence( $seed );
    
    $options = getHmmbuildOptionsFromPfam( $seed );

    $options =~ s/-f//;
    $options =~ s/-g//;
    $options =~ s/-s//;
    $options =~ s/-F//;
    $options =~ s/-A//;
    $options =~ s/-o\s+\S+//;
    $options =~ s/(\s|^)[^-]\S+/ /g;

    if ( $options =~ /--prior/ ) {
        my $basename = basename( $seed );
        $basename .= ".PRIOR";
        $options =~ s/--prior/--prior $PRIOR_FILE_DIR$basename/;
    }

    # Remove for versions of HMMER lower than 2.2.
    if ( $options =~ /--informat\s+\S+/ ) {
        $options =~ s/--informat\s+\S+/--informat SELEX/; 
    }
    else {
        $options = "--informat SELEX ".$options;
    }

    system( "$HMMBUILD $options $outfile $seed" )
    && &dieWithUnexpectedError( "Could not execute \"$HMMBUILD $options $outfile $seed\"" ); 
    return $options;

} ## executeHmmbuild.




# One argument:
# Pfam align name.
# Last modified: 02/26/01
sub getHmmbuildOptionsFromPfam {

    my $infile      = $_[ 0 ];
    my $return_line = "";
    my $result      = "";

    &testForTextFilePresence( $infile );

    open( GHO, $infile ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
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




# Purpose. Aligns a FASTA file to a Pfam alignment using an HMM profile.
# Five arguemnts:
# 1. Pfam flat file name
# 2. Name of FASTA file to append
# 3. HMM profile file name
# 4. outputfile name
# 5. 1 use --mapali, --withali otherwise (in hmmalign) 
# Returns 1 if successful, -1 if no alignment was made because
# E value of HMMSEARCH output was larger than $E_VALUE_THRESHOLD.
# Last modified: 07/11/01
sub alignWithHmmalign {
    my $alignment           = $_[ 0 ];
    my $query               = $_[ 1 ];
    my $hmm                 = $_[ 2 ];
    my $outfile             = $_[ 3 ];
    my $use_mapali          = $_[ 4 ];
    my $E                   = 2000;
    my $ali = "--withali";
    
    if ( $use_mapali == 1 ) {
        $ali = "--mapali";
    }

    &testForTextFilePresence( $alignment );
    &testForTextFilePresence( $query );
    &testForTextFilePresence( $hmm );

    system( "$HMMSEARCH $hmm $query > $temp_dir/HMMSEARCHOUT" )
    &&  &dieWithUnexpectedError( "Could not execute \"$HMMSEARCH $hmm $query > $temp_dir/HMMSEARCHOUT\"" );

    
   
    $E = &getEvalue( "$temp_dir/HMMSEARCHOUT" );
    if ( $E == 2000 ) {
        &dieWithUnexpectedError( "No E-value found in \"$temp_dir/HMMSEARCHOUT\"" );
    }
    elsif ( $E > $E_VALUE_THRESHOLD ) {
        unlink( "$temp_dir/HMMSEARCHOUT" );
        return ( -1 );
    }

    system( "$P7EXTRACT -d $temp_dir/HMMSEARCHOUT > $temp_dir/GDF" )
    && &dieWithUnexpectedError( "Could not execute \"$P7EXTRACT -d $temp_dir/HMMSEARCHOUT > $temp_dir/GDF\"" );

    
    system( "$MULTIFETCH  -d -g $query $temp_dir/GDF > $temp_dir/MULTIFETCHOUT" )
    && &dieWithUnexpectedError( "Could not execute \"$MULTIFETCH  -d -g $query $temp_dir/GDF > $temp_dir/MULTIFETCHOUT\"" );

    # Checks if score was too low to have made a reasonable alignment. 
    unless ( -s "$temp_dir/MULTIFETCHOUT" ) {
        unlink( "$temp_dir/HMMSEARCHOUT", "$temp_dir/GDF", "$temp_dir/MULTIFETCHOUT" );
        return ( -1 );
    }

    system( "$HMMALIGN -o $outfile $ali $alignment $hmm $temp_dir/MULTIFETCHOUT >/dev/null 2>&1" )
    && &dieWithUnexpectedError( "Could not execute \"$HMMALIGN -o $outfile $ali $alignment $hmm $temp_dir/MULTIFETCHOUT\"" );

    if ( unlink( "$temp_dir/HMMSEARCHOUT", "$temp_dir/GDF","$temp_dir/MULTIFETCHOUT" ) != 3 ) {
        &dieWithUnexpectedError( "Could not delete (a) file(s)" );
    }

    return 1; 
} ## alignWithHmmalign




# Gets the E value for complete sequences (score includes all domains)
# from a HMMSEARCH output file.
# One argument: the HMMSEARCH output file name
# Returns the E value, 2000 if no E value found
# Last modified: 07/11/01
sub getEvalue {

    my $infile      = $_[ 0 ];
    my $return_line = "";
    my $flag        = 0;
    my $E           = 2000;

    &testForTextFilePresence( $infile );

    open( E, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    while ( $return_line = <E> ) {
        
        # "Sequence    Description                                 Score    E-value  N"  
        if ( $return_line =~ /Sequence.+Description.+Score.+E.+value.+N/ ) { 
            $flag = 1; 
        }
        # "QUERY_HUMAN                                             657.4   1.3e-198   1"
        elsif ( $flag == 1 && $return_line =~ /\s+(\S+)\s+\d+\s*$/ ) { 
            $E = $1;
            close( E );
            return $E;
        }
        
    }
    close( E );
    return $E;

} ## getEvalue



# Four/Five arguments:
# 1. Number of bootstraps
# 2. bsp (bootstrap positions) file
# 3. Infile (alignment)
# 4. Outfile (bootstrapped according to bsp file)
# 5. Number of processors
# Last modified: 01/30/02
sub executeBootstrap_cz {
    my $boots       = $_[ 0 ];
    my $bsp_file    = $_[ 1 ];
    my $infile      = $_[ 2 ];
    my $outfile     = $_[ 3 ];
    my $processors  = $_[ 4 ];

    if ( defined( $processors ) && ( $processors > 1 ) ) {
        system( "$BOOTSTRAP_CZ $boots $infile $bsp_file $outfile $processors" )
        && &dieWithUnexpectedError( "Could not execute \"$BOOTSTRAP_CZ $boots $infile $bsp_file $outfile $processors\"" );

    }
    else {
        system( "$BOOTSTRAP_CZ $boots $infile $bsp_file $outfile" )
        && &dieWithUnexpectedError( "Could not execute \"$BOOTSTRAP_CZ $boots $infile $bsp_file $outfile\"" );
    }

} ## executeBootstrap_cz





# One argument:
# options for DoRIO.main.
# Last modified: 02/26/01
sub executeDoRIO {

    my $options = $_[ 0 ]; 
   
    system( "$DORIO $options >/dev/null 2>&1" )
    && &dieWithUnexpectedError( "Could not execute \"$DORIO $options\"" );

    return;

} ## executeDoRIO











# -----------------------------------------------------------
#                               These deal with the alignment
# -----------------------------------------------------------




# Counts sequences from a Pfam flat file or
# in a PHYLIP interleaved aligment.
# One arguments: Pfam flat file name.
# Returns the number of sequences.
# Last modified: 07/10/01
sub countSeqsInPfamAlign {
    my $infile             = $_[ 0 ];
    my $return_line        = "";
    my $saw_sequence_line  = 0;
    my $number_of_seqs     = 0;
  
    &testForTextFilePresence( $infile );

    open( C, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    while ( $return_line = <C> ) {

        if ( $saw_sequence_line == 1
        && !&containsPfamNamedSequence( $return_line )
        && !&isPfamCommentLine( $return_line ) ) {
            last;
        }
        if ( &isPfamSequenceLine( $return_line )
        && $return_line !~ /^\s*\d+\s+\d+/ ) { 
            if ( $saw_sequence_line == 0 ) {
                $saw_sequence_line = 1;
            }
            $number_of_seqs++;
        }  
    }
    close( C );
    return $number_of_seqs;

} ## countSeqsInPfamAlign




# This gets the complete name(s) of a sequence from a Pfam alignment.
# I.e. it adds "/xxx-xxx".
# 2 arguments:
# 1. Infile (alignment)
# 2. Name of query
# Returns a String-array of all the complete names found.
# Last modified: 03/04/01
sub getCompleteName {
    
    my $infile         = $_[ 0 ];
    my $query_name     = $_[ 1 ];
    my $return_line    = "";
    my @complete_names = ();
    my $complete_name  = "";
    my $i              = 0;
    
    &testForTextFilePresence( $infile );

    $query_name =~ s/\/.*//;

    open( INGCN, $infile ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    while ( $return_line = <INGCN> ) {
        if ( $return_line =~ /^\s*$query_name(\S+)\s+.+/ ) {
            $complete_name = $query_name.$1;
            if ( $i > 0 && $complete_names[ 0 ] eq $complete_name ) {
                # Now, we saw of all of them.
                last;
            }
            $complete_names[ $i++ ] = $complete_name;
        }
    }
    
    close( INGCN );
    return @complete_names;
} ## getCompleteName




# Removes sequences from a Pfam flat file.
# It can remove all sequences not from species listed in a species names file.
# It can remove all sequences which do not have a SWISS-PROT name (XXXX_XXXXX)
# It can remove all sequences which are "TrEMBL" fragments.
# Six arguments:
# 1. Pfam flat file name
# 2. outfile name
# 3. Name of the query - not to be removed
#    (use " " to not use this functionality)
# 4. species names file (will be ignored if " ")
# 5. 1 to NOT remove non-SWISS_PROT seqs.
# 6. 1 to remove TrEMBL seqs with "(FRAGMENT)" in their DE line.
#    (Only used if non SWISS_PROT seqswill not be removed)
# Returns the number of sequences in the resulting alignment.
# If a query name is given, it returns -1 if query is not found in alignment,
# -10 if the name is not unique.
# Last modified: 05/11/01
sub removeSeqsFromPfamAlign {
    my $infile                   = $_[ 0 ];
    my $outfile                  = $_[ 1 ];
    my $query                    = $_[ 2 ];
    my $species_names_file       = $_[ 3 ];
    my $keep_non_sp              = $_[ 4 ];
    my $remove_frags             = $_[ 5 ];
    my $return_line              = "";
    my $name                     = "";
    my $seq                      = "";
    my $saw_sequence_line        = 0;
    my $number_of_seqs           = 0;
    my $saw_query                = 0;
    my $query_given              = 0;
    my $species_names_file_given = 0;
    my $length_of_name           = 0;
    my %AC_OS                    = (); # AC -> species name (TrEMBL)
    my %AC_DE                    = (); # AC -> description (TrEMBL)
    my $AC                       = "";
    my $DE                       = "";
    my $OS                       = "";

    &testForTextFilePresence( $infile );

    if ( $query =~ /\S/ ) {
        $query_given = 1;
    }
    if ( $species_names_file =~ /\S/ ) {
        $species_names_file_given = 1;
        &readSpeciesNamesFile( $species_names_file );
    }
    
    if ( $keep_non_sp == 1 
    || ( $query_given == 1 && !&startsWithSWISS_PROTname( $query ) ) ) {

        &testForTextFilePresence( $TREMBL_ACDEOS_FILE );
       
        # Fill up hash $AC_OS and $AC_DE.
        open( HH, "$TREMBL_ACDEOS_FILE" ) || &dieWithUnexpectedError( "Cannot open file \"$TREMBL_ACDEOS_FILE\"" );
        while ( $return_line = <HH> ) {
            if ( $return_line =~ /(\S+);([^;]*);(\S+)/ ) {
                $AC_OS{ $1 } = $3;
                if ( $remove_frags == 1 ) {
                    $AC_DE{ $1 } = $2;
                }
            }
        }
        close( HH ); 
    }

    open( OUT_RNSP, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
    open( IN_RNSP, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
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
            $return_line =~ /(\S+)\s+(\S+)/;
            $name = $1;
            $seq  = $2;
            if ( $query_given == 1 && $name eq $query ) {
                $saw_query++;
            }
            if ( ( $query_given == 1 && $name ne $query )
            || $query_given != 1 ) {
                if ( !&startsWithSWISS_PROTname( $name ) ) {
                    if ( $keep_non_sp != 1 ) {
                        next;
                    }
                    else {
                        $name =~ /(\S+)\//;
                        $AC = $1;
                        unless( exists( $AC_OS{ $AC } ) ) {
                            #ACs not present in "ACDEOS" file.
                            next;
                        }
                        $OS = $AC_OS{ $AC };
                        if ( !$OS || $OS eq "" ) {
                            &dieWithUnexpectedError( "species for \"$AC\" not found" );
                        }
                        if ( $species_names_file_given == 1 ) { 
                            unless( exists( $Species_names_hash{ $OS } ) ) {
                                next;
                            }
                        }
                        if ( $remove_frags == 1 ) {
                            $DE = $AC_DE{ $AC };
                            if ( $DE && $DE =~ /\(FRAGMENT\)/ ) {
                                next;
                            }
                        }
                        $name =~ s/\//_$OS\//;
                    }
                }
                else {
                    if ( $species_names_file_given == 1 ) {   
                        if ( $name =~ /_([A-Z0-9]{1,5})/ ) {
                            unless( exists( $Species_names_hash{ $1 } ) ) {
                                next;
                            }
                        }
                        # remove everything whose species cannot be determined.
                        else {
                            next;
                        }
                    }
                }
            }
            elsif ( $query_given == 1 && $name eq $query
            && !&startsWithSWISS_PROTname( $query ) ) {
                # Adding species to non SWISS-PROT query
                $name =~ /(\S+)\//;
                $AC   = $1;
                unless( exists( $AC_OS{ $AC } ) ) {
                #ACs not present in "ACDEOS" file.
                    &userError( "Could not establish species of query.\n Check file \"$TREMBL_ACDEOS_FILE\"." );
                }
                $OS = $AC_OS{ $AC };
                if ( !$OS || $OS eq "" ) {
                    &dieWithUnexpectedError( "species for \"$AC\" not found" );
                }
                $name =~ s/\//_$OS\//;
            }

            $length_of_name = length( $name );

            if ( $length_of_name > ( $LENGTH_OF_NAME - 1 ) ) {
                &userError( "Name \"$name\" is too long." );
            }

            for ( my $j = 0; $j <= ( $LENGTH_OF_NAME - $length_of_name - 1 ); ++$j ) {
	            $name .= " ";
            }
            
            $return_line = $name.$seq."\n";
        }
        
        print OUT_RNSP $return_line;
        if ( $saw_sequence_line == 1 ) {
            $number_of_seqs++;
        }
    }
    close( IN_RNSP );
    close( OUT_RNSP );
    if ( $query_given == 1 ) {
        if ( $saw_query < 1 ) {
            return -1;
        }
        elsif ( $saw_query > 1 ) {
            return -10;
        }
    }
    return $number_of_seqs;

} ## removeSeqsFromPfamAlign




# One argument:
# 1. PWD file
# "Returns" a Hash of Strings (=keys) containing all the names found in PWD file
# Last modified: 05/29/01
sub getNamesFromPWDFile {
    my $infile        = $_[ 0 ];
    my $return_line   = "";
    my $i             = 0;
    my $saw_dist_line = 0;

    &testForTextFilePresence( $infile );

    open( GN_IN, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );

    while ( $return_line = <GN_IN> ) {
        if ( $saw_dist_line == 1 && $return_line =~ /^\s*(\d+)\s*$/ ) {
            if ( $1 != $i ) {
                &dieWithUnexpectedError( "Failed sanity check" );
            } 
            last;
        }
        elsif ( $return_line =~ /^\s*(\S+)\s+\S+/ ) {
            $names_in_pwd_file{ $1 } = 0;
            $i++;
            $saw_dist_line = 1;
        }
    }
    close( GN_IN );
    return;
} ## getNamesFromPWDFile




# Moves sequences which start with query name (argument 1) 
# to the last positions in pfam alignment sepecified by argument 2.
# Removes seqs present in argument 4, unless for query name.
# Four arguments:
# 1. Query name
# 2. Infile (alignment)
# 3. Outfile (=infile with query seq moved to the bottom)
# 4. Array of seq names to remove, unless for query name
# Last modified: 06/25/01
sub moveToLast {
    my $query       = $_[ 0 ];
    my $infile      = $_[ 1 ];
    my $outfile     = $_[ 2 ];
    my @to_remove   = @{ $_[ 3 ] };  # @{} tells Perl that this is a list.
    my $return_line = "";
    my $query_line  = "";
    my $n           = "";    

    &testForTextFilePresence( $infile );

    open( MTL_IN, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    open( MTL_OUT, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );

    W: while ( $return_line = <MTL_IN> ) {
        if ( &isPfamCommentLine( $return_line ) 
        && ( !isRFline( $return_line ) || $mode != 1 ) ) {
            next W;
        }
        if ( @to_remove > 1 ) {
            foreach $n ( @to_remove ) {
                if ( $n ne $query && $return_line =~ /^\s*$n\s+/ ) {
                    next W;
                }
            }
        }
        if ( $return_line =~ /^\s*$query\s+/ ) {
            $query_line = $return_line;
        }
        elsif ( $query_line ne "" 
        && ( $return_line !~ /\S+/ || isRFline( $return_line ) ) ) {
            print MTL_OUT $query_line;
            print MTL_OUT $return_line;
            $query_line = "";
        }
        else {
            print MTL_OUT $return_line;
        }
    }
    if ( $query_line ne "" ) {
        print MTL_OUT $query_line;
    }

    close( MTL_IN );
    close( MTL_OUT );

    return;

} ## moveToLast










# -----------------------------------------------------------
#                                                      Others
# -----------------------------------------------------------




# This gets the complete name of a TrEMBL sequence from a Pfam alignment.
# I.e. it adds the species between "_" and "/XXX-XXX".
# 2 arguments:
# 1. Infile (alignment)
# 2. Name of query
# Returns the complete name found.
# Last modified: 04/25/01
sub getCompleteNameForTrEMBLquerySeq {
    
    my $infile         = $_[ 0 ];
    my $query_name     = $_[ 1 ];
    my $return_line    = "";
    my $complete_name  = "";
    my $before_slash   = "";
    my $after_slash    = "";
    
    &testForTextFilePresence( $infile );

    $query_name =~ /(.+)\/.+/;
    $before_slash = $1;

    $query_name =~ /.+\/(.+)/;
    $after_slash  = $1;

    open( INGCN, $infile ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    while ( $return_line = <INGCN> ) {
        if ( $return_line =~ /^\s*($before_slash.+\/$after_slash)/ ) {
            $complete_name = $1;
            last;
        }
    }
    close( INGCN );
    if ( $complete_name eq "" ) { 
        &userError( "Could not find \"$query_name\" in \"$alignment\"." );
    }
    return $complete_name;
} ## getCompleteNameForTrEMBLquerySeq




# One argument:
# Pfam align name.
# Last modified: 02/26/01
sub getDescriptionFromPfam {

    my $infile      = $_[ 0 ];
    my $return_line = "";
    my $result      = "";

    &testForTextFilePresence( $infile );

    open( INGDPF, $infile ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    while ( $return_line = <INGDPF> ) {
        if ( $return_line =~ /^\s*#=DE\s+(.+)/ ) {
            $result = $1;
            close( INGDPF );
            return $result;
        }
    }
    close( INGDPF );
    return $result;

} ## getDescriptionFromPfam



# Reads in (SWISS-PROT) species names from a file.
# Names must be separated by newlines.
# Lines beginning with "#" are ignored.
# A possible "=" and everything after is ignored.
# One argument: species-names-file name
# Last modified: 04/24/01
sub readSpeciesNamesFile {
    my $infile = $_[ 0 ];
    my $return_line = "";
    my $species     = "";

    &testForTextFilePresence( $infile );

    open( IN_RSNF, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
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




# This reads a raw sequence file or a FASTA sequence file
# and saves it as a "cleaned up" FASTA sequence file.
# If no > line is in the file, it creates one with new sequence name.
# If a > line is in the file, it modifes it:
# white space -> _, ";" ":" "," or "|" -> "~", deletes everything after ( or [;
# length is limited to 40 characters.
# Error if $new_seq_name is "" and  no > line in the file.
# Two/three arguments:
# 1. infile name
# 2. outfile name
# 3. new sequence name for > line(will be ignored if "")
# If new sequence name is "":
# returns the contents of the ">" line after modification.
# If new sequence name is specified:
# return new sequence name.
# Last modified: 03/04/01
sub seqFile2CleanedUpFastaFile {
    my $infile        = $_[ 0 ];
    my $outfile       = $_[ 1 ];
    my $new_seq_name  = $_[ 2 ];
    my $return_line   = "";
    my $mod_desc      = "";
    my $saw_desc_line = 0;

    &testForTextFilePresence( $infile );

    open( IN_CUFF, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    open( OUT_CUFF, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
   
    while ( $return_line = <IN_CUFF> ) {
	    if ( $return_line =~ /\w/ && $return_line !~ /^\s*#/ ) {
	        if ( $return_line =~ /^\s*>/ ) {
                    if ( $new_seq_name eq "" && $return_line !~ /_/ ) {
                        &userError( "Description line of query file appears not to\n contain any species information. Use \"N=\" option." );
                    }
                    elsif ( $new_seq_name eq "" ) {
                        $return_line =~ s/^\s*>\s*(.*?)\s*/>$1/; # Removes spaces before and after >.
                        $return_line = substr( $return_line, 0, $LENGTH_OF_NAME - 1 );
                        $return_line =~ s/[\(\[].*//;            # Removes "(" or "[" and everything after.
                        $return_line =~ s/\s+$//;                # Removes spaces at end.
                        $return_line =~ s/\s+/_/g;               # Replaces all white spaces with "_".
                        $return_line =~ s/[;:,\|]/~/g;           # Replaces all ";", ":", ",", or "|" with "~".                
                        $return_line =~ />\s*(\S+)/;
                        $mod_desc = $1;
                        $return_line .= "\n";
                    }
                    else {
                        $return_line = ">".$new_seq_name."\n";
                        $mod_desc = $new_seq_name;
                    }
                    $saw_desc_line = 1;
	        }
                else {
                    if ( $saw_desc_line != 1 ) {
                        if ( $new_seq_name ne "" ) {
                            print OUT_CUFF ( ">".$new_seq_name."\n" );
                            $mod_desc = $new_seq_name;
                        }
                        else {
                            &userError( "Query file is not a FASTA file\n and option \"N=\" has not been used." );
                        }
                        $saw_desc_line = 1;
                    }
                    $return_line =~ s/[^a-zA-Z\r\n\f]//g;    # Removes non-letters from sequence.
		}
	       
	        if ( $return_line =~ /\w/ ) {
        	    print OUT_CUFF $return_line;
	        }
	    }
    }
    close( IN_CUFF );
    close( OUT_CUFF );

    return $mod_desc;
} ## seqFile2CleanedUpFastaFile




# Purpose. Gets description for TrEMBL seqs,
# from a file which contains the AC, DE, and OS
# and which has to be generated from a TrEMBL flat file db
# using "extractTrembl.pl".
# The same file is used in method "addSpeciesToNonSPseqs".
# Two arguments:
# 1. "ACDEOS" file (AC, DE, OS from TrEMBL db)
# 2. AC ("_species/..." is removed)
#    Format: AC;DE;OS\n
# Last modified: 02/14/02
sub getDescriptionFromTrEMBL_ACDEOSfile {
    my $ACDEOS = $_[ 0 ];
    my $AC     = $_[ 1 ];
    my $DE     = "";
    
    # Fill up (huge) hash, if not already done.
    unless ( %AC_DE ) {
        &testForTextFilePresence( $ACDEOS );
        open( ACDEOS, "$ACDEOS" ) || &dieWithUnexpectedError( "Cannot open file \"$ACDEOS\"" );
        while ( $return_line = <ACDEOS> ) {
            if ( $return_line =~ /(\S+);([^;]+);/ ) {
                $AC_DE{ $1 } = $2;   
            }
        }
        close( ACDEOS ); 
    }

    $AC =~ s/_.+//;

    unless( exists( $AC_DE{ $AC } ) ) {
        #AC not present in "ACDEOS" file.
        return "-";
    }

    $DE = $AC_DE{ $AC };
   
    if ( !$DE || $DE eq "" ) {
        $DE = "-";
    }

    return $DE;

} ## getDescriptionFromTrEMBL_ACDEOSfile



# Purpose. Gets description for SP seqs,
# from a file which contains the AC, DE, and OS
# and which has to be generated from a sprot.dat flat file db
# using "extractSWISS-PROT.pl".
# Two arguments:
# 1. "ACDEOS" file (AC, DE, OS from SWISS-PROT db)
# 2. SWISS-PROT AC (XXXX_XXXX)
#    Format: AC;DE;OS\n
# Last modified: 02/12/02
sub getDescriptionFromSWISSPROT_ACDEOSfile {
    my $SPACDEOS = $_[ 0 ];
    my $AC       = $_[ 1 ];
    my $DE       = "";
    
    # Fill up (huge) hash, if not already done.
    unless ( %SP_AC_DE ) {
        &testForTextFilePresence( $SPACDEOS );
        open( ACDEOS, "$SPACDEOS" ) || &dieWithUnexpectedError( "Cannot open file \"$SPACDEOS\"" );
        while ( $return_line = <ACDEOS> ) {
            if ( $return_line =~ /(\S+);([^;]+);/ ) {
                $SP_AC_DE{ $1 } = $2;   
            }
        }
        close( ACDEOS ); 
    }

    $AC =~ s/\/.+//;

    unless( exists( $SP_AC_DE{ $AC } ) ) {
        #AC not present in "ACDEOS" file.
        return "-";
    }

    $DE = $SP_AC_DE{ $AC };
   
    if ( !$DE || $DE eq "" ) {
        $DE = "-";
    }

    return $DE;

} ## getDescriptionFromSWISSPROT_ACDEOSfile









# -----------------------------------------------------------
#                                                     Helpers
# -----------------------------------------------------------



# One argument:
# Numeric value to be rounded to int.
# Last modified: 10/17/01
sub roundToInt {
    my $x = $_[ 0 ];
    unless ( $x eq "-" ) {
        $x = int ( $x + 0.5 );
    }
    return $x;
} ## roundToInt

        

# Removes files.
# Last modified: 03/10/01
sub cleanUpTempDir {
    unlink( $temp_dir."/MAKETREEOUT".$TREE_FILE_SUFFIX, $temp_dir."/MAKETREEOUT".$LOG_FILE_SUFFIX,
            $temp_dir."/MAKETREEOUT".$ALIGN_FILE_SUFFIX, $temp_dir."/MAKETREEOUT".$MULTIPLE_TREES_FILE_SUFFIX,
            $temp_dir."/MAKETREEOUT".$SUFFIX_PWD_NOT_BOOTS, $temp_dir."/".$DO_RIO_TEMP_OUTFILE,
            $temp_dir."/ALIGN1",  $temp_dir."/ALIGN2", $temp_dir."/QUERY_SEQ", $temp_dir."/NBD_NJ_TREE",
            $temp_dir."/ALIGN2_BOOTSTRAPPED", $temp_dir."/ALIGN2_PROCESSED", $temp_dir."/DIST_TO_QUERY",
            $temp_dir."/DISTs_TO_QUERY", $temp_dir."/HMMALIGNOUT", $temp_dir."/NBD_INC_QUERY", $temp_dir."/PWD_INC_QUERY",
            $temp_dir."/HMMFILE", $temp_dir."/MOVETOLASTOUT" );
    rmdir( $temp_dir );
} ## cleanUpTempDir












# -----------------------------------------------------------
#                          Command line and arguments, Errors
# -----------------------------------------------------------



# One argument:
# the command line.
# Last modified: 03/08/01
sub analyzeCommandLine {

    my $args = "";
    my $arg  = "";
    my $char = "";
   
   

    $mode = shift( @_ );

    if ( $mode != 1 && $mode != 2 && $mode != 3 && $mode != 4 ) {
        &errorInCommandLine( "Mode can only be: 1, 2, 3, or 4." );
    } 

    
    foreach $args ( @_ ) {

        $args =~ s/\s//g;

        $char = substr( $args, 0, 1 );
       
        
        if ( length( $args ) > 1 ) {
            $arg = substr( $args, 2 );
        }

        if ( $char =~ /A/ ) {
            if (  $alignment ne "" ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            if ( $mode == 3 || $mode == 4 ) { 
                &userErrorCheckForTextFileExistence( $arg );
            }
            $alignment = $arg;
        }
        elsif ( $char =~ /B/ ) {
            if ( $t_sn != $THRESHOLD_SN_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $t_sn = $arg;
        }
        elsif ( $char =~ /C/ ) {
            if ( $description == 1 || $complete_description == 1 ) {
                &errorInCommandLine( "Entered same argument twice or conflicting arguments: \"D\" and \"C\"." );
            }
            $complete_description = 1;
        }
        elsif ( $char =~ /D/ ) {
            if ( $description == 1 || $complete_description == 1 ) {
                &errorInCommandLine( "Entered same argument twice or conflicting arguments: \"D\" and \"C\"." );
            }
            $description = 1;
        }
        elsif ( $char =~ /E/ ) {
            if ( $long_output != 0 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $long_output = 1;
        }
        elsif ( $char =~ /F/ ) {
            if ( $hmm_file ne "" || $hmm_name ne "" || $seed_aln_for_hmmbuild ne "") {
                &errorInCommandLine( "Entered same argument twice or conflicting arguments: \"F=\", \"H=\" and \"b=\"." );
            }
            if ( $mode == 1 || $mode == 2 ) {
                &errorInCommandLine( "Can not use \"F=\" in modes 1 or 2." );
            }
            &userErrorCheckForTextFileExistence( $arg );
            $hmm_file = $arg;
        }
        elsif ( $char =~ /G/ ) {
            if ( $species_names_file ne " " ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            &userErrorCheckForTextFileExistence( $arg );
            $species_names_file = $arg;
        }
        elsif ( $char =~ /H/ ) {
            if ( $hmm_name ne "" || $hmm_file ne "" || $seed_aln_for_hmmbuild ne "" ) {
                &errorInCommandLine( "Entered same argument twice or conflicting arguments: \"F=\", \"H=\" and \"b=\"." );
            }
            if ( $mode == 1 || $mode == 2 ) {
                &errorInCommandLine( "Can not use \"H=\" in modes 1 or 2." );
            }
            $hmm_name = $arg;
        }
        elsif ( $char =~ /I/ ) {
            if ( $safe_nhx != 0 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $safe_nhx = 1;
        }
        elsif ( $char =~ /K/ ) {
            if ( $keep != 0 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $keep = 1;
        }
        elsif ( $char =~ /L/ ) {
            if ( $t_orthologs != $THRESHOLD_ORTHOLOGS_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $t_orthologs = $arg;
        }
        elsif ( $char =~ /N/ ) {
            if ( $query_name ne "" ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $query_name = $arg;
        }
        elsif ( $char =~ /O/ ) {
            if ( $outfile ne "" ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $outfile = $arg;
        }
        elsif ( $char =~ /P/ ) {
            if ( $sort != $SORT_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $sort = $arg;
        }
        elsif ( $char =~ /Q/ ) {
            if ( $seqX_file ne "" ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            &userErrorCheckForTextFileExistence( $arg );
            $seqX_file = $arg;
        }
        elsif ( $char =~ /S/ ) {
            if ( $species_tree_file ne "" ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            &userErrorCheckForTextFileExistence( $arg );
            $species_tree_file = $arg;
        }
        elsif ( $char =~ /T/ ) {
            if ( $mode == 1 || $mode == 2 ) {
                &errorInCommandLine( "Matrix cannot be changed in modes 1 and 2 (is dictated by \"\$MATRIX_FOR_PWD\" for mode 1)." );
            }
            if ( $arg eq "J" ) {
                $matrix_n = 0;
            }
            elsif ( $arg eq "P" ) {
                $matrix_n = 1;
            }
            elsif ( $arg eq "B" ) {
                $matrix_n = 2;
            }
            elsif ( $arg eq "M" ) {
                $matrix_n = 3;
            }
            elsif ( $arg eq "V" ) {
                $matrix_n = 5;
            }
            elsif ( $arg eq "W" ) {
                $matrix_n = 6;
            }
            else {
                &errorInCommandLine( "Use T=J for JTT, P for PAM, B for BLOSUM62, M for mtREV24, V for VT, W for WAG." );
            }
        }
        elsif ( $char =~ /U/ ) {
            if ( $t_orthologs_dc != $THRESHOLD_ORTHOLOGS_DEFAULT_DC ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $t_orthologs_dc = $arg;
        }
        elsif ( $char =~ /X/ ) {
            if ( $warn_more_than_one_ortho
                 != $WARN_MORE_THAN_ONE_ORTHO_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $warn_more_than_one_ortho = $arg;
        }
        elsif ( $char =~ /Y/ ) {
            if ( $warn_no_orthos != $WARN_NO_ORTHOS_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $warn_no_orthos = $arg;
        }
        elsif ( $char =~ /Z/ ) {
            if ( $warn_one_ortho != $WARN_ONE_ORTHO_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $warn_one_ortho = $arg;
        }
        elsif ( $char =~ /a/ ) {
            if ( $boostraps_for_makeTree != $BOOSTRAPS_FOR_MAKETREE_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            if ( $mode == 1 || $mode == 2 ) {
                &errorInCommandLine( "Modes 1 and 2: Cannot change bootstrap value. Do not use \"a=\"." );
            }
            $boostraps_for_makeTree = $arg;
            if ( $boostraps_for_makeTree < 10 ) {
                &errorInCommandLine( "Bootsraps cannot be smaller than 10." );
            }
        }
        elsif ( $char =~ /b/ ) {
            if ( $hmm_name ne "" || $hmm_file ne "" || $seed_aln_for_hmmbuild ne "" ) {
                &errorInCommandLine( "Entered same argument twice or conflicting arguments: \"F=\", \"H=\" and \"b=\"." );
            }
            if ( $mode == 1 || $mode == 2 ) {
                &errorInCommandLine( "Can not use \"b=\" in modes 1 or 2." );
            }
            &userErrorCheckForTextFileExistence( $arg );
            $seed_aln_for_hmmbuild = $arg;
        }
        elsif ( $char =~ /f/ ) {
            if ( $no_frags ne 0 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $no_frags = 1;
        }
        elsif ( $char =~ /j/ ) {
            if ( $temp_dir ne "" ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $temp_dir = $arg;
        }
        elsif ( $char =~ /p/ ) {
            if ( $output_ultraparalogs != 0 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $output_ultraparalogs = 1;
        }
        elsif ( $char =~ /s/ ) {
            if ( $non_sp != 1 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $non_sp = 0;
        }
        elsif ( $char =~ /v/ ) {
            $t_ultra_paralogs = $arg;
        }
        elsif ( $char =~ /x/ ) {
            if ( $output_HTML == 1 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $output_HTML = 1;
        }
        elsif ( $char =~ /y/ ) {
            if ( $seed_for_makeTree != $SEED_FOR_MAKETREE_DEFAULT ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $seed_for_makeTree = $arg;
        }
        elsif ( $char =~ /\+/ ) {
            if ( $parallel != 0 ) {
                &errorInCommandLine( "Entered same argument twice." );
            }
            $parallel = 1;
        }
        else {
            &errorInCommandLine( "Unknown option: \"$args\"." );
        }
    }
} ## analyzeCommandLine




# Last modified: 03/08/01
sub CheckArguments {
    
    if ( $outfile eq "" ) {
        &errorInCommandLine( "Outfile not specified. Use \"O=\"." );
    }
    if ( $alignment eq "" ) {
        &errorInCommandLine( "Need to specify a Pfam alignment file. Use \"A=\"." );
    }
    if ( -e $outfile ) {
        &userError( "\"$outfile\" already exists." );
    } 

    if ( $sort < 0 || $sort > 17 ) {
        &errorInCommandLine( "Sort priority (\"P=\") must be between 0 and 15." );
    }
 
    if ( $parallel == 1 && $mode != 1 ) {
        &errorInCommandLine( "Parallelization only implemented for mode 1." );
    } 

    if ( $mode == 1 || $mode == 2 ) {
        
        if ( $species_names_file =~ /\S/ ) {
            &errorInCommandLine( "Modes 1 and 2: Cannot use species names file. Do not use \"G=\"." );
        }
        if ( $non_sp == 0 ) {
            &errorInCommandLine( "Can not use \"s\" in modes 1 or 2." );
        }
        if ( $no_frags == 1 ) {
            &errorInCommandLine( "Can not use \"f\" in modes 1 or 2." );
        }
    }

    if ( $mode == 1 || $mode == 3 ) {
        if ( $seqX_file eq "" ) {
            &errorInCommandLine( "Modes 1 and 3: Need to specify a query file. Use \"Q=\"." );
        }
    }

    if ( $mode == 3 ) {
        if ( $hmm_name eq "" && $hmm_file eq "" && $seed_aln_for_hmmbuild eq "" ) {
            &errorInCommandLine( "Mode 3: Need to specify either a HMM name (\"H=\"), a HMM file (\"F=\") or build a HMM (\"b=\")." );
        }
    }

    if ( $mode == 1 ) {
        if ( $hmm_name ne "" || $hmm_file ne "" || $seed_aln_for_hmmbuild ne "" ) {
            &errorInCommandLine( "Mode 1: Must not specify a HMM name (\"H=\"), a HMM file (\"F=\") or build a HMM (\"b=\")." );
        }
    }

    if ( $mode == 2 || $mode == 4 ) {
        if ( $seqX_file ne "" ) {
            &errorInCommandLine( "Modes 2 and 4: Must not specify a query file. Do not use \"Q=\".\n" );
        }
        if ( $query_name eq "" ) {
            &errorInCommandLine( "Modes 2 and 4: Must specify a query name. Use \"N=\"." );
        }
        if ( $hmm_name ne "" || $hmm_file ne "" || $seed_aln_for_hmmbuild ne "" ) {
            &errorInCommandLine( "Modes 2 and 4: Cannot specify a HMM name (\"H=\"), a HMM file (\"F=\") or build a HMM (\"b=\")." );
        }

    }
    
    if ( $non_sp != 1 && $no_frags == 1 ) {
        &errorInCommandLine( "\"Fragments\" are assumed to be only found in non SWISS-PROT seqs.\n Do not use \"f\" together with \"s\"." );
    }

    if ( $output_HTML == 1 ) {
       if ( $mode != 1 ) {
           &errorInCommandLine( "Output in HTML (for web server) only for mode 1." );
       }
    }
    
    if ( $output_ultraparalogs == 0 && $t_ultra_paralogs != $T_ULTRA_PARALOGS_DEFAULT ) {
        &errorInCommandLine( "Use \"p\" to output ultra paralogs (cannot use \"v=\" without \"p\")." );
    }

    if ( $non_sp == 1 &&  ( $mode == 3 || $mode == 4 ) ) {
        unless ( ( -s $TREMBL_ACDEOS_FILE  ) && ( -f $TREMBL_ACDEOS_FILE  ) && ( -T $TREMBL_ACDEOS_FILE ) ) {
            my $message = "AC, DE, and OS-file not found.\n";
            $message .= " If non SWISS-PROT sequences are not to be removed from the\n";
            $message .= " Pfam alignment (\"s\" option), variable \"\$TREMBL_ACDEOS_FILE\" needs\n";
            $message .= " to point to a file containing AC, DE, and OS from TrEMBL. Such a\n";
            $message .= " file can be generated with \"extractTrembl.pl\".\n";
            $message .= " Currently, \"TREMBL_ACDEOS_FILE\" points to:\n";
            $message .= " $TREMBL_ACDEOS_FILE";
            &userError( $message );
        }
    }

    unless ( ( -s $species_tree_file  ) && ( -f $species_tree_file  ) && ( -T $species_tree_file ) ) {
        my $message = "Species tree file not found.\n";
        $message .= " A valid species tree must be specified.\n";
        $message .= " Either, use \"S=\" option, or set variable\n";
        $message .= " \"\$SPECIES_TREE_FILE_DEFAULT\".\n";
        $message .= " Currently, this program looks for a species tree at:\n";
        $message .= " $species_tree_file";
        &userError( $message );
    }

    if ( $hmm_name ne "" ) {
        unless ( ( -s $PFAM_HMM_DB  ) && ( -f $PFAM_HMM_DB ) ) {
            my $message = "HMMER model db file not found.\n";
            $message .= " If \"H=\" option is used, a valid HMMER model db needs\n";
            $message .= " to be specified with variable \"\$PFAM_HMM_DB\".\n";
            $message .= " Currently, \"\$PFAM_HMM_DB\" points to:\n";
            $message .= " $PFAM_HMM_DB";
            &userError( $message );
        }
    }
} ## CheckArguments



# Last modfied: 06/25/01
sub userErrorCheckForTextFileExistence {
    my $file = $_[ 0 ];
    unless ( ( -s $file ) && ( -f $file ) && ( -T $file ) ) {
        &userError( "\"$file\" does not exist or is not a plain text file." );
    }
} ## checkForFileExistence



# One argument: the error message.
# Last modified: 04/26/01
sub errorInCommandLine {
    
    my $error = $_[ 0 ];

    print " \n";
    print " rio.pl  version: $VERSION\n";
    print " ------\n";
    print " \n";
    print " Error in command line:\n";
    if ( $error ne "" ) {
        print " $error";
    }
    print " \n\n";
    print " Type \"rio.pl\" (no arguments) for more information.\n";  
    print " \n";
    exit( -1 );
} ## errorInCommandLine




# One argument: the error message.
# Last modified: 04/26/01
sub userError {
    
    my $error = $_[ 0 ];

    print " \n";
    print " rio.pl  version: $VERSION\n";
    print " ------\n";
    print " \n";
    print " Error:\n";
    if ( $error ne "" ) {
        print " $error";
    }
    print " \n\n";
    print " Type \"rio.pl\" (no arguments) for more information.\n";  
    print " \n";
    &cleanUpTempDir();
    exit( -1 );
} ## UserError






# Last modified: 04/26/01
sub printHelp {

    print " \n";
    print " rio.pl  version: $VERSION\n";
    print " ------\n\n";

    print <<END;
 Copyright (C) 2000-2002 Washington University School of Medicine
 and Howard Hughes Medical Institute
 All rights reserved

 Created: 11/25/00
 Author: Christian M. Zmasek
 zmasek\@genetics.wustl.edu
 http://www.genetics.wustl.edu/eddy/people/zmasek/

 Last modified 05/26/02

 Available at : http://www.genetics.wustl.edu/eddy/forester/
 RIO webserver: http://www.rio.wustl.edu/

 Reference:
 Zmasek C.M. and Eddy S.R. (2002)
 RIO: Analyzing proteomes by automated phylogenomics using
 resampled inference of orthologs.
 BMC Bioinformatics 3:14
 http://www.biomedcentral.com/1471-2105/3/14/

 It is highly recommended that you read this paper before
 installing and/or using RIO. (Included in the RIO 
 distribution as PDF: "RIO.pdf".)
 

 Before rio.pl can be used, some variables in rio_module.pm need to be set, 
 as described in RIO_INSTALL.



 Usage: rio.pl <Mode: 1, 2, 3, or 4> <tagged arguments, single letter arguments>
 -----


 Examples:
 ---------

 % RIO1.1/perl/rio.pl 1 A=aconitase Q=RIO1.1/LEU2_HAEIN N=QUERY_HAEIN O=out1 p I C E
 
 % RIO1.1/perl/rio.pl 2 A=aconitase N=LEU2_LACLA/5-449 O=out2 p I C E

 % RIO1.1/perl/rio.pl 3 A=/path/to/my/pfam/Full/aconitase H=aconitase Q=RIO1.1/LEU2_HAEIN N=QUERY_HAEIN O=out3 p I C E

 % RIO1.1/perl/rio.pl 4 A=/path/to/my/pfam/Full/aconitase N=LEU2_LACLA/5-449 O=out4 p I C E

 % RIO1.1/perl/rio.pl 3 A=/path/to/my/pfam/Full/aconitase b=/path/to/my/pfam/Seed/aconitase Q=RIO1.1/LEU2_HAEIN N=QUERY_HAEIN O=out5 p I C E



 Modes:
 ------

 1: RIO analysis based on precalculated pairwise distances
    alignment does not contain query sequence

 2: RIO analysis based on precalculated pairwise distances
    alignment does contain query sequence

 3: RIO analysis based on Pfam alignments,
    alignment does not contain query sequence

 4: RIO analysis based on Pfam alignments,
    alignment does contain query sequence



 Tagged arguments:
 -----------------

 No "G=", "H=", "F=", "T=", "a=", "b=", "s", "f" in modes 1 and 2.


 A=<String> Pfam alignment name (mandatory). This specifies the alignment
            against which the RIO analysis is to be performed.
            In modes 1 and 2: Pfam model (alignment) name
                              (e.g. "A=aconitase").
            In modes 3 and 4: Pfam alignment path/name
                              (e.g. "A=/path/to/your/pfam/Full/aconitase").

 Q=<String> Path/name of file containing the query sequence
            (in FASTA format or raw sequence) (mandatory in modes 1 and 3).

 N=<String> Query name (mandatory). This must include the SWISS-PROT code
            for the species of the query after a "_" (e.g. "N=QUERY_HAEIN").
            If the query sequence is already in the alignment (modes 2 and 4)
            the complete name needs to be specified -- including "/xxx-xxx".

 O=<String> Output file path/name (mandatory).

 T=<char>   Model for pairwaise distance calculation:
            J=JTT, B=BLOSUM 62, M=mtREV24, V=VT, W=WAG, P=PAM.
            BLOSUM 62 is default.
            (Not in modes 1 and 2; these modes use \$MATRIX_FOR_PWD instead.)

            In modes 1 and 3, a HMM is needed to align the query sequence to
            the alignment and either one of the following options must be
            employed:
 H=<String> HMM name: This uses hmmfetch to retrieve a HMM from
            \$PFAM_HMM_DB.
 F=<String> HMM file: This directly reads the HMM from a file.

 S=<String> Species tree file path/name (in NHX format) (optional).
            If not specified, \$SPECIES_TREE_FILE_DEFAULT is used.

 G=<String> Species names file (optional). Only sequences associated with
            species found in this file are used.
            In the species names file, individual species names must be
            separated by newlines and lines starting with "#" are ignored.    
            While only sequences associated with species found in the species
            tree ("S=") are used for the actual RIO analysis, this allows to 
            remove sequences prior to tree calculation (which is the most
            time consuming step).  

 P=<int>    Sort priority (default is 12):
            0 : Ortholog  
            1 : Ortholog,         Super ortholog
            2 : Super ortholog,   Ortholog
            3 : Ortholog,         Distance
            4 : Distance,         Ortholog
            5 : Ortholog,         Super ortholog,   Distance  
            6 : Ortholog,         Distance,         Super ortholog
            7 : Super ortholog,   Ortholog,         Distance
            8 : Super ortholog,   Distance,         Ortholog
            9 : Distance,         Ortholog,         Super ortholog
           10 : Distance,         Super ortholog,   Ortholog
           11 : Ortholog,         Subtree neighbor, Distance 
           12 : Ortholog,         Subtree neighbor, Super ortholog,   Distance (default)
           13 : Ortholog,         Super ortholog,   Subtree neighbor, Distance
           14 : Subtree neighbor, Ortholog,         Super ortholog,   Distance
           15 : Subtree neighbor, Distance,         Ortholog,         Super ortholog   
           16 : Ortholog,         Distance,         Subtree neighbor, Super ortholog 
           17 : Ortholog,         Subtree neighbor, Distance,         Super ortholog 

 a=<int>    Bootstraps for tree construction (not in modes 1 and 2).
            Default is 100.   

 L=<int>    Threshold for orthologies for output. Default is 0.
 v=<int>    Threshold for ultra-paralogies for output. Default is 50.

 U=<int>    Threshold for orthologies for distance calculation. Default is 60.

 X=<int>    In case of more than one putative orthologs:
            number of sd the distance query - LCA has to differ
            from the mean to generate a warning. Default is 2.

 Y=<int>    In case of no putative orthologs:
            number of sd the distance query - root has to differ
            from mean to generate a warning. Default is 2.

 Z=<double> In case of one putative ortholog:
            threshold for factor between the two distances to their
            LCA (larger/smaller) to generate a warning. Default is 2.

 B=<int>    Threshold for subtree-neighborings. Default is 0.

 b=<String> Build HMM from seed alignment with "hmmbuild -s" (optional).
            This is to prevent from finding multiple domains per sequence
            (i.e. prevents "cutting" the query sequence). Give path/name to
            Seed with this.

 j=<String> Name for temporary directory (optional).

 y=<int>    Seed for random number generator. Default is 41.

 I          Create and save a rooted, with duplication vs speciation,
            and orthology information annotated gene tree.
            If precalculated distances are used (modes 1 and 2): this gene
            tree is a NJ tree calculated based on the non-bootstrap resampled
            (original) pairwise distances.
            If precalculated distances are not used (modes 3 and 4): this gene
            is a consenus tree with ML branch length values and is also
            annotated with bootstrap values for each node.

            Options for output:
 p          Output ultra-paralogs.
 D          Description from SWISS-PROT and TrEMBL.
 C          Complete description from SWISS-PROT and TrEMBL.
 E          118 character output instead of 78 character output.

 K          Keep intermediate files (they will go into the same directory 
            as the output file, their names are the same as of the output
            file, with various suffixes added).

 s          Ignore non SWISS-PROT sequences (i.e. sequences from TrEMBL)
            in the Pfam alignment.

 f          Try to ignore TrEMBL "fragments" (sequences with "fragment" in
            their description).

 +          Parallel, use machines listed in file \$NODE_LIST.

 x          RIO used as web server -- HTML output. 
 

END
    exit( 0 );
    
} ## printHelp


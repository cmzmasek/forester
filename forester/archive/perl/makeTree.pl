#!/usr/bin/perl -W

# makeTree.pl
# -----------
# Copyright (C) 1999-2003 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek 
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Last modified 04/06/04
#
#
#  Requirements  makeTree is part of the RIO/FORESTER suite of programs.
#  ------------  Many of its global variables are set via rio_module.pm.
#
#
#  Note. Use xt.pl (for Pfam alignments) or mt.pl (for other alignments)
#  to run makeTree.pl on whole directories of alignments files.       
#
#
#
#  Usage
#  -----
#
#  Tree calculation based on a Pfam/Clustal W alignment
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#     makeTree.pl [-options] <input alignment in SELEX (Pfam), PHYLIP
#     sequential format, or Clustal W output> <outputfile>
#     [path/name for temporary directory to be created]
#
#     Example:
#     "% makeTree.pl -UTB1000S41NDXV /DB/PFAM/Full/IL5 IL5_tree"
#
#
#  Tree calculation based on precalculated pairwise distances
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Consensus tree will have no branch length values.
#  Precalculated pairwise distances are the output of "pfam2pwd.pl",
#  number of bootstraps needs to match the one used for the pwds.
#
#     makeTree.pl <-options, includes "F"> <pwdfile: boostrapped pairwise
#     distances> <outputfile> [path/name for temporary directory
#     to be created]
#
#     Example:
#     "% makeTree.pl -FB100S21XV /pfam2pwd_out/IL5.pwd IL5_tree"
#
#
#  Tree calculation based on precalculated pairwise distances
#  and matching alignment
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Consensus tree will have branch length values.
#  Precalculated pairwise distances and the matching (processed)
#  alignment are the output of "pfam2pwd.pl", number of bootstraps
#  needs to match the one used for the pwds, matrix needs to match
#  the one used for the pwds.
#
#     makeTree.pl <-options, includes "UF"> <pwdfile: boostrapped pairwise
#     distances> <alnfile: corresponding alignment>
#     <outputfile> [path/name for temporary directory to be created]
#         
#     Example:
#     "% makeTree.pl -UFLB100S21XV /pfam2pwd_out/IL5.pwd /pfam2pwd_out/IL5.aln IL5_tree"
#
#
#  Options
#  -------
#
#  N  : Suggestion to remove columns in the alignment which contain gaps.
#       Gaps are not removed, if, after removal of gaps, the resulting alignment would
#       be shorter than $MIN_NUMBER_OF_AA. Default is not to remove gaps.
#  Bx : Number of bootstrapps. B0: do not bootstrap. Default is 100 bootstrapps.
#       The number of bootstrapps should be divisible by 10.
#  U  : Use TREE-PUZZLE to calculate ML branchlengths for consesus tree, in case of 
#       bootstrapped analysis.
#  J  : Use JTT matrix (Jones et al. 1992) in TREE-PUZZLE, default: PAM.
#  L  : Use BLOSUM 62 matrix (Henikoff-Henikoff 92) in TREE-PUZZLE, default: PAM.
#  M  : Use mtREV24 matrix (Adachi-Hasegawa 1996) inTREE-PUZZLE, default: PAM.
#  W  : Use WAG matrix (Whelan-Goldman 2000) in TREE-PUZZLE, default: PAM.
#  T  : Use VT matrix (Mueller-Vingron 2000) in TREE-PUZZLE, default: PAM.
#  P  : Let TREE-PUZZLE choose which matrix to use, default: PAM.
#  E  : Exact parameter estimates in TREE-PUZZLE, default: Approximate.
#       Model of rate heterogeneity in TREE-PUZZLE (default: Uniform rate)
#  g  : 8 Gamma distributed rates
#  t  : Two rates (1 invariable + 1 variable)
#  m  : Mixed (1 invariable + 8 Gamma rates)
#  R  : Randomize input order in PHYLIP NEIGHBOR.
#  A  : Use PHYLIP PROTPARS instead of NEIGHBOR (and no pairwise distance calculation).
#  jx : Number of jumbles when using PHYLIP PROTPARS (random seed set with Sx).
#  Sx : Seed for random number generator(s). Must be 4n+1. Default is 9.
#  X  : To keep multiple tree file (=trees from bootstrap resampled alignments).
#  D  : To keep (and create in case of bootstrap analysis) pairwise distance matrix file.
#       This is created form the not resampled (original) alignment.
#  C  : Calculate pairwise distances only (no tree). Bootstrap is always 1.
#       No other files are generated.  
#  F  : Pairwise distance (pwd) file as input (instead of alignment).
#       No -D, -C, and -N options available in this case.
#  V  : Verbose.
#  #  : Only for rio.pl: Do not calculate consensus tree ("I" option in rio.pl).
#
#
#
# History:
# -------
#
# 09/06/03: Added "#" option (to be used only for rio.pl).
# 03/24/04: Do not replace "?" with "-" in method pfam2phylip.
       

use strict;

use FindBin;
use lib $FindBin::Bin;
use rio_module2;

my $VERSION                = "4.210";

my $TEMP_DIR_DEFAULT       = "/tmp/maketree"; # Where all the infiles, outfiles, etc will be created.
   
my $remove_gaps            = 0;   # 0: do not remove gaps;  1: remove gaps
my $bootstraps             = 100; # 0,1: do not bootstrap. Default: 100
my $puzzle_consensus_tree  = 0;   # 0: no; 1: yes. No is default.
my $matrix                 = 1;   # 0 = JTT
                                  # 1 = PAM - default 
                                  # 2 = BLOSUM 62
                                  # 3 = mtREV24
                                  # 5 = VT
                                  # 6 = WAG 
                                  # 7 = auto
my $rate_heterogeneity     = 0;   # 0 = Uniform rate (default)
                                  # 1 = 8 Gamma distributed rates
                                  # 2 = Two rates (1 invariable + 1 variable)
                                  # 3 = Mixed (1 invariable + 8 Gamma rates)
my $randomize_input_order  = 0;   # 0: do not randomize input order; 1 jumble
my $seed                   = 9;   # Seed for random number generators. Default: 9
my $keep_multiple_trees    = 0;   # 0: delete multiple tree file
                                  # 1: do not delete multiple tree file
my $keep_distance_matrix   = 0;   # 1: (create and) keep; 0: do not (create and) keep
my $verbose                = 0;   # 0: no; 1: yes
my $pairwise_dist_only     = 0;   # 0: no; 1: yes
my $start_with_pwd         = 0;   # 0: no; 1: yes
my $start_with_pwd_and_aln = 0;   # 0: no; 1: yes
my $no_consenus_tree       = 0;   # 0: no; 1: yes
my $exact_parameter_est    = 0;   # 0: no; 1: yes
my $use_protpars           = 0;   # 0: no; 1: yes
my $protpars_jumbles       = 0;

my %seqnames       = ();           # number   =>  seqname 
my %numbers        = ();           # seqname  =>  number
my $options        = "";
my $infile         = "";
my $pwdfile        = "";
my $outfile        = "";
my $outfilenhx     = "";
my $logfile        = "";
my $alignfile      = "";
my $multitreefile  = "";
my $distancefile   = "";
my $log            = "";
my $number_of_aa   = 0;
my $orig_length    = 0;
my $ii             = 0;
my $temp_dir       = "";
my $current_dir    = "";
my @out            = ();
my $number_of_seqs = 0;



unless ( @ARGV == 2 || @ARGV == 3 || @ARGV == 4 || @ARGV == 5 ) {
    &printUsage();
    exit ( -1 ); 
}



# Analyzes the options:
# ---------------------

if ( $ARGV[ 0 ] =~ /^-.+/ ) {
    
    unless ( @ARGV > 2 ) {
        &printUsage();
        exit ( -1 ); 
    }
    $options = $ARGV[ 0 ];
    
    if ( $options =~ /F/ && $options !~ /U/ ) {
        if ( @ARGV != 3 && @ARGV != 4 ) {
            &printUsage();
            exit ( -1 ); 
        
        }
        $start_with_pwd = 1;
        $infile         = "";
        $pwdfile        = $ARGV[ 1 ];
        
        $outfile = $ARGV[ 2 ];
        if ( @ARGV == 4 ) {
            $temp_dir = $ARGV[ 3 ];
        }
        
    }
    elsif ( $options =~ /F/ && $options =~ /U/ ) {
        if ( @ARGV != 4 && @ARGV != 5 ) {
            &printUsage();
            exit ( -1 ); 
        }
        $start_with_pwd         = 1;
        $start_with_pwd_and_aln = 1;
        $pwdfile        = $ARGV[ 1 ];
        $infile         = $ARGV[ 2 ];
        $outfile = $ARGV[ 3 ];
        if ( @ARGV == 5 ) {
            $temp_dir = $ARGV[ 4 ];
        }
        
    }
    else {
        if ( @ARGV != 3 && @ARGV != 4 ) {
            &printUsage();
            exit ( -1 ); 
        }
        $infile  = $ARGV[ 1 ];
        $outfile = $ARGV[ 2 ];
        if ( @ARGV == 4 ) {
            $temp_dir = $ARGV[ 3 ];
        }
    }
    
    if ( $options =~ /N/ && $start_with_pwd != 1 ) {
        $remove_gaps    = 1; # do remove gaps 
    }
    if ( $options =~ /B(\d+)/ ) {
        $bootstraps = $1;
        if ( $bootstraps <= 1 ) {
            $bootstraps = 0;
        }
        elsif ( $bootstraps <= 9 ) {
            $bootstraps = 0;
            print "\n\nMAKETREE: WARNING: Bootstrap number must be devisable by 10,\nno bootstrapping.\n\n";
        }   
        elsif ( $bootstraps % 10 != 0 ) {
            $bootstraps = $bootstraps - $bootstraps % 10; # to ensure $bootstraps % 10 == 0
            print "\n\nMAKETREE: WARNING: Bootstrap number must be devisable by 10,\nhas been set to $bootstraps.\n\n";
        }    
    }
    if ( $options =~ /A/ ) {
        $use_protpars = 1 # PROTPARS
    }
    if ( $options =~ /j(\d+)/ ) {
        $protpars_jumbles = $1;
        if ( $protpars_jumbles < 0 ) {
            $protpars_jumbles = 0;
        }
    }
    if ( $options =~ /J/ ) {
        $matrix = 0;      # JTT
    }
    if ( $options =~ /L/ ) {
        $matrix = 2;      # Blossum
    }
    if ( $options =~ /M/ ) {
        $matrix = 3;      # mtREV24
    }
    if ( $options =~ /T/ ) {
        $matrix = 5;      # VT
    }
    if ( $options =~ /W/ ) {
        $matrix = 6;      # WAG
    }
    if ( $options =~ /P/ ) {
        $matrix = 7;      # auto
    }
    if ( $options =~ /R/ ) {
        $randomize_input_order = 1;
    }
    if ( $options =~ /S(\d+)/ ) {
        $seed = $1; 
    }
    if ( $options =~ /U/ ) {
        $puzzle_consensus_tree = 1;
    }
    if ( $options =~ /X/ ) {
        $keep_multiple_trees = 1; 
    }
    if ( $options =~ /D/ && $start_with_pwd != 1 ) {
        $keep_distance_matrix = 1; 
    }
    if ( $options =~ /V/ ) {
        $verbose = 1; 
    }
    if ( $options =~ /C/ && $start_with_pwd != 1 ) {
        $pairwise_dist_only = 1; 
    }
    if ( $options =~ /E/ ) {
        $exact_parameter_est = 1; 
    }
    if ( $options =~ /g/ ) {
        $rate_heterogeneity = 1; 
    }
    if ( $options =~ /t/ ) {
        $rate_heterogeneity = 2; 
    }
    if ( $options =~ /m/ ) {
        $rate_heterogeneity = 3; 
    }
    if ( $options =~ /#/ ) {
        $no_consenus_tree = 1; 
    }
    if ( $protpars_jumbles > 0 && $use_protpars != 1 ) {
        &printUsage();
        exit ( -1 ); 
    }
    if ( $use_protpars == 1 ) {
        if ( $randomize_input_order >= 1
        || $start_with_pwd == 1 
        || $keep_distance_matrix == 1
        || $pairwise_dist_only == 1 ) {
            &printUsage();
            exit ( -1 );
        }
        if ( $bootstraps > 1 && $protpars_jumbles < 1 ) {
            $protpars_jumbles = 1;
        }
    }
    
}

else {
    unless ( @ARGV == 2 || @ARGV == 3 ) {
        &printUsage();
        exit ( -1 );  
    }
    $infile  = $ARGV[ 0 ];
    $outfile = $ARGV[ 1 ];
    if ( @ARGV == 3 ) {
        $temp_dir = $ARGV[ 2 ];
    } 
}




$current_dir = `pwd`;
$current_dir =~ s/\s//;

if ( $outfile !~ /^\// ) {
    # outfile is not absolute path.
    $outfile = $current_dir."/".$outfile;
}



if ( $pairwise_dist_only == 1 ) {
    $bootstraps            = 0;
    $keep_multiple_trees   = 0;
    $puzzle_consensus_tree = 0;
    $randomize_input_order = 0;
    $start_with_pwd        = 0;
    $keep_distance_matrix  = 1; 
}

if ( $bootstraps < 2 ) {
    $no_consenus_tree = 0;
}

# TREE-PUZZLE sets the option in this way:
# If two rates or mixed, exact parameter estimates are used.
if ( $rate_heterogeneity == 2
||   $rate_heterogeneity == 3 ) {
    $exact_parameter_est = 1
}

$logfile       = $outfile.$LOG_FILE_SUFFIX;
$alignfile     = $outfile.$ALIGN_FILE_SUFFIX;
$multitreefile = $outfile.$MULTIPLE_TREES_FILE_SUFFIX;
$distancefile  = $outfile.$SUFFIX_PWD_NOT_BOOTS;

if ( $outfile =~ /\.nhx$/i ) {
    $outfilenhx    = $outfile;
    $logfile       =~ s/\.nhx//i;
    $alignfile     =~ s/\.nhx//i;
    $outfile       =~ s/\.nhx//i;
    $multitreefile =~ s/\.nhx//i;
    $distancefile  =~ s/\.nhx//i;
}  
else {
    $outfilenhx    = $outfile.".nhx";
}

if ( -e $outfilenhx ) {
    die "\n\nmakeTree: \"$outfilenhx\" already exists.\n\n";
}
if ( $infile ne "" ) {
    unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
        die "\n\nmakeTree: Input alignment file \"$infile\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
}
if ( $start_with_pwd == 1 ) {
    unless ( ( -s $pwdfile ) && ( -f $pwdfile ) && ( -T $pwdfile ) ) {
        die "\n\nmakeTree: Pairwise distance file \"$pwdfile\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
}



# Prints out the options:
# -----------------------


$log = "\n$0 logfile:\n";
$log = $log."Version: $VERSION\n\n";


if ( $start_with_pwd == 1 ) {
    $log = $log."Input pairwise distance file (bootstrapped): $pwdfile\n";
}
if ( $infile ne ""  ) {
    $log = $log."Input alignment                     : $infile\n";
}

if ( $no_consenus_tree != 1 ) {
    $log = $log."Output tree file                    : $outfilenhx\n";
}

if ( $keep_multiple_trees == 1 && $bootstraps >= 2 ) {
    $log = $log."Output  multiple trees file         : $multitreefile\n";
}
if ( $keep_distance_matrix ) {
    $log = $log."Output pairwise distance file       : $distancefile\n";
}

$log = $log."Bootstraps                          : $bootstraps\n";

if ( $start_with_pwd != 1 && $use_protpars != 1 ) {    
    $log = $log."Prgrm to calculate pairwise dist.   : TREE-PUZZLE\n";
}

if ( $use_protpars == 1 ) {
    $log = $log."Program to calculate tree           : PHYLIP PROTPARS\n";
    $log = $log."Number of jumbles in PROTPARS       : $protpars_jumbles\n";
}
else {
    $log = $log."Program to calculate tree           : PHYLIP NEIGHBOR (NJ)\n";
}

if ( $puzzle_consensus_tree == 1 ) {
    $log = $log."Prgrm to calculate ML branch lenghts: TREE-PUZZLE\n";
}
if ( $puzzle_consensus_tree == 1 || $start_with_pwd != 1 ) {
    $log = $log."Model                               : ";
    if ( $matrix == 0 ) { 
        $log = $log."JTT (Jones et al. 1992)\n";
    }
    elsif ( $matrix == 2 ) {
        $log = $log."BLOSUM 62 (Henikoff-Henikoff 92)\n";
    }
    elsif ( $matrix == 3 ) {
        $log = $log."mtREV24 (Adachi-Hasegawa 1996)\n";
    }
    elsif ( $matrix == 5 ) {
        $log = $log."VT (Mueller-Vingron 2000)\n";
    }
    elsif ( $matrix == 6 ) {
        $log = $log."WAG (Whelan-Goldman 2000)\n";
    }
    elsif ( $matrix == 7 ) {
        $log = $log."auto\n";
    }
    else {
        $log = $log."PAM (Dayhoff et al. 1978)\n";
    }
}
$log = $log."Model of rate heterogeneity         : ";
if ( $rate_heterogeneity == 1 ) { 
    $log = $log."8 Gamma distributed rates\n";
}
elsif ( $rate_heterogeneity == 2 ) {
    $log = $log."Two rates (1 invariable + 1 variable)\n";
}
elsif ( $rate_heterogeneity == 3 ) {
    $log = $log."Mixed (1 invariable + 8 Gamma rates)\n";
}
else {
    $log = $log."Uniform rate\n";
}
if ( $randomize_input_order >= 1 ) {
    $log = $log."Randomize input order in NEIGHBOR   : yes\n";
}
$log = $log."Seed for random number generators   : $seed\n";
if ( $exact_parameter_est == 1 ) {
    $log = $log."Exact parameter estimates in TREE-PUZZLE\n";
}

$log = $log."Start time/date                     : ".`date`;




# That's where the mischief starts....
# ------------------------------------

$ii = 0;

my $time_st = time;

if ( $temp_dir eq "" ) {
    $temp_dir = $TEMP_DIR_DEFAULT;
}

$temp_dir = $temp_dir.$time_st.$ii;

while ( -e $temp_dir ) {
    $ii++;
    $temp_dir = $temp_dir.$time_st.$ii;
}

mkdir(  $temp_dir, 0700 )
|| die "\n\n$0: Unexpected error: Could not create <<$temp_dir>>: $!\n\n";

unless ( ( -e $temp_dir ) && ( -d $temp_dir ) ) {
    die "\n\n$0: Unexpected error: <<$temp_dir>> does not exist, or is not a directory.\n\n";
}


if ( $start_with_pwd != 1 ) {
    system( "cp", $infile, $temp_dir."/INFILE" );
    unless ( chmod ( 0600, $temp_dir."/INFILE" ) ) {
        warn "\n\n$0: Could not chmod. $!\n\n";
    }
    $infile = "INFILE";
}


chdir ( $temp_dir ) 
|| die "\n\n$0: Unexpected error: Could not chdir to <<$temp_dir>>: $!\n\n";


if ( $start_with_pwd != 1 ) {

    @out = &DoPfam2phylip( $infile, $alignfile, $remove_gaps );
    $number_of_aa   = $out[ 0 ];
    $orig_length    = $out[ 1 ];
    $number_of_seqs = $out[ 2 ];

    system( "cp", $alignfile, "infile" );
    
    if ( $use_protpars != 1 ) {
        # Calculating the pairwise distances (saved in file "infile"): "puzzle"

        system( "cp", $alignfile, "align" ); 

        if ( $bootstraps > 1 ) {

            &executeSeqboot( $seed, $bootstraps );

            if ( $keep_distance_matrix ) {
                system( "mv", "outfile", "outfile___" );
                system( "cp", "align", "infile" );
                &executePuzzle( "infile",
                                $matrix,
                                $exact_parameter_est,
                                $rate_heterogeneity );
                system( "mv", "infile.dist", $distancefile );
                system( "mv", "outfile___", "outfile" );
            }
            unlink( "infile" ); # Necessary, since "infile" is puzzle's default input. 
            system( "mv", "outfile", "IN" );

            &executePuzzleBootstrapped( "IN",
                                        $matrix,
                                        $exact_parameter_est,
                                        $rate_heterogeneity );

            $pwdfile = "IN.dist";

        }
        else {

            &executePuzzle( "infile",
                            $matrix,
                            $exact_parameter_est,
                            $rate_heterogeneity );

            if ( $keep_distance_matrix ) {
                system( "cp outdist $distancefile" );
            } 
            $pwdfile = "infile.dist";
        }

        unlink( "infile.tree" );

        if ( $pairwise_dist_only == 1 ) {
            unlink( "infile", "align", "INFILE", "outdist", $alignfile );
            chdir( $current_dir ) 
            || die "\n\n$0: Unexpected error: Could not chdir to <<$current_dir>>: $!\n\n";

            rmdir( $temp_dir )
            || die "\n\n$0: Unexpected error: Could not remove <<$temp_dir>>: $!\n\n";

            print "\n\n$0 finished.\n\n";
            print "Output pairwise distance file written as: $distancefile\n\n";
            print "\n\nmakeTree successfully terminated.\n\n";
            exit( 0 );
        }
    
    } ## if ( $use_protpars != 1 )

} ## if ( $start_with_pwd != 1 )

 
# Calculating the tree (saved in file "infile"):

if ( $use_protpars != 1 ) {
    unlink( "infile" );
    &executeNeighbor( $pwdfile, $bootstraps, $randomize_input_order, $seed, 1 );
}
else {
    if ( $bootstraps > 1 ) {
        &executeSeqboot( $seed, $bootstraps );
        unlink( "infile" );
        system( "mv", "outfile", "infile" );
    }    
    &executeProtpars( "infile", $bootstraps, $protpars_jumbles, $seed );
}

unlink( "outfile" );

if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
   
    system( "cp", "outtree", $multitreefile );
}


system( "mv", "outtree", "intree" );

if ( $bootstraps > 1 ) {
    if ( $no_consenus_tree != 1 ) {

        # Consense:
        &executeConsense( "intree" );

        if ( $puzzle_consensus_tree == 1 ) { 

            system( "cp", "outtree", "treefile_consense" );
            system( "mv", "outtree", "intree" );

            # Puzzle for ML branch lenghts:
            # The alignment is read from infile by default.
            # The tree is read from intree by default.

            if ( $start_with_pwd_and_aln == 1 ) {
                &pfam2phylipMatchOnly( $infile,
                                       "infile",
                                       1 );
            }
            elsif ( $use_protpars != 1 ) {
                system( "mv", "align", "infile" ); # align = original alignment in phylip interleaved.
            } 

            &executePuzzleToCalculateBranchLenghts( $matrix,
                                                    $exact_parameter_est,
                                                    $rate_heterogeneity );

            unlink( "outfile", "outdist" );
            system( "mv", "outtree", "outree_puzzle" );

            # Transfer
            &executeTransfersBranchLenghts( "outree_puzzle", "treefile_consense", $outfilenhx );

        }
        else {
            unlink( "outfile", "align" );
            system( "mv", "outtree", $outfilenhx );
        }
    }
    else {
        unlink( "outfile", "align" );
    
    }
}
else {
    unlink( "align", "infile.dist" );
    if ( $start_with_pwd != 1 ) {
        system( "mv intree $outfilenhx" );
    }

}


unlink( "treefile_consense", "outtree", "outree_puzzle",
        "infile", "intree", "align", "INFILE", "IN", "IN.dist", "outdist"  );


$log = $log."Finish time/date                    : ".`date`;

if ( $start_with_pwd != 1 ) {
    $log = $log."Removed gaps                        : ";
    if ( $remove_gaps == 1 ) {
        $log = $log."yes\n";
    }
    else {
        $log = $log."no\n";
    }
    $log = $log."Columns in alignment used           : $number_of_aa\n";
    $log = $log."Columns in original alignment       : $orig_length\n";
    $log = $log."Number of sequences in alignment    : $number_of_seqs\n";
}


open( OUT, ">$logfile" ) || die "\n$0: Cannot create file <<$logfile>>: $!\n";
print OUT $log;
close( OUT );


chdir( $current_dir ) 
|| die "\n\n$0:Unexpected error: Could not chdir to <<$current_dir>>: $!\n\n";


rmdir( $temp_dir )
|| die "\n\n$0:Unexpected error: Could not remove <<$temp_dir>>: $!\n\n";

if ( $verbose == 1 ) {
    print "\n\n$0 finished.\n";
    if ( $no_consenus_tree != 1 ) {
        print "Output tree written as    : $outfilenhx\n";
    }
    print "Log written as            : $logfile\n";
    if ( $start_with_pwd != 1 ) {
        print "Alignment written as      : $alignfile\n";
    }
    if ( $keep_multiple_trees == 1 && $bootstraps >= 2 ) {
        print "Multiple trees written as : $multitreefile\n";
    }
    if ( $keep_distance_matrix ) {
        print "Distance matrix written as: $distancefile\n";
    }
}


exit( 0 ); 
    




# Methods:
# --------




# Executes pfam2phylip.
# If resulting alignment is too short due to the removal
# of gaps, is does not remove gaps.
# Three arguments:
# 1. infile
# 2. outfile
# 3. remove gaps: 1 to remove gaps; 0: do not remove gaps 
# Last modified: 06/04/01
sub DoPfam2phylip {
    my $in         = $_[ 0 ]; 
    my $out        = $_[ 1 ];
    my $option     = $_[ 2 ];
    my $aa         = 0;
    my @output     = ();

    if ( $option == 1 ) {
        @output = &pfam2phylip( $in, $out, 1 );
        $aa = $output[ 0 ];
        if ( $aa < 0 ) {
            die "\n\n$0: DoPfam2phylip: Unexpected error.\n\n";
        }
        if ( $aa < $MIN_NUMBER_OF_AA ) {
            unlink( $out );
            $option      = 0;
            $remove_gaps = 0;
        }
    }
    if ( $option == 0 ) {    # Must be another "if" (no elsif of else)!
        @output = &pfam2phylip( $in, $out, 2 );
        # 2 is to substitute non-letters with "-" in the sequence.
        $aa = $output[ 0 ];
        if ( $aa <= 0 ) {
             die "\n\n$0: DoPfam2phylip: Unexpected error.\n\n";
        }
    }
    return @output;
}



# Two arguments:
# 1. seed for random number generator
# 2. number of bootstraps
# Reads in "infile" by default.
sub executeSeqboot {

    my $s    = $_[ 0 ];
    my $bs   = $_[ 1 ];
    my $verb = "";
    
    &testForTextFilePresence( $infile );

    if ( $verbose != 1 ) {
        $verb = "
2";
    }

   
    system( "$SEQBOOT << !
r
$bs$verb
Y
$s
!" )
    && die "$0: Could not execute \"$SEQBOOT\"";
    
    return;
    
}




# One/two/three argument(s):
# Reads in tree from "intree" by default. (Presence of "intree" automatically 
# switches into "User defined trees" mode.)
# 1. matrix option: 0 = JTT; 2 = BLOSUM 62; 3 = mtREV24;
#    5 = VT; 6 = WAG; 7 = auto; PAM otherwise
# 2. Parameter estimates: 1 for "Exact (slow)"; "Approximate (faster)" otherwise
# 3. Model of rate heterogeneity:
#    1 for "8 Gamma distributed rates"
#    2 for "Two rates (1 invariable + 1 variable)"
#    3 for "Mixed (1 invariable + 8 Gamma rates)"
#    otherwise: Uniform rate
# Last modified: 09/08/03 (added 2nd and 3rd parameter)
sub executePuzzleToCalculateBranchLenghts {
    my $matrix_option              = $_[ 0 ];
    my $parameter_estimates_option = $_[ 1 ];
    my $rate_heterogeneity_option  = $_[ 2 ];
    my $i             = 0;
    my $mat           = "";
    my $est           = "";
    my $rate          = "";
    
    unless ( ( -s "infile" ) && ( -f "infile" ) && ( -T "infile" ) ) {
        die "\n$0: executePuzzleToCalculateBranchLenghts: <<infile>> does not exist, is empty, or is not a plain textfile.\n";
    }
    unless ( ( -s "intree" ) && ( -f "intree" ) && ( -T "intree" ) ) {
        die "\n$0: executePuzzleToCalculateBranchLenghts: <<intree>> does not exist, is empty, or is not a plain textfile.\n";
    }

    $mat = setModelForPuzzle( $matrix_option );
    if ( $parameter_estimates_option ) {
        $est = &setParameterEstimatesOptionForPuzzle( $parameter_estimates_option );
    }
    if ( $rate_heterogeneity_option ) {
        $rate = &setRateHeterogeneityOptionForPuzzle( $rate_heterogeneity_option );
    }
    
    system( "$PUZZLE << !
$mat$est$rate
x
y
!" )
    && die "$0: Could not execute \"$PUZZLE\"";
    
    return;
    
}







# Three/four arguments:
# 1. Name of file containing tree with correct branch lengths
# 2. Name of file containing tree with correct bootstraps
# 3. Outputfilename
# 4. R to reroot both trees in the same manner (use for FITCH,
#    since this changes to rooting.
sub executeTransfersBranchLenghts {
    my $tree_with_bl = $_[ 0 ];
    my $tree_with_bs = $_[ 1 ];
    my $out          = $_[ 2 ];
    my $reroot       = $_[ 3 ];
    my $R            = "";

    if ( $reroot && $reroot eq "R" ) {
        $R = "R";
    }

    &testForTextFilePresence( $tree_with_bl );
    &testForTextFilePresence( $tree_with_bs );
    
    system( "$TRANSFERSBRANCHLENGHTS $tree_with_bl $tree_with_bs $out $R" )
    && die "$0: Could not execute \"$TRANSFERSBRANCHLENGHTS $tree_with_bl $tree_with_bs $out $R\"";
    
    
    return;
}



# Called by method DoPfam2phylip.
# This reads a multiple sequence alignment file in Pfam format,
# Phylip's sequential format, or ClustalW (".aln")output and saves them
# in Phylip's sequential or interleaved format.
# (Those two are the same in this case, since all the seqs will be
# one line in length (no returns)).
# It returns (1st) the number of aa (columns) in the resulting
# alignment and the (2nd) number of aa (columns) in the original
# alignment.
#
# Reads a file containing a sequence alignment in the following format
# (as used in Pfam):
#  #comments      <- empty lines and lines begining with # (not mandatory)
#  name1 kal
#  name2 kal
#                 <- at least one empty line between blocks
#  name1 kale
#  name2 k.le
#
# Saves it in the "sequential" format of phylip:
#  number of OTUs length of aa seqs
#  name1     kalkale
#  name2     kalk-le
#
# Three arguments:
# 1. infile name
# 2. outfile name
# 3. 1  : Removes colums with a gap (non-letter character)
#    2  : Substitutes non-letter characters (except "?") in the sequence with "-".
#
# Last modified: 03/24/04
# Changes:
# 03/24/04: Do not replace "?" with "-"
#
sub pfam2phylip { 

    my $infile              = $_[ 0 ];
    my $outfile             = $_[ 1 ];
    my $options             = $_[ 2 ]; # 1: remove gaps; 2: non-letters (except "?") -> "-"
    my $return_line         = "";
    my $x                   = 0;
    my $y                   = 0;
    my $x_offset            = 0;
    my $original_length     = 0;
    my @seq_name            = ();
    my @seq_array           = ();
    my $seq                 = "";
    my $max_x               = 0;
    my $max_y               = 0;
    my $m                   = 0;
    my $n                   = 0;
    my $i                   = 0;
    my $move                = 0;
    my $saw_a_sequence_line = 0;

    if ( -e $outfile ) {
        die "\n$0: pfam2phylip: <<$outfile>> already exists.\n";
    }

    &testForTextFilePresence( $infile );

    open( INPP, "$infile" ) || die "\n$0: pfam2phylip: Cannot open file <<$infile>>: $!\n";

    until ( $return_line !~ /^\s*\S+\s+\S+/ && $saw_a_sequence_line == 1 ) {
        if ( $return_line =~ /^\s*\S+\s+\S+/ 
        && $return_line !~ /^\s*#/ 
        && $return_line !~ /^\s*\d+\s+\d+/
        && $return_line !~ /^\s*CLUSTAL/ ) {
            $saw_a_sequence_line = 1;
            $return_line =~ /^\s*(\S+)\s+(\S+)/;
            $seq_name[ $y ] = $1;
            $seq = $2;
            $seq_name[ $y ] = substr( $seq_name[ $y ], 0, $LENGTH_OF_NAME );
           
            for ( $x = 0; $x <= length( $seq ) - 1; $x++ ) {
                $seq_array[ $x ][ $y ] = substr( $seq, $x, 1 );
            }
            if ( $x_offset < length( $seq ) ) {
                $x_offset = length( $seq );
            }
            $y++;
        }
        $return_line = <INPP>;
        if ( !$return_line ) {
            last;
        }
    }

    $max_y = $y;
    $y     = 0;
    $max_x = 0;

    while ( $return_line = <INPP> ) {
        if ( $return_line =~ /^\s*(\S+)\s+(\S+)/
        && $return_line !~ /^\s*#/ 
        && $return_line !~ /^\s*\d+\s+\d+/ ) {
            $return_line =~ /^\s*\S+\s+(\S+)/;
            $seq = $1;
            for ( $x = 0; $x <= length( $seq ) - 1; $x++ ) {
                $seq_array[ $x + $x_offset ][ $y % $max_y ] = substr( $seq, $x, 1 );
            }
            if ( $max_x < length( $seq ) ) {
                $max_x = length( $seq );
            }
            $y++;
            if ( ( $y % $max_y ) == 0 ) {
                $y = 0;
                $x_offset = $x_offset + $max_x;
                $max_x = 0;
            }
        }
    }
    $original_length = $x_offset;
    close( INPP );


    # Removes "gap-columns" (gaps = everything except a-z characters):
    if ( $options == 1 ) {
        $move = 0;

        COLUMN: for ( $x = 0; $x <= $x_offset - 1; $x++ ) {  # goes through all aa positions (columns)

            for ( $y = 0; $y <= $max_y - 1; $y++ ) { # goes through all aas in a particular position

                unless ( $seq_array[ $x ][ $y ] && $seq_array[ $x ][ $y ] =~ /[a-z]/i ) {
                    $move++;
                    next COLUMN;
                }
            }

            # If this point is reached, column must be OK = no gaps.
            if ( $move >= 1 ) {
                for ( $m = 0; $m <= $max_y; $m++ ) {       
                    for ( $n = $x; $n <= $x_offset; $n++ ) {
                        $seq_array[ $n - $move ][ $m ] = $seq_array[ $n ][ $m ];
                    }
                }  
                $x_offset = $x_offset - $move;
                $x = $x - $move;
                $move = 0;
            }
        }
        if ( $move >= 1 ) {
            for ( $m = 0; $m <= $max_y; $m++ ) {       
                for ( $n = $x; $n <= $x_offset; $n++ ) {
                    $seq_array[ $n - $move ][ $m ] = $seq_array[ $n ][ $m ];
                }
            }   
            $x_offset = $x_offset - $move;
            $x = $x - $move;
            $move = 0;
        }
    }


    # Writes the file:

    open( OUTPP, ">$outfile" ) || die "\n$0: pfam2phylip: Cannot create file <<$outfile>>: $!\n";
    print OUTPP "$max_y $x_offset\n";
    for ( $y = 0; $y < $max_y; $y++ ) {
        print OUTPP "$seq_name[ $y ]";
        for ( $i = 0; $i <= ( $LENGTH_OF_NAME - length( $seq_name[ $y ] ) - 1 ); $i++ ) {
	        print OUTPP " ";
        }
        for ( $x = 0; $x <= $x_offset - 1; $x++ ) {
            if ( $options == 2 ) {
                if ( $seq_array[ $x ][ $y ] ) {
                    $seq_array[ $x ][ $y ] =~s /[^a-zA-Z\?]/-/;
                }
                else {
                    $seq_array[ $x ][ $y ] = "-";
                }
            }
            print OUTPP "$seq_array[ $x ][ $y ]";
        }
        print OUTPP "\n";
    }  
    close( OUTPP );

    return $x_offset, $original_length, $max_y;

} ## pfam2phylip




sub printUsage {

    print "\n";
    print " makeTree.pl  version $VERSION\n";
    print " -----------\n";

    print <<END;

 Copyright (C) 1999-2003 Washington University School of Medicine
 and Howard Hughes Medical Institute
 All rights reserved

 Author: Christian M. Zmasek 
         zmasek\@genetics.wustl.edu
         http://www.genetics.wustl.edu/eddy/people/zmasek/

 Last modified 09/06/03


  Requirements  makeTree is part of the RIO/FORESTER suite of programs.
  ------------  Many of its global variables are set via rio_module.pm.



  Note. Use xt.pl (for Pfam alignments) or mt.pl (for other alignments) 
  to run makeTree.pl on whole directories of alignments files. 



  Usage
  -----

  Tree calculation based on a Pfam/Clustal W alignment
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      makeTree.pl [-options] <input alignment in SELEX (Pfam), PHYLIP
      sequential format, or Clustal W output> <outputfile>
      [path/name for temporary directory to be created]

     Example:
     "% makeTree.pl -UTB1000S41NDXV /DB/PFAM/Full/IL5 IL5_tree"


  Tree calculation based on precalculated pairwise distances
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Consensus tree will have no branch length values.
  Precalculated pairwise distances are the output of "pfam2pwd.pl",
  number of bootstraps needs to match the one used for the pwds.

     makeTree.pl <-options, includes "F"> <pwdfile: boostrapped pairwise
     distances> <outputfile> [path/name for temporary directory
     to be created]

     Example:
     "% makeTree.pl -FB100S21XV /pfam2pwd_out/IL5.pwd IL5_tree"


  Tree calculation based on precalculated pairwise distances
  and matching alignment
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Consensus tree will have branch length values.
  Precalculated pairwise distances and the matching (processed)
  alignment are the output of "pfam2pwd.pl", number of bootstraps
  needs to match the one used for the pwds, matrix needs to match
  the one used for the pwds.

     makeTree.pl <-options, includes "UF"> <pwdfile: boostrapped pairwise
     distances> <alnfile: corresponding alignment>
     <outputfile> [path/name for temporary directory to be created]
         
     Example:
     "% makeTree.pl -UFLB100S21XV /pfam2pwd_out/IL5.pwd /pfam2pwd_out/IL5.aln IL5_tree"


  Options
  -------

  N  : Suggestion to remove columns in the alignment which contain gaps.
       Gaps are not removed, if, after removal of gaps, the resulting alignment would
       be shorter than $MIN_NUMBER_OF_AA. Default is not to remove gaps.
  Bx : Number of bootstrapps. B0: do not bootstrap. Default is 100 bootstrapps.
       The number of bootstrapps should be divisible by 10.
  U  : Use TREE-PUZZLE to calculate ML branchlengths for consesus tree, in case of 
       bootstrapped analysis.
  J  : Use JTT matrix (Jones et al. 1992) in TREE-PUZZLE, default: PAM.
  L  : Use BLOSUM 62 matrix (Henikoff-Henikoff 92) in TREE-PUZZLE, default: PAM.
  M  : Use mtREV24 matrix (Adachi-Hasegawa 1996) inTREE-PUZZLE, default: PAM.
  W  : Use WAG matrix (Whelan-Goldman 2000) in TREE-PUZZLE, default: PAM.
  T  : Use VT matrix (Mueller-Vingron 2000) in TREE-PUZZLE, default: PAM.
  P  : Let TREE-PUZZLE choose which matrix to use, default: PAM.
  E  : Exact parameter estimates in TREE-PUZZLE, default: Approximate.
       Model of rate heterogeneity in TREE-PUZZLE (default: Uniform rate)
  g  : 8 Gamma distributed rates
  t  : Two rates (1 invariable + 1 variable)
  m  : Mixed (1 invariable + 8 Gamma rates)
  R  : Randomize input order in PHYLIP NEIGHBOR.
  A  : Use PHYLIP PROTPARS instead of NEIGHBOR (and no pairwise distance calculation).
  jx : Number of jumbles when using PHYLIP PROTPARS (random seed set with Sx).
  Sx : Seed for random number generator(s). Must be 4n+1. Default is 9.
  X  : To keep multiple tree file (=trees from bootstrap resampled alignments).
  D  : To keep (and create in case of bootstrap analysis) pairwise distance matrix file.
       This is created form the not resampled (original) alignment.
  C  : Calculate pairwise distances only (no tree). Bootstrap is always 1.
       No other files are generated.  
  F  : Pairwise distance (pwd) file as input (instead of alignment).
       No -D, -C, and -N options available in this case.
  V  : Verbose.
  #  : Only for rio.pl: Do not calculate consensus tree ("I" option in rio.pl).


END
  
} ## printUsage

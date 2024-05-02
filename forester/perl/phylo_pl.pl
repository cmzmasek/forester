#!/usr/bin/perl -W
#
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
#
#
#
#  Requirements  phylo_pl is part of the FORESTER libraries.
#  ------------  Many of its global variables are set via forester.pm.
#
#
#  Note. Use xt.pl (for Pfam alignments) or mt.pl (for other alignments)
#  to run phylo_pl.pl on whole directories of alignments files.       
#
#
#
#
# =========================
#
# METHOD ORDER (IMPORTANT!)
# 1. FastME
# 2. phylip NJ
# 3. phylip fitch FM
# 4. phylip fitch ME
# 5. Raxml
# 6. phylip proml
# 7. phylip protpars
# 8. all
#==========================

use strict;

use FindBin;
use lib $FindBin::Bin;
use forester;

my $VERSION                = "1.0.3";
my $LAST_MODIFIED          = "2020/05/23";

my $RAXML_MODEL_BASE_PROT  = "PROTCAT";
my $RAXML_MODEL_BASE_NUC   = "GTRCAT";
my $RAXML_ALGORITHM        = "a";
my $RAXML_T_OPTION         = "-T 4";


my $TEMP_DIR_DEFAULT       = "/tmp/phylo_pl_"; # Where all the infiles, outfiles, etc will be created.
   
my $bootstraps             = 100; # 0,1: do not bootstrap. Default: 100
my $matrix                 = 5;   # 0 = JTT
                                  # 1 = PAM 
                                  # 2 = BLOSUM 62
                                  # 3 = mtREV24
                                  # 5 = VT 
                                  # 6 = WAG 
                                  # 7 = auto in puzzle
                                  # 8 = DCMut in PHYML, VT in TREE-PUZZLE
                                  # 9 = HKY [na]
                                  # 10 = TN [na]
                                  # 11 = GTR [na]
                                  # 12 = SH [na]
                                  # 13 = LG
my $rate_heterogeneity     = 0;   # 0 = Uniform rate (default)
                                  # 1 = 8 Gamma distributed rates
                                  # 2 = Two rates (1 invariable + 1 variable)
                                  # 3 = Mixed (1 invariable + 8 Gamma rates)
my $seed                   = 9;   # Seed for random number generators. Default: 9
my $keep_multiple_trees    = 0;   # 0: delete multiple tree file
                                  # 1: do not delete multiple tree file
my $exact_parameter_est    = 0;   # 0: no; 1: yes

my $jumbles                = 2;
my $use_fastme             = 0;   # 0: no; 1: yes
my $use_phylip_nj          = 0;   # 0: no; 1: yes
my $use_phylip_fitch_fm    = 0;   # 0: no; 1: yes
my $use_phylip_fitch_me    = 0;   # 0: no; 1: yes
my $use_raxml              = 0;   # 0: no; 1: yes
my $use_proml              = 0;   # 0: no; 1: yes
my $use_protpars           = 0;   # 0: no; 1: yes
my $use_global_rearr       = 0;   # 0: no; 1: yes
my $estimate_invar_sites   = 0;   # 0: no; 1: yes

my $fastme_init_tree_opt   = "NJ";

my %seqnames       = ();           # number   =>  seqname 
my %numbers        = ();           # seqname  =>  number
my $options        = "";
my $infile         = "";
my $pwdfile        = "";
my $outfile        = "";
my $logfile        = "";
my $multipwdfile   = "";
my $distancefile   = "";
my $log            = "";
my $ii             = 0;
my $temp_dir       = "";
my $current_dir    = "";
my @out            = ();
my $number_of_seqs = 0;
my $number_of_aa   = 0;

my $use_pwd_based_methods = 0;

print( "\n");
print( "phylo_pl $VERSION ($LAST_MODIFIED)\n" );
print( "__________________________________\n");
print( "\n\n");



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
   
    if ( @ARGV != 3 && @ARGV != 4 ) {
        &printUsage();
        exit ( -1 ); 
    }
    $infile  = $ARGV[ 1 ];
    $outfile = $ARGV[ 2 ];
    if ( @ARGV == 4 ) {
        $temp_dir = $ARGV[ 3 ];
    }
    if ( $options =~ /B(\d+)/ ) {
        $bootstraps = $1;
        if ( $bootstraps <= 1 ) {
            $bootstraps = 0;
        }
        elsif ( $bootstraps <= 9 ) {
            $bootstraps = 0;
            print "\n\nphylo_pl: WARNING: Bootstrap number must be devisable by 10,\nno bootstrapping.\n\n";
        }   
        elsif ( $bootstraps % 10 != 0 ) {
            $bootstraps = $bootstraps - $bootstraps % 10; # to ensure $bootstraps % 10 == 0
            print "\n\nphylo_pl: WARNING: Bootstrap number must be devisable by 10,\nset to $bootstraps.\n\n";
        }    
    }
    if ( $options =~ /n/ ) {
        $use_phylip_nj = 1; 
    }
    if ( $options =~ /q/ ) {
        $use_fastme = 1;
    }
    if ( $options =~ /f/ ) {
        $use_phylip_fitch_fm = 1; 
    }
    if ( $options =~ /e/ ) {
        $use_phylip_fitch_me = 1; 
    }
    if ( $options =~ /x/ ) {
        $use_raxml = 1;
    }
    if ( $options =~ /o/ ) {
        $use_proml = 1;
    }
    if ( $options =~ /p/ ) {
        $use_protpars = 1;
    }
    if ( $options =~ /G/ ) {
        $use_global_rearr = 1;
    } 
    if ( $options =~ /I/ ) {
        $estimate_invar_sites = 1;
    }
    if ( $options =~ /j(\d+)/ ) {
        $jumbles = $1;
        if ( $jumbles < 1 ) {
            $jumbles = 0;
        }
    }
    if ( $options =~ /J/ ) {
        $matrix = 0;      # JTT
    }
    if ( $options =~ /P/ ) {
        $matrix = 1;      # PAM
    }
    if ( $options =~ /O/ ) {
        $matrix = 2;      # Blosum 62
    }
    if ( $options =~ /M/ ) {
        $matrix = 3;      # mtREV24
    }
    if ( $options =~ /W/ ) {
        $matrix = 6;      # WAG
    }
    if ( $options =~ /A/ ) {
        $matrix = 7;      # auto
    }
    if ( $options =~ /D/ ) {
        $matrix = 8;      # DCMut in RAXML, VT in PUZZLE
    }
    if ( $options =~ /H/ ) {
        $matrix = 9;      # HKY
    }
    if ( $options =~ /T/ ) {
        $matrix = 10;      # TN
    }
    if ( $options =~ /Z/ ) {
        $matrix = 11;      # GTR 
    }
    if ( $options =~ /C/ ) {
        $matrix = 12;      # SH
    }
    if ( $options =~ /L/ ) {
        $matrix = 13;      # LG
    }
    if ( $options =~ /S(\d+)/ ) {
        $seed = $1; 
    }
    if ( $options =~ /X/ ) {
        $keep_multiple_trees = 1; 
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

if ( $use_fastme    != 1 &&
     $use_phylip_nj != 1 &&
     $use_phylip_fitch_fm != 1 &&
     $use_phylip_fitch_me != 1 &&
     $use_raxml != 1 &&
     $use_proml != 1 &&         
     $use_protpars != 1 ) {
    
     $use_fastme = 1;
     $use_phylip_nj = 1;
     $use_phylip_fitch_fm = 1;
     $use_phylip_fitch_me = 1;
     $use_raxml = 1;
     $use_proml = 1;         
     $use_protpars = 1; 
}


if ( $use_fastme    == 1 ||
     $use_phylip_nj == 1 ||
     $use_phylip_fitch_fm == 1 ||
     $use_phylip_fitch_me == 1 ) { 
    $use_pwd_based_methods = 1;
}
else {
    $use_pwd_based_methods = 0;
}    

$current_dir = `pwd`;
$current_dir =~ s/\s//;

if ( $outfile !~ /^\// ) {
    # outfile is not absolute path.
    $outfile = $current_dir."/".$outfile;
}



# TREE-PUZZLE sets the option in this way:
# If two rates or mixed, exact parameter estimates are used.
if ( $rate_heterogeneity == 2
||   $rate_heterogeneity == 3 ) {
    $exact_parameter_est = 1
}


if ( $outfile =~ /\.xml$/i ) {
    $outfile =~ s/\.xml//i;
} 
elsif ( $outfile =~ /\.aln$/i ) {
    $outfile =~ s/\.aln//i;
}
elsif ( $outfile =~ /\.fasta$/i ) {
    $outfile =~ s/\.fasta//i;
}
elsif ( $outfile =~ /\.fas$/i ) {
    $outfile =~ s/\.fas//i;
}
elsif ( $outfile =~ /\.seqs$/i ) {
    $outfile =~ s/\.seqs//i;
}  


$logfile          = $outfile.$LOG_FILE_SUFFIX;
$multipwdfile     = $outfile.$MULTIPLE_PWD_FILE_SUFFIX;
$distancefile     = $outfile.$SUFFIX_PWD_NOT_BOOTS;

&dieIfFileExists( $logfile );
&dieIfFileExists( $multipwdfile );
&dieIfFileExists( $distancefile );


my $fastme_outtree    = $outfile."_fme.xml";
my $fastme_stats      = $outfile."_fme_stats.txt";
my $phylip_nj_outtree = $outfile."_pnj.xml";
my $phylip_fm_outtree = $outfile."_pfm.xml";
my $phylip_me_outtree = $outfile."_pme.xml";
my $raxml_outtree     = $outfile."_raxml.xml";
my $proml_outtree     = $outfile."_proml.xml";
my $protpars_outtree  = $outfile."_ppp.xml";
my $all_outtree       = $outfile."_comb.xml";

my $multitreefile_fastme    = $outfile."_fme".$MULTIPLE_TREES_FILE_SUFFIX;
my $multitreefile_phylip_nj = $outfile."_pnj".$MULTIPLE_TREES_FILE_SUFFIX;
my $multitreefile_phylip_fm = $outfile."_pfm".$MULTIPLE_TREES_FILE_SUFFIX;
my $multitreefile_phylip_me = $outfile."_pme".$MULTIPLE_TREES_FILE_SUFFIX;
my $multitreefile_raxml     = $outfile."_raxml".$MULTIPLE_TREES_FILE_SUFFIX;
my $multitreefile_proml     = $outfile."_proml".$MULTIPLE_TREES_FILE_SUFFIX;
my $multitreefile_protpars  = $outfile."_ppp".$MULTIPLE_TREES_FILE_SUFFIX;

if ( $use_fastme == 1 ) {
    &dieIfFileExists( $fastme_outtree );
    if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
        &dieIfFileExists( $multitreefile_fastme );
    }
}
if( $use_phylip_nj == 1 ) {
    &dieIfFileExists( $phylip_nj_outtree );
    if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
        &dieIfFileExists( $multitreefile_phylip_nj );
    }
}
if( $use_phylip_fitch_fm == 1 ) {
    &dieIfFileExists( $phylip_fm_outtree );
    if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
        &dieIfFileExists( $multitreefile_phylip_fm );
    }
}
if( $use_phylip_fitch_me == 1 ) {
    &dieIfFileExists( $phylip_me_outtree );
    if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
        &dieIfFileExists( $multitreefile_phylip_me );
    }
}
if( $use_raxml == 1 ) {
    &dieIfFileExists( $raxml_outtree );
    if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
        &dieIfFileExists( $multitreefile_raxml );
    }
}
if( $use_proml == 1 ) {
    &dieIfFileExists( $proml_outtree );
    if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
        &dieIfFileExists( $multitreefile_proml );
    }
} 
if( $use_protpars == 1 ) {
    &dieIfFileExists( $protpars_outtree );
    if ( $keep_multiple_trees == 1 && $bootstraps > 1 ) {
        &dieIfFileExists( $multitreefile_protpars );
    }
} 
if ( $bootstraps > 1 ) {
     &dieIfFileExists( $all_outtree );
}
if ( $infile ne "" ) {
    unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
        die "\n\nphylo_pl: Input alignment file \"$infile\" does not exist, is empty, or is not a plain textfile.\n\n";
    }
}




# Prints out the options:
# -----------------------

$log = "\n$0 logfile:\n";
$log = $log."Version: $VERSION\n\n";


if ( $infile ne ""  ) {
    $log = $log."Input                               : $infile\n";
}

if ( $keep_multiple_trees == 1 && $bootstraps >= 2 ) {
    $log = $log."Multiple distance matrices          : $multipwdfile\n";
}


$log = $log."Bootstraps                          : $bootstraps\n";

if ( $use_pwd_based_methods == 1 ) {    
    $log = $log."Prgrm to calculate pairwise dist.   : TREE-PUZZLE (version: $PUZZLE_VERSION)\n";
}


if ( $use_fastme == 1 ) {
    $log = $log."Program to calculate tree           : FastME (version: $FASTME_VERSION)\n";   
}
if ( $use_phylip_nj == 1 ) {
    $log = $log."Program to calculate tree           : PHYLIP NEIGHBOR NJ (version: $PHYLIP_VERSION)\n";
}
if ( $use_phylip_fitch_fm == 1 ) {
    $log = $log."Program to calculate tree           : PHYLIP FITCH Fitch-Margoliash (version: $PHYLIP_VERSION)\n";
}
if ( $use_phylip_fitch_me == 1 ) {
    $log = $log."Program to calculate tree           : PHYLIP FITCH Minimal Evolution (version: $PHYLIP_VERSION)\n";
}
if ( $use_raxml == 1 ) {
    $log = $log."Program to calculate tree           : RAxML-NG (version: $RAXMLNG_VERSION)\n";
}
if ( $use_proml == 1 ) {
    $log = $log."Program to calculate tree           : PHYLIP PROML (uses PAM unless JTT selected) (version: $PHYLIP_VERSION)\n";
}
if ( $use_protpars == 1 ) {
    $log = $log."Program to calculate tree           : PHYLIP PROTPARS (with global rearrangements) (version: $PHYLIP_VERSION)\n";
}
if ( $use_phylip_fitch_fm == 1 || $use_phylip_fitch_me == 1 || $use_protpars == 1 || $use_proml ) {
    $log = $log."Number of jumbles (input order rand): $jumbles\n"; 
    
}
if ( $use_phylip_fitch_fm == 1 || $use_phylip_fitch_me == 1 || $use_proml ) {
    if ( $use_global_rearr == 1 ) {
        $log = $log."Global rearrangements               : true\n";
    }
    else {
        $log = $log."Global rearrangements               : false\n";
        
    }
}

if ( $bootstraps > 0 ) {
    $log = $log."Prgrm to calculate ML branch lenghts: TREE-PUZZLE (version: $PUZZLE_VERSION)\n";
}

$log = $log."Model                               : ";
if ( $matrix == 0 ) { 
    $log = $log."JTT (Jones et al. 1992)\n";
}
elsif ( $matrix == 1 ) {
    $log = $log."PAM (Dayhoff et al. 1978)\n";
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
    $log = $log."auto in TREE-PUZZLE\n";
}
elsif ( $matrix == 8 ) {
    $log = $log."DCMut (Kosial and Goldman, 2005) in RAxML; VT in TREE-PUZZLE\n";
}
elsif ( $matrix == 9 ) {
    $log = $log."HKY (Hasegawa et al. 1985) in TREE-PUZZLE, GTR in RAxML\n";
}
elsif ( $matrix == 10 ) {
    $log = $log."TN (Tamura-Nei 1993) in TREE-PUZZLE, GTR in RAxML\n";
}
elsif ( $matrix == 11 ) {
    $log = $log."GTR (e.g. Lanave et al. 1980)in TREE-PUZZLE and RAxM\n";
}
elsif ( $matrix == 12 ) {
    $log = $log."SH (Schoeniger-von Haeseler 1994) in TREE-PUZZLE, GTR in RAxML\n";
}
elsif ( $matrix == 13 ) {
    $log = $log."LG model (Le and Gascuel, 2008) in RAxML; WAG in TREE-PUZZLE\n";
}
else {
    &dieWithUnexpectedError( "Unknown model: matrix=$matrix" );
}
if ( $use_raxml == 1 ) {
    $log = $log."Model of rate heterogeneity in RAxML: G8\n";
}
$log = $log."Model of rate heterogeneity (PUZZLE): ";
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
$log = $log."Seed for random number generators   : $seed\n";
if ( $exact_parameter_est == 1 ) {
    $log = $log."Exact parameter estimates in TREE-PUZZLE\n";
}

$log = $log."Start time/date                     : ".`date`;




# That's where the mischief starts....
# ------------------------------------

$ii = 0;

srand();
my $range = 1000000;
my $random_number = int( rand( $range ) );

if ( $temp_dir eq "" ) {
    $temp_dir = $TEMP_DIR_DEFAULT;
}

$temp_dir = $temp_dir.$random_number.$ii;

while ( -e $temp_dir ) {
    $ii++;
    $temp_dir = $temp_dir.$random_number.$ii;
}

mkdir( $temp_dir, 0700 )
|| die "\n\n$0: Could not create <<$temp_dir>>: $!\n\n";

unless ( ( -e $temp_dir ) && ( -d $temp_dir ) ) {
    die "\n\n$0: <<$temp_dir>> does not exist, or is not a directory.\n\n";
}


&cp( $infile, $temp_dir."/INFILE__" );
unless ( chmod ( 0600, $temp_dir."/INFILE__" ) ) {
    warn "\n\n$0: Could not chmod. $!\n\n";
}
$infile = "INFILE__";

chdir ( $temp_dir ) 
|| die "\n\n$0: Could not chdir to <<$temp_dir>>: $!\n\n";
&cp( $infile, "infile" );
@out = &getNumberOfSeqsAndAas( $infile );
$number_of_seqs = $out[ 0 ];
$number_of_aa   = $out[ 1 ];

my $SEQBOOT_OUTFILE = "seqboot_outfile"; 

if (  $bootstraps > 1 && ( $use_pwd_based_methods == 1
                                    || $use_proml == 1
                                    || $use_protpars == 1 ) ) {
    &executeSeqboot( $seed, $bootstraps );
    &mv( "outfile", $SEQBOOT_OUTFILE );
    &rm( "infile" );
}

&cp( $infile, "align" ); 

if ( $use_pwd_based_methods == 1 ) {
    # Calculating the pairwise distances (saved in file "infile"): "puzzle"
    if ( $bootstraps > 1 ) {
        &executePuzzleBootstrapped( $SEQBOOT_OUTFILE,
                                    $matrix,
                                    $number_of_seqs,
                                    $exact_parameter_est,
                                    $rate_heterogeneity );

        $pwdfile = $SEQBOOT_OUTFILE.".dist";
    }
    else {
        &executePuzzle( "infile",
                        $matrix,
                        $number_of_seqs,
                        $exact_parameter_est,
                        $rate_heterogeneity );
        $pwdfile = "infile.dist";
    }
} 

&rm( "infile" );

# Methods based on alignment
# --------------------------
my $OUTTREE_RAXML    = "outtree_rax";
my $OUTTREE_PROML    = "outtree_proml";
my $OUTTREE_PROTPARS = "outtree_protpars";

my $CONSENSUS_RAXML    = "consensus_raxml";
my $CONSENSUS_PROML    = "consensus_proml";
my $CONSENSUS_PROTPARS = "consensus_protpars";

my $OUTTREES_ALL = "outtrees_all";
my $all_count = 0;

if ( $use_raxml == 1 ) {
    my $model = "---";
    if ( $matrix == 0 ) { 
        $model = "JTT";
    }
    elsif ( $matrix == 1 ) {
        $model = "DAYHOFF";
    }
    elsif ( $matrix == 2 ) {
        $model = "BLOSUM62";
    }
    elsif ( $matrix == 3 ) {
        $model = "MTREV";
    }
    elsif ( $matrix == 5 ) {
        $model = "VT";
    }
    elsif ( $matrix == 6 ) {
        $model = "WAG";
    }
    elsif ( $matrix == 7 ) {
        $model = "VT";
    }
    elsif ( $matrix == 8 ) {
        $model = "DCMUT";
    }
    elsif ( $matrix == 13 ) {
        $model = "LG";
    }
    elsif ( $matrix == 9 ) {
        $model = "GTR";
    }
    elsif ( $matrix == 10 ) {
        $model = "TN";
    }
    elsif ( $matrix == 11 ) {
        $model = "GTR";
    }
    elsif ( $matrix == 12 ) {
        $model = "SH";
    }
    else {
        &dieWithUnexpectedError( "Unknown model: matrix=$matrix" );
    }
    print( "\n========== RAxML begin =========\n\n" );    

    # NOTE. RaxML does its own bootstrapping.
    &executeRaxmlNG( "align", $model."+G8+F", $bootstraps );


    print( "\n========== RAxML end =========\n\n" );
    
    &rm( "align.raxml.log" );

    &rm( "align.raxml.bestModel" );
    &rm( "align.raxml.mlTrees" );
    if ( $bootstraps > 1 ) {
        &rm( "align.raxml.bestTree" );
    }
    &rm( "align.raxml.bestTreeCollapsed" );
    &rm( "align.raxml.startTree" );
    &rm( "align.raxml.rba" );
    &rm( "align.raxml.reduced.phy" );
    &rm( "align.raxml.*.TMP" );
    &rm( "align.raxml.log" );
    &rm( "align.raxml.ckp" );

 print( "\n========== a=========\n\n" );
    #&rm( "RAxML_parsimonyTree.phylopl" );
    #&rm( $outfile."_raxml_info" );
    #&mv( "RAxML_info.phylopl", $outfile."_raxml_info" );

    if ($bootstraps > 1 ) {
    print( "\n========== b=========\n\n" );
        &mv( "align.raxml.support", $CONSENSUS_RAXML );
        #&rm( "RAxML_bipartitionsBranchLabels.phylopl" );
        &append( "align.raxml.bootstraps", $OUTTREES_ALL );
        #if ( $keep_multiple_trees == 1 ) {
        #    &mv( "RAxML_bootstrap.phylopl", $multitreefile_raxml );
        #}
        #else {
        #    &rm( "RAxML_bootstrap.phylopl" );
        #}
    }
    &rm( "align.raxml.bootstraps" );
  
    $all_count++;
  
}

if ( $use_proml == 1 ) {
    my $input = "";
    if ( $bootstraps > 1 ) {
        $input = $SEQBOOT_OUTFILE;
    }
    else {
        $input = "align";
    }    
    print( "\n========== PHYLIP PROML begin =========\n\n" );    
    # Five arguments:
    # 1. name of alignment file (in correct format!)
    # 2. number of bootstraps
    # 3. jumbles: 0: do not jumble; >=1 number of jumbles
    # 4. seed for random number generator
    # 5. 1 for PAM instead of JTT
    my $use_pam = 1;
    if ( $matrix == 0 ) {
        $use_pam = 0;
    }
    &executeProml( $input, $bootstraps, $jumbles, $seed, $use_pam, $use_global_rearr );
    print( "\n========== PHYLIP PROML end =========\n\n" );
    &mv( "outtree", $OUTTREE_PROML );
    &rm( "outfile" ); 
    if ( $bootstraps > 1 ) {
        &append( $OUTTREE_PROML, $OUTTREES_ALL );
        $all_count++;
    }
}


if ( $use_protpars == 1 ) {
    my $input = "";
    if ( $bootstraps > 1 ) {
        $input = $SEQBOOT_OUTFILE;
    }
    else {
       $input = "align";
    }    
    print( "\n========== PHYLIP PROTPARS begin =========\n\n" );    
    &executeProtpars( $input, $bootstraps, $jumbles, $seed );
    print( "\n========== PHYLIP PROTPARS end =========\n\n" );
    &mv( "outtree", $OUTTREE_PROTPARS );
    &rm( $outfile."_protpars_outfile" ); 
    &mv( "outfile", $outfile."_protpars_outfile" );
    if ( $bootstraps > 1 ) {
        &append( $OUTTREE_PROTPARS, $OUTTREES_ALL );
        $all_count++;
    }
}



# Methods based on PWD
# --------------------
my $OUTTREE_FASTME    = "outtree_fastme";
my $OUTTREE_PHYLIP_NJ = "outtree_phylip_nj";
my $OUTTREE_PHYLIP_FM = "outtree_phylip_fm";
my $OUTTREE_PHYLIP_ME = "outtree_phylip_me";


my $CONSENSUS_FASTME    = "consensus_fastme";
my $CONSENSUS_PHYLIP_NJ = "consensus_phylip_nj";
my $CONSENSUS_PHYLIP_FM = "consensus_phylip_fm";
my $CONSENSUS_PHYLIP_ME = "consensus_phylip_me";
my $CONSENSUS_ALL       = "consensus_all";


if ( $use_fastme == 1 ) {
    print( "\n========== FASTME begin =========\n\n" );
    &executeFastme( $pwdfile, $bootstraps, "output.tre" );
    print( "\n========== FASTME end ===========\n\n" );
    if ( $bootstraps > 1 ) {
        &mv( "seqboot_outfile.dist__fastme_stat.txt", $fastme_stats );
    }
    else {
        &mv( "infile.dist__fastme_stat.txt", $fastme_stats );
    }
    &rm( "seqboot_outfile.dist_" );
    &mv( "output.tre", $OUTTREE_FASTME );
    if ( $bootstraps > 1 ) {
        &append( $OUTTREE_FASTME, $OUTTREES_ALL ); 
        $all_count++;
    }
}
if ( $use_phylip_nj ) {
    print( "\n========== PHYLIP NEIGHBOR begin =========\n\n" );
    &executeNeighbor( $pwdfile, $bootstraps, $seed, 0 );
    print( "\n========== PHYLIP NEIGHBOR end =========\n\n" );
    &mv( "outtree", $OUTTREE_PHYLIP_NJ );
    &rm( "outfile" );
    if ( $bootstraps > 1 ) {
        &append( $OUTTREE_PHYLIP_NJ, $OUTTREES_ALL );
        $all_count++;
    }
}
if ( $use_phylip_fitch_fm ) {
    print( "\n========== PHYLIP FITCH FM begin =========\n\n" );
    &executeFitch( $pwdfile, $bootstraps, $seed, $jumbles, 0, "FM", $use_global_rearr );
    print( "\n========== PHYLIP FITCH FM end =========\n\n" );
    &mv( "outtree", $OUTTREE_PHYLIP_FM );
    &rm( "outfile" );
    if ( $bootstraps > 1 ) {
        &append( $OUTTREE_PHYLIP_FM, $OUTTREES_ALL );
        $all_count++;
    }    
}
if (  $use_phylip_fitch_me ) {
    print( "\n========== PHYLIP FITCH ME begin =========\n\n" );
    &executeFitch( $pwdfile, $bootstraps, $seed, $jumbles, 0, "ME", $use_global_rearr );
    print( "\n========== PHYLIP FITCH ME end =========\n\n" );
    &mv( "outtree", $OUTTREE_PHYLIP_ME );
    &rm( "outfile" );
    if ( $bootstraps > 1 ) {
        &append( $OUTTREE_PHYLIP_ME, $OUTTREES_ALL );
        $all_count++;
    } 
}


if ( $bootstraps > 1 ) {
    # Consense:
    if ( $use_fastme == 1 ) {
        &consense( $OUTTREE_FASTME, $CONSENSUS_FASTME );
    }
    if ( $use_phylip_nj == 1 ) {
        &consense( $OUTTREE_PHYLIP_NJ, $CONSENSUS_PHYLIP_NJ );
    }
    if ( $use_phylip_fitch_fm == 1 ) {
        &consense( $OUTTREE_PHYLIP_FM, $CONSENSUS_PHYLIP_FM );
    }
    if ( $use_phylip_fitch_me == 1 ) {
        &consense( $OUTTREE_PHYLIP_ME, $CONSENSUS_PHYLIP_ME );
    }
    if ( $use_proml == 1 ) {
        &consense( $OUTTREE_PROML, $CONSENSUS_PROML );
    } 
    if ( $use_protpars == 1 ) {
        &consense( $OUTTREE_PROTPARS, $CONSENSUS_PROTPARS );
    } 
    if ( $all_count > 1 ) {
        &consense( $OUTTREES_ALL, $CONSENSUS_ALL );
    }
    else {
        &rm( $OUTTREES_ALL );
    }
   
    my $INTREE_FOR_PUZZLE = "intree"; #why so serious?
    &rm( $INTREE_FOR_PUZZLE );
    system( "touch", $INTREE_FOR_PUZZLE )
    && die("\n\n$0: could not \"touch $INTREE_FOR_PUZZLE\": $!\n\n");
     
    if ( $use_fastme == 1 ) {
        &append( $CONSENSUS_FASTME, $INTREE_FOR_PUZZLE );
    }
    if ( $use_phylip_nj == 1 ) {
        &append( $CONSENSUS_PHYLIP_NJ, $INTREE_FOR_PUZZLE );
    }
    if ( $use_phylip_fitch_fm == 1 ) {
        &append( $CONSENSUS_PHYLIP_FM, $INTREE_FOR_PUZZLE );
    }
    if ( $use_phylip_fitch_me == 1 ) {
        &append( $CONSENSUS_PHYLIP_ME, $INTREE_FOR_PUZZLE );
    }
    if ( $use_raxml == 1 ) {
        # Needed, because TREE-PUZZLE adds internal labels for all subsequent trees
        # when evaluating given trees (this seems a strange behaviour).
        removeSupportValues( $CONSENSUS_RAXML, $CONSENSUS_RAXML."_support_removed" );
        &append( $CONSENSUS_RAXML."_support_removed", $INTREE_FOR_PUZZLE );
        &rm( $CONSENSUS_RAXML."_support_removed" );
    }
    if ( $use_proml == 1 ) {
        &append( $CONSENSUS_PROML, $INTREE_FOR_PUZZLE );
    }
    if ( $use_protpars == 1 ) {
        &append( $CONSENSUS_PROTPARS, $INTREE_FOR_PUZZLE );
    }
    if ( $all_count > 1 ) {
        &append( $CONSENSUS_ALL, $INTREE_FOR_PUZZLE );
    }
    

    # Puzzle for ML branch lenghts:
    # The alignment is read from infile by default.
    # The tree is read from intree by default.
    &rm( "infile" );
    &mv( "align", "infile" ); # align = original alignment in phylip interleaved.
    
    &executePuzzleToCalculateBranchLenghts( $matrix,
                                            $exact_parameter_est,
                                            $rate_heterogeneity );

    my $OUTTREE_PUZZLE = "outtree_puzzle";
 
    &rm( $outfile."_puzzle_outfile" ); 
   
    &mv( "outfile", $outfile."_puzzle_outfile" );
    &mv( "outtree", $OUTTREE_PUZZLE );
    &rm( "outdist" );
    &rm( "intree" );


    # Transfer
    # --------
    my $counter = 0;
    if ( $use_fastme == 1 ) {
        &executeSupportTransfer( $OUTTREE_PUZZLE, $CONSENSUS_FASTME, $fastme_outtree, $counter++ );
        &rm( $CONSENSUS_FASTME );
    }
    if ( $use_phylip_nj == 1 ) {
        &executeSupportTransfer( $OUTTREE_PUZZLE, $CONSENSUS_PHYLIP_NJ, $phylip_nj_outtree, $counter++ );
        &rm( $CONSENSUS_PHYLIP_NJ );
    }
    if ( $use_phylip_fitch_fm == 1 ) {
        &executeSupportTransfer( $OUTTREE_PUZZLE, $CONSENSUS_PHYLIP_FM, $phylip_fm_outtree, $counter++ );
        &rm( $CONSENSUS_PHYLIP_FM );
    }
    if ( $use_phylip_fitch_me == 1 ) {
        &executeSupportTransfer( $OUTTREE_PUZZLE, $CONSENSUS_PHYLIP_ME, $phylip_me_outtree, $counter++ );
        &rm( $CONSENSUS_PHYLIP_ME );
    }
    if ( $use_raxml == 1 ) {
        &to_phyloxml( $CONSENSUS_RAXML, $raxml_outtree, 1, 1 );
        $counter++;
    }
    if ( $use_proml == 1 ) {
        &executeSupportTransfer( $OUTTREE_PUZZLE, $CONSENSUS_PROML, $proml_outtree, $counter++ );
        &rm( $CONSENSUS_PROML );
    } 
    if ( $use_protpars == 1 ) {
        &executeSupportTransfer( $OUTTREE_PUZZLE, $CONSENSUS_PROTPARS, $protpars_outtree, $counter++ );
        &rm( $CONSENSUS_PROTPARS );
    } 
    if ( $all_count > 1 ) {
        &executeSupportTransfer( $OUTTREE_PUZZLE, $CONSENSUS_ALL, $all_outtree, $counter++ );
        &rm( $CONSENSUS_ALL );
    }
    
    # Clean up
    # --------
    &rm( $OUTTREE_PUZZLE );
    &rm( $SEQBOOT_OUTFILE );
    if ( $keep_multiple_trees == 1 ) {
        if ( $use_fastme == 1 ) {
            &mv( $OUTTREE_FASTME, $multitreefile_fastme );
        }
        if ( $use_phylip_nj == 1 ) {
            &mv( $OUTTREE_PHYLIP_NJ, $multitreefile_phylip_nj );
        }
        if ( $use_phylip_fitch_fm == 1 ) {
            &mv( $OUTTREE_PHYLIP_FM, $multitreefile_phylip_fm );
        }
        if ( $use_phylip_fitch_me == 1 ) {
            &mv( $OUTTREE_PHYLIP_ME, $multitreefile_phylip_me );
        }
        if ( $use_proml == 1 ) {
            &mv( $OUTTREE_PROML, $multitreefile_proml );
        }
        if ( $use_protpars == 1 ) {
            &mv( $OUTTREE_PROTPARS, $multitreefile_protpars );
        }
        &mv( $pwdfile, $multipwdfile );
    }
    else {
        if ( $use_fastme == 1 ) {
            &rm( $OUTTREE_FASTME );
        }
        if ( $use_phylip_nj == 1 ) {
            &rm( $OUTTREE_PHYLIP_NJ );
        }
        if ( $use_phylip_fitch_fm == 1 ) {
            &rm( $OUTTREE_PHYLIP_FM );
        }
        if ( $use_phylip_fitch_me == 1 ) {
            &rm( $OUTTREE_PHYLIP_ME );
        }
        if ( $use_proml == 1 ) {
            &rm( $OUTTREE_PROML );
        }
        if ( $use_protpars == 1 ) {
            &rm( $OUTTREE_PROTPARS );
        }
        &rm( $pwdfile );
    }
    if ( $all_count > 1 ) {
        &rm( $OUTTREES_ALL );
    }
} # if ( $bootstraps > 1 )
else {
    &rm( "infile.dist" );
    &rm( "infile.dist_" );
    &rm( "infile.puzzle" );
    if ( $use_fastme == 1 ) {
        &to_phyloxml( $OUTTREE_FASTME, $fastme_outtree, 0, 1 );
    }
    if ( $use_phylip_nj == 1 ) {
        &to_phyloxml( $OUTTREE_PHYLIP_NJ, $phylip_nj_outtree, 0, 1);
    }
    if ( $use_phylip_fitch_fm == 1 ) {
        &to_phyloxml( $OUTTREE_PHYLIP_FM, $phylip_fm_outtree, 0, 1 );
    }
    if ( $use_phylip_fitch_me == 1 ) {
        &to_phyloxml( $OUTTREE_PHYLIP_ME, $phylip_me_outtree, 0, 1 );
    }
    if ( $use_raxml == 1 ) {
         &to_phyloxml( "align.raxml.bestTree", $raxml_outtree, 1, 1 );
         &rm("align.raxml.bestTree");
    }
    if ( $use_proml == 1 ) {
        &to_phyloxml( $OUTTREE_PROML, $proml_outtree, 0, 1 );
    }
    if ( $use_protpars == 1 ) {
        &to_phyloxml( $OUTTREE_PROTPARS, $protpars_outtree, 0, 1 );
    }
} # if ( $bootstraps > 1 )

&rm( $infile );
&rm( "infile" );
&rm( "align" );
&rm( "align.reduced" );


$log = $log."Finish time/date                    : ".`date`;

if ( $bootstraps > 1 ) {
    $log = $log."Puzzle output file                  : ".$outfile."_puzzle_outfile\n";
}
$log = $log."Columns in alignment                : $number_of_aa\n";
$log = $log."Number of sequences in alignment    : $number_of_seqs\n";
if ( $all_count > 1 ) {
    $log = $log."Combined consensus                  : $all_outtree\n";
} 


if ( $bootstraps > 1 ) {
    $log = $log."\n\n";
    $log = $log."Simple support value statistics (trees are numbered the same as for TREE PUZZLE output)\n";
    $log = $log."------------------------------- \n";
    $log = $log."\n";
}    

open( OUT, ">$logfile" ) || die "\n$0: Cannot create file <<$logfile>>: $!\n";
print OUT $log;
close( OUT );

if ( $bootstraps > 1 ) {
    # Simple support statistics
    # -------------------------
    my $SS_OUT = $temp_dir."/ss_out";
    my @phylos = ();
    my $ounter = 0;
    if ( $use_fastme == 1 ) {
        $phylos[ $ounter++ ] = $fastme_outtree;
    }
    if ( $use_phylip_nj == 1 ) {
        $phylos[ $ounter++ ] = $phylip_nj_outtree;
    }
    if ( $use_phylip_fitch_fm == 1 ) {
        $phylos[ $ounter++ ] = $phylip_fm_outtree;
    }
    if ( $use_phylip_fitch_me == 1 ) {
        $phylos[ $ounter++ ] = $phylip_me_outtree;
    }
    if ( $use_raxml == 1 ) {
        $phylos[ $ounter++ ] = $raxml_outtree;
    }
    if ( $use_proml == 1 ) {
        $phylos[ $ounter++ ] = $proml_outtree;
    }
    if ( $use_protpars == 1 ) {
        $phylos[ $ounter++ ] = $protpars_outtree;
    }
    if ( $all_count > 1) {
        $phylos[ $ounter++ ] = $all_outtree;
    }
    &executeSupportStatistics( $SS_OUT, @phylos );
    &append( $SS_OUT, $logfile );
    &rm( $SS_OUT );
    
    # Append parts of puzzle output file
    # ----------------------------------
    if ( $all_count > 1 ) {
        &parsePuzzleOutfile( $outfile."_puzzle_outfile", $logfile );
    }
}

chdir( $current_dir ) 
|| die "\n\n$0: Could not chdir to <<$current_dir>>: $!\n\n";

rmdir( $temp_dir )
|| print "\n\n$0: Warning: Could not remove <<$temp_dir>>: $!\n\n";

print "\n\n\n$0 successfully completed.\n\n";

exit( 0 ); 
    


# Methods:
# --------

# Six arguments:
# 1. DNA or Amino-Acids sequence filename (PHYLIP format)
# 2. Model
# 3. Replicates (bootstrap)
# NOTE. RaxML does its own bootstrapping.
sub executeRaxmlNG {
    my $msa            = $_[ 0 ];
    my $model          = $_[ 1 ];
    my $replicates     = $_[ 2 ];

    &testForTextFilePresence( $msa );
    my $command = "$RAXMLNG --tree rand{10},pars{10} --msa-format PHYLIP --model $model --msa $msa";

    if ( $replicates > 1 ) {
        $command = $command . " --all --bs-trees $replicates";
    }

    print( "\n$command\n");

    system( $command )
    && &dieWithUnexpectedError( $command );

}


sub to_phyloxml {
    my $from = $_[ 0 ];
    my $to   = $_[ 1 ];
    my $internal_names_are_boots = $_[ 2 ];
    my $extract_taxonomy = $_[ 3 ];
    &dieIfFileExists( $to );
    &dieIfFileNotExists( $from );
    my $command = "$NEWICK_TO_PHYLOXML -f=nn $from $to";
    if ( $internal_names_are_boots == 1 ) {
        $command = $command . " -i";
    }
    if ( $extract_taxonomy  == 1 ) {
        $command = $command . " -xt";
    }
    system( $command  )
    && die "$0: Could not execute \"$command \"";
    &rm( $from );
}


sub mv {
    my $from = $_[ 0 ];
    my $to   = $_[ 1 ];
    &dieIfFileExists( $to );
    &dieIfFileNotExists( $from );
    system( "mv", $from, $to )
    && die "\n\n$0: could not move \"$from\" to \"$to\": $!\n\n";
}

sub cp {
    my $from = $_[ 0 ];
    my $to   = $_[ 1 ];
    &dieIfFileExists( $to );
    &dieIfFileNotExists( $from );
   
    system( "cp", $from, $to )
    && die "\n\n$0: could not copy \"$from\" to \"$to\": $!\n\n";
}

sub rm {
    my $f = $_[ 0 ];
    unlink( $f );
}

sub consense {
    my $multi_in     = $_[ 0 ]; 
    my $consense_out = $_[ 1 ];
    &executeConsense( $multi_in );
    &mv( "outtree", $consense_out );
    &rm( "outfile" );
    
}    



# 1. file to be appended
# 2. file to append to
sub append {
    my $to_be_appended = $_[ 0 ];
    my $append_to      = $_[ 1 ];
    &dieIfFileNotExists( $to_be_appended );
    system( "cat $to_be_appended >> $append_to" )
    && die "\n\n$0: could not execute \"cat $to_be_appended >> $append_to\": $!\n\n";
    
}

sub dieIfFileExists {
    my $file = $_[ 0 ]; 
    if ( -e $file ) {
        die "\n\n$0: \"$file\" already exists\n\n";
    }
} 

sub dieIfFileNotExists {
    my $file = $_[ 0 ]; 
    unless ( ( -s $file ) && ( -f $file ) ) {
       die( "\n\n$0: \"$file\" does not exist or is empty" );
    }
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

   
    $verb = "
2";
    
    system( "$SEQBOOT << !
r
$bs$verb
Y
$s
!" )
    && die "$0: Could not execute \"$SEQBOOT\"";
  
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
    && die "$0: Could not execute \"$PUZZLE\" (mat=$mat est=$est rate=$rate)";
    
}

# Two arguments:
# 1. puzzle outfile
# 2. file to append to
sub parsePuzzleOutfile {
    my $puzzle_outfile    = $_[ 0 ];
    my $file_to_append_to = $_[ 1 ];
    &testForTextFilePresence( $puzzle_outfile );
    open( OUT, ">>$file_to_append_to" ) || &dieWithUnexpectedError( "Cannot open \"$file_to_append_to\"" );
    open( IN, "$puzzle_outfile" ) || &dieWithUnexpectedError( "Cannot open file \"$puzzle_outfile\"" );
    my $return_line;
    my $read = 0;
    print OUT "\nTREE PUZZLE output\n";
    print OUT "------------------\n";
    while ( $return_line = <IN> ) {
        if ( $return_line =~/COMPARISON OF USER TREES/ ) {
            $read = 1;
        }                         
        elsif( $return_line =~/TIME STAMP/  ) {
            $read = 0;
        }
        elsif( $read ) {
            print OUT $return_line;
        }
    }
    close( IN );
    close( OUT );
}   

# Three/four arguments:
# 1. Name of file containing tree with correct branch lengths
# 2. Name of file containing tree with correct bootstraps
# 3. Outputfilename
# 4. Index of tree with correct branch lengths, in case more than one in file
# Last modified: 2007.11.27
sub executeSupportTransfer {
    my $tree_with_bl = $_[ 0 ];
    my $tree_with_bs = $_[ 1 ];
    my $out          = $_[ 2 ];
    my $index        = $_[ 3 ];
  
    &testForTextFilePresence( $tree_with_bl );
    &testForTextFilePresence( $tree_with_bs );
    my $command = "$SUPPORT_TRANSFER $tree_with_bl $tree_with_bs $out $index";
    system( $command )
    && die "$0: Could not execute \"$command\"";
}


# Two or more arguments:
# 1. outfile
# 2. phylogeny 1 with support values 
# 3. phylogeny 2 with support values 
# 4. ...
sub executeSupportStatistics {
    my $outfile      = $_[ 0 ];
    &dieIfFileExists( $outfile );
    my $phylos = "";
    for( my $i = 1; $i < scalar(@_); ++$i ) {
        &testForTextFilePresence( $_[ $i ] );
        $phylos .= $_[ $i ]." ";
    }    
    my $command = "$SUPPORT_STATISTICS -o=$outfile $phylos";
    system( "$command" )
    && die "$0: Could not execute \"$command\"";
}


sub getNumberOfSeqsAndAas { 
    my $infile = $_[ 0 ];
    my $seqs = 0;
    my $aa   = 0;
    open( IN, "$infile" ) || die "\n$0: Cannot open file <<$infile>>: $!\n";
    while( <IN> ) { 
        if ( $_ =~ /^\s*(\d+)\s+(\d+)\s*$/ ) { 
            $seqs = $1;
            $aa   = $2;
        } 
    }
    close( IN );
    
    if (  $seqs == 0 ||  $aa  == 0 ) {
        die( "\n$0: Could not get number of seqs and aa from: $infile" );
    }
    return $seqs, $aa;
}



sub removeSupportValues {
    my $infile  = $_[ 0 ];
    my $outfile = $_[ 1 ];
    &testForTextFilePresence( $infile );
    open( OUT, ">$outfile" ) || &dieWithUnexpectedError( "Cannot create file \"$outfile\"" );
    open( IN, "$infile" ) || &dieWithUnexpectedError( "Cannot open file \"$infile\"" );
    while ( my $line = <IN> ) {
        $line =~ s/\)\d+\.?\d*:/\):/g;
        print OUT "$line";
    }
    close( OUT );
    close( IN );   
}




# Six arguments:
# 1. name of alignment file (in correct format!)
# 2. number of bootstraps
# 3. jumbles: 0: do not jumble; >=1 number of jumbles
# 4. seed for random number generator
# 5. 1 for PAM instead of JTT
# 6. 1 to use globale rearragements
sub executeProml {
    my $align            = $_[ 0 ];
    my $bs               = $_[ 1 ];
    my $rand             = $_[ 2 ];
    my $s                = $_[ 3 ];
    my $use_pam          = $_[ 4 ];
    my $use_global_rearr = $_[ 5 ];
    my $jumble = "";
    my $multi  = "";
    my $pam    = ""; 
   
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
   
    
   if ( $use_pam == 1 ) {
        $pam = "
P
P";
    }
    
    my $global = "";
    if ( $use_global_rearr == 1 ) { 
        $global = "
G"; 
    }

    system( "$PROML  2>&1 << !
$align$jumble$multi$pam$global
3
Y
!" )
    && &dieWithUnexpectedError( "Could not execute \"$PROML $align$jumble$multi$pam$global\"" );
    # 3: Do NOT print out tree
      
    return;

} ## executeProml


sub printUsage {

    print <<END;
Copyright (C) 2017 Christian M Zmasek
All rights reserved

Author: Christian M Zmasek
cmzmasek at yahoo dot com
https://sites.google.com/site/cmzmasek/home/software/forester

  Requirements  phylo_pl is part of the FORESTER collection of programs.
  ------------  Many of its global variables are set via forester.pm.

  Note. Use xt.pl (for Pfam alignments) or mt.pl (for other alignments) 
  to run phylo_pl.pl on whole directories of alignments files. 

  Usage
  -----

      phylo_pl.pl [-options] <input alignment in SELEX (Pfam), PHYLIP
      sequential format, or Clustal W output> <outputfile>
      [path/name for temporary directory to be created]

     Example:
     "% phylo_pl.pl -B100q\@1nbS9X IL5.aln IL5_tree"

  Options
  -------
 
  Bx : Number of bootstraps. B0: do not bootstrap. Default is 100 bootstrapps.
       The number of bootstrapps should be divisible by 10.
  J  : Use JTT model (Jones et al. 1992) in TREE-PUZZLE and/or PHYML, RAXML, default: VT (Mueller-Vingron 2000).
  O  : Use BLOSUM 62 model (Henikoff-Henikoff 92) in TREE-PUZZLE and/or PHYML, RAXML, default: VT.
  M  : Use mtREV24 model (Adachi-Hasegawa 1996) in TREE-PUZZLE and/or PHYML, default: VT.
  W  : Use WAG model (Whelan-Goldman 2000) in TREE-PUZZLE and/or PHYML, RAXML, default: VT.
  P  : Use PAM model (Dayhoff et al. 1978) in TREE-PUZZLE and/or PHYML, RAXML, default: VT.
  L  : Use LG model (Le and Gascuel, 2008) in PHYML, RAXML; WAG in TREE-PUZZLE.
  D  : Use DCMut model (Kosial and Goldman, 2005) in PHYML, RAXML; VT in TREE-PUZZLE.
  A  : Let TREE-PUZZLE choose which model to use, default: VT.
  H  : Use HKY (Hasegawa et al. 1985) in TREE-PUZZLE [for nucleic acids]
  T  : Use TN (Tamura-Nei 1993) in TREE-PUZZLE [for nucleic acids]
  Z  : Use GTR (e.g. Lanave et al. 1980) in TREE-PUZZLE [for nucleic acids]
  C  : Use SH (Schoeniger-von Haeseler 1994) in TREE-PUZZLE [for nucleic acids]
  E  : Exact parameter estimates in TREE-PUZZLE, default: Approximate.
       Model of rate heterogeneity in TREE-PUZZLE (default: Uniform rate):
  g  : 8 Gamma distributed rates
  t  : Two rates (1 invariable + 1 variable)
  m  : Mixed (1 invariable + 8 Gamma rates)
  q  : Use FastME
  n  : Use PHYLIP Neighbor (NJ).                    
  f  : Use PHYLIP Fitch.
  e  : Use PHYLIP Minimal Evolution.
  x  : Use RAxML.
  o  : Use PHYLIP proml. 
  p  : Use PHYLIP protpars.
  rx : Number of relative substitution rate categories in PHYML (default is 4).
  jx : Number of jumbles (input order randomization) for PHYLIP FM, ME, PROTPARS, and PROML (default is 2) (random seed set with Sx).
  I  : Estimate proportion of invariable sites in RAXML and/or PHYML (otherwise, proportion "0.0" is used in PHYML)
  G  : to turn on global rearrangements in PHYLIP FM, ME, and PROML
  Sx : Seed for random number generator(s). Must be 4n+1. Default is 9.
  X  : To keep multiple tree file (=trees from bootstrap resampled alignments) and 
       pairwise distance matrix file (in case of bootstrap analysis).
  
END
  
} ## printUsage

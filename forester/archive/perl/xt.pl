#!/usr/bin/perl -W

# xt.pl
# -----
#
# Copyright (C) 2003 Christian M. Zmasek
# All rights reserved
#
# Author: Christian M. Zmasek
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Version: 1.010
# Last modified 03/25/03
#
#
#
# Calculates trees based on Pfam alignments or precalculated distances using
# makeTree.pl.


use strict;
use FindBin;
use lib $FindBin::Bin;
use rio_module;


# To use _your_ species list make $MY_SPECIES_NAMES_FILE point to it.
# To use _your_ TrEMBL ACDEOS make $MY_TREMBL_ACDEOS_FILE point to it.

my $MY_SPECIES_NAMES_FILE = $SPECIES_NAMES_FILE; # $SPECIES_NAMES_FILE is inherited 
                                                 # from rio_module.pm

my $MY_TREMBL_ACDEOS_FILE = $TREMBL_ACDEOS_FILE; # $TREMBL_ACDEOS_FILE is inherited 
                                                 # from rio_module.pm
                                                 
my $MY_TEMP_DIR           = $TEMP_DIR_DEFAULT;   # $TEMP_DIR_DEFAULT is inherited 
                                                 # from rio_module.pm 

my $LOGFILE               = "00_xt_logfile";
my $PWD_SUFFIX            = ".pwd";
my $ALN_SUFFIX            = ".aln";

my $use_precalc_pwd       = 0;  # 0: input is Pfam aligments ($input_dir must point to "/Pfam/Full/").
                                # 1: input is precalculated pairwise distancs ($input_dir must point to "").
my $use_precalc_pwd_and_aln = 0;# 0: otherwise
                                # 1: input is precalculated pairwise distancs
                                # _and_ alns,$use_precalc_pwd = 1 ($input_dir must point to alns).
my $add_species           = 0;  # "I": 0: do nothing with species information.
                                # "S": 1: add species code to TrEMBL sequences and ignore sequences from 
                                #         species not in $MY_SPECIES_NAMES_FILE (only if input is Pfam aligments).
my $options               = ""; # Options for makeTree.pl, see makeTree.pl.
                                # Do not use F [Pairwise distance (pwd) file as input (instead of alignment)]
                                # since this is determined with $USE_PRECALC_PWD
my $min_seqs              = 0;  # Minimal number of sequences (TREE-PUZZLE needs at least four seqs).
                                # Ignored if $USE_PRECALC_PWD = 1
my $max_seqs              = 0;  # Maximal number of sequences.
                                # Ignored if $USE_PRECALC_PWD = 1
my $input_dir             = "";
my $input_dir_aln         = ""; # for .aln files
my $output_dir            = "";

my $i                     = 0;
my $seqs                  = 0;
my $filename              = "";
my @filenames             = ();
my %AC_OS                 = (); # AC -> species name
my %Species_names_hash    = ();
my $too_small             = 0;
my $too_large             = 0;
my $already_present       = 0;
my @too_small_names       = ();
my @too_large_names       = ();
my @already_present_names = ();


# Analyzes the options:
# ---------------------

unless ( @ARGV == 3 || @ARGV == 4 || @ARGV == 6 ) {
    &printUsage();
}

if ( @ARGV == 3 ) {
    $use_precalc_pwd         = 1;
    $use_precalc_pwd_and_aln = 0;
    $options         = $ARGV[ 0 ]; 
    $input_dir       = $ARGV[ 1 ];
    $output_dir      = $ARGV[ 2 ];
    $add_species     = 0; 
}
elsif ( @ARGV == 4 ) {
    $use_precalc_pwd         = 1;
    $use_precalc_pwd_and_aln = 1;
    $options         = $ARGV[ 0 ]; 
    $input_dir       = $ARGV[ 1 ];
    $input_dir_aln   = $ARGV[ 2 ];
    $output_dir      = $ARGV[ 3 ];
    $add_species     = 0;
    $input_dir_aln = &addSlashAtEndIfNotPresent( $input_dir_aln ); 
}
else {
    $use_precalc_pwd         = 0;
    $use_precalc_pwd_and_aln = 0;
    $add_species     = $ARGV[ 0 ];
    $options         = $ARGV[ 1 ]; 
    $min_seqs        = $ARGV[ 2 ];
    $max_seqs        = $ARGV[ 3 ];
    $input_dir       = $ARGV[ 4 ];
    $output_dir      = $ARGV[ 5 ];
    if ( $min_seqs < 4 ) {
        $min_seqs = 4;
    }
    if ( $add_species eq "I" ) {
        $add_species = 0;
    }
    elsif ( $add_species eq "S" ) {
        $add_species = 1;
    }
    else {
        print( "\nFirst must be either \"I\" [Ignore species] or\n\"S\" [add Species code to TrEMBL sequences and ignore sequences from species not in $MY_SPECIES_NAMES_FILE].\n\n" );
        &printUsage(); 
    }    
}



$input_dir   = &addSlashAtEndIfNotPresent( $input_dir ); 
$output_dir  = &addSlashAtEndIfNotPresent( $output_dir );
$MY_TEMP_DIR = &addSlashAtEndIfNotPresent( $MY_TEMP_DIR );




# This adds a "-" before the options for makeTree:
# ------------------------------------------------
unless ( $options =~ /^-/ ) {
    $options = "-".$options;
}



# If based on pwd, species are "fixed" and certain options for makeTree 
# are not applicable and option "F" is mandatory:
# ---------------------------------------------------------------------
if ( $use_precalc_pwd == 1 ) {
    $options =~ s/D//g;
    $options =~ s/C//g;
    $options =~ s/N//g;
    unless ( $options =~ /F/ ) { 
        $options = $options."F";
    }
}
else {
    $options =~ s/F//g;
}

if ( $use_precalc_pwd_and_aln == 1 ) {
    unless ( $options =~ /U/ ) { 
        $options = $options."U";
    }
}
if ( $use_precalc_pwd_and_aln == 0 && $use_precalc_pwd == 1 ) {
    $options =~ s/U//g;
}    




# If species are to be considered, speices names file and TrEMBL ACDEOS
# files need to be read in:
# ---------------------------------------------------------------------
if ( $add_species == 1 ) {
    print "\nXT.PL: Reading species names file...\n";  
    &readSpeciesNamesFile( $MY_SPECIES_NAMES_FILE );
    print "\nXT.PL: Reading TrEMBL ACDEOS file...\n";
    &readTrEMBL_ACDEOS_FILE( $MY_TREMBL_ACDEOS_FILE );
}



# This creates the temp file:
# --------------------------

my $time = time;
my $ii   = 0;

my $temp_file = $MY_TEMP_DIR."xt".$time.$ii;

while ( -e $temp_file ) {
    $ii++;
    $temp_file = $MY_TEMP_DIR."xt".$time.$ii;
}



&startLogfile();

opendir( DIR, $input_dir ) || error( "Cannot open directory \"$input_dir\": $!" );

$i = 0;

while( defined( $filename = readdir( DIR ) ) ) {
    if ( $filename =~ /^\.\.?$/ ) {
        next;
    }
    if ( $use_precalc_pwd == 1 && $filename !~ /$PWD_SUFFIX$/ ) {
        next
    }
    $filenames[ $i ] = $filename;
    $i++;
}

close( DIR );

$i = 0;

FOREACH: foreach $filename ( @filenames ) {

    # If the corresponding tree seems to already exists, do next one.
    my $fn = $filename; 
    if ( $use_precalc_pwd == 1 ) {
        $fn =~ s/$PWD_SUFFIX$//;
    }
    if ( -e "$output_dir$fn.nhx" ) {
        $already_present_names[ $already_present++ ] = $fn;
        next FOREACH;
    }

    if ( $use_precalc_pwd != 1 ) {
    
        if ( $add_species == 1 ) {
            
            # 1. Pfam flat file name
            # 2. outfile name
            # Returns the number of sequences in the resulting alignment.
            $seqs = &removeSeqsFromPfamAlign( $input_dir.$filename, $temp_file );
        
        }
        else {   
            # Gets the number of seqs in the alignment.
            open( F, "$input_dir"."$filename" );
            while( <F> ) {
                if ( $_ =~/^#.+SQ\s+(\d+)\s*$/ ) {
                    $seqs = $1;
                    last;
                }
            }
            close( F );
        }

        if ( $seqs < $min_seqs ) {
            $too_small_names[ $too_small++ ] = $filename;
            next FOREACH;
        }
        if ( $seqs > $max_seqs ) {
            $too_large_names [ $too_large++ ] = $filename;
            next FOREACH;
        }
    }

    print "\n\n\n\n";
    print "XT.PL\n";
    if ( $use_precalc_pwd == 1 ) {
        print "working on: $filename\n";
    }
    else {
        print "working on: $filename [$seqs seqs]\n";
    }
    print "[tree calculation $i]\n";
    print "=====================================================================\n\n\n";


    unlink( "$output_dir$filename.aln", "$output_dir$filename.log" );

    print( "XT.PL: executing:\n" );
    
    my $inputfile = "";
    
    if ( $add_species == 1 ) {
        $inputfile = $temp_file;
    }
    else {
        $inputfile = $input_dir.$filename;
    }
    
    if ( $use_precalc_pwd == 1 ) {
        $filename =~ s/$PWD_SUFFIX$//;
    }
    
    if ( $use_precalc_pwd_and_aln == 1 ) {
        $inputfile = $inputfile." ".$input_dir_aln.$filename.$ALN_SUFFIX;
    }
   
    my $command = "$MAKETREE $options $inputfile $output_dir$filename.nhx";
    
    print( "$command\n" );
    system( $command ) && &error( "Could not execute \"$command\"" );
   
    if ( $add_species == 1 ) {
        if ( unlink( $temp_file ) != 1 ) {
            &error( "Unexpected: Could not delete \"$temp_file\"" );
        }
    }

    $i++;

}

&finishLogfile();

print( "\n\n\nXT.PL: Done!\n" );
print( "Wrote \"$LOGFILE\".\n\n" );

exit( 0 );






sub error{

    my $text = $_[ 0 ];

    print( "\nxt.pl: ERROR:\n" );
    print( "$text\n\n" );

    exit( -1 );

} ## dieWithUnexpectedError

# Similar to the method with the same name in "rio.pl".
# Removes sequences from a Pfam flat file.
# Adds species to TrEMBL seqs.
# It can remove all sequences not from species listed in a species names file.
# Two arguments:
# 1. Pfam flat file name
# 2. outfile name
# Returns the number of sequences in the resulting alignment.
# Last modified: 02/22/03
sub removeSeqsFromPfamAlign {
    my $infile                   = $_[ 0 ];
    my $outfile                  = $_[ 1 ];
    my $return_line              = "";
    my $saw_sequence_line        = 0;
    my $number_of_seqs           = 0;
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



# Last modified: 05/18/01
sub readTrEMBL_ACDEOS_FILE {
    my $infile      = $_[ 0 ];
    my $return_line = "";
     
    unless ( ( -s $infile  ) && ( -f $infile  ) && ( -T $infile ) ) {
        &error( "\"$infile\" does not exist, is empty, or is not a plain textfile" );
    }
    # Fill up (huge) hashs.
    open( HH, "$infile" ) || &error(  "Unexpected error: Cannot open file \"$infile\"" );
    while ( $return_line = <HH> ) {
       
        if ( $return_line =~ /(\S+);[^;]*;(\S+)/ ) {
             $AC_OS{ $1 } = $2;
        }
    }
    close( HH ); 
} ## readTrEMBL_ACDEOS_FILE



# Last modified: 05/17/01
sub startLogfile {
    if ( -e "$LOGFILE" ) {
        &error( "logfile \"$LOGFILE\" already exists, rename it or place it in another directory" );
    } 

    open( L, ">$LOGFILE" ) || &error( "Cannot create logfile: $!" );
    if ( $use_precalc_pwd != 1 ) {
        print L "Trees are based directly on Pfam alignments\n";
        if ( $add_species == 1 ) {
            print L "Add species code to TrEMBL sequences and ignore sequences\nfrom species not in $MY_SPECIES_NAMES_FILE\n";
        }
        else {
            print L "Do nothing with species information\n";
        }
    }
    else {
        print L "Trees are based on precalculated pairwise distances\n";
    }
    if ( $use_precalc_pwd_and_aln == 1 ) {
        print L "and the matching alignments\n";
    }
    print L "Options for makeTree: $options\n";
    if ( $use_precalc_pwd != 1 ) {
        print L "Min seqs            : $min_seqs\n";
        print L "Max seqs            : $max_seqs\n";
    }
    if ( $add_species == 1 ) {
        print L "TrEMBL ACDEOS file  : $MY_TREMBL_ACDEOS_FILE\n";
        print L "Species names file  : $MY_SPECIES_NAMES_FILE\n";
    }
    print L "Input directory     : $input_dir\n";
    if ( $use_precalc_pwd_and_aln == 1 ) {
    print L "Input directory aln : $input_dir_aln\n";
    }
    print L "Output directory    : $output_dir\n";
    print L "Start date          : ".`date`;
    
} ## startLogfile



# Last modified: 05/17/01
sub finishLogfile {
    my $j = 0;
    print L "\n\n";
    print L "Successfully calculated $i trees.\n";
    if ( $use_precalc_pwd != 1 ) {
        print L "Too large alignments (>$max_seqs): $too_large\n";
        print L "Too small alignments (<$min_seqs): $too_small\n";
    }
    print L "Alignments for which a tree appears to already exist: $already_present\n";
    print L "Finish date        : ".`date`."\n\n";
    if ( $use_precalc_pwd != 1 ) {
        print L "List of the $too_large alignments which were ignored because they\n";
        print L "contained too many sequences (>$max_seqs) [after pruning]:\n";
        for ( $j = 0; $j < $too_large; ++$j ) { 
            print L "$too_large_names[ $j ]\n";
        } 
        print L "\n\n";
        print L "List of the $too_small alignments which were ignored because they\n";
        print L "contained not enough sequences (<$min_seqs) [after pruning]:\n";
        for ( $j = 0; $j < $too_small; ++$j ) { 
            print L "$too_small_names[ $j ]\n";
        }
    }
    print L "\n\n";
    print L "List of the $already_present alignments which were ignored because\n";
    print L "a tree appears to already exist:\n";
    for ( $j = 0; $j < $already_present; ++$j ) { 
        print L "$already_present_names[ $j ]\n";
    }  
    print L "\n";
    close( L );
} ## finishLogfile


sub printUsage {
    print "\n";
    print " xt.pl\n";
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
    print " Tree construction using makeTree.pl based on directories\n"; 
    print " of Pfam alignments or precalculated pairwise distances.\n"; 
    print "\n";
    print "\n";
    print " Usage\n";
    print " -----\n";
    print "\n"; 
    print " Input is Pfam aligments:\n"; 
    print "   xt.pl <I or S> <options for makeTree.pl> <minimal number of seqs>\n";
    print "   <maximal number of seqs> <input directory: Pfam aligments> <output directory>\n";
    print "\n";
    print " Input is precalculated pairwise distancs:\n"; 
    print "   xt.pl <options for makeTree.pl> <input directory: pairwise distances \"$PWD_SUFFIX\"> <output directory>\n";
    print "\n";
    print " Input is precalculated pairwise distancs and corresponding alignment files:\n"; 
    print "   xt.pl <options for makeTree.pl> <input directory: pairwise distances \"$PWD_SUFFIX\">\n";
    print "   <input directory: corresponding alignment files \"$ALN_SUFFIX\"> <output directory>\n";
    print "\n";
    print "\n";
    print " Examples\n";
    print " --------\n";
    print "\n";
    print "   \"xt.pl S NS21UTRB100DX 4 200 DB/PFAM/Full/ trees/\"\n";
    print "\n";
    print "   \"xt.pl FLB100R /pfam2pwd_out/ trees/\"\n";
    print "\n";
    print "   \"xt.pl FULB100R /pfam2pwd_out/ /pfam2pwd_out/ trees/\"\n";
    print "\n";
    print "\n";
    print " Options\n";
    print " -------\n";
    print "\n";
    print " I: ignore species information (use all sequences)\n";
    print " S: add species codes to TrEMBL sequences and ignore sequences\n";
    print "    from species not in $MY_SPECIES_NAMES_FILE,\n";
    print "    species codes are extracted from $MY_TREMBL_ACDEOS_FILE\n";
    print "\n";
    print "\n";
    print " Options for makeTree\n";
    print " --------------------\n";
    print "\n";
    print " N  : Suggestion to remove columns in the alignment which contain gaps.\n"; 
    print "      Gaps are not removed, if, after removal of gaps, the resulting\n";
    print "      alignment would be shorter than $MIN_NUMBER_OF_AA aa (\$MIN_NUMBER_OF_AA).\n";
    print "      Default is not to remove gaps.\n";
    print " Bx : Number of bootstrapps. B0: do not bootstrap. Default: 100 bootstrapps.\n";
    print "      The number of bootstrapps should be divisible by 10.\n";
    print " U  : Use TREE-PUZZLE to calculate ML branchlengths for consesus tree, in case of\n";
    print "      bootstrapped analysis.\n";
    print " J  : Use JTT matrix (Jones et al. 1992) in TREE-PUZZLE, default: PAM.\n";
    print " L  : Use BLOSUM 62 matrix (Henikoff-Henikoff 92) in TREE-PUZZLE, default: PAM.\n";
    print " M  : Use mtREV24 matrix (Adachi-Hasegawa 1996) in TREE-PUZZLE, default: PAM.\n";
    print " W  : Use WAG matrix (Whelan-Goldman 2000) in TREE-PUZZLE, default: PAM.\n";
    print " T  : Use VT matrix (Mueller-Vingron 2000) in TREE-PUZZLE, default: PAM.\n";
    print " P  : Let TREE-PUZZLE choose which matrix to use, default: PAM.\n";
    print " R  : Randomize input order in PHYLIP NEIGHBOR.\n";
    print " Sx : Seed for random number generator(s). Must be 4n+1. Default is 9.\n";
    print " X  : To keep multiple tree file (=trees from bootstrap resampled alignments).\n";
    print " D  : To keep (and create, in case of bootstrap analysis) pairwise distance\n";
    print "      matrix file. This is created form the not resampled aligment.\n";
    print " C  : Calculate pairwise distances only (no tree). Bootstrap is always 1.\n";
    print "      No other files are generated.\n";
    print " F  : Pairwise distance (pwd) file as input (instead of alignment).\n";
    print "      No -D, -C, and -N options available in this case.\n";
    print " V  : Verbose\n";
    print "\n";
    exit( -1 );

}

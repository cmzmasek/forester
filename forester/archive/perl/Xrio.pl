#!/usr/bin/perl -w
#
# Xrio.pl
# -------
# Copyright (C) 1999-2001 Washington University School of Medicine
# and Howard Hughes Medical Institute
# All rights reserved
#
# Author: Christian M. Zmasek 
#         zmasek@genetics.wustl.edu
#         http://www.genetics.wustl.edu/eddy/people/zmasek/
#
# Created: 03/01/01
#
# Last modified 06/22/01


# Objective. Runs "rio.pl" for each Pfam assignment in "infile".
#
# Usage. rio.pl <infile> <species names file> <output directory> <outfile> <logfile> 
#
# species names file: list of species to use for analysis.
#
# This version uses the CE number as identifier for output files.\n";
#
# Format for infile: 
#  
#  >>3R5.2 CE19648    (CAMBRIDGE) TR:Q9XWB1 protein_id:CAA21778.1
#  //
#  
#  >>4R79.1 CE19649   Zinc-binding metalloprotease domain (CAMBRIDGE) protein_id:CAB63429.1
#  =Astacin  Astacin (Peptidase family M12A)                296.3    3.8e-85   1
#  //
#  
#  >>4R79.2 CE19650   Ras family (CAMBRIDGE) TR:Q9XXA4 protein_id:CAA20282.1
#  =ras           Ras family                                208.8    8.1e-59   1
#  =FA_desaturase Fatty acid desaturase                       4.5        1.5   1
#  =UPF0117       Domain of unknown function DUF36            3.1        3.5   1
#  =arf           ADP-ribosylation factor family            -46.0    1.5e-05   1
#  //                                                                           
# 
#
#

# Xrio.pl /nfs/wol2/people/zmasek/wormpep43_hmmpfam6.2/wormpep43_Hmmpfam_6.2 /nfs/wol2/people/zmasek/species_trees/tree_of_life_bin_1-4_species_list /nfs/wol2/people/zmasek/XrioTEST3 /nfs/wol2/people/zmasek/XrioTEST3/OUTFILE1 /nfs/wol2/people/zmasek/XrioTEST3/LOG1





use strict;

use FindBin;
use lib $FindBin::Bin;
use rio_module;

   $RIO_PL                  = "rio.pl";
my $VERSION                 = "3.000";

my $FASTA_DB                = "/nfs/wol2/people/zmasek/DB/wormpep/wormpep43";
my $QUERY_SPECIES           = "CAEEL";
my $SPECIES_TREE            = $SPECIES_TREE_FILE_DEFAULT;

my $RIOPL_OPTIONS           = "T=B P=6 L=0 R=0 U=80 V=0 X=2 Y=2 Z=2 C E I";

my $TEMP_DIR                = "/tmp/Xriopl";  # Where all the temp files, etc will be created.

my %Species_names_hash      = ();

my $infile                  = "";
my $outfile                 = ""; # Huge file of all rio outputs.
my $logfile                 = ""; # Lists all sequences which have been analyzed successfully.
my $output_directory        = "";
my $species_names_file      = "";


my $return_line             = "";
my $ID                      = "";
my $pfam_name               = "";
my $E_value                 = 0;
my $score                   = 0;
my $GA                      = 0;
my $temp_dir                = "";
my $outname                 = "";
my %outnames                = ();
my $seqs                    = 0;
my $ii                      = 0;
my $time                    = 0;
my $successful              = 0;
my $query_not_aligned       = 0;
my $pwd_not_present         = 0;
my $already_done            = 0;
my $start_date              = "";
my $old_fh                  = "";
my %AC_OS                   = (); # AC -> species name for TrEMBL seqs
my %AC_DE                   = (); # AC -> description for TrEMBL seqs

my $description_line        = "";
my $message1                = "";
my $message2                = "";

$start_date = `date`;

if ( @ARGV != 5 ) {
    &errorInCommandLine();
    exit ( -1 ); 
}

$infile             = $ARGV[ 0 ];
$species_names_file = $ARGV[ 1 ];
$output_directory   = $ARGV[ 2 ];
$outfile            = $ARGV[ 3 ];
$logfile            = $ARGV[ 4 ];


if ( -e $outfile ) {
    die "\n\n$0: <<$outfile>> already exists.\n\n";
}
unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
    die "\n\n$0: <<$infile>> does not exist, is empty, or is not a plain textfile.\n\n";
}
unless ( ( -s $species_names_file ) && ( -f $species_names_file ) && ( -T $species_names_file ) ) {
    die "\n\n$0: <<$species_names_file>> does not exist, is empty, or is not a plain textfile.\n\n";
}
unless ( ( -s $TREMBL_ACDEOS_FILE  ) && ( -f $TREMBL_ACDEOS_FILE  ) && ( -T $TREMBL_ACDEOS_FILE ) ) {
    die "\n\n$0: <<$TREMBL_ACDEOS_FILE>> does not exist, is empty, or is not a plain textfile.\n\n";
}
unless ( ( -e $output_directory ) && ( -d $output_directory ) ) {
    die "\n\n$0: <<$output_directory>> does not exist, or is not a directory.\n\n";
}



# Reads in the species file:
# -------------------------- 
&readSpeciesNamesFile( $species_names_file );



# Reads in the file containing AC, DE and OS for TrEMBL seqs:
# ----------------------------------------------------------- 
open( HH, "$TREMBL_ACDEOS_FILE" ) || die "\n\n$0: Unexpected error: Cannot open file <<$TREMBL_ACDEOS_FILE>>: $!\n\n";
while ( $return_line = <HH> ) {
    if ( $return_line =~ /(\S+);([^;]*);(\S+)/ ) {
        $AC_OS{ $1 } = $3;
        $AC_DE{ $1 } = $2;
    }
}
close( HH ); 
    


# Reads in outnames in logfile, if present:
# -----------------------------------------
if ( ( -s $logfile  ) ) {
    open( L, "$logfile" ) || die "\n\n$0: Unexpected error: Cannot open file <<$logfile>>: $!\n\n";
    while ( $return_line = <L> ) {
        if ( $return_line =~ /\s*(\S+)/ ) {
            $outnames{ $1 } = 0;
        }
    }
    close( L ); 
}



# Creates the temp directory:
# ---------------------------

$ii = 0;

$time = time;

$temp_dir = $TEMP_DIR.$time.$ii;

while ( -e $temp_dir ) {
    $ii++;
    $temp_dir = $TEMP_DIR.$time.$ii;
}

mkdir(  $temp_dir, 0777 )
|| die "\n\n$0:Unexpected error: Could not create <<$temp_dir>>: $!\n\n";

unless ( ( -e $temp_dir ) && ( -d $temp_dir ) ) {
    die "\n\n$0:Unexpected error: <<$temp_dir>> does not exist, or is not a directory: $!\n\n";
}



$message1 = "# $0\n".
            "# Version                : $VERSION\n".
            "# Date started           : $start_date".
            "# Infile                 : $infile\n".
            "# Species names file     : $species_names_file\n".
            "# Output directory       : $output_directory\n".
            "# Outfile                : $outfile\n".
            "# RIO PWD directory      : $RIO_PWD_DIRECTORY\n".
            "# RIO BSP directory      : $RIO_BSP_DIRECTORY\n".
            "# RIO NBD directory      : $RIO_NBD_DIRECTORY\n".
            "# RIO ALN directory      : $RIO_ALN_DIRECTORY\n".
            "# RIO HMM directory      : $RIO_HMM_DIRECTORY\n".
            "# Fasta db               : $FASTA_DB\n".    
            "# Species of query       : $QUERY_SPECIES\n".
            "# Species tree           : $SPECIES_TREE\n".     
            "# rio.pl options         : $RIOPL_OPTIONS\n\n\n"; 

open( IN, "$infile"  ) || die "\n\n$0: Cannot open file <<$infile>>: $!\n\n";
open( LOG, ">> $logfile" ) || die "\n\n$0: Cannot open file <<$logfile>>: $!\n\n";


# Turns off buffering for LOG.
$old_fh = select( LOG );
$| = 1;
select( $old_fh );


$ID = "";

W: while ( $return_line = <IN> ) {

    if ( $return_line =~ /^\s*>>.*(CE\d+)/ ) {
        $ID = $1;
        $return_line =~ /^\s*>>(.+)/;
        $description_line = $1;
    }
    elsif ( $return_line =~ /^\s*\/\// ) {
        $ID = "";
    }
    elsif ( $return_line =~ /^\s*=(\S+)\s+.+\s+(\S+)\s+(\S+)\s+\S+\s*$/
    && $ID ne "" ) {
 
        $pfam_name = $1;
        $score     = $2;          
        $E_value   = $3;

        $outname = $ID.".".$pfam_name;                          
                    
        # Checks if already done.
        if ( %outnames && exists( $outnames{ $outname } ) ) {  
            $already_done++;
            next W; 
        }
        
        &executeHmmfetch( $PFAM_HMM_DB, $pfam_name, $temp_dir."/HMMFILE" );
        
        $GA = &getGA1cutoff( $temp_dir."/HMMFILE" );
        unlink( $temp_dir."/HMMFILE" );       
 
        if ( $GA == 2000 ) {
            die "\n\n$0: Unexpected error: No GA cutoff found for \"$pfam_name\".\n\n";
        }
        elsif ( $score < $GA ) {
            next W;
        }
         
        if ( -s $output_directory."/".$outname ) {  
            unlink( $output_directory."/".$outname );
        }


        $message1 .= "\n\n".
                     "# ############################################################################\n".
                     "# Annotation: $description_line\n".
                     "# HMM       : $pfam_name\n".
                     "# score     : $score\n".
                     "# E-value   : $E_value\n";



        unless ( ( -s $RIO_PWD_DIRECTORY.$pfam_name.$SUFFIX_PWD ) ) {
            $pwd_not_present++;
            $message1 .= "# No PWD file for this family.\n".
                         "# ############################################################################\n";
            next W;
        }

       
        unless ( ( -s $PFAM_SEED_DIRECTORY."/".$pfam_name ) && ( -f $PFAM_SEED_DIRECTORY."/".$pfam_name ) && ( -T $PFAM_SEED_DIRECTORY."/".$pfam_name ) ) {
            die "\n\n$0: Error: Pfam seed alignment <<$PFAM_SEED_DIRECTORY"."/"."$pfam_name>> not present.\n\n";
        }
      

        &getSequenceFromFastaFile( $FASTA_DB,
                                   $temp_dir."/QUERY",
                                   $ID );

        &performRIO( $pfam_name,                               # A= 
                     $temp_dir."/QUERY",                       # Q=
                     $output_directory."/".$outname,           # O=
                     $ID."_".$QUERY_SPECIES,                   # N=
                     $SPECIES_TREE,                            # S= 
                     $RIOPL_OPTIONS,                           # L=0 R=0 U=70 V=0 X=2 Y=2 Z=2 C E K I x
                     $temp_dir."/riopltempdir" );              # j=



        if ( -s $output_directory."/".$outname ) {  
            $successful++;                  
        }
        else {
            $message1 .= "# Query has not been aligned (E value too low).\n".
                         "# ############################################################################\n";
            $query_not_aligned++;
        }

        if ( unlink( $temp_dir."/QUERY" ) != 1 ) {
            die "\n$0: Unexpected error: File(s) could not be deleted.\n";
        }



        if ( -s $output_directory."/".$outname ) {
            open( OUT_MESG_ONE, ">$temp_dir/_message1_" ) || die "\n\n$0: Cannot create file \"$temp_dir/_message1_\": $!\n\n";
            print OUT_MESG_ONE ( $message1 );
            close( OUT_MESG_ONE );

            $message1 = "";

            open( OUT_MESG_TWO, ">$temp_dir/_message2_" ) || die "\n\n$0: Cannot create file \"$temp_dir/_message2_\": $!\n\n";
            print OUT_MESG_TWO ( "# Successful calculations                  : $successful\n" );
            print OUT_MESG_TWO ( "# No calculation due to absence of PWD file: $pwd_not_present\n" );
            print OUT_MESG_TWO ( "# Calculation already performed            : $already_done\n" );
            print OUT_MESG_TWO ( "# ############################################################################\n" );
            close( OUT_MESG_TWO );

            if ( -s $outfile ) {
                system( "cat $outfile $temp_dir/_message1_ $output_directory/$outname $temp_dir/_message2_ > $outfile"."___" )
                && die "\n\n$0: Could not execute \"cat $outfile $temp_dir/_message1_ $output_directory/$outname $temp_dir/_message2_ > $outfile"."___\": $!\n\n";
                system( "mv", $outfile."___", $outfile )
                && die "\n\n$0: Could not execute \"mv $outfile"."___ $outfile\": $!\n\n";
            }
            else {
                system( "cat $temp_dir/_message1_ $output_directory/$outname $temp_dir/_message2_ > $outfile" )
                && die "\n\n$0: Could not execute \"cat $temp_dir/_message1_ $output_directory/$outname $temp_dir/_message2_ > $outfile\": $!\n\n";

            }

            print LOG "$outname\n";

            unlink( "$temp_dir/_message1_", "$temp_dir/_message2_" ); 

        }

        

    } ## End of elsif ( $return_line =~ /^\s*=(\S+)\s+.+\s+(\S+)\s+(\S+)\s+\S+$/ && $ID ne "" )

} ## End of while ( $return_line = <IN> )

close( IN );
close( LOG );     


open( OUT_MESG_TWO, ">$temp_dir/_message2_" ) || die "\n$0: Cannot create file \"$temp_dir/_message2_\": $!\n";
print OUT_MESG_TWO ( "\n\n# Xrio.pl successfully terminated.\n" );
print OUT_MESG_TWO ( "# Started   : $start_date" );
print OUT_MESG_TWO ( "# Terminated: ".`date`."\n" );
print OUT_MESG_TWO ( "# Successful calculations                  : $successful\n" );
print OUT_MESG_TWO ( "# No calculation due to absence of PWD file: $pwd_not_present\n" );
print OUT_MESG_TWO ( "# Calculation already performed            : $already_done\n\n" );
close( OUT_MESG_TWO ); 

if ( -s $outfile ) {
    if ( $message1 ne "" ) { 
        open( OUT_MESG_ONE, ">$temp_dir/_message1_" ) || die "\n$0: Cannot create file \"$temp_dir/_message1_\": $!\n";
        print OUT_MESG_ONE ( $message1 );
        close( OUT_MESG_ONE );
        system( "cat $outfile $temp_dir/_message1_ $temp_dir/_message2_ > $outfile"."___" )
        && die "$0: Could not execute \"cat $outfile $temp_dir/_message1_ $temp_dir/_message2_ > $outfile"."___\": $!";
    } 
    else {
        system( "cat $outfile $temp_dir/_message2_ > $outfile"."___" )
        && die "$0: Could not execute \"cat $outfile $temp_dir/_message2_ > $outfile"."___\": $!";
    }
    system( "mv", $outfile."___", $outfile )
    && die "$0: Could not execute \"mv $outfile"."___ $outfile\": $!";     
}
else {
    open( OUT_MESG_ONE, ">$temp_dir/_message1_" ) || die "\n$0: Cannot create file \"$temp_dir/_message1_\": $!\n";
    print OUT_MESG_ONE ( $message1 );
    close( OUT_MESG_ONE );
    system( "cat $temp_dir/_message1_ $temp_dir/_message2_ > $outfile" )
    && die "$0: Could not execute \"cat $temp_dir/_message1_ $temp_dir/_message2_ > $outfile\": $!";
}

unlink( "$temp_dir/_message1_", "$temp_dir/_message2_" ); 

rmdir( $temp_dir ) || die "\n$0: Unexpected failure (could not remove: $temp_dir): $!\n";

print( "\n\nXrio.pl successfully terminated.\n" );
print( "Successful calculations                  : $successful\n" );
print( "No calculation due to absence of PWD file: $pwd_not_present\n" );
print( "Calculation already performed            : $already_done\n" );
print( "Started   : $start_date" );
print( "Terminated: ".`date`."\n" );
print( "\n" );

exit( 0 );



# Methods
# -------



# Gets the gathering cutoff per sequence from a HMM file.
# 
# One argument: the HMM file name
# Returns the gathering cutoff per sequence, 2000 upon failure
# Last modified: 07/11/01
sub getGA1cutoff {

    my $infile      = $_[ 0 ];
    my $return_line = "";
    my $GA          = 2000;

    &testForTextFilePresence( $infile );

    open( H, "$infile" ) || die "\n\n$0: Unexpected error: Cannot open file <<$infile>>: $!";
    while ( $return_line = <H> ) {
        
        if ( $return_line =~ /^GA\s+(\S+)/ ) { 
            $GA = $1;
            close( H );
            return $GA;
        }
        
    }
    close( H );
    return $GA;

} ## getGA1cutoff





#  1.  A= Name of Pfam family
#  2.  Q= Query file
#  3.  O= Output
#  4.  N= query Name
#  5.  S= Species tree file
#  6.  more options, such I K m 
#  7.  j= Name for temporary directory
sub performRIO {
    my $pfam_name              = $_[ 0 ];
    my $query_file             = $_[ 1 ];
    my $output_file            = $_[ 2 ];
    my $name_for_query         = $_[ 3 ];
    my $species_tree_file      = $_[ 4 ];
    my $more_options           = $_[ 5 ];
    my $tmp_file_rio           = $_[ 6 ];

    my $options_for_rio = "";

    $options_for_rio .= ( " A=".$pfam_name );
    $options_for_rio .= ( " Q=".$query_file );
    $options_for_rio .= ( " O=".$output_file );
    $options_for_rio .= ( " N=".$name_for_query );
    $options_for_rio .= ( " S=".$species_tree_file );
    $options_for_rio .= ( " j=".$tmp_file_rio );
    $options_for_rio .= ( " ".$more_options );

    system( "$RIO_PL 1 $options_for_rio" )
    && die "$0: performRIO: Could not execute \"$RIO_PL 1 $options_for_rio\": $!\n"; 

} ## performRIO



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

    unless ( ( -s $infile ) && ( -f $infile ) && ( -T $infile ) ) {
        &Error( "\"$infile\" does not exist,\n is empty, or is not a plain textfile." );
    }

    open( IN_RSNF, "$infile" ) || die "\n\n$0: Unexpected error: Cannot open file <<$infile>>: $!";
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


 
# Searches the > line of a multiple seq file for a
# query, saves the found entries.
# Three arguments:
# 1. multi Fasta file to search through
# 2. outputfile name
# 3. query
# Last modified: 03/05/01
sub getSequenceFromFastaFile {
 
    my $inputfile  = $_[ 0 ];
    my $outputfile = $_[ 1 ];
    my $query      = $_[ 2 ];
    my $hits       = 0;    

    open( IN_GSFF,  "$inputfile"   )
    || die "\n$0: getSequenceFromFastaFile: Cannot open file <<$inputfile>>: $!\n";
    open( OUT_GSFF, ">$outputfile" )
    || die "\n$0: getSequenceFromFastaFile: Cannot create file <<$outputfile>>: $!\n";
   

    while ( $return_line = <IN_GSFF> ) {
        if ( $return_line =~ /^\s*>.*$query\s+/ ) {
            $hits++;
            print $return_line;
            print OUT_GSFF $return_line;
            $return_line = <IN_GSFF>;
            while ( $return_line && $return_line =~ /^\s*[^>]/ ) {
                print OUT_GSFF $return_line;
                $return_line = <IN_GSFF>;
            }
            last; # In Wormpep there _are_ ambigous CE numbers.
        }
        
    }
    
    close( IN_GSFF );
    close( OUT_GSFF );
    if ( $hits < 1 ) {
        die "\n$0: getSequenceFromFastaFile: Unexpected error: <<$query>> not found.\n";
    }
    if ( $hits > 1 ) {
        die "\n$0: getSequenceFromFastaFile: Unexpected error: <<$query>> is ambigous.\n";
    }

} ## getSequenceFromFastaFile




# Last modified: 03/08/01
sub errorInCommandLine {

    print "\n";
    print " Xrio.pl  $VERSION\n";
    print " -------\n";
    print "\n";
    print " Christian Zmasek (zmasek\@genetics.wustl.edu)\n";
    print "\n";
    print " Purpose. Runs \"rio.pl\" for each Pfam assignment in \"infile\".\n";
    print "\n";
    print " Usage. rio.pl <infile> <species names file> <output directory> <outfile> <logfile>\n";
    print "\n";
    print " infile: has the following format (defined per example):\n";
    print "   >>4R79.1 CE19649   Zinc-binding metalloprotease domain (CAMBRIDGE) protein_id:CAB63429.1\n";
    print "   =Astacin  Astacin (Peptidase family M12A)                296.3    3.8e-85   1\n";
    print "   //\n";
    print "\n";
    print "   >>4R79.2 CE19650   Ras family (CAMBRIDGE) TR:Q9XXA4 protein_id:CAA20282.1\n";
    print "   =ras           Ras family                                208.8    8.1e-59   1\n";
    print "   =FA_desaturase Fatty acid desaturase                       4.5        1.5   1\n";
    print "   =UPF0117       Domain of unknown function DUF36            3.1        3.5   1\n";
    print "   =arf           ADP-ribosylation factor family            -46.0    1.5e-05   1\n";
    print "   //\n";
    print "\n";
    print " species names file: list of species to use for analysis.\n";
    print "\n";
    print " This version uses the CE number as identifier for output files.\n";
    print "\n";

    exit( -1 );
    
} ## errorInCommandLine



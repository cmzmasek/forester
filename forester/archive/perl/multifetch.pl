#!/usr/bin/perl

# multifetch.pl [options] <list of seqs>
#
# Fetch all the seqs on the list. The list is a file with one line
# per sequence; the first field is the key.
#
# Options:
#    -d            : domain fetch - list is in GDF format
#    -n <extra 5'> : include this many extra residues upstream (-d only)
#    -c <extra 3'> : include this many extra residues downstream (-d only)
#    -f            : fetch in FASTA instead of native format
#    -g <file>     : use getseq from <file>, not fetch from main databases.
#                    This always gives FASTA output.
#    -D <database> : specify a source database, same usage as getseq:
#                           -Dsw  SwissProt
#                           -Dpir PIR
#                           -Dem  EMBL
#                           -Dgb  GenBank
#                           -Dwp  WormPep
#                           -Dowl OWL 
 

use FindBin;
use lib $FindBin::Bin;
use rio_module;
require "getopts.pl";


&Getopts('c:n:dfg:D:');
if ($opt_c) { $extra_c  = $opt_c; }
if ($opt_n) { $extra_n  = $opt_n; }
if ($opt_d) { $domains  = 1;      }
if ($opt_f) { $fmtarg = "-Ffasta";} else {$fmtarg = ""; }
if ($opt_g) { $filearg = "-d$opt_g ";} else {$filearg = ""; }
if ($opt_D) { $dbarg = "-D$opt_D "; } else {$dbarg = ""; }


while (<>) {
    if ($domains) {
	if (($name, $from, $to, $source) = /^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/){
	    if ($from < $to) {
		$from -= $opt_n;
		$to   += $opt_c;
	    }
	    else {
		$from += $opt_n;
		$to   -= $opt_c;
	    }
            
	    system("$SFE $filearg $dbarg $fmtarg -r \"$name\" -f $from -t $to \"$source\"")
	    && die "\n\n$0: Unexpected error: Could not execute \"$SFE $filearg $dbarg $fmtarg -r \"$name\" -f $from -t $to \"$source\"\": $!";
        }
    } else {
	if (/^\s*(\S+)/) {
	    $key = $1;
	     
	    system("$SFE $filearg $dbarg $fmtarg \"$key\"")
            && die "\n\n$0: Unexpected error: Could not execute \"$SFE $filearg $dbarg $fmtarg \"$key\"\": $!";
	}
    }
}

# 01/30/02
# CZ
# Added usage of rio_module.pm, $SFE for sfetch.

# Thu Apr 10 18:27:40 1997
#   - added -D option
#   - simplified from six different getseq calls to two


#! /usr/bin/perl

# Usage: p7extract.pl <hmmsearch output file>
#
# Converts hmmsearch output to GLF or GDF. GLF is the default.
# Order is sorted by bit score
#
# Options:
#    -C             : extract Pfam coverage statistics (NAR paper)
#    -d             : extract domains in GDF format 
#    -t <x>         : report only hits better than evalue of <x>
#    -s             : include scores in output
#    -e             : include evalues in output  
#    -l             : include negative log evalues in output for easy sorting
#
# Note: p7extract.pl -sel gives the extended GLF format expected by
#       the PROFMARK benchmark scripts

require "getopts.pl";

$ethresh = 0;

&Getopts('Cdt:sel');
if ($opt_C) { $coverage_mode = 1; $gdfmode = 1;}
if ($opt_d) { $gdfmode   = 1; }   
if ($opt_t) { $ethresh   = $opt_t; }
if ($opt_s) { $do_scores = 1; }
if ($opt_e) { $do_eval   = 1; }
if ($opt_l) { $do_log    = 1; }

$status = 1;			# -C only: assume we will fail, 'til proven otherwise

while (<>)
{
    if (/^Query HMM:\s+(\S+)/)            {$hmmname = $1;}        
    if (/^Scores for complete sequences/) {$indom = 0; $inseq = 1;}
    if (/^Parsed for domains/)            {$indom = 1; $inseq = 0;}
    if (/^Histogram of all scores/)       {$indom = 0; $inseq = 0;}
    if (/^Total sequences searched/)      {$status = 0;} # looks like we've seen output
    
    if ( $inseq &&
	(($id, $sc, $ev, $nd) = /(\S+).+\s(\S+)\s+(\S+)\s+(\d+)\s*$/))
    {
	if (($ethresh == 0 || $ev < $ethresh) && $show_key{$id} == 0) 
	{
	    if (! $gdfmode) {
		$show_key{$id} = 1;	   # remember this name
		$show_sc{$id}  = $sc;
		$show_ev{$id}  = $ev;
	    }
	    $numseqs++;
	}
    }

    if ($gdfmode && $indom &&
	(($id, $sqfrom, $sqto, $sc, $ev) =
	 /(\S+)\s+\S+\s+(\d+)\s+(\d+).+\s(\S+)\s+(\S+)\s*$/))
    {
	if (($ethresh == 0 || $ev < $ethresh) && $show_key{$id} == 0) 
	{
	    $key               = "$id/$sqfrom-$sqto";
	    $show_key{$key}    = 1;
	    $show_id{$key}     = $id;
	    $show_sqfrom{$key} = $sqfrom;
	    $show_sqto{$key}   = $sqto;
	    $show_sc{$key}     = $sc;
	    $show_ev{$key}     = $ev;

	    $numdomains++;

	    $domsize = $sqto - $sqfrom + 1;
	    if ($domsize < 0) { $domsize *= -1; }
	    $numresidues += $domsize;
	}
    }

}
 
if ($coverage_mode)
{
    if ($status == 0) {
	printf "%-20s %6d %6d %6d\n", $hmmname, $numseqs, $numdomains, $numresidues;
	exit 0;
    } else {
	printf "%-20s [FAILED]\n", $hmmname;
	exit 1;
    }

}
	

	
foreach $key (sort byscore keys(%show_key))
{
    if ($gdfmode)
    {
	printf("%-24s\t%6d\t%6d\t%15s",
	       $key, $show_sqfrom{$key}, $show_sqto{$key}, $show_id{$key})
    } else {			 
	printf("%-24s", $key);
    }
				# Optional extensions to GDF/GLF
    if ($do_scores) { printf("\t%8s",    $show_sc{$key}); } 
    if ($do_eval)   { printf("\t%12s",   $show_ev{$key}); }
    if ($do_log)    { printf("\t%12.1f", -log($show_ev{$key})); }
    print "\n";
}
    
sub byscore { 
    $show_sc{$b} <=> $show_sc{$a};
}

sub byevalue { 
    $show_ev{$a} <=> $show_ev{$b};
}


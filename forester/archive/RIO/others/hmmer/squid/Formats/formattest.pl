#! /usr/local/bin/perl

$binpath = shift;		

# Suck in the regression data on our file format test suite.
#

print "Format test suite...\t";

open(DAT,"regression.dat") || die "failed to open regression.dat";
$nfiles = 0;
while (<DAT>) {
    if (/^\#/) { next; }
    if (/^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
	$filename[$nfiles] = $1;
	$format[$nfiles]   = $2;
	$seqtype[$nfiles]  = $3;
	$nseq[$nfiles]     = $4;
	$nres[$nfiles]     = $5;
	$shortest[$nfiles] = $6;
	$longest[$nfiles]  = $7;
	if ($8  eq "yes") { $autodetect[$nfiles]   = 1; } else { $autodetect[$nfiles]   = 0; }
	if ($9  eq "yes") { $is_alignment[$nfiles] = 1; } else { $is_alignment[$nfiles] = 0; }
	if ($10 eq "yes") { $is_singleseq[$nfiles] = 1; } else { $is_singleseq[$nfiles] = 0; }
	$nfiles++;
    }
}
close(DAT);

# Test 1.
# Run seqstat on every file in two modes;
# autodetecting (if format allows it), then forcing a format with --informat.
#
for ($i = 0; $i < $nfiles; $i++) {
    if ($autodetect[$i]) {
	$output = `$binpath/seqstat $filename[$i]`;
	if ($? != 0) { die "seqstat failed, autodetecting, on $filename[$i]"; }
	($ns, $nr, $fr, $to) = &parse_seqstat($output);
	if ($ns != $nseq[$i] ||
	    $nr != $nres[$i] ||
	    $fr != $shortest[$i] ||
	    $to != $longest[$i])
	{ die "seqstat regression failed, autodetecting, on $filename[$i]"; }
    }
    $output = `$binpath/seqstat --informat $format[$i] $filename[$i]`;
	if ($? != 0) { die "seqstat failed, using --informat, on $filename[$i]"; }
	($ns, $nr, $fr, $to) = &parse_seqstat($output);
	if ($ns != $nseq[$i] ||
	    $nr != $nres[$i] ||
	    $fr != $shortest[$i] ||
	    $to != $longest[$i])
	{ die "seqstat regression failed, using --informat, on $filename[$i]"; }
}

# Test 2.
# Reformatting tests.
#
for ($i = 0; $i < $nfiles; $i++) {
    for ($j = 0; $j < $nfiles; $j++) {
	if (! $is_alignment[$i] && $is_alignment[$j]) { next; } # can't convert unaligned to aligned
	if (! $is_singleseq[$i] && $is_singleseq[$j]) { next; } # can't convert multiple seqs to single seq format

	`$binpath/sreformat --informat $format[$i] $format[$j] $filename[$i] > formattest.tmp`;
	if ($? != 0) { die "sreformat failed ($format[$i] to $format[$j]) on $filename[$i]"; }
	$output = `$binpath/seqstat --informat $format[$j] formattest.tmp`;
	if ($? != 0) { die "seqstat failed after sreformat ($format[$i] to $format[$j]) on $filename[$i]"; }
	($ns, $nr, $fr, $to) = &parse_seqstat($output);
	if ($ns != $nseq[$i] ||
	    $nr != $nres[$i] ||
	    $fr != $shortest[$i] ||
	    $to != $longest[$i])
	{ die "seqstat regression failed after sreformat ($format[$i] to $format[$j]) on $filename[$i]"; }
    }
}

print "passed.\n";
unlink "formattest.tmp";


# Function: parse_seqstat(file)
#
# Returns the number of sequences in the file,
# and their maximum and minimum length, and their avg. len.
# Dies if 'seqstat' fails.
#
sub parse_seqstat {
    local($output) = shift;
    my ($nseq, $nres, $fromlen, $tolen);	

    if ($output =~ /Number of sequences:\s+(\d+)/)      {$nseq    = $1; }
    if ($output =~ /Total # residues:\s+(\d+)/)         {$nres    = $1; }
    if ($output =~ /Smallest:\s+(\d+)/)                 {$fromlen = $1; }
    if ($output =~ /Largest:\s+(\d+)/)                  {$tolen   = $1; }
    ($nseq, $nres, $fromlen, $tolen);
}



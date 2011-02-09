#! /usr/bin/perl

# Unpack a pfam flatfile, containing many alignments,
# into separate SELEX-format alignment files.
#
# Assumes that ID is the first line in a record, 
# that SQ is the last line before the alignment starts,
# and that there is one aligned sequence per line.
#


################################################################
# PFAMSERVER - The Washington University/St. Louis Pfam web server
# Copyright (C) 1995-1999 Washington University School of Medicine
# Copyright (C) 1995-1999 Sanger Centre/Genome Research Ltd.
# Copyright (C) 1998-1999 Karolinska Institutet Center for Genomics Research
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYRIGHT and LICENSE
#     for details.
# 
################################################################

$cpl = 50;			# 50 sequence characters per line
$/ = "\n//";		# paragraph mode on // separators

while (<>) {
    $in_alignment = 0;
    $nseq = 0;
    @lines = split(/^/);
    while ($line = shift(@lines)) {
	if ($in_alignment) {
	    if    ($line =~ /^\#/) { next; }
	    elsif ($line =~ /^(\S+)\s+(\S+)/) {
		$name[$nseq] = $1;
		$aseq[$nseq] = $2;
		$nseq++;
	    }
	}
	elsif ($line =~ /^\#=GF ID   (\S+)\s*$/) { 
	    $root = $1;
	    print "working on $root\n";
	    if (-e "$root") {
		system ("mv $root $root.orig");
		print "$root exists -- moved to $root.orig\n";
	    }
	    open(SELEX,">$root") || die;
	    print SELEX "#=ID $root\n";
	}
	elsif ($line =~ /^\#=GF AC   (.+)$/) { print SELEX "#=AC $1\n"; }
	elsif ($line =~ /^\#=GF DE   (.+)$/) { print SELEX "#=DE $1\n"; }

	elsif ($line =~ /^\#=GF GA   (\S+)\s+(\S+)/) 
	{ print SELEX "#=GA $1 $2\n"; }
	
	elsif ($line =~ /^\#=GF TC   (\S+) (\S+)/) 
	{ print SELEX "#=TC $1 $2\n"; }
	    
	elsif ($line =~ /^\#=GF NC   (\S+) (\S+)/) 
	{ print SELEX "#=NC $1 $2\n"; }
	
	elsif ($line =~ /^\#=GF SQ   \d+/) {
	    print SELEX "# $line";
	    $in_alignment = 1;
	}
	elsif ($line =~ /^\/\//) {
	    last;
	}
	else {
	    print SELEX "# $line";
	}
    }
    
				# figure out maximum name length
    $maxnamelen = 0;
    for ($idx = 0; $idx < $nseq; $idx++) {
	if (length($name[$idx]) > $maxnamelen) {
	    $maxnamelen = length($name[$idx]);
	}
    }
				# break the alignment across
				# multiple lines
    $alen = length($aseq[0]);
    for ($pos = 0; $pos < $alen; $pos += $cpl) {
	for ($idx = 0; $idx < $nseq; $idx++) {
	    printf(SELEX "%-${maxnamelen}s %s\n", 
		   $name[$idx], substr($aseq[$idx], $pos, $cpl));
	}
	print SELEX "\n";
    }
    close SELEX;
}
	    

#!/usr/local/bin/ruby -w


input = ARGV[ 0 ]

system( "hmmscan --nobias --domtblout " + input + "_hmmscan_260_10 -E 10 Pfam-A.hmm " + input + ".ni.fasta" )

system( "hsp " + input + "_hmmscan_260_10 " + input + "_hmmscan_260_10_domain_table" )

system( "d2f " + input + "_hmmscan_260_10_domain_table " + input + ".ni.fasta " + input + "_hmmscan_260_10.dff" )
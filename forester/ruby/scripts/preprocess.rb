#!/usr/local/bin/ruby -w


if ( ARGV == nil || ARGV.length != 2 )
  puts( "usage: preprocess.rb <input> <path to Pfam A HMMs>" )         
  exit( -1 )
end 

input = ARGV[ 0 ]

pfam = ARGV[ 1 ]

system( "hmmscan --nobias --domtblout " + input + "_hmmscan_260_10 -E 10 " + pfam + " "  + input + ".ni.fasta" )

system( "hsp " + input + "_hmmscan_260_10 " + input + "_hmmscan_260_10_domain_table" )

system( "d2f " + input + "_hmmscan_260_10_domain_table " + input + ".ni.fasta " + input + "_hmmscan_260_10.dff" )
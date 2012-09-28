#!/usr/local/bin/ruby -w
#
# = run_phylo_pipeline
#
# Copyright::  Copyright (C) 2010 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id Exp $
#
#


module Evoruby

  class RunPhyloPipeline

    PFAM      = "/home/czmasek/DATA/PFAM/PFAM260X/"
    HMMSCAN  = "/home/czmasek/SOFTWARE/HMMER/hmmer-3.0/src/hmmscan"
    HSP       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/hsp.rb"
    D2F       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/d2f.rb"
    DSX       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/dsx.rb"


    def run
      unless ARGV.length == 4
        error "arguments are: <fasta formatted inputfile> <hmm-name> <min-length> <neg e-value exponent>"
      end




      input       = ARGV[ 0 ]
      hmm         = ARGV[ 1 ]
      length      = ARGV[ 2 ].to_i
      e_value_exp = ARGV[ 3 ].to_i

      hmmscan_option = "--nobias"
      e_for_hmmscan = 100

      if e_value_exp < 0
        error "e-value exponent cannot be negative"
      end
      if length <= 1
        error "length exponent cannot be smaller than or equal to 1"
      end

      base_name = nil
      if input.downcase.end_with?( ".fasta" )
        base_name = input[ 0 .. input.length - 7 ]
      elsif input.downcase.end_with?( ".fsa" )
        base_name = input[ 0 .. input.length - 5 ]
      else
        base_name = input
      end

      puts "1. hmmscan:"
      cmd = "#{HMMSCAN} #{hmmscan_option} --domtblout #{base_name}_hmmscan_#{e_for_hmmscan.to_s} -E #{e_for_hmmscan.to_s} #{PFAM}Pfam-A.hmm #{input}"
      run_command( cmd )
      puts

      puts "2. hmmscan to simple domain table:"
      cmd = "#{HSP} #{base_name}_hmmscan_#{e_for_hmmscan.to_s} #{base_name}_hmmscan_#{e_for_hmmscan.to_s}_domain_table"
      run_command( cmd )
      puts

      puts "3. domain table to forester format:"
      cmd = "#{D2F} -e=10 #{base_name}_hmmscan_#{e_for_hmmscan.to_s}_domain_table #{input} #{base_name}_hmmscan_#{e_for_hmmscan.to_s}.dff"
      run_command( cmd )
      puts

      puts "4. dsx:"
      cmd = "#{DSX} -d -e=1e-#{e_value_exp.to_s} -l=#{length} #{hmm} #{base_name}.hmmsscan_#{e_for_hmmscan.to_s} #{input} #{base_name}_#{hmm}_e#{e_value_exp.to_s}_#{length}"
      run_command( cmd )
      puts

    end

    def run_command( cmd )
      puts cmd
      `#{cmd}`
    end

    def error msg
      puts
      puts msg
      puts
      exit
    end

  end

  p = RunPhyloPipeline.new()

  p.run()

end




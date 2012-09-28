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



    def run
      unless ARGV.length == 4
        error "arguments are: <fasta formatted inputfile> <hmm-name> <min-length> <neg e-value exponent>"
      end

      hmmscan   = "/home/czmasek/SOFTWARE/HMMER/hmmer-3.0/src/hmmscan"
      hmmsearch = "/home/czmasek/SOFTWARE/HMMER/hmmer-3.0/src/hmmsearch"
      hsp       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/hsp.rb"
      d2f       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/d2f.rb"
      dsx       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/dsx.rb"

      input       = ARGV[ 0 ]
      hmm         = ARGV[ 1 ]
      length      = ARGV[ 2 ]
      e_value_exp = ARGV[ 3 ]
      do_domain_combination_analysis = true

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

      if do_domain_combination_analysis

        puts "hmmscan:"
        cmd = "#{hmmscan} --nobias --domtblout #{base_name}_hmmscan_10 -E 10 #{PFAM}Pfam-A.hmm #{input}"
        run_command( cmd )
        puts

        puts "hmmscan to simple domain table:"
        cmd = "#{hsp} #{base_name}_hmmscan_10 #{base_name}_hmmscan_10_domain_table"
        run_command( cmd )
        puts

        puts "domain table to forester format:"
        cmd = "#{d2f} -e=10 #{base_name}_hmmscan_10_domain_table #{input} #{base_name}_hmmscan_10.dff"
        run_command( cmd )
        puts

      end

      puts "hmmsearch:"
      cmd = "#{hmmsearch} --nobias -E 1000 --domtblout #{base_name}.hmmsearch_#{hmm}  #{PFAM}PFAM_A_HMMs/#{hmm}.hmm #{input}"
      run_command( cmd )
      puts

      puts "dsx:"
      cmd = "#{dsx} -d -e=1e-#{e_value_exp.to_s} -l=#{length} #{hmm} #{base_name}.hmmsearch_#{hmm} #{input} #{base_name}_#{hmm}_e#{e_value_exp.to_s}_#{length}"
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




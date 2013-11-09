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

    PFAM      = "/home/czmasek/DATA/PFAM/PFAM270X/"
    HMMSCAN  = "/home/czmasek/SOFTWARE/HMMER/hmmer-3.0/src/hmmscan"
    HSP       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/hsp.rb"
    D2F       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/d2f.rb"
    DSX       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/dsx.rb"
    TAP       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/tap.rb"

    def run
      unless ARGV.length >= 2 && ARGV.length <= 4
        error "arguments are:  <min-length> " +
         "<neg E-value exponent for domain extraction> [E-value for hmmscan, default is 10] [hmmscan option, default is --nobias, --max for no heuristics]"
      end

      length      = ARGV[ 0 ].to_i
      e_value_exp = ARGV[ 1 ].to_i

      e_for_hmmscan = 10
      hmmscan_option = "--nobias"

      if ARGV.length == 4
        hmmscan_option = ARGV[ 3 ]
      end
      if ARGV.length == 3 || ARGV.length == 4
        e_for_hmmscan = ARGV[ 2 ].to_i
      end

      if e_value_exp < 0
        error "E-value exponent for domain extraction cannot be negative"
      end
      if length <= 1
        error "length cannot be smaller than or equal to 1"
      end
      if e_for_hmmscan < 1
        error "E-value for hmmscan cannot be smaller than 1"
      end

      input_files = Dir.entries(".").select { |f| !File.directory?( f ) && f.downcase.end_with?( ".fasta" ) }

      puts "Input files:"
      input_files.each do | input |
        puts input
      end
      puts

      counter = 1
      input_files.each do | input |

        puts counter.to_s + "/" +  input_files.size.to_s + " " + input + ": "

        counter += 1

        hmm_name = ""

        if input.downcase.end_with?( "_ni.fasta" )
          hmm_name = input[ 0 .. input.length - 10 ]
        elsif input.downcase.end_with?( ".fasta" )
          hmm_name = input[ 0 .. input.length - 7 ]
          puts
          puts "a. identifier normalization:"
          cmd = "#{TAP} #{input} #{hmm_name}_ni.fasta #{hmm_name}.nim"
          run_command( cmd )
          input = hmm_name + "_ni.fasta"
        else
          error "illegal name: " + input
        end

        Dir.mkdir( hmm_name )

        puts
        puts "b. hmmscan:"
        cmd = "#{HMMSCAN} #{hmmscan_option} --domtblout #{hmm_name}/#{hmm_name}_hmmscan_#{e_for_hmmscan.to_s} -E #{e_for_hmmscan.to_s} #{PFAM}Pfam-A.hmm #{input}"
        run_command( cmd )
        puts

        puts "c. hmmscan to simple domain table:"
        cmd = "#{HSP} #{hmm_name}/#{hmm_name}_hmmscan_#{e_for_hmmscan.to_s} #{hmm_name}/#{hmm_name}_hmmscan_#{e_for_hmmscan.to_s}_domain_table"
        run_command( cmd )
        puts

        puts "d. domain table to forester format:"
        cmd = "#{D2F} -e=10 #{hmm_name}/#{hmm_name}_hmmscan_#{e_for_hmmscan.to_s}_domain_table #{input} #{hmm_name}/#{hmm_name}_hmmscan_#{e_for_hmmscan.to_s}.dff"
        run_command( cmd )
        puts

        puts "e. dsx:"
        cmd = "#{DSX} -d -e=1e-#{e_value_exp.to_s} -l=#{length} #{hmm_name} #{hmm_name}/#{hmm_name}_hmmscan_#{e_for_hmmscan.to_s} #{input} #{hmm_name}/#{hmm_name}__#{hmm_name}__ee#{e_value_exp.to_s}_#{length}"
        run_command( cmd )
        puts

      end

    end

    def run_command cmd
      puts cmd
      `#{cmd}`
    end

    def get_base_name n
      if n.downcase.end_with?( "_ni.fasta" )
        n[ 0 .. n.length - 10 ]
      elsif n.downcase.end_with?( ".fasta" )
        n[ 0 .. n.length - 7 ]
      else
        error "illegal name: " + n
      end
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




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

require 'fileutils'

module Evoruby

  class RunPhyloPipeline

    LAUNCH_ANALYSIS = true
    HOME          = "/home/czmasek/"
    FORESTER_RUBY = "#{HOME}SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/"
    PFAM          = "#{HOME}DATA/PFAM/PFAM270X/"
    HMMSCAN       = "#{HOME}SOFTWARE/HMMER/hmmer-3.0/src/hmmscan"
    HSP           = "#{FORESTER_RUBY}hsp.rb"
    D2F           = "#{FORESTER_RUBY}d2f.rb"
    DSX           = "#{FORESTER_RUBY}dsx.rb"
    TAP           = "#{FORESTER_RUBY}tap.rb"
    PF            = "#{FORESTER_RUBY}phylogeny_factory.rb"
    TEMPLATE_FILE = '00_phylogeny_factory.template'

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
        id_norm = false
        orig_input = input

        if input.downcase.end_with?( "_ni.fasta" )
          hmm_name = input[ 0 .. input.length - 10 ]
        elsif input.downcase.end_with?( ".fasta" )
          hmm_name = input[ 0 .. input.length - 7 ]
          unless File.exist? hmm_name
            id_norm = true
            puts
            puts "a. identifier normalization:"
            cmd = "#{TAP} #{input} #{hmm_name}_ni.fasta #{hmm_name}.nim"
            run_command( cmd )
            input = hmm_name + "_ni.fasta"
          else
            input = hmm_name + "/" + hmm_name + "_ni.fasta"
            unless File.exist? input
              error "expected to already exist: " + input
            end
            puts "a. identifier normalization already done:" + input
          end
        else
          error "illegal name: " + input
        end

        unless File.exist? hmm_name
          Dir.mkdir( hmm_name )
        end

        puts
        hmmscan_output = hmm_name + "/" + hmm_name + "_hmmscan_" + e_for_hmmscan.to_s
        unless File.exist? hmmscan_output
          puts "b. hmmscan:"
          cmd = "#{HMMSCAN} #{hmmscan_option} --domtblout #{hmmscan_output} -E #{e_for_hmmscan.to_s} #{PFAM}Pfam-A.hmm #{input}"
          run_command( cmd )
        else
          puts "b. hmmscan output already exists: " + hmmscan_output
        end
        puts


        hsp_output = hmm_name + "/" + hmm_name + "_hmmscan_#{e_for_hmmscan.to_s}_domain_table"
        unless File.exist? hsp_output
          puts "c. hmmscan to simple domain table:"
          cmd = "#{HSP} #{hmmscan_output} #{hsp_output}"
          run_command( cmd )
        else
          puts "c. hmmscan to simple domain table output already exists: " + hsp_output
        end
        puts

        d2f_output = "#{hmm_name}/#{hmm_name}_hmmscan_#{e_for_hmmscan.to_s}.dff"
        unless File.exist? d2f_output
          puts "d. domain table to forester format:"
          cmd = "#{D2F} -e=10 #{hsp_output} #{input} #{d2f_output}"
          run_command( cmd )
        else
          puts "d. domain table to forester format output already exists: " + d2f_output
        end
        puts


        dsx_output = "#{hmm_name}/#{hmm_name}__#{hmm_name}__ee#{e_value_exp.to_s}_#{length}"
        unless File.exist? d2f_output
          puts "e. dsx:"
          cmd = "#{DSX} -d -e=1e-#{e_value_exp.to_s} -l=#{length} #{hmm_name} #{hmmscan_output} #{input} #{dsx_output}"
          run_command( cmd )
        else
          puts "e. dsx output already exists: " + dsx_output
        end
        puts

        if id_norm
          FileUtils.mv "#{hmm_name}_ni.fasta", "#{hmm_name}/#{hmm_name}_ni.fasta"
          FileUtils.mv "#{hmm_name}.nim", "#{hmm_name}/#{hmm_name}.nim"
          FileUtils.cp orig_input, "#{hmm_name}/#{orig_input}"
        end

        msa_dir = hmm_name + "/msa_ee#{e_value_exp.to_s}_#{length}"
        msa_100_dir =hmm_name + "/msa100_ee#{e_value_exp.to_s}_#{length}"

        unless File.exist? msa_dir
          Dir.mkdir( msa_dir )
        end
        unless File.exist? msa_100_dir
          Dir.mkdir( msa_100_dir )
        end

        run_1 = false
        run_100 = false

        unless File.exist? "#{msa_dir}/#{dsx_output}"
          run_1 = true
          FileUtils.cp "#{dsx_output}.fasta", "#{msa_dir}/#{dsx_output}"
        end

        unless File.exist? "#{msa_100_dir}/#{dsx_output}"
          run_100 = true
          FileUtils.cp "#{dsx_output}.fasta", "#{msa_100_dir}/#{dsx_output}"
        end

        if File.exist?( TEMPLATE_FILE )
          if run_1
            FileUtils.cp TEMPLATE_FILE, msa_dir
          end
          if run_100
            FileUtils.cp TEMPLATE_FILE, msa_100_dir
          end

          if LAUNCH_ANALYSIS
            puts "f. analysis:"
            if run_1
              Dir.chdir msa_dir
              run_command "#{PF} -b=1 -s"
              Dir.chdir "../.."
            end
            if run_100
              Dir.chdir msa_100_dir
              run_command "#{PF} -b=100 -s"
              Dir.chdir "../.."
            end
            puts
          end
        end

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




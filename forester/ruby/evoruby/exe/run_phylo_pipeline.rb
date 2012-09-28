#!/usr/local/bin/ruby -w
#
# = run_phylo_pipeline
#
# Copyright::  Copyright (C) 2010 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: run_phylo_pipeline.rb,v 1.15 2010/10/09 02:35:42 cmzmasek Exp $
#
#


#  hmmscan --nobias --domtblout <BACTH_CHIPI>_hmmscan_250_10 -E 10 /home/czmasek/DATA/PFAM/PFAM250/Pfam-A.hmm <BACTH_CHIPI>.fasta

#  hsp <BACTH_CHIPI>_hmmscan_250_10 <BACTH_CHIPI>_hmmscan_250_10_domain_table

#  d2f -e=10 <BACTH_CHIPI>_hmmscan_250_10_domain_table <BACTH_CHIPI>.fasta <BACTH_CHIPI>_hmmscan_250_10.dff

# hmmsearch --nobias -E 1000 --domtblout <BACTH_CHIPI>.hmmsearch_SusD  <~/DATA/PFAM/PFAM250/PFAM_A_HMMs/SusD.hmm> BACTH_CHIPI.fasta

# dsx -dd -e=<1e-2> -l=<200> <BACTH_CHIPI>.hmmsearch_SusD <BACTH_CHIPI>.fasta BACTH_CHIPI_e2_200


module Evoruby

  class RunPhyloPipeline

    def run
      unless ARGV.length == 4
        puts
        puts "arguments are: [inputfile].fasta [hmm-name] [min-length] [neg e-value exponent]"
        puts
        exit
      end

      pfam      = "/home/czmasek/DATA/PFAM/PFAM260X/"
      hmmscan   = "/home/czmasek/SOFTWARE/HMMER/hmmer-3.0/src/hmmscan"
      hmmsearch = "/home/czmasek/SOFTWARE/HMMER/hmmer-3.0/src/hmmsearch"
      hsp       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/hsp.rb"
      d2f       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/d2f.rb"
      dsx       = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/dsx.rb"

      base_name   = ARGV[ 0 ]
      hmm         = ARGV[ 1 ]
      length      = ARGV[ 2 ]
      e_value_exp = ARGV[ 3 ]
      do_domain_combination_analysis = true

      if do_domain_combination_analysis

        cmd = "#{hmmscan} --nobias --domtblout #{base_name}_hmmscan_10 -E 10 #{pfam}Pfam-A.hmm #{base_name}.fasta"
        run_command( cmd )

        cmd = "#{hsp} #{base_name}_hmmscan_10 #{base_name}_hmmscan_10_domain_table"
        run_command( cmd )

        cmd = "#{d2f} -e=10 #{base_name}_hmmscan_10_domain_table #{base_name}.fasta #{base_name}_hmmscan_10.dff"
        run_command( cmd )

      end

      cmd = "#{hmmsearch} --nobias -E 1000 --domtblout #{base_name}.hmmsearch_#{hmm}  #{pfam}PFAM_A_HMMs/#{hmm}.hmm #{base_name}.fasta"
      run_command( cmd )

      cmd = "#{dsx} -dd -e=1e-#{e_value_exp.to_s} -l=#{length} #{base_name}.hmmsearch_#{hmm} #{base_name}.fasta #{base_name}_#{hmm}_e#{e_value_exp.to_s}_#{length}"
      run_command( cmd )

    end

    def run_command( cmd )
      puts cmd
      `#{cmd}`
    end

  end

  p = RunPhyloPipeline.new()

  p.run()

end




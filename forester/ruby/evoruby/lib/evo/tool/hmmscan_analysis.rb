#
# = lib/evo/tool/hmmscan_summary.rb - HmmscanSummary class
#
# Copyright::  Copyright (C) 2012 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: hmmscan_parser.rb,v 1.5 2010/12/13 19:00:11 cmzmasek Exp $
#

require 'set'

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/io/parser/hmmscan_parser'
require 'lib/evo/msa/msa'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/writer/fasta_writer'

module Evoruby

  class HmmscanAnalysis

    PRG_NAME       = "hsp"
    PRG_VERSION    = "2.001"
    PRG_DESC       = "hmmscan summary"
    PRG_DATE       = "2013.10.23"
    COPYRIGHT      = "2013 Christian M Zmasek"
    CONTACT        = "phyloxml@gmail.com"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    DELIMITER_OPTION              = "d"
    SPECIES_OPTION                = "s"
    I_E_VALUE_THRESHOLD_OPTION    = "ie"
    FS_E_VALUE_THRESHOLD_OPTION   = "pe"
    HMM_FOR_PROTEIN_OUTPUT        = "m"
    IGNORE_DUF_OPTION             = "i"
    PARSE_OUT_DESCRIPITION_OPTION = "a"
    HELP_OPTION_1                 = "help"
    HELP_OPTION_2                 = "h"

    USE_AVOID_HMMS = false
    AVOID_HHMS = [ "RRM_1", "RRM_2", "RRM_3", "RRM_4", "RRM_5", "RRM_6" ]
    LIMIT_FOR_CLOSE_DOMAINS = 20

    def initialize
      @domain_counts = Hash.new
    end

    def run

      #   Util.print_program_information( PRG_NAME,
      #     PRG_VERSION,
      #     PRG_DESC,
      #     PRG_DATE,
      #     COPYRIGHT,
      #     CONTACT,
      #     WWW,
      #     STDOUT )

      begin
        cla = CommandLineArguments.new( ARGV )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      end

      if ( cla.is_option_set?( HELP_OPTION_1 ) ||
           cla.is_option_set?( HELP_OPTION_2 ) )
        print_help
        exit( 0 )
      end

      if ( cla.get_number_of_files != 1 &&  cla.get_number_of_files != 2 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push( DELIMITER_OPTION )
      allowed_opts.push( I_E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( FS_E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( IGNORE_DUF_OPTION )
      allowed_opts.push( PARSE_OUT_DESCRIPITION_OPTION )
      allowed_opts.push( HMM_FOR_PROTEIN_OUTPUT )
      allowed_opts.push( SPECIES_OPTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
          "unknown option(s): " + disallowed,
          STDOUT )
      end

      inpath = cla.get_file_name( 0 )
      seq_file_path = nil

      if ( cla.get_number_of_files == 2 )
        seq_file_path = cla.get_file_name( 1 )
      end

      msa = nil
      if seq_file_path != nil
        msa = read_fasta_file(seq_file_path )
      end

      column_delimiter = "\t"
      if ( cla.is_option_set?( DELIMITER_OPTION ) )
        begin
          column_delimiter = cla.get_option_value( DELIMITER_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      i_e_value_threshold = -1
      if ( cla.is_option_set?( I_E_VALUE_THRESHOLD_OPTION ) )
        begin
          i_e_value_threshold = cla.get_option_value_as_float( I_E_VALUE_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( i_e_value_threshold < 0.0 )
          Util.fatal_error( PRG_NAME, "attempt to use a negative i-E-value threshold", STDOUT )
        end
      end

      fs_e_value_threshold = -1
      if ( cla.is_option_set?( FS_E_VALUE_THRESHOLD_OPTION ) )
        begin
          fs_e_value_threshold = cla.get_option_value_as_float( FS_E_VALUE_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( fs_e_value_threshold < 0.0 )
          Util.fatal_error( PRG_NAME, "attempt to use a negative E-value threshold", STDOUT )
        end
      end

      hmm_for_protein_outputs = []
      if ( cla.is_option_set?( HMM_FOR_PROTEIN_OUTPUT ) )
        begin
          hmm_for_protein_output = cla.get_option_value( HMM_FOR_PROTEIN_OUTPUT )
          hmm_for_protein_outputs = hmm_for_protein_output.split( "~" );
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      species = "HUMAN"
      if ( cla.is_option_set?( SPECIES_OPTION ) )
        begin
          species = cla.get_option_value( SPECIES_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      ignore_dufs = false
      if ( cla.is_option_set?( IGNORE_DUF_OPTION ) )
        ignore_dufs = true
      end

      parse_descriptions = false
      if ( cla.is_option_set?( PARSE_OUT_DESCRIPITION_OPTION ) )
        parse_descriptions = true
      end


      begin
        parse( inpath,
          column_delimiter,
          i_e_value_threshold,
          ignore_dufs,
          parse_descriptions,
          fs_e_value_threshold,
          hmm_for_protein_outputs,
          species,
          msa        )
      rescue IOError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      end

    end # def run

    private

    def read_fasta_file( input )
      f = MsaFactory.new()
      msa = nil
      begin
        msa = f.create_msa_from_file( input, FastaParser.new() )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end
      msa
    end


    # raises ArgumentError, IOError
    def parse( inpath,
        column_delimiter,
        i_e_value_threshold,
        ignore_dufs,
        get_descriptions,
        fs_e_value_threshold,
        hmm_for_protein_outputs,
        species,
        msa )

      Util.check_file_for_readability( inpath )

      hmmscan_parser = HmmscanParser.new( inpath )
      results = hmmscan_parser.parse


      query     = ""
      desc      = ""
      model     = ""
      env_from  = ""
      env_to    = ""
      i_e_value = ""

      hmmscan_results_per_protein = []

      prev_query = ""

      results.each do | r |
        model     = r.model
        query     = r.query
        i_e_value = r.i_e_value
        env_from  = r.env_from
        env_to    = r.env_to


        if !prev_query.empty? && prev_query != query
          if !hmmscan_results_per_protein.empty?
            process_hmmscan_results_per_protein( hmmscan_results_per_protein,
              fs_e_value_threshold,
              hmm_for_protein_outputs,
              i_e_value_threshold,
              species,
              msa            )
          end
          hmmscan_results_per_protein.clear
        end
        prev_query = query

        if USE_AVOID_HMMS
          if !AVOID_HHMS.include? r.model
            hmmscan_results_per_protein << r
          end
        else
          hmmscan_results_per_protein << r
        end

      end

      if !hmm_for_protein_outputs.empty? && !hmmscan_results_per_protein.empty?
        process_hmmscan_results_per_protein( hmmscan_results_per_protein,
          fs_e_value_threshold,
          hmm_for_protein_outputs,
          i_e_value_threshold,
          species,
          msa        )
      end


    end # def parse

    def process_id( id )
      if  id =~ /(sp|tr)\|\S+\|(\S+)/
        id = $2
      end
      id
    end



    def process_hmmscan_results_per_protein( hmmscan_results_per_protein,
        fs_e_value_threshold,
        target_hmms,
        i_e_value_threshold,
        species,
        msa )

      raise StandardError, "target hmms is empty" if target_hmms.length < 1
      raise StandardError, "results is empty" if hmmscan_results_per_protein.length < 1

      # filter according to i-Evalue threshold
      # abort if fs Evalue too high


      if fs_e_value_threshold >= 0.0
        hmmscan_results_per_protein.each do | r |
          target_hmms.each do | hmm |
            if r.model == hmm && r.fs_e_value > fs_e_value_threshold
              return
            end
          end
        end
      end

      #  dcs = []
      hmmscan_results_per_protein_filtered = []
      matched = Set.new

      hmmscan_results_per_protein.each do | r |
        if i_e_value_threshold < 0 || r.i_e_value <= i_e_value_threshold
          hmmscan_results_per_protein_filtered << r
          target_hmms.each do | hmm |
            if r.model == hmm

              matched << hmm
              break
            end
          end

        end
      end

      if matched.length < target_hmms.length
        return
      end
      if  hmmscan_results_per_protein_filtered.length < 1
        return
      end

      hmmscan_results_per_protein_filtered.sort! { |r1,r2| r1.env_from <=> r2.env_from }

      owns = []
      target_hmms.each do | hmm |
        hmmscan_results_per_protein_filtered.each do | r |
          if r.model == hmm
            owns << r
            break
          end
        end
      end

      s = ""
      query = nil
      owns.each do | own |
        s << own.query + "\t"
        query = own.query
      end
      s << species + "\t"
      owns.each do | own |
        s << own.fs_e_value.to_s + "\t"
      end
      owns.each do | own |
        s << own.qlen.to_s + "\t" #TODO !
      end

      #   dcs.each do | dc |
      #     s << dc.to_s + "\t"
      #   end
      s << hmmscan_results_per_protein_filtered.length.to_s + "\t"
      hmmscan_results_per_protein_filtered.each do | r |
        s << r.model + " "
      end
      s << "\t"



      # overview = make_overview( hmmscan_results_per_protein_filtered, hmm_for_protein_output )

      #   s << overview  + "\t"

      #  s << calc_linkers(  hmmscan_results_per_protein_filtered, hmm_for_protein_output )  + "\t"

      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |

        if  prev_r != nil
          s << make_interdomain_sequence( r.env_from - prev_r.env_to - 1 )

          if ( target_hmms.length == 2 && prev_r.model == target_hmms[ 0 ] && r.model == target_hmms[ 1 ] )
            puts "xxx"
            linker( prev_r.env_to, r.env_from, query, msa )
          end


        else
          s << make_interdomain_sequence( r.env_from, false )

        end
        s << r.model
        s << "["
        s << r.env_from.to_s << "-" << r.env_to.to_s
        s << " " << r.i_e_value.to_s
        s << "]"
        prev_r = r
      end
      # s << make_interdomain_sequence( own.qlen - prev_r.env_from, false )
      puts s
    end


    def linker( first, last , query , msa)
      puts first.to_s + "-" + last.to_s
      if ( last - first >= 1 )
        seq = msa.get_by_name( query, true, false )
        linker = seq.get_subsequence( first -1 , last - 1 )
        puts linker.get_sequence_as_string
      end

    end


    def calc_linkers(  hmmscan_results_per_protein_filtered, hmm_for_protein_output )
      linkers = ""
      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if r.model == hmm_for_protein_output
          if  prev_r != nil
            linkers << ( r.env_from - prev_r.env_to - 1 ).to_s + " "
          end
          prev_r = r
        end
      end
      linkers
    end



    def make_overview( hmmscan_results_per_protein_filtered, hmm_for_protein_output )
      overview = ""
      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if r.model == hmm_for_protein_output
          if prev_r == nil
            overview << hmm_for_protein_output
          else
            if  ( r.env_from - prev_r.env_to - 1 ) <= LIMIT_FOR_CLOSE_DOMAINS
              overview << "~" << hmm_for_protein_output
            else
              overview << "----" << hmm_for_protein_output
            end
          end
          prev_r = r
        end
      end
      overview
    end

    def make_interdomain_sequence( d, mark_short = true )
      s = ""
      d /= 20
      if d >= 10
        s << "----//----"
      elsif d >= 1
        d.times do
          s << "-"
        end
      elsif mark_short
        s << "~"
      end
      s
    end


    def print_help()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb [options] <hmmscan outputfile> <outputfile>" )
      puts()
      puts( "  options: -" + DELIMITER_OPTION + ": column delimiter for outputfile, default is TAB" )
      puts( "           -" + I_E_VALUE_THRESHOLD_OPTION  + ": i-E-value threshold, default is no threshold" )
      puts( "           -" + PARSE_OUT_DESCRIPITION_OPTION  + ": parse query description (in addition to query name)" )
      puts( "           -" + IGNORE_DUF_OPTION  + ": ignore DUFs" )
      puts( "           -" + FS_E_VALUE_THRESHOLD_OPTION  + ": E-value threshold for full protein sequences, only for protein summary" )
      puts( "           -" + HMM_FOR_PROTEIN_OUTPUT + ": HMM for protein summary" )
      puts( "           -" + SPECIES_OPTION + ": species for protein summary" )
      puts()
    end

  end # class

end # module Evoruby
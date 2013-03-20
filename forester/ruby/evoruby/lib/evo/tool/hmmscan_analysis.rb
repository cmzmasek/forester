#
# = lib/evo/tool/hmmscan_summary.rb - HmmscanSummary class
#
# Copyright::  Copyright (C) 2012 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: hmmscan_parser.rb,v 1.5 2010/12/13 19:00:11 cmzmasek Exp $
#

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/io/parser/hmmscan_parser'
require 'lib/evo/msa/msa'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/parser/fasta_parser'

module Evoruby

  class HmmscanAnalysis

    PRG_NAME       = "hsa"
    PRG_VERSION    = "1.000"
    PRG_DESC       = "hmmscan analysis"
    PRG_DATE       = "130319"
    COPYRIGHT      = "2013 Christian M Zmasek"
    CONTACT        = "phyloxml@gmail.com"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    I_E_VALUE_THRESHOLD_OPTION    = "ie"
    FS_E_VALUE_THRESHOLD_OPTION   = "pe"
    TARGET_MODELS                 = "m"
    EXTRACTION                    = "x"
    HELP_OPTION_1                 = "help"
    HELP_OPTION_2                 = "h"

    USE_AVOID_HMMS = true
    AVOID_HHMS = [ "RRM_1", "RRM_2", "RRM_3", "RRM_4", "RRM_5", "RRM_6" ]
    LIMIT_FOR_CLOSE_DOMAINS = 20


    def run

      Util.print_program_information( PRG_NAME,
        PRG_VERSION,
        PRG_DESC,
        PRG_DATE,
        COPYRIGHT,
        CONTACT,
        WWW,
        STDOUT )

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

      if ( cla.get_number_of_files != 1 && cla.get_number_of_files != 3 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push( I_E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( FS_E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( TARGET_MODELS )
      allowed_opts.push( EXTRACTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
          "unknown option(s): " + disallowed,
          STDOUT )
      end

      inpath = cla.get_file_name( 0 )
      Util.check_file_for_readability( inpath )
      seq_file_path = nil
      extraction_output = nil

      if ( cla.get_number_of_files == 3 )
        seq_file_path = cla.get_file_name( 1 )
        Util.check_file_for_readability(  seq_file_path  )
        extraction_output = cla.get_file_name( 2 )
        if File.exist?( extraction_output )
          Util.fatal_error( PRG_NAME, "error: [#{extraction_output}] already exists" )
        end

      end

      msa = nil
      if seq_file_path != nil
        msa = read_fasta_file( seq_file_path )
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

      target_models = []
      if ( cla.is_option_set?( TARGET_MODELS ) )
        begin
          hmm_for_protein_output = cla.get_option_value( TARGET_MODELS )
          target_models = hmm_for_protein_output.split( "/" );
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      x_models = []
      if ( cla.is_option_set?( EXTRACTION ) )
        begin
          hmm_for_protein_output = cla.get_option_value( EXTRACTION )
          x_models = hmm_for_protein_output.split( "~" );
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end

      begin
        parse( inpath,
          i_e_value_threshold,
          fs_e_value_threshold,
          target_models,
          x_models,
          msa,
          extraction_output )
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
        i_e_value_threshold,
        fs_e_value_threshold,
        target_models,
        x_models,
        msa,
        extraction_output )

      extraction_output_file = nil
      if extraction_output != nil
        extraction_output_file = File.open( extraction_output, "a" )
      end

      hmmscan_parser = HmmscanParser.new( inpath )

      results = hmmscan_parser.parse

      hmmscan_results_per_protein = []

      query     = ""
      prev_query = ""
      results.each do | r |
        query = r.query
        if !prev_query.empty? && prev_query != query
          if !hmmscan_results_per_protein.empty?
            process_hmmscan_results_per_protein( hmmscan_results_per_protein,
              fs_e_value_threshold,
              target_models,
              x_models,
              i_e_value_threshold,
              msa,
              extraction_output_file )
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

      if !hmmscan_results_per_protein.empty?
        process_hmmscan_results_per_protein( hmmscan_results_per_protein,
          fs_e_value_threshold,
          target_models,
          x_models,
          i_e_value_threshold,
          msa,
          extraction_output_file )
      end
      if extraction_output_file != nil
        extraction_output_file.close
      end
    end # def parse


    def process_hmmscan_results_per_protein( hmmscan_results_per_protein,
        fs_e_value_threshold,
        target_hmms,
        x_models,
        i_e_value_threshold,
        msa,
        extraction_output_file )

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

      query = nil
      qlen = nil
      owns.each do | own |
        if query == nil
          query = own.query
          qlen = own.qlen
        else
          raise StandardError, "failed sanity check" if query != own.query || qlen != own.qlen
          raise StandardError, "failed sanity check: qlen != own.qlen" if qlen != own.qlen
        end
      end

      species = nil
      if msa != nil
        seq = get_sequence( query, msa  )
        species = get_species( seq )
        raise StandardError, "could not get species" if species == nil || species.empty?
        if x_models != nil &&  x_models.length > 0
          extract_linkers( hmmscan_results_per_protein_filtered, x_models, seq,  extraction_output_file )
        end
      end

      s = query + "\t"
      s << species + "\t"
      owns.each do | own |
        s << own.fs_e_value.to_s + "\t"
      end

      s << qlen.to_s + "\t"

      s << hmmscan_results_per_protein_filtered.length.to_s + "\t"
      hmmscan_results_per_protein_filtered.each do | r |
        s << r.model + " "
      end
      s << "\t"
      s <<  make_overview_da( hmmscan_results_per_protein_filtered )
      s << "\t"
      s << make_detailed_da( hmmscan_results_per_protein_filtered, qlen )
      puts s

    end


    def extract_linkers( hmmscan_results_per_protein_filtered, x_models, seq, extraction_output_file )
      raise StandardError, "extraction output file is nil" if extraction_output_file == nil
      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if  prev_r != nil
          if ( x_models.length == 2 && prev_r.model == x_models[ 0 ] && r.model == x_models[ 1 ] )
            extract_linker( prev_r.env_to, r.env_from, seq, extraction_output_file )
          end
        end
        prev_r = r
      end
    end


    def get_sequence( query, msa  )
      seq = nil
      indices = msa.find_by_name_start( query , true )
      if indices.length != 1
        if query[ -1, 1 ] == "|"
          query.chop!
        end
        seq = msa.get_by_name_pattern( /\b#{Regexp.quote(query)}\b/ )
      else
        seq = msa.get_sequence( indices[ 0 ] )
      end
      seq
    end


    def get_species( seq )
      species = nil
      if seq.get_name =~ /\[([A-Z0-9]{3,5})\]/
        species = $1
      end
      species
    end


    def extract_linker( first, last , seq,  output_file )
      if ( last - first >= 1 )
        output_file.print( ">" + seq.get_name + " [" +  first.to_s + "-" + last.to_s +  "]" + "\n")
        output_file.print( seq.get_subsequence( first - 1 , last - 1 ).get_sequence_as_string + "\n" )
      end
    end


    def make_detailed_da( hmmscan_results_per_protein_filtered , qlen )
      s = ""
      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if  prev_r != nil
          s << make_interdomain_sequence( r.env_from - prev_r.env_to - 1 )
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
      s << make_interdomain_sequence( qlen - prev_r.env_to, false )
      s
    end


    def make_overview_da( hmmscan_results_per_protein_filtered )
      overview = ""
      prev_r = nil
      hmmscan_results_per_protein_filtered.each do | r |
        if prev_r == nil
          overview << r.model
        else
          if  ( r.env_from - prev_r.env_to - 1 ) <= LIMIT_FOR_CLOSE_DOMAINS
            overview << "~" << r.model
          else
            overview << "----" << r.model
          end
        end
        prev_r = r

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
      puts( "  options: -" + I_E_VALUE_THRESHOLD_OPTION  + ": i-E-value threshold, default is no threshold" )
      puts( "           -" + FS_E_VALUE_THRESHOLD_OPTION  + ": E-value threshold for full protein sequences, only for protein summary" )
      puts( "           -" + TARGET_MODELS + ": target HMMs" )
      puts()
    end

  end # class

end # module Evoruby
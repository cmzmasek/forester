#
# = lib/evo/apps/domain_sequence_extractor.rb - DomainSequenceExtractor class
#
# Copyright::  Copyright (C) 2012 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id:Exp $


require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/io/parser/hmmscan_domain_extractor'

module Evoruby

  class DomainSequenceExtractor

    PRG_NAME       = "dsx"
    PRG_VERSION    = "2.000"
    PRG_DESC       = "extraction of domain sequences from hmmscan output"
    PRG_DATE       = "20121001"
    COPYRIGHT      = "2012 Christian M Zmasek"
    CONTACT        = "phylosoft@gmail.com"
    WWW            = "www.phylosoft.org"

    E_VALUE_THRESHOLD_OPTION           = 'e'
    LENGTH_THRESHOLD_OPTION            = 'l'
    ADD_POSITION_OPTION                = 'p'
    ADD_DOMAIN_NUMBER_OPTION           = 'd'
    ADD_SPECIES                        = 's'
    MIN_LINKER_OPT                     = 'ml'
    LOG_FILE_SUFFIX                    = '_domain_seq_extr.log'
    PASSED_SEQS_SUFFIX                 = '_with_passing_domains.fasta'
    FAILED_SEQS_SUFFIX                 = '_with_no_passing_domains.fasta'
    HELP_OPTION_1                      = 'help'
    HELP_OPTION_2                      = 'h'

    def run()

      Util.print_program_information( PRG_NAME,
        PRG_VERSION,
        PRG_DESC ,
        PRG_DATE,
        COPYRIGHT,
        CONTACT,
        WWW,
        STDOUT )

      ld = Constants::LINE_DELIMITER

      begin
        cla = CommandLineArguments.new( ARGV )
      rescue ArgumentError
        Util.fatal_error( PRG_NAME, "error: " + $!, STDOUT )
      end

      if ( cla.is_option_set?( HELP_OPTION_1 ) ||
           cla.is_option_set?( HELP_OPTION_2 ) )
        print_help
        exit( 0 )
      end

      if ( cla.get_number_of_files != 4 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push( E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( ADD_POSITION_OPTION )
      allowed_opts.push( ADD_DOMAIN_NUMBER_OPTION )
      allowed_opts.push( LENGTH_THRESHOLD_OPTION )
      allowed_opts.push( ADD_SPECIES )
      allowed_opts.push( MIN_LINKER_OPT )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
          "unknown option(s): " + disallowed,
          STDOUT )
      end

      domain_id           = cla.get_file_name( 0 )
      hmmsearch_output    = cla.get_file_name( 1 )
      fasta_sequence_file = cla.get_file_name( 2 )
      outfile             = cla.get_file_name( 3 )

      if outfile.downcase.end_with?( ".fasta" )
        outfile = outfile[ 0 .. outfile.length - 7 ]
      elsif outfile.downcase.end_with?( ".fsa" )
        outfile = outfile[ 0 .. outfile.length - 5 ]
      end


      add_position = false
      if ( cla.is_option_set?( ADD_POSITION_OPTION ) )
        add_position = true
      end

      add_domain_number           = false
      if ( cla.is_option_set?( ADD_DOMAIN_NUMBER_OPTION ) )
        add_domain_number = true
      end

      add_species = false
      if cla.is_option_set? ADD_SPECIES
        add_species = true
      end

      e_value_threshold = -1.0
      if ( cla.is_option_set?( E_VALUE_THRESHOLD_OPTION ) )
        begin
          e_value_threshold = cla.get_option_value_as_float( E_VALUE_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Forester::Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( e_value_threshold < 0.0 )
          Forester::Util.fatal_error( PRG_NAME, "attempt to use a negative E-value threshold", STDOUT )
        end
      end

      length_threshold = -1
      if ( cla.is_option_set?( LENGTH_THRESHOLD_OPTION ) )
        begin
          length_threshold = cla.get_option_value_as_int( LENGTH_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Forester::Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( length_threshold < 0)
          Forester::Util.fatal_error( PRG_NAME, "attempt to use a negative length threshold", STDOUT )
        end
      end


      min_linker = nil
      if ( cla.is_option_set?( MIN_LINKER_OPT ) )
        begin
          min_linker = cla.get_option_value_as_int( MIN_LINKER_OPT )
        rescue ArgumentError => e
          Forester::Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( !min_linker || min_linker > 100 || min_linker < -100 )
          Forester::Util.fatal_error( PRG_NAME, "unexpected value for min linker " + min_linker.to_s, STDOUT )
        end
      end


      log = String.new

      puts()
      puts( "Domain                                 : " + domain_id )
      log << "Domain                                 : " + domain_id + ld
      puts( "Hmmscan outputfile                     : " + hmmsearch_output )
      log << "Hmmscan outputfile                     : " + hmmsearch_output + ld
      puts( "Fasta sequencefile (complete sequences): " + fasta_sequence_file )
      log << "Fasta sequencefile (complete sequences): " + fasta_sequence_file + ld
      puts( "Outputfile                             : " + outfile + ".fasta" )
      log << "Outputfile                             : " + outfile + ld
      puts( "Passed sequences outfile (fasta)       : " + outfile + PASSED_SEQS_SUFFIX )
      log << "Passed sequences outfile (fasta)       : " + outfile + PASSED_SEQS_SUFFIX + ld
      puts( "Failed sequences outfile (fasta)       : " + outfile + FAILED_SEQS_SUFFIX )
      log << "Failed sequences outfile (fasta)       : " + outfile + FAILED_SEQS_SUFFIX + ld
      puts( "Logfile                                : " + outfile + LOG_FILE_SUFFIX )
      log <<  "Logfile                                : " + outfile + LOG_FILE_SUFFIX + ld
      if ( e_value_threshold >= 0.0 )
        puts( "E-value threshold : " + e_value_threshold.to_s )
        log << "E-value threshold : " + e_value_threshold.to_s + ld
      else
        puts( "E-value threshold : no threshold" )
        log << "E-value threshold : no threshold" + ld
      end
      if ( length_threshold > 0 )
        puts( "Length threshold  : " + length_threshold.to_s )
        log << "Length threshold  : " + length_threshold.to_s + ld
      else
        puts( "Length threshold  : no threshold" )
        log << "Length threshold  : no threshold" + ld
      end
      if ( min_linker )
        puts( "Min linker        : " + min_linker.to_s )
        log << "Min linker        :  " + min_linker.to_s +  ld

      end


      if ( add_position )
        puts( "Add positions (rel to complete seq) to extracted domains: true" )
        log << "Add positions (rel to complete seq) to extracted domains: true" + ld
      else
        puts( "Add positions (rel to complete seq) to extracted domains: false" )
        log << "Add positions (rel to complete seq) to extracted domains: false" + ld
      end

      if ( add_domain_number )
        puts( "Add numbers to extracted domains (in case of more than one domain per complete seq): true" )
        log << "Add numbers to extracted domains (in case of more than one domain per complete seq): true" + ld
      else
        puts( "Add numbers to extracted domains (in case of more than one domain per complete seq): false" )
        log << "Add numbers to extracted domains (in case of more than one domain per complete seq): false" + ld
      end

      puts

      domain_count = 0
      begin
        parser = HmmscanDomainExtractor.new()
        domain_count = parser.parse( domain_id,
          hmmsearch_output,
          fasta_sequence_file,
          outfile,
          outfile + PASSED_SEQS_SUFFIX,
          outfile + FAILED_SEQS_SUFFIX,
          e_value_threshold,
          length_threshold,
          add_position,
          add_domain_number,
          add_species,
          min_linker,
          log )
      rescue ArgumentError, IOError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )

      rescue Exception => e
        puts e.backtrace
        Util.fatal_error( PRG_NAME, "unexpected exception: " + e.to_s, STDOUT )

      end

      puts
      Util.print_message( PRG_NAME, "extracted a total of " + domain_count.to_s + " domains" )
      Util.print_message( PRG_NAME, "wrote;               " + outfile + ".fasta")
      Util.print_message( PRG_NAME, "wrote:               " + outfile + LOG_FILE_SUFFIX )
      Util.print_message( PRG_NAME, "wrote:               " + outfile + PASSED_SEQS_SUFFIX )
      Util.print_message( PRG_NAME, "wrote:               " + outfile + FAILED_SEQS_SUFFIX )

      begin
        f = File.open( outfile + LOG_FILE_SUFFIX, 'a' )
        f.print( log )
        f.close
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      puts
      Util.print_message( PRG_NAME, "OK" )
      puts

    end

    def print_help()
      puts()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb [options] <domain> <hmmscan outputfile> <file containing complete sequences in fasta format> <outputfile>" )
      puts()
      puts( "  options: -" + E_VALUE_THRESHOLD_OPTION  + "=<f>: iE-value threshold, default is no threshold" )
      puts( "           -" + LENGTH_THRESHOLD_OPTION   + "=<i>: length threshold, default is no threshold" )
      puts( "           -" + ADD_POSITION_OPTION  + ": to add positions (rel to complete seq) to extracted domains" )
      puts( "           -" + ADD_DOMAIN_NUMBER_OPTION  + ": to add numbers to extracted domains (in case of more than one domain per complete seq) (example \"domain~2-3\")" )
      puts( "           -" + ADD_SPECIES  + ": to add species [in brackets]" )
      puts( "           -" + MIN_LINKER_OPT  + "=<i>: to extract pairs of same domains with a distance inbetween shorter than a given value" )
      puts()
    end

  end # class DomainSequenceExtractor

end # module Evoruby

#
# = lib/evo/tool/multi_domain_seq_extractor - MultiDomainSeqExtractor class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/io/parser/hmmscan_multi_domain_extractor'

module Evoruby
  class MultiDomainSeqExtractor

    PRG_NAME       = "mdsx"
    PRG_VERSION    = "1.000"
    PRG_DESC       = "Extraction of multi domain sequences from hmmscan output"
    PRG_DATE       = "20170220"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    E_VALUE_THRESHOLD_DEFAULT = 0.1
    LENGTH_THRESHOLD_DEFAULT  = 50

    E_VALUE_THRESHOLD_OPTION           = 'e'
    LENGTH_THRESHOLD_OPTION            = 'l'
    ADD_POSITION_OPTION                = 'p'
    ADD_DOMAIN_NUMBER_OPTION           = 'd'
    ADD_SPECIES                        = 's'
    LOG_FILE_SUFFIX                    = '_multi_domain_seq_extr.log'
    PASSED_SEQS_SUFFIX                 = '_with_passing_domains.fasta'
    FAILED_SEQS_SUFFIX                 = '_with_no_passing_domains.fasta'
    HELP_OPTION_1                      = 'help'
    HELP_OPTION_2                      = 'h'
    def run()

      Util.print_program_information( PRG_NAME,
      PRG_VERSION,
      PRG_DESC ,
      PRG_DATE,
      WWW,
      STDOUT )

      if ( ARGV == nil || ( ARGV.length < 1 )  )
        print_help
        exit( -1 )
      end

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

      unless ( cla.get_number_of_files == 2 || cla.get_number_of_files == 3 || cla.get_number_of_files == 4 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push( E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( ADD_POSITION_OPTION )
      allowed_opts.push( ADD_DOMAIN_NUMBER_OPTION )
      allowed_opts.push( LENGTH_THRESHOLD_OPTION )
      allowed_opts.push( ADD_SPECIES )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
        "unknown option(s): " + disallowed,
        STDOUT )
      end

      e_value_threshold = E_VALUE_THRESHOLD_DEFAULT
      if ( cla.is_option_set?( E_VALUE_THRESHOLD_OPTION ) )
        begin
          e_value_threshold = cla.get_option_value_as_float( E_VALUE_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( e_value_threshold < 0.0 )
          Util.fatal_error( PRG_NAME, "attempt to use a negative E-value threshold", STDOUT )
        end
      end

      length_threshold = LENGTH_THRESHOLD_DEFAULT
      if ( cla.is_option_set?( LENGTH_THRESHOLD_OPTION ) )
        begin
          length_threshold = cla.get_option_value_as_int( LENGTH_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( length_threshold < 0)
          Util.fatal_error( PRG_NAME, "attempt to use a negative length threshold", STDOUT )
        end
      end

      domain_id           = cla.get_file_name( 0 )
      hmmscan_output      = cla.get_file_name( 1 )
      fasta_sequence_file = ""
      outfile             = ""

      if (cla.get_number_of_files == 4)
        fasta_sequence_file = cla.get_file_name( 2 )
        outfile             = cla.get_file_name( 3 )
      elsif (cla.get_number_of_files == 2 || cla.get_number_of_files == 3)
        if (cla.get_number_of_files == 3)
          fasta_sequence_file = cla.get_file_name( 2 )
        else
          hmmscan_index = hmmscan_output.index(Constants::HMMSCAN)
          if ( hmmscan_index != nil )
            prefix = hmmscan_output[0 .. hmmscan_index-1 ]
            suffix = Constants::ID_NORMALIZED_FASTA_FILE_SUFFIX
            files = Dir.entries( "." )
            matching_files = Util.get_matching_files( files, prefix, suffix)
            if matching_files.length < 1
              Util.fatal_error( PRG_NAME, 'no file matching [' + prefix +
              '...' + suffix + '] present in current directory: need to indicate <file containing complete sequences in fasta format> as second argument' )
            end
            if matching_files.length > 1
              Util.fatal_error( PRG_NAME, 'more than one file matching [' +
              prefix  + '...' + suffix + '] present in current directory: need to indicate <file containing complete sequences in fasta format> as second argument' )
            end
            fasta_sequence_file = matching_files[ 0 ]
          else
            Util.fatal_error( PRG_NAME, 'input files do not seem in format for standard analysis pipeline, need to explicitly indicate all' )
          end
        end
        hmmscan_index = hmmscan_output.index(Constants::HMMSCAN)
        if ( hmmscan_index != nil )
          outfile = hmmscan_output.sub(Constants::HMMSCAN, "_") + "_" + domain_id
          e = e_value_threshold >= 1 ? e_value_threshold.to_i : e_value_threshold.to_s.sub!('.','_')
          outfile += "_E" + e.to_s
          outfile += "_L" + length_threshold.to_i.to_s
        else
          Util.fatal_error( PRG_NAME, 'input files do not seem in format for standard analysis pipeline, need to explicitly indicate all' )
        end
      end

      if outfile.downcase.end_with?( ".fasta" )
        outfile = outfile[ 0 .. outfile.length - 7 ]
      elsif outfile.downcase.end_with?( ".fsa" )
        outfile = outfile[ 0 .. outfile.length - 5 ]
      end

      add_position = false
      if ( cla.is_option_set?( ADD_POSITION_OPTION ) )
        add_position = true
      end

      add_domain_number = false
      if ( cla.is_option_set?( ADD_DOMAIN_NUMBER_OPTION ) )
        add_domain_number = true
      end

      add_species = false
      if cla.is_option_set?( ADD_SPECIES )
        add_species = true
      end

      log = String.new
      ld = Constants::LINE_DELIMITER

      puts()
      puts( "Domain                                                                             : " + domain_id )
      log << "Domain                                                                             : " + domain_id + ld
      puts( "Hmmscan outputfile                                                                 : " + hmmscan_output )
      log << "Hmmscan outputfile                                                                 : " + hmmscan_output + ld
      puts( "Fasta sequencefile (complete sequences)                                            : " + fasta_sequence_file )
      log << "Fasta sequencefile (complete sequences)                                            : " + fasta_sequence_file + ld
      puts( "Outputfile                                                                         : " + outfile + ".fasta" )
      log << "Outputfile                                                                         : " + outfile + ld
      puts( "Passed sequences outfile (fasta)                                                   : " + outfile + PASSED_SEQS_SUFFIX )
      log << "Passed sequences outfile (fasta)                                                   : " + outfile + PASSED_SEQS_SUFFIX + ld
      puts( "Failed sequences outfile (fasta)                                                   : " + outfile + FAILED_SEQS_SUFFIX )
      log << "Failed sequences outfile (fasta)                                                   : " + outfile + FAILED_SEQS_SUFFIX + ld
      puts( "Logfile                                                                            : " + outfile + LOG_FILE_SUFFIX )
      log << "Logfile                                                                            : " + outfile + LOG_FILE_SUFFIX + ld
      if ( e_value_threshold >= 0.0 )
        puts( "iE-value threshold                                                                 : " + e_value_threshold.to_s )
        log << "iE-value threshold                                                                 : " + e_value_threshold.to_s + ld
      else
        puts( "iE-value threshold                                                                 : no threshold" )
        log << "iE-value threshold                                                                 : no threshold" + ld
      end
      if ( length_threshold > 0 )
        puts( "Length threshold (env)                                                             : " + length_threshold.to_s )
        log << "Length threshold (env)                                                             : " + length_threshold.to_s + ld
      else
        puts( "Length threshold (env)                                                             : no threshold" )
        log << "Length threshold  (env)                                                            : no threshold" + ld
      end
      if ( add_position )
        puts( "Add positions (rel to complete seq) to extracted domains                           : true" )
        log << "Add positions (rel to complete seq) to extracted domains                           : true" + ld
      else
        puts( "Add positions (rel to complete seq) to extracted domains                           : false" )
        log << "Add positions (rel to complete seq) to extracted domains                           : false" + ld
      end

      if ( add_domain_number )
        puts( "Add numbers to extracted domains (in case of more than one domain per complete seq): true" )
        log << "Add numbers to extracted domains (in case of more than one domain per complete seq): true" + ld
      else
        puts( "Add numbers to extracted domains (in case of more than one domain per complete seq): false" )
        log << "Add numbers to extracted domains (in case of more than one domain per complete seq): false" + ld
      end

      puts
      log <<  ld

      domain_count = 0
      begin
        parser = HmmscanMultiDomainExtractor.new()
        domain_count = parser.parse( domain_id,
        hmmscan_output,
        fasta_sequence_file,
        outfile,
        outfile + PASSED_SEQS_SUFFIX,
        outfile + FAILED_SEQS_SUFFIX,
        e_value_threshold,
        length_threshold,
        add_position,
        add_domain_number,
        add_species,
        log )
      rescue ArgumentError, IOError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )

      rescue Exception => e
        puts e.backtrace
        Util.fatal_error( PRG_NAME, "unexpected exception: " + e.to_s, STDOUT )

      end

      puts
      Util.print_message( PRG_NAME, "extracted a total of " + domain_count.to_s + " domains" )
      Util.print_message( PRG_NAME, "wrote: " + outfile + ".fasta")
      Util.print_message( PRG_NAME, "wrote: " + outfile + LOG_FILE_SUFFIX )
      Util.print_message( PRG_NAME, "wrote: " + outfile + PASSED_SEQS_SUFFIX )
      Util.print_message( PRG_NAME, "wrote: " + outfile + FAILED_SEQS_SUFFIX )

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
      puts( "  " + PRG_NAME + ".rb [options] <target domain> <hmmscan outputfile> [file containing complete sequences in fasta format] [outputfile]" )
      puts()
      puts( "  options: -" + E_VALUE_THRESHOLD_OPTION  + "=<f> : iE-value threshold for target domain, default is " + E_VALUE_THRESHOLD_DEFAULT.to_s )
      puts( "           -" + LENGTH_THRESHOLD_OPTION   + "=<i> : length threshold target domain (env), default is " + LENGTH_THRESHOLD_DEFAULT.to_s )
      puts( "           -" + ADD_DOMAIN_NUMBER_OPTION  + "     : to add numbers to extracted domains (in case of more than one domain per complete seq) (example \"domain~2-3\")" )
      puts( "           -" + ADD_POSITION_OPTION  + "     : to add positions (rel to complete seq) to extracted domains" )
      puts( "           -" + ADD_SPECIES  + "     : to add species [in brackets]" )
      puts()
      puts( "Examples:" )
      puts
      puts( "  " + PRG_NAME + ".rb -d -e=1e-6 -l=50 Pkinase P53_hmmscan_#{Constants::PFAM_V_FOR_EX}_10 P53_ni.fasta P53_#{Constants::PFAM_V_FOR_EX}_10_Pkinase_E1_0e-06_L50" )
      puts
      puts( "  " + PRG_NAME + ".rb -d -e=1e-6 -l=50 Pkinase P53_hmmscan_#{Constants::PFAM_V_FOR_EX}_10 P53_ni.fasta" )
      puts
      puts( "  " + PRG_NAME + ".rb -d -e=1e-6 -l=50 Pkinase P53_hmmscan_#{Constants::PFAM_V_FOR_EX}_10" )
      puts()
    end

  end # class DomainSequenceExtractor

end # module Evoruby

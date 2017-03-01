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

    LOG_FILE_SUFFIX                    = '_MDSX.log'
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

      unless ( cla.get_number_of_files == 2 || cla.get_number_of_files == 3 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
        "unknown option(s): " + disallowed,
        STDOUT )
      end

      domain_id           = cla.get_file_name( 0 )
      hmmscan_output      = cla.get_file_name( 1 )
      fasta_sequence_file = ""
      outfile             = ""

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
        outfile = hmmscan_output.sub(Constants::HMMSCAN, "_") + "_MDSX"
      else
        Util.fatal_error( PRG_NAME, 'input files do not seem in format for standard analysis pipeline, need to explicitly indicate all' )
      end

      log = String.new
      ld = Constants::LINE_DELIMITER

      #      puts()
      #      puts( "Domain                                                                             : " + domain_id )
      #      log << "Domain                                                                             : " + domain_id + ld
      #      puts( "Hmmscan outputfile                                                                 : " + hmmscan_output )
      #      log << "Hmmscan outputfile                                                                 : " + hmmscan_output + ld
      #      puts( "Fasta sequencefile (complete sequences)                                            : " + fasta_sequence_file )
      #      log << "Fasta sequencefile (complete sequences)                                            : " + fasta_sequence_file + ld
      #      puts( "Outputfile                                                                         : " + outfile + ".fasta" )
      #      log << "Outputfile                                                                         : " + outfile + ld
      #      puts( "Passed sequences outfile (fasta)                                                   : " + outfile + PASSED_SEQS_SUFFIX )
      #      log << "Passed sequences outfile (fasta)                                                   : " + outfile + PASSED_SEQS_SUFFIX + ld
      #      puts( "Failed sequences outfile (fasta)                                                   : " + outfile + FAILED_SEQS_SUFFIX )
      #      log << "Failed sequences outfile (fasta)                                                   : " + outfile + FAILED_SEQS_SUFFIX + ld
      #      puts( "Logfile                                                                            : " + outfile + LOG_FILE_SUFFIX )
      #      log << "Logfile                                                                            : " + outfile + LOG_FILE_SUFFIX + ld

      puts
      log <<  ld

      domain_count = 0
      begin
        parser = HmmscanMultiDomainExtractor.new()
        domain_count = parser.parse( domain_id,
        hmmscan_output,
        fasta_sequence_file,
        outfile,
        log )
      rescue ArgumentError, IOError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )

      rescue Exception => e
        puts e.backtrace
        Util.fatal_error( PRG_NAME, "unexpected exception: " + e.to_s, STDOUT )

      end

      puts
      Util.print_message( PRG_NAME, "extracted a total of " + domain_count.to_s + " domains" )
      #  Util.print_message( PRG_NAME, "wrote: " + outfile + ".fasta")
      #  Util.print_message( PRG_NAME, "wrote: " + outfile + LOG_FILE_SUFFIX )
      #  Util.print_message( PRG_NAME, "wrote: " + outfile + PASSED_SEQS_SUFFIX )
      #  Util.print_message( PRG_NAME, "wrote: " + outfile + FAILED_SEQS_SUFFIX )

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
      puts( "  " + PRG_NAME + ".rb <da> <hmmscan outputfile> [file containing complete sequences in fasta format]" )
      puts()
      puts( "  options: -"  )
      puts()
      puts( "Examples:" )
      puts
      puts( "  " + PRG_NAME + ".rb " )
      puts

      puts()
    end

  end # class DomainSequenceExtractor

end # module Evoruby

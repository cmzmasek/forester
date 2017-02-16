#
# = lib/evo/apps/fasta_extractor.rb - FastaExtractor class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'

module Evoruby
  class FastaExtractor

    PRG_NAME       = "fastx"
    PRG_VERSION    = "1.000"
    PRG_DESC       = "extraction of molecular sequences from a fasta file"
    PRG_DATE       = "170215"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

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
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      if ( cla.is_option_set?( HELP_OPTION_1 ) ||
      cla.is_option_set?( HELP_OPTION_2 ) )
        print_help
        exit( 0 )
      end

      if ( cla.get_number_of_files != 3 )
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

      input_file  = cla.get_file_name( 0 )
      query       = cla.get_file_name( 1 )
      output_file = cla.get_file_name( 2 )

      if  !File.exist?( input_file )
        Util.fatal_error( PRG_NAME, "error: input file [#{input_file}] does not exist" )
      end
      if File.exist?( output_file   )
        Util.fatal_error( PRG_NAME, "error: [#{output_file}] already exists" )
      end

      results = extract_sequences( query, input_file, output_file )

      Util.print_message( PRG_NAME, "matched: " + results )
      Util.print_message( PRG_NAME, "wrote:   " + output_file )
      Util.print_message( PRG_NAME, "OK" )

    end

    def extract_sequences( query, fasta_file, output_file )
      output = File.open( output_file, "a" )
      matching_state = false
      matches = 0
      total = 0
      File.open( fasta_file ) do | file |
        while line = file.gets
          if !Util.is_string_empty?( line )
            if line =~ /^\s*>/
              total += 1
              if total % 10000 == 0
                STDOUT.write "\r#{matches}/#{total}"
                STDOUT.flush
              end
              if line =~ /#{query}/
                matching_state = true
                matches += 1
                output.print( line )
              else
                matching_state = false
              end
            elsif matching_state
              output.print( line )
            end
          end
        end
      end
      output.close()
      matches.to_s + "/" + total.to_s
    end

    def print_help()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb <input fasta file> <query> <output file>" )
      puts()
      puts( "Examples:" )
      puts
      puts( "  " + PRG_NAME + ".rb Pfam-A.fasta kinase kinases" )
      puts()
    end

  end # class FastaExtractor
end
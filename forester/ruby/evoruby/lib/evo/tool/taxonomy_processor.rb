#
# = lib/evo/apps/taxonomy_processor - TaxonomyProcessor class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/msa/msa'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/parser/general_msa_parser'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/writer/phylip_sequential_writer'
require 'lib/evo/util/command_line_arguments'

module Evoruby
  class TaxonomyProcessor

    PRG_NAME       = "tap"
    PRG_DATE       = "170214"
    PRG_DESC       = "Replacement of labels in multiple sequence files"
    PRG_VERSION    = "2.004"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    EXTRACT_TAXONOMY_OPTION = "t"
    ANNOTATION_OPTION       = "a"
    HELP_OPTION_1           = "help"
    HELP_OPTION_2           = "h"
    def run()

      Util.print_program_information( PRG_NAME,
      PRG_VERSION,
      PRG_DESC,
      PRG_DATE,
      WWW,
      STDOUT )

      if ( ARGV == nil || ( ARGV.length < 1 ) )
        print_help()
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

      input      = nil
      output     = nil
      list_file  = nil

      if cla.get_number_of_files == 3
        input     = cla.get_file_name( 0 )
        output    = cla.get_file_name( 1 )
        list_file = cla.get_file_name( 2 )
      elsif cla.get_number_of_files == 1
        input     = cla.get_file_name( 0 )
        i = nil
        if input.downcase.end_with?( ".fasta" )
          i = input[ 0 .. input.length - 7 ]
        elsif input.downcase.end_with?( ".fsa" )
          i = input[ 0 .. input.length - 5 ]
        else
          i = input
        end
        output    = i + Constants::ID_NORMALIZED_FASTA_FILE_SUFFIX
        list_file = i + Constants::ID_MAP_FILE_SUFFIX
      else
        print_help()
        exit(-1)
      end

      allowed_opts = Array.new
      allowed_opts.push( EXTRACT_TAXONOMY_OPTION )
      allowed_opts.push( ANNOTATION_OPTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME, "unknown option(s): " + disallowed )
      end

      extract_taxonomy = false
      if ( cla.is_option_set?( EXTRACT_TAXONOMY_OPTION ) )
        extract_taxonomy = true
      end

      annotation = nil
      if ( cla.is_option_set?( ANNOTATION_OPTION ) )
        annotation = cla.get_option_value( ANNOTATION_OPTION )
      end

      if ( File.exist?( output ) )
        Util.fatal_error( PRG_NAME, "outfile [" + output + "] already exists" )
      end
      if ( File.exist?( list_file ) )
        Util.fatal_error( PRG_NAME, "list file [" + list_file + "] already exists" )
      end
      if ( !File.exist?( input) )
        Util.fatal_error( PRG_NAME, "infile [" + input + "] does not exist" )
      end

      fasta_like = Util.looks_like_fasta?( input )

      puts()
      puts( "Input alignment : " + input )
      puts( "Output alignment: " + output )
      puts( "Name list       : " + list_file )
      if ( fasta_like )
        puts( "Format          : Fasta"  )
      else
        puts( "Format          : Phylip like" )
      end
      if ( extract_taxonomy )
        puts( "Extract taxonomy: true"  )
      end
      if ( annotation != nil )
        puts( "Annotation      : " + annotation )
      end
      puts()

      f = MsaFactory.new()
      begin
        if ( fasta_like )
          msa = f.create_msa_from_file( input, FastaParser.new() )
        else
          msa = f.create_msa_from_file( input, GeneralMsaParser.new() )
        end
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "failed to read file: " + e.to_s )
      end

      if ( msa == nil || msa.get_number_of_seqs() < 1 )
        Util.fatal_error( PRG_NAME, "failed to read MSA" )
      end
      begin
        Util.check_file_for_writability( list_file )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_, STDOUT )
      end

      lf = File.open( list_file, "a" )
      for i in 0 ... msa.get_number_of_seqs
        seq = msa.get_sequence( i )
        seq.set_name( modify_name( seq.get_name(), i, lf, extract_taxonomy, annotation ) )
      end
      io = MsaIO.new()
      w = nil
      if ( fasta_like )
        w = FastaWriter.new()
      else
        w = PhylipSequentialWriter.new()
      end
      w.set_max_name_length( 9 )
      w.clean( true )
      begin
        io.write_to_file( msa, output, w )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "failed to write file: " + e.to_s )
      end
      lf.close()
      Util.print_message( PRG_NAME, "wrote: " + list_file )
      Util.print_message( PRG_NAME, "wrote: " + output )
      Util.print_message( PRG_NAME, "next steps in standard analysis pipeline: hmmscan followed by hsp.rb")
      Util.print_message( PRG_NAME, "OK" )
    end

    private

    def modify_name( desc, counter, file, extract_taxonomy, annotation )
      new_desc = nil
      desc.gsub!( /\s+/, ' ' )
      if extract_taxonomy
        if desc =~/\s\[(([A-Z9][A-Z]{2}[A-Z0-9]{2})|RAT|PIG|PEA|CAP)\]/
          new_desc = counter.to_s( 16 ) + "_" + $1
        else
          Util.fatal_error( PRG_NAME, "could not get taxonomy from: " + desc )
        end
      else
        new_desc = counter.to_s( 16 )
      end
      if (annotation != nil)
        new_desc = new_desc + annotation
        file.print( new_desc + "\t" + desc + " " + annotation + "\n" )
      else
        file.print( new_desc + "\t" + desc + "\n" )
      end
      if ( new_desc.length > 9)
        Util.fatal_error( PRG_NAME, "shortened identifier [" +
        new_desc + "] is too long (" + new_desc.length.to_s + " characters)" )
      end
      new_desc
    end

    def print_help()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb [options] <input sequences> [output sequences] [output id list]" )
      puts()
      puts( "  options: -" + EXTRACT_TAXONOMY_OPTION + "    : to extract taxonomy information from bracketed expressions" )
      puts( "           -" + ANNOTATION_OPTION + "=<s>: to add an annotation to all entries" )
      puts()
      puts( "  [next steps in standard analysis pipeline: hmmscan followed by hsp.rb]")
      puts()
      puts( "Example:" )
      puts()
      puts( "  " + PRG_NAME + ".rb P53.fasta" )
      puts()
    end

  end # class TaxonomyProcessor

end # module Evoruby

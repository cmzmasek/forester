#
# = lib/evo/apps/taxonomy_processor - TaxonomyProcessor class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: taxonomy_processor.rb,v 1.26 2010/12/13 19:00:11 cmzmasek Exp $


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
    PRG_DATE       = "2012.09.27"
    PRG_DESC       = "replacement of species names in multiple sequence files"
    PRG_VERSION    = "1.02"
    COPYRIGHT      = "2012 Christian M Zmasek"
    CONTACT        = "phylosoft@gmail.com"
    WWW            = "www.phylosoft.org"

    SIMPLE = true

    EXTRACT_TAXONOMY_OPTION = "t"

    def initialize()
      @taxonomies = Hash.new()
    end

    def run()

      Util.print_program_information( PRG_NAME,
        PRG_VERSION,
        PRG_DESC,
        PRG_DATE,
        COPYRIGHT,
        CONTACT,
        WWW,
        STDOUT )

      if ( ARGV == nil || ( ARGV.length != 1 && ARGV.length != 2 && ARGV.length != 3 && ARGV.length != 4 && ARGV.length != 5 && ARGV.length != 6 ) )
        puts( "Usage: #{PRG_NAME}.rb [options] [input map file] <input sequences> [output sequences] [output id list]" )
        puts()
        puts( "  options: -" + EXTRACT_TAXONOMY_OPTION + ": to extract taxonomy information from bracketed expression" )
        puts()
        exit( -1 )
      end

      begin
        cla = CommandLineArguments.new( ARGV )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end


      mapfile   = nil
      input     = nil
      output    = nil
      list_file = nil



      if cla.get_number_of_files == 4
        mapfile   = cla.get_file_name( 0 )
        input     = cla.get_file_name( 1 )
        output    = cla.get_file_name( 2 )
        list_file = cla.get_file_name( 3 )
      elsif cla.get_number_of_files == 3
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
        output    = i + "_ni.fasta"
        list_file = i + ".nim"
      end


      allowed_opts = Array.new
      allowed_opts.push( EXTRACT_TAXONOMY_OPTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME, "unknown option(s): " + disallowed )
      end

      extract_taxonomy = false
      if ( cla.is_option_set?( EXTRACT_TAXONOMY_OPTION ) )
        extract_taxonomy = true
      end

      if ( File.exists?( output ) )
        Util.fatal_error( PRG_NAME, "outfile [" + output + "] already exists" )
      end
      if ( File.exists?( list_file ) )
        Util.fatal_error( PRG_NAME, "list file [" + list_file + "] already exists" )
      end
      if ( !File.exists?( input) )
        Util.fatal_error( PRG_NAME, "infile [" + input + "] does not exist" )
      end
      if ( mapfile != nil && !File.exists?( mapfile ) )
        Util.fatal_error( PRG_NAME, "mapfile [" + mapfile + "] does not exist" )
      end

      fasta_like = Util.looks_like_fasta?( input )

      puts()
      if mapfile != nil
        puts( "Map file        : " + mapfile )
      end
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
      puts()

      species_map = Hash.new
      if mapfile != nil
        File.open( mapfile ) do | file |
          while line = file.gets
            if ( line =~/(.+)#(.+)/ || line =~/(.+)\s+(.+)/ )
              species_map[ $1 ] = $2
              Util.print_message( PRG_NAME, "mapping: " + $1 + ' => ' + $2 )
            end
          end
        end
      end

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

      #removed = msa.remove_redundant_sequences!( true )
      #if removed.size > 0
      #  Util.print_message( PRG_NAME, "going to ignore the following " + removed.size.to_s + " redundant sequences:" )
      #  removed.each { | seq_name |
      #    puts seq_name
      #  }
      #  Util.print_message( PRG_NAME, "will process " + msa.get_number_of_seqs.to_s + " non redundant sequences" )
      #end

      lf = File.open( list_file, "a" )
      for i in 0 ... msa.get_number_of_seqs
        seq  = msa.get_sequence( i )
        seq.set_name( Util::normalize_seq_name( modify_name( seq.get_name(), i, lf, species_map, extract_taxonomy ), 10 ) )
      end

      io = MsaIO.new()
      w = nil
      if ( fasta_like )
        w = FastaWriter.new()
      else
        w = PhylipSequentialWriter.new()
      end
      w.set_max_name_length( 10 )
      w.clean( true )
      begin
        io.write_to_file( msa, output, w )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "failed to write file: " + e.to_s )
      end
      lf.close()
      if ( @taxonomies.length > 0 )
        Util.print_message( PRG_NAME, "number of unique taxonomies: " + @taxonomies.length.to_s )
      end
      Util.print_message( PRG_NAME, "wrote: " + list_file )
      Util.print_message( PRG_NAME, "wrote: " + output )
      Util.print_message( PRG_NAME, "OK" )
    end

    private

    def modify_name( desc, counter, file, species_map, extract_taxonomy )
      new_desc = nil
      my_species = nil
      # if desc =~ /^>?\s*\S{1,10}_([0-9A-Z]{3,5})/
      if desc =~ /^>?\s*\S{1,10}_([A-Z]{3,5})/
        new_desc = counter.to_s( 16 ) + "_" + $1
      elsif SIMPLE
        new_desc = counter.to_s( 16 )
      elsif extract_taxonomy
        if ( desc.count( "[" ) != desc.count( "]" ) )
          Util.fatal_error( PRG_NAME, "illegal bracket count in: " + desc )
        end
        species = nil
        species_map.each_key do | key |
          if desc =~ /[\b|_]#{key}\b/  # Added boundaries to prevent e.g. RAT matching ARATH.
            species = species_map[ key ]
            new_desc = counter.to_s( 16 ) + "_" + species
            break
          end
        end
        if species == nil
          if desc =~/.*\[(\S{3,}?)\]/
            species = $1
            species.strip!
            species.upcase!
            species.gsub!( /\s+/, " " )
            species.gsub!( /-/, "" )
            species.gsub!( /\)/, "" )
            species.gsub!( /\(/, "" )
            species.gsub!( /\'/, "" )
            if species =~ /\S+\s\S+/ || species =~ /\S{3,5}/
              if species =~ /(\S+)\s(\S+)/
                code = $1[ 0..2 ] + $2[ 0..1 ]
              elsif  species =~ /\S{3,5}/
                code = species
              elsif species.count( " " ) > 2
                species =~ /(\S+)\s+(\S+)\s+(\S+)$/
                third_last = $1
                second_last = $2
                last = $3
                code = code[ 0 ] + third_last[ 0 ] + second_last[ 0 ] + last[ 0 ] + last[ last.size - 1 ]
              elsif species.count( " " ) > 1
                species =~ /(\S+)\s+(\S+)$/
                second_last = $1
                last = $2
                code = code[ 0..1 ] + second_last[ 0 ] + last[ 0 ] + last[ last.size - 1 ]
              end
              new_desc = counter.to_s( 16 ) + "_" + code
              if @taxonomies.has_key?( code )
                if ( !@taxonomies.has_value?( species ) )
                  Util.fatal_error( PRG_NAME, "code [#{code}] is not unique in [#{desc}]" )
                end
              else
                if ( @taxonomies.has_value?( species ) )
                  Util.fatal_error( PRG_NAME, "genome [#{species}] is not unique in [#{desc}]" )
                else
                  @taxonomies[ code ] = species
                end
              end
            else
              Util.fatal_error( PRG_NAME, "illegal format [#{species}] in: " + desc )
            end
          else
            Util.fatal_error( PRG_NAME, "illegal format in: " + desc )
          end
        end
      else
        species = nil
        my_species = nil
        species_map.each_key do | key |
          if desc =~ /#{key}/
            species = species_map[ key ]
            species = species.gsub( /\s+/, "" )
            species = species.gsub( /_/, " " )
            my_species = species
            if species =~ /(\S+)\s+(\S+)/
              species = $1[0..2] + $2[0..1]
            end
            species = species.gsub( /\s+/, "" )
            species = species.slice(0, 5)
            species.upcase!
            break
          end
        end
        if species == nil
          Util.fatal_error( PRG_NAME, "species not found in: " + desc  )
        else
          new_desc = counter.to_s( 16 ) + "_" + species
        end
      end
      if new_desc == nil
        Util.fatal_error( PRG_NAME, "failed to extract species from: " + desc  )
      end
      if my_species != nil
        file.print( new_desc + ": " + desc + " [" + my_species + "]" + "\n" )
      else
        file.print( new_desc + ": " + desc + "\n" )
      end
      new_desc
    end

  end # class TaxonomyProcessor

end # module Evoruby

#
# = lib/evo/apps/multi_sequence_extractor.rb - MultiSequenceExtractor class
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: multi_sequence_extractor.rb,v 1.10 2010/12/13 19:00:11 cmzmasek Exp $


require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/msa/msa'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/util/command_line_arguments'



module Evoruby

  class MultiSequenceExtractor

    PRG_NAME                           = "mse"
    PRG_VERSION                        = "1.03"
    PRG_DESC                           = "extraction of sequences by name from multiple multi-sequence ('fasta') files"
    PRG_DATE                           = "131127"
    COPYRIGHT                          = "2008-2013 Christian M Zmasek"
    CONTACT                            = "phylosoft@gmail.com"
    WWW                                = "https://sites.google.com/site/cmzmasek/home/software/forester"
    HELP_OPTION_1                      = 'help'
    HELP_OPTION_2                      = 'h'

    EXT_OPTION                          = 'e'
    EXTRACT_LINKERS_OPTION              = 'l'
    LOG_SUFFIX                          = ".mse_log"
    LINKERS_SUFFIX                      = ".linkers"
    FASTA_SUFFIX                        = ".fasta"
    FASTA_WITH_NORMALIZED_IDS_SUFFIX    = ".ni.fasta"
    NORMALIZED_IDS_MAP_SUFFIX           = ".nim"
    PROTEINS_LIST_FILE_SEPARATOR        = "\t"

    def initialize()
      @file_to_msa = Hash.new
      @seqs = 0
    end

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
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      if ( cla.is_option_set?( HELP_OPTION_1 ) ||
           cla.is_option_set?( HELP_OPTION_2 ) )
        print_help
        exit( 0 )
      end

      if ( cla.get_number_of_files != 4 && cla.get_number_of_files != 5 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push(EXT_OPTION)
      allowed_opts.push(EXTRACT_LINKERS_OPTION)

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
          "unknown option(s): " + disallowed,
          STDOUT )
      end

      seq_names_files_suffix = cla.get_file_name( 0 )
      input_dir              = cla.get_file_name( 1 )
      out_dir                = cla.get_file_name( 2 )
      out_dir_doms           = cla.get_file_name( 3 )
      mapping_file            = nil

      if ( cla.get_number_of_files == 5 )
        mapping_file = cla.get_file_name( 4 )
        begin
          Util.check_file_for_readability( mapping_file )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s )
        end
      end

      extension = 0
      if cla.is_option_set?(EXT_OPTION)
        extension = cla.get_option_value_as_int(EXT_OPTION)
        if extension < 0
          extension = 0
        end
      end

      extract_linkers = false
      if cla.is_option_set?(EXTRACT_LINKERS_OPTION)
        extract_linkers = true
      end

      if  !File.exist?( input_dir )
        Util.fatal_error( PRG_NAME, "error: input directory [#{input_dir}] does not exist" )
      end
      if  !File.exist?( out_dir )
        Util.fatal_error( PRG_NAME, "error: output directory [#{out_dir}] does not exist" )
      end
      if  !File.exist?( out_dir_doms )
        Util.fatal_error( PRG_NAME, "error: output directory [#{out_dir_doms}] does not exist" )
      end
      if !File.directory?( input_dir )
        Util.fatal_error( PRG_NAME, "error: [#{input_dir}] is not a directory" )
      end
      if !File.directory?( out_dir )
        Util.fatal_error( PRG_NAME, "error:  [#{out_dir}] is not a directory" )
      end
      if !File.directory?( out_dir_doms )
        Util.fatal_error( PRG_NAME, "error:  [#{out_dir_doms}] is not a directory" )
      end


      log = String.new

      log << "Program            : " + PRG_NAME + ld
      log << "Version            : " + PRG_VERSION + ld
      log << "Program date       : " + PRG_DATE + ld

      puts()
      puts( "Sequence names files suffix: " + seq_names_files_suffix )
      log << "Sequence names files suffix: " + seq_names_files_suffix + ld
      puts( "Input dir                  : " + input_dir )
      log << "Input dir                  : " + input_dir + ld
      puts( "Output dir                 : " + out_dir )
      log << "Output dir                 : " + out_dir + ld
      puts( "Output dir domains         : " + out_dir_doms )
      log << "Output dir domains         : " + out_dir_doms + ld
      if ( mapping_file != nil )
        puts( "Mapping file               : " + mapping_file )
        log << "Mapping file               : " + mapping_file + ld
      end
      if extension > 0
        puts( "Extension                  : " + extension.to_s )
        log << "Extension                  : " + extension.to_s + ld
      end
      if extract_linkers
        puts( "Extract linkers            : true" )
        log << "Extract linkers            : true" + ld
      end
      log << "Date                       : " + Time.now.to_s + ld
      puts

      if ( mapping_file != nil )
        species_codes_to_paths = extract_mappings( mapping_file )
      end

      input_files = obtain_inputfiles( input_dir, seq_names_files_suffix )

      counter = 0

      input_files.each { |input_file|
        counter += 1
        puts
        puts
        puts counter.to_s + "/" + input_files.size.to_s
        read_seq_family_file( input_file,
          seq_names_files_suffix,
          input_dir,
          species_codes_to_paths,
          log,
          out_dir,
          out_dir_doms,
          mapping_file,
          extension,
          extract_linkers )
      }
      puts
      Util.print_message( PRG_NAME, "OK" )
      puts

    end


    def read_seq_family_file( input_file,
        seq_names_files_suffix,
        input_dir,
        species_codes_to_paths,
        log,
        out_dir,
        out_dir_doms,
        mapping_file,
        extension,
        extract_linkers )

      begin
        Util.check_file_for_readability( input_file )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end
      basename = File.basename( input_file, seq_names_files_suffix )
      out_file_path_fasta_file                = out_dir + Constants::FILE_SEPARATOR + basename + FASTA_SUFFIX
      out_file_path_normalized_ids_fasta_file = out_dir + Constants::FILE_SEPARATOR + basename + FASTA_WITH_NORMALIZED_IDS_SUFFIX
      out_file_path_ids_map                   = out_dir + Constants::FILE_SEPARATOR + basename + NORMALIZED_IDS_MAP_SUFFIX
      doms_out_file_path_fasta_file           = out_dir_doms + Constants::FILE_SEPARATOR + basename + "_domains" + FASTA_SUFFIX
      doms_ext_out_file_path_fasta_file           = nil
      if extension > 0
        doms_ext_out_file_path_fasta_file = out_dir_doms + Constants::FILE_SEPARATOR + basename + "_domains_ext_" + extension.to_s + FASTA_SUFFIX
      end
      begin
        Util.check_file_for_writability( out_file_path_fasta_file )
        Util.check_file_for_writability( out_file_path_normalized_ids_fasta_file )
        Util.check_file_for_writability( out_file_path_ids_map  )
        Util.check_file_for_writability( doms_out_file_path_fasta_file  )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      ids_map_writer = nil
      begin
        ids_map_writer = File.open( out_file_path_ids_map, 'a' )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      current_species         = ""
      current_msa            = nil
      new_msa                = Msa.new
      new_msa_normalized_ids = Msa.new
      new_msa_domains = Msa.new
      new_msa_domains_extended = Msa.new
      per_species_counter = 0
      linkers = ""

      puts basename

      File.open( input_file ) do | file |
        species_counter = 1
        while line = file.gets
          line.strip!
          if !Util.is_string_empty?( line ) && !(line =~ /\s*#/ )
            values = line.split( PROTEINS_LIST_FILE_SEPARATOR )
            mod_line = nil
            if ( values.length < 2 )
              Util.fatal_error( PRG_NAME, "unexpected format: " + line )
            end
            species = values[ 0 ]
            #if species == "BRADI" || species == "ASPNG" || species == "SCLSC" || species == "PTEVA"  || species == "EIMTE"
            #  next
            #end
            seq_name = values[ 1 ]
            domain_ranges = nil
            if ( values.length > 3 )
              domain_ranges_block = values[ 3 ]
              domain_ranges = domain_ranges_block.split( "/" )
            end
            if ( species != current_species )
              current_species = species
              my_file = input_dir + Constants::FILE_SEPARATOR + current_species

              if ( !File.exist?( my_file ) )
                if species_codes_to_paths == nil
                  Util.fatal_error( PRG_NAME, "error: [#{my_file}] not found and no mapping file provided" )
                elsif ( !species_codes_to_paths.has_key?( current_species ) )
                  Util.fatal_error( PRG_NAME, "error: species [#{current_species}] not found in mapping file [#{mapping_file}]" )
                end
                my_file = species_codes_to_paths[ current_species ]
              end
              my_path = File.expand_path( my_file )
              my_readlink = my_path
              if ( File.symlink?( my_path ) )
                my_readlink = File.readlink( my_path )
              end
              current_msa = read_fasta_file( my_file )

              if ( per_species_counter > 0 )
                print_counts( per_species_counter, log, Constants::LINE_DELIMITER )
                per_species_counter = 0
              end
              puts " " + species_counter.to_s +  ":" + current_species + " [" + my_readlink + "]"
              log << species_counter.to_s <<  ": " << current_species << " [" + my_readlink + "]" << Constants::LINE_DELIMITER
              species_counter += 1
            end
            log << "   " << seq_name << Constants::LINE_DELIMITER
            per_species_counter = per_species_counter + 1
            seq = nil

            indices = current_msa.find_by_name_start( seq_name, true )
            if indices.size == 1
              seq =  current_msa.get_sequence( indices[ 0 ] )
            elsif indices.size == 0
              # Not found, try finding by partial match.
              begin
                seq = current_msa.get_by_name( seq_name, true, true )
              rescue ArgumentError => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s )
              end
            else
              Util.fatal_error( PRG_NAME, "error: seq name \"" + seq_name + "\" not unique"  )
            end

            normalized_id = per_species_counter.to_s( 16 ).upcase +
             "_" + current_species

            per_species_counter.to_i

            ids_map_writer.write( normalized_id + "\t" + seq.get_name + Constants::LINE_DELIMITER )

            orig_name = nil
            if seq != nil
              orig_name = seq.get_name
              seq.set_name( seq.get_name + " [" + current_species + "]" )
              new_msa.add_sequence( seq )
            else
              Util.fatal_error( PRG_NAME, "unexected error: seq is nil" )
            end

            if domain_ranges != nil
              first = true
              prev_to = -1

              domain_ranges.each { |range|
                if range != nil && range.length > 0
                  s = range.split("-")
                  from = s[ 0 ].to_i
                  to = s[ 1 ].to_i
                  new_msa_domains.add_sequence( Sequence.new( orig_name +
                       " [" + from.to_s +
                       "-" + to.to_s +
                       "] [" + basename + "] [" +
                       current_species + "]",
                      seq.get_sequence_as_string[from..to] ) )
                  if extension > 0
                    from_e = from - extension
                    if from_e < 0
                      from_e = 0
                    end
                    to_e = to + extension
                    if to_e > seq.get_sequence_as_string.length - 1
                      to_e = seq.get_sequence_as_string.length - 1
                    end
                    new_msa_domains_extended.add_sequence( Sequence.new( orig_name +
                         " [" + from.to_s +
                         "-" + to.to_s  +
                         "] [extended by " +
                         extension.to_s +
                         "] [" + basename + "] [" +
                         current_species + "]",
                        seq.get_sequence_as_string[ from_e..to_e ] ) )
                  end # extension > 0
                  if  extract_linkers
                    if first
                      first = false
                      f = 0
                      t = from - 1
                      if extension > 0
                        f = from - extension
                      end
                      mod_line = line + "\t[" + get_linker_sequence( f, t, seq ) + "|"
                    else
                      mod_line += get_linker_sequence( prev_to + 1, from - 1, seq ) + "|"
                    end
                    prev_to = to
                  end
                end # range != nil && range.length > 0
              }
              if extract_linkers && prev_to > 0
                f = prev_to + 1
                t = seq.get_sequence_as_string.length - 1
                if extension > 0
                  t = prev_to + extension
                end
                mod_line += get_linker_sequence( f, t, seq ) + "]"
              end
            end

            new_msa_normalized_ids.add_sequence( Sequence.new( normalized_id, seq.get_sequence_as_string ) )
            if mod_line
              linkers << mod_line + Constants::LINE_DELIMITER
            end
          end # !Util.is_string_empty?( line ) && !(line =~ /\s*#/ )
        end # while line = file.gets

      end

      ids_map_writer.close

      if ( per_species_counter > 0 )
        print_counts( per_species_counter, log, Constants::LINE_DELIMITER )
      end

      io = MsaIO.new()

      fasta_writer = FastaWriter.new()
      fasta_writer.remove_gap_chars
      fasta_writer.clean

      begin
        io.write_to_file( new_msa, out_file_path_fasta_file, fasta_writer )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      if new_msa_domains != nil
        begin
          io.write_to_file( new_msa_domains, doms_out_file_path_fasta_file, fasta_writer )
        rescue Exception => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s )
        end
      end

      if extension > 0 && new_msa_domains_extended != nil
        begin
          io.write_to_file( new_msa_domains_extended, doms_ext_out_file_path_fasta_file, fasta_writer )
        rescue Exception => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s )
        end
      end

      begin
        io.write_to_file( new_msa_normalized_ids, out_file_path_normalized_ids_fasta_file, fasta_writer )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      begin
        f = File.open( out_dir + Constants::FILE_SEPARATOR + basename +  LOG_SUFFIX , 'a' )
        f.print( log )
        f.close
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      if extract_linkers && linkers != nil
        begin
          f = File.open( out_dir + Constants::FILE_SEPARATOR + basename +  LINKERS_SUFFIX , 'a' )
          f.print( linkers )
          f.close
        rescue Exception => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s )
        end
      end
    end


    def get_linker_sequence( from, to, seq )
      if from < 0
        from = 0
      end
      if to > seq.get_sequence_as_string.length - 1
        to = seq.get_sequence_as_string.length - 1
      end
      if from > to
        return ""
      else
        return from.to_s + "-" + to.to_s + ":" + seq.get_sequence_as_string[ from..to ]
      end
    end

    def obtain_inputfiles( input_dir, seq_names_files_suffix )
      input_files = Array.new()
      Dir.foreach( input_dir ) { |file_name|
        if file_name.index( seq_names_files_suffix ) == ( file_name.size - seq_names_files_suffix.size )
          input_files.push( input_dir + Constants::FILE_SEPARATOR + file_name )
        end
      }
      input_files
    end

    def extract_mappings( mapping_file )
      species_code_to_path = Hash.new()
      File.open( mapping_file ) do | file |
        while line = file.gets
          if ( !Util.is_string_empty?( line ) && !(line =~ /\s*#/ ) )
            if ( line =~ /(\S+)\s+(\S+)/ )
              species = $1
              path = $2
              if ( species_code_to_path.has_key?( species ) )
                Util.fatal_error( PRG_NAME, "error: species code [#{species}] is not unique" )
              end
              if ( species_code_to_path.has_value?( path ) )
                Util.fatal_error( PRG_NAME, "error: path [#{path}] is not unique" )
              end
              if ( !File.exist?( path ) )
                Util.fatal_error( PRG_NAME, "error: file [#{path}] does not exist" )
              end
              if ( !File.file?( path ) )
                Util.fatal_error( PRG_NAME, "error: [#{path}] is not a regular file" )
              end
              if ( !File.readable?( path ) )
                Util.fatal_error( PRG_NAME, "error: file [#{path}] is not readable" )
              end
              if ( File.size( path ) < 10000 )
                Util.fatal_error( PRG_NAME, "error: file [#{path}] appears too small" )
              end
              if ( !Util.looks_like_fasta?( path ) )
                Util.fatal_error( PRG_NAME, "error: file [#{path}] does not appear to be a fasta file" )
              end
              species_code_to_path[ species ] = path
              puts species + " -> " + path
            end
          end
        end
      end
      species_code_to_path
    end

    def print_counts( per_species_counter, log, ld )
      puts "   sum: " + per_species_counter.to_s
      log << "   sum: " + per_species_counter.to_s + ld
    end

    def read_fasta_file( input )
      if @file_to_msa.has_key?( input )
        return @file_to_msa[ input ]
      end

      f = MsaFactory.new()
      msa = nil
      begin
        msa = f.create_msa_from_file( input, FastaParser.new() )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end
      if @seqs <= 400000
        @file_to_msa[ input ] = msa
        @seqs += msa.get_number_of_seqs
        puts "   total seqs in memory: " + @seqs.to_s
      end
      msa
    end

    def print_help()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb <sequence names file suffix> <input dir containing sequence names files " +
         "and possibly genome multiple-sequence ('fasta') files> <output directory for sequences> <output directory for domains> [mapping file for " +
         "genome multiple-sequence ('fasta') files not in input dir]" )
      puts()
      puts( "  option: -" + EXT_OPTION  + "=<int>: to extend extracted domains" )
      puts( "          -" + EXTRACT_LINKERS_OPTION  + "      : to extract linkers" )
      puts()
      puts( "  " + "Example: \"mse.rb .prot . protein_seqs domain_seqs ../genome_locations.txt\"" )
      puts()
    end

  end # class MultiSequenceExtractor
end
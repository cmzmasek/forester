#
# = lib/evo/apps/msa_processor.rb - MsaProcessor class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: msa_processor.rb,v 1.33 2010/12/13 19:00:10 cmzmasek Exp $
#

require 'date'
require 'set'

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/writer/phylip_sequential_writer'
require 'lib/evo/io/writer/nexus_writer'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/parser/general_msa_parser'
require 'lib/evo/io/writer/msa_writer'

module Evoruby

  class MsaProcessor

    PRG_NAME       = "msa_pro"
    PRG_DATE       = "131112"
    PRG_DESC       = "processing of multiple sequence alignments"
    PRG_VERSION    = "1.08"
    COPYRIGHT      = "2008-2010 Christian M Zmasek"
    CONTACT        = "phylosoft@gmail.com"
    WWW            = "www.phylosoft.org"


    NAME_LENGTH_DEFAULT                = 10
    WIDTH_DEFAULT_FASTA                = 60
    INPUT_TYPE_OPTION                  = "i"
    OUTPUT_TYPE_OPTION                 = "o"
    MAXIMAL_NAME_LENGTH_OPTION         = "n"
    DIE_IF_NAME_TOO_LONG               = "d"
    WIDTH_OPTION                       = "w"
    CLEAN_UP_SEQ_OPTION                = "c"
    REM_RED_OPTION                     = "rem_red"
    REMOVE_GAP_COLUMNS_OPTION          = "rgc"
    REMOVE_GAP_ONLY_COLUMNS            = "rgoc"
    REMOVE_COLUMNS_GAP_RATIO_OPTION    = "rr"
    REMOVE_ALL_GAP_CHARACTERS_OPTION   = "rg"
    REMOVE_ALL_SEQUENCES_LISTED_OPTION = "r"
    KEEP_ONLY_SEQUENCES_LISTED_OPTION  = "k"

    KEEP_MATCHING_SEQUENCES_OPTION     = "mk"
    REMOVE_MATCHING_SEQUENCES_OPTION   = "mr"

    TRIM_OPTION                        = "t"
    REMOVE_SEQS_GAP_RATIO_OPTION       = "rsgr"
    REMOVE_SEQS_NON_GAP_LENGTH_OPTION  = "rsl"
    SPLIT                              = "split"
    LOG_SUFFIX                         = "_msa_pro.log"
    HELP_OPTION_1                      = "help"
    HELP_OPTION_2                      = "h"


    def initialize()
      @input_format_set = false
      @output_format_set = false
      @fasta_input      = false
      @phylip_input     = true
      @name_length      = NAME_LENGTH_DEFAULT
      @name_length_set  = false
      @width            = WIDTH_DEFAULT_FASTA     # fasta only
      @pi_output        = true
      @fasta_output     = false
      @nexus_output     = false
      @clean            = false  # phylip only
      @rgc              = false
      @rgoc             = false
      @rg               = false  # fasta only
      @rem_red          = false
      @die_if_name_too_long  = false
      @rgr              = -1
      @rsgr             = -1
      @rsl              = -1
      @remove_matching  = nil
      @keep_matching    = nil

      @seqs_name_file   = nil
      @remove_seqs      = false
      @keep_seqs        = false
      @trim             = false
      @split            = -1
      @first            = -1
      @last             = -1
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

      if ( ARGV == nil || ARGV.length < 1 )
        Util.print_message( PRG_NAME, "Illegal number of arguments" )
        print_help
        exit( -1 )
      end

      begin
        cla = CommandLineArguments.new( ARGV )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "Error: " + e.to_s, STDOUT )
      end

      if ( cla.is_option_set?( HELP_OPTION_1 ) ||
           cla.is_option_set?( HELP_OPTION_2 ) )
        print_help
        exit( 0 )
      end

      if ( cla.get_number_of_files != 2 || ARGV.length < 2 )
        Util.print_message( PRG_NAME, "Illegal number of arguments" )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push( INPUT_TYPE_OPTION )
      allowed_opts.push( OUTPUT_TYPE_OPTION )
      allowed_opts.push( MAXIMAL_NAME_LENGTH_OPTION )
      allowed_opts.push( WIDTH_OPTION )
      allowed_opts.push( CLEAN_UP_SEQ_OPTION )
      allowed_opts.push( REMOVE_GAP_COLUMNS_OPTION )
      allowed_opts.push( REMOVE_GAP_ONLY_COLUMNS )
      allowed_opts.push( REMOVE_COLUMNS_GAP_RATIO_OPTION )
      allowed_opts.push( REMOVE_ALL_GAP_CHARACTERS_OPTION )
      allowed_opts.push( REMOVE_ALL_SEQUENCES_LISTED_OPTION )
      allowed_opts.push( KEEP_ONLY_SEQUENCES_LISTED_OPTION )
      allowed_opts.push( TRIM_OPTION )
      allowed_opts.push( REMOVE_SEQS_GAP_RATIO_OPTION )
      allowed_opts.push( REMOVE_SEQS_NON_GAP_LENGTH_OPTION )
      allowed_opts.push( SPLIT )
      allowed_opts.push( REM_RED_OPTION )
      allowed_opts.push( KEEP_MATCHING_SEQUENCES_OPTION )
      allowed_opts.push( REMOVE_MATCHING_SEQUENCES_OPTION )
      allowed_opts.push( DIE_IF_NAME_TOO_LONG )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
          "unknown option(s): " + disallowed )
      end

      input = cla.get_file_name( 0 )
      output = cla.get_file_name( 1 )

      analyze_command_line( cla )

      begin
        Util.check_file_for_readability( input )
      rescue IOError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      begin
        Util.check_file_for_writability( output )
      rescue IOError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      if ( @rg )
        set_pi_output( false )
        set_fasta_output( true )
        set_nexus_output( false )
      end

      if ( !@input_format_set )
        fasta_like = false
        begin
          fasta_like = Util.looks_like_fasta?( input )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s )
        end
        @fasta_input = fasta_like
        @phylip_input = !fasta_like
        if ( !@output_format_set )
          @fasta_output = fasta_like
          @pi_output = !fasta_like
          @nexus_output = false
        end
      end

      ld = Constants::LINE_DELIMITER
      log = PRG_NAME + " " + PRG_VERSION + " [" + PRG_DATE + "]" + " LOG" + ld
      now = DateTime.now
      log << "Date/time: " + now.to_s + ld

      puts()
      puts( "Input alignment  : " + input )
      log << "Input alignment  : " + input + ld
      puts( "Output alignment : " + output )
      log << "Output alignment : " + output + ld
      if ( @phylip_input )
        puts( "Input is         : Phylip, or something like it" )
        log << "Input is         : Phylip, or something like it" + ld
      elsif ( @fasta_input )
        puts( "Input is         : Fasta" )
        log << "Input is         : Fasta" + ld
      end
      if( @rgr >= 0 )
        puts( "Max col gap ratio: " + @rgr.to_s )
        log << "Max col gap ratio: " + @rgr.to_s + ld
      elsif ( @rgc )
        puts( "Remove gap colums" )
        log << "Remove gap colums" + ld
      elsif( @rgoc )
        puts( "Remove gap only colums" )
        log << "Remove gap only colums" + ld
      end
      if ( @clean )
        puts( "Clean up         : true" )
        log << "Clean up         : true" + ld
      end

      if ( @pi_output )
        puts( "Output is        : Phylip interleaved" )
        log << "Output is        : Phylip interleaved" + ld
      elsif ( @fasta_output )
        puts( "Output is        : Fasta" )
        log << "Output is        : Fasta" + ld
        if ( @width )
          puts( "Width            : " + @width.to_s )
          log << "Width            : " + @width.to_s + ld
        end
        if ( @rg )
          puts( "Remove all gap characters (alignment is destroyed)" )
          log << "Remove all gap characters (alignment is destroyed)" + ld
        end
      elsif ( @nexus_output )
        puts( "Output is        : Nexus" )
        log << "Output is        : Nexus" + ld
      end
      if ( @name_length_set || !@fasta_output )
        puts( "Max name length  : " + @name_length.to_s )
        log << "Max name length  : " + @name_length.to_s + ld
      end
      if( @rsgr >= 0 )
        puts( "Remove sequences for which the gap ratio > " + @rsgr.to_s )
        log << "Remove sequences for which the gap ratio > " + @rsgr.to_s + ld
      end
      if( @rsl >= 0 )
        puts( "Remove sequences with less than "  + @rsl.to_s + " non-gap characters" )
        log << "Remove sequences with less than "  + @rsl.to_s + " non-gap characters" + ld
      end
      if ( @remove_seqs )
        puts( "Remove sequences listed in: " + @seqs_name_file )
        log << "Remove sequences listed in: " + @seqs_name_file + ld
      elsif ( @keep_seqs )
        puts( "Keep only sequences listed in: " + @seqs_name_file )
        log << "Keep only sequences listed in: " + @seqs_name_file + ld
      end
      if ( @trim )
        puts( "Keep only columns from: "+ @first.to_s + " to " + @last.to_s )
        log << "Keep only columns from: "+ @first.to_s + " to " + @last.to_s + ld
      end
      if ( @rem_red )
        puts( "Remove redundant sequences: true" )
        log << "Remove redundant sequences: true" + ld
      end
      if ( @split > 0 )
        puts( "Split            : " + @split.to_s )
        log << "Split            : " + @split.to_s + ld
      end
      puts()

      f = MsaFactory.new()

      msa = nil

      begin
        if ( @phylip_input )
          msa = f.create_msa_from_file( input, GeneralMsaParser.new() )
        elsif ( @fasta_input )
          msa = f.create_msa_from_file( input, FastaParser.new() )
        end
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      end

      if ( msa.is_aligned() )
        Util.print_message( PRG_NAME, "Length of original alignment         : " + msa.get_length.to_s )
        log << "Length of original alignment         : " + msa.get_length.to_s + ld
        gp = msa.calculate_gap_proportion
        Util.print_message( PRG_NAME, "Gap-proportion of original alignment : " + gp.to_s )
        log << "Gap-proportion of original alignment : " +  gp.to_s + ld
      else
        Util.print_message( PRG_NAME, "Input is not aligned" )
        log << "Input is not aligned" + ld
      end

      all_names = Set.new()
      for i in 0 ... msa.get_number_of_seqs()
        current_name = msa.get_sequence( i ).get_name
        if all_names.include?( current_name )
          Util.print_warning_message( PRG_NAME, "sequence name [" + current_name + "] is not unique" )
        else
          all_names.add( current_name )
        end
      end

      begin

        if ( @remove_seqs || @keep_seqs )
          names = Util.file2array( @seqs_name_file, true )
          if ( names == nil ||  names.length() < 1 )
            error_msg = "file \"" + @seqs_name_file.to_s + "\" appears empty"
            Util.fatal_error( PRG_NAME, error_msg )
          end

          if ( @remove_seqs )
            c = 0
            for i in 0 ... names.length()
              to_delete = msa.find_by_name( names[ i ], true, false )
              if ( to_delete.length() < 1 )
                error_msg = "sequence name \"" + names[ i ] + "\" not found"
                Util.fatal_error( PRG_NAME, error_msg )
              elsif ( to_delete.length() > 1 )
                error_msg = "sequence name \"" + names[ i ] + "\" is not unique"
                Util.fatal_error( PRG_NAME, error_msg )
              else
                msa.remove_sequence!( to_delete[ 0 ] )
                c += 1
              end
            end
            Util.print_message( PRG_NAME, "Removed " + c.to_s + " sequences" )
            log <<  "Removed " + c.to_s + " sequences" + ld
          elsif ( @keep_seqs )
            msa_new = Msa.new()
            r = 0
            k = 0
            for j in 0 ... msa.get_number_of_seqs()
              if ( names.include?( msa.get_sequence( j ).get_name() ) )
                msa_new.add_sequence( msa.get_sequence( j ) )
                k += 1
              else
                r += 1
              end
            end
            msa = msa_new
            Util.print_message( PRG_NAME, "Kept    " + k.to_s + " sequences" )
            log << "Kept    " + k.to_s + " sequences" + ld
            Util.print_message( PRG_NAME, "Removed " + r.to_s + " sequences" )
            log << "removed " + r.to_s + " sequences" + ld
          end
        end

        if ( @trim )
          msa.trim!( @first, @last )
        end
        if( @rgr >= 0 )
          msa.remove_gap_columns_w_gap_ratio!( @rgr )
        elsif ( @rgc )
          msa.remove_gap_columns!()
        elsif( @rgoc )
          msa.remove_gap_only_columns!()
        end
        if( @rsgr >= 0 )
          n = msa.get_number_of_seqs()
          removed = msa.remove_sequences_by_gap_ratio!( @rsgr )
          k = msa.get_number_of_seqs()
          r = n - k
          Util.print_message( PRG_NAME, "Kept    " + k.to_s + " sequences" )
          log << "Kept    " + k.to_s + " sequences" + ld
          Util.print_message( PRG_NAME, "Removed " + r.to_s + " sequences"  )
          log << "Removed " + r.to_s + " sequences:" + ld
          removed.each { | seq_name |
            log << "         " + seq_name  + ld
          }
        end
        if( @rsl >= 0 )
          n = msa.get_number_of_seqs()
          removed = msa.remove_sequences_by_non_gap_length!( @rsl )
          k = msa.get_number_of_seqs()
          r = n - k
          Util.print_message( PRG_NAME, "Kept    " + k.to_s + " sequences" )
          log << "Kept    " + k.to_s + " sequences" + ld
          Util.print_message( PRG_NAME, "Removed " + r.to_s + " sequences" )
          log << "Removed " + r.to_s + " sequences:" + ld
          removed.each { | seq_name |
            log << "         " + seq_name  + ld
          }
        end
        if ( @keep_matching )
          n = msa.get_number_of_seqs
          to_be_removed = Set.new
          for ii in 0 ...  n
            seq = msa.get_sequence( ii )
            if !seq.get_name.downcase.index( @keep_matching.downcase )
              to_be_removed.add( ii )
            end
          end
          to_be_removed_ary = to_be_removed.to_a.sort.reverse
          to_be_removed_ary.each { | index |
            msa.remove_sequence!( index )
          }
          # msa = sort( msa )
        end
        if ( @remove_matching )
          n = msa.get_number_of_seqs
          to_be_removed = Set.new
          for iii in 0 ... n

            seq = msa.get_sequence( iii )

            if seq.get_name.downcase.index( @remove_matching.downcase )
              to_be_removed.add( iii )
            end
          end
          to_be_removed_ary = to_be_removed.to_a.sort.reverse
          to_be_removed_ary.each { | index |
            msa.remove_sequence!( index )
          }
          msa = sort( msa )
        end



        if ( @split > 0 )
          begin
            msas = msa.split( @split, true )
            io = MsaIO.new()
            w = MsaWriter
            if ( @pi_output )
              w = PhylipSequentialWriter.new()
              w.clean( @clean )
              w.set_max_name_length( @name_length )
            elsif( @fasta_output )
              w = FastaWriter.new()
              w.set_line_width( @width )
              if ( @rg )
                w.remove_gap_chars( true )
                Util.print_warning_message( PRG_NAME, "removing gap character, the output is likely to become unaligned" )
                log << "removing gap character, the output is likely to become unaligned" + ld
              end
              w.clean( @clean )
              if ( @name_length_set )
                w.set_max_name_length( @name_length )
              end
            elsif( @nexus_output )
              w = NexusWriter.new()
              w.clean( @clean )
              w.set_max_name_length( @name_length )
            end
            i = 0
            for m in msas
              i = i + 1
              io.write_to_file( m, output + "_" + i.to_s, w )
            end
            Util.print_message( PRG_NAME, "wrote " + msas.length.to_s + " files"  )
            log << "wrote " + msas.length.to_s + " files" + ld
          rescue Exception => e
            Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
          end

        end
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      end

      if ( @split <= 0 )

        unless ( @rg )
          if ( msa.is_aligned() )
            Util.print_message( PRG_NAME, "Length of processed alignment        : " + msa.get_length.to_s )
            log <<  "Length of processed alignment        : " + msa.get_length.to_s + ld
            gp = msa.calculate_gap_proportion
            Util.print_message( PRG_NAME, "Gap-proportion of processed alignment: " + gp.to_s )
            log << "Gap-proportion of processed alignment: " +  gp.to_s + ld
          else
            min = 0
            max = 0
            sum = 0
            first = true
            for s in 0 ... msa.get_number_of_seqs
              seq = msa.get_sequence( s )
              l = seq.get_length
              sum += l
              if l > max
                max = l
              end
              if first || l < min
                min = l
              end
              first = false
            end
            avg = sum / msa.get_number_of_seqs
            Util.print_message( PRG_NAME, "Output is not aligned" )
            log << "Output is not aligned" + ld
            Util.print_message( PRG_NAME, "Shortest sequence                    : " + min.to_s )
            log <<  "Shortest sequence                    : " + min.to_s + ld
            Util.print_message( PRG_NAME, "Longest sequence                     : " + max.to_s )
            log <<  "Longest sequence                     : " + max.to_s + ld
            Util.print_message( PRG_NAME, "Average length                       : " + avg.to_s )
            log <<  "Average length                       : " + avg.to_s + ld

          end
        end

        if @rem_red
          removed = msa.remove_redundant_sequences!( true, false )
          if removed.size > 0
            identicals = msa.get_identical_seqs_detected
            log << "the following " + identicals.size.to_s + " sequences are identical:" + ld
            identicals.each { | s |
              log << s + ld
            }
            log << "ignoring the following " + removed.size.to_s + " redundant sequences:" + ld
            removed.each { | seq_name |
              log << seq_name + ld
            }
            Util.print_message( PRG_NAME, "will store " + msa.get_number_of_seqs.to_s + " non-redundant sequences" )
            log << "will store " + msa.get_number_of_seqs.to_s + " non-redundant sequences" + ld
          end
        end

        io = MsaIO.new()

        w = MsaWriter

        if ( @pi_output )
          w = PhylipSequentialWriter.new()
          w.clean( @clean )
          w.set_max_name_length( @name_length )
          w.set_exception_if_name_too_long( @die_if_name_too_long )
        elsif( @fasta_output )
          w = FastaWriter.new()
          w.set_line_width( @width )
          if ( @rg )
            w.remove_gap_chars( true )
            Util.print_warning_message( PRG_NAME, "removing gap characters, the output is likely to become unaligned"  )
            log << "removing gap character, the output is likely to become unaligned" + ld
          end
          w.clean( @clean )
          if ( @name_length_set )
            w.set_max_name_length( @name_length )
            w.set_exception_if_name_too_long( @die_if_name_too_long )
          end
        elsif( @nexus_output )
          w = NexusWriter.new()
          w.clean( @clean )
          w.set_max_name_length( @name_length )
          w.set_exception_if_name_too_long( @die_if_name_too_long )
        end

        begin
          io.write_to_file( msa, output, w )
        rescue Exception => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s )
        end

        Util.print_message( PRG_NAME, "Number of sequences in output        : " + msa.get_number_of_seqs.to_s )
        log << "Number of sequences in output        : " + msa.get_number_of_seqs.to_s + ld

        begin
          f = File.open( output + LOG_SUFFIX, 'a' )
          f.print( log )
          f.close
        rescue Exception => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s )
        end


      end
      Util.print_message( PRG_NAME, "OK" )
      puts
    end


    private

    def sort( msa )
      names = Set.new
      for i in 0 ... msa.get_number_of_seqs
        name = msa.get_sequence( i ).get_name
        names.add( name )
      end
      sorted_ary = names.to_a.sort
      new_msa = Msa.new
      sorted_ary.each { | seq_name |
        seq = msa.get_sequence( msa.find_by_name( seq_name, true, false )[ 0 ] )
        new_msa.add_sequence( seq )
      }
      new_msa
    end

    def set_fasta_input( fi = true )
      @fasta_input = fi
      @input_format_set = true
    end
    def set_phylip_input( pi = true )
      @phylip_input = pi
      @input_format_set = true
    end
    def set_name_length( i )
      @name_length = i
      @name_length_set = true
    end
    def set_width( i )
      @width = i
    end
    def set_fasta_output( fo = true )
      @fasta_output = fo
      @output_format_set = true
    end
    def set_pi_output( pso = true )
      @pi_output = pso
      @output_format_set = true
    end
    def set_nexus_output( nexus = true )
      @nexus_output = nexus
      @output_format_set = true
    end
    def set_clean( c = true )
      @clean = c
    end
    def set_remove_gap_columns( rgc = true )
      @rgc = rgc
    end
    def set_remove_gap_only_columns( rgoc = true )
      @rgoc = rgoc
    end
    def set_remove_gaps( rg = true )
      @rg = rg
    end
    def set_remove_gap_ratio( rgr )
      @rgr = rgr
    end
    def set_remove_seqs_gap_ratio( rsgr )
      @rsgr = rsgr
    end
    def set_remove_seqs_min_non_gap_length( rsl )
      @rsl = rsl
    end
    def set_remove_seqs( file )
      @seqs_name_file = file
      @remove_seqs    = true
      @keep_seqs      = false
    end
    def set_keep_seqs( file )
      @seqs_name_file = file
      @keep_seqs      = true
      @remove_seqs    = false
    end
    def set_trim( first, last )
      @trim            = true
      @first           = first
      @last            = last
    end
    def set_remove_matching( remove )
      @remove_matching  = remove
    end
    def set_keep_matching( keep )
      @keep_matching = keep
    end
    def set_rem_red( rr )
      @rem_red = rr
    end



    def set_split( s )
      if ( s > 0 )
        @split            = s
        @clean            = false  # phylip only
        @rgc              = false
        @rgoc             = false
        @rg               = false  # fasta only
        @rgr              = -1
        @rsgr             = -1
        @rsl              = -1
        @seqs_name_file   = nil
        @remove_seqs      = false
        @keep_seqs        = false
        @trim             = false
        @first            = -1
        @last             = -1
      end
    end

    def analyze_command_line( cla )
      if ( cla.is_option_set?( INPUT_TYPE_OPTION ) )
        begin
          type = cla.get_option_value( INPUT_TYPE_OPTION )
          if ( type == "p" )
            set_phylip_input( true )
            set_fasta_input( false )
          elsif ( type == "f" )
            set_fasta_input( true )
            set_phylip_input( false )
          end
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( OUTPUT_TYPE_OPTION ) )
        begin
          type = cla.get_option_value( OUTPUT_TYPE_OPTION )
          if ( type == "p" )
            set_pi_output( true )
            set_fasta_output( false )
            set_nexus_output( false )
          elsif ( type == "f" )
            set_pi_output( false )
            set_fasta_output( true )
            set_nexus_output( false )
          elsif ( type == "n" )
            set_pi_output( false )
            set_fasta_output( false )
            set_nexus_output( true )
          end
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( MAXIMAL_NAME_LENGTH_OPTION ) )
        begin
          l = cla.get_option_value_as_int( MAXIMAL_NAME_LENGTH_OPTION )
          set_name_length( l )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( WIDTH_OPTION ) )
        begin
          w = cla.get_option_value_as_int( WIDTH_OPTION )
          set_width( w )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( CLEAN_UP_SEQ_OPTION ) )
        set_clean( true )
      end
      if ( cla.is_option_set?( REMOVE_GAP_COLUMNS_OPTION ) )
        set_remove_gap_columns( true )
      end
      if ( cla.is_option_set?( REM_RED_OPTION ) )
        set_rem_red( true )
      end
      if ( cla.is_option_set?( REMOVE_GAP_ONLY_COLUMNS ) )
        set_remove_gap_only_columns( true )
      end
      if ( cla.is_option_set?( REMOVE_ALL_GAP_CHARACTERS_OPTION ) )
        set_remove_gaps( true )
      end
      if ( cla.is_option_set?( REMOVE_COLUMNS_GAP_RATIO_OPTION ) )
        begin
          f = cla.get_option_value_as_float( REMOVE_COLUMNS_GAP_RATIO_OPTION )
          set_remove_gap_ratio( f )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( REMOVE_ALL_SEQUENCES_LISTED_OPTION ) )
        begin
          s = cla.get_option_value( REMOVE_ALL_SEQUENCES_LISTED_OPTION )
          set_remove_seqs( s )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( KEEP_ONLY_SEQUENCES_LISTED_OPTION ) )
        begin
          s = cla.get_option_value( KEEP_ONLY_SEQUENCES_LISTED_OPTION )
          set_keep_seqs( s )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( TRIM_OPTION ) )
        begin
          s = cla.get_option_value( TRIM_OPTION )
          if ( s =~ /(\d+)-(\d+)/ )
            set_trim( $1.to_i(), $2.to_i() )
          else
            puts( "illegal argument" )
            print_help
            exit( -1 )
          end
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( REMOVE_SEQS_GAP_RATIO_OPTION ) )
        begin
          f = cla.get_option_value_as_float( REMOVE_SEQS_GAP_RATIO_OPTION )
          set_remove_seqs_gap_ratio( f )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( REMOVE_SEQS_NON_GAP_LENGTH_OPTION ) )
        begin
          f = cla.get_option_value_as_int( REMOVE_SEQS_NON_GAP_LENGTH_OPTION )
          set_remove_seqs_min_non_gap_length( f )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( SPLIT ) )
        begin
          s = cla.get_option_value_as_int( SPLIT )
          set_split( s )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end

      end
      if ( cla.is_option_set?( REMOVE_MATCHING_SEQUENCES_OPTION ) )
        begin
          s = cla.get_option_value( REMOVE_MATCHING_SEQUENCES_OPTION )
          set_remove_matching( s )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( KEEP_MATCHING_SEQUENCES_OPTION ) )
        begin
          s = cla.get_option_value( KEEP_MATCHING_SEQUENCES_OPTION )
          set_keep_matching( s )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
      end
      if ( cla.is_option_set?( DIE_IF_NAME_TOO_LONG ) )
        @die_if_name_too_long = true
      end


    end

    def print_help()
      puts()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb [options] <input alignment> <output>" )
      puts()
      puts( "  options: -" + INPUT_TYPE_OPTION + "=<input type>: f for fasta, p for phylip selex type" )
      puts( "           -" + OUTPUT_TYPE_OPTION + "=<output type>: f for fasta, n for nexus, p for phylip sequential (default)" )
      puts( "           -" + MAXIMAL_NAME_LENGTH_OPTION + "=<n>: n=maximal name length (default for phylip 10, for fasta: unlimited )" )
      puts( "           -" + DIE_IF_NAME_TOO_LONG + ": die if sequence name too long" )
      puts( "           -" + WIDTH_OPTION + "=<n>: n=width (fasta output only, default is 60)" )
      puts( "           -" + CLEAN_UP_SEQ_OPTION + ": clean up sequences" )
      puts( "           -" + REMOVE_GAP_COLUMNS_OPTION + ": remove gap columns" )
      puts( "           -" + REMOVE_GAP_ONLY_COLUMNS + ": remove gap-only columns" )
      puts( "           -" + REMOVE_COLUMNS_GAP_RATIO_OPTION + "=<n>: remove columns for which ( seqs with gap / number of sequences > n )" )
      puts( "           -" + REMOVE_ALL_GAP_CHARACTERS_OPTION + ": remove all gap characters (destroys alignment, fasta output only)" )
      puts( "           -" + REMOVE_ALL_SEQUENCES_LISTED_OPTION + "=<file>: remove all sequences listed in file" )
      puts( "           -" + KEEP_ONLY_SEQUENCES_LISTED_OPTION + "=<file>: keep only sequences listed in file" )
      puts( "           -" + TRIM_OPTION + "=<first>-<last>: remove columns before first and after last" )
      puts( "           -" + REMOVE_SEQS_GAP_RATIO_OPTION + "=<n>: remove sequences for which the gap ratio > n (after column operations)" )
      puts( "           -" + REMOVE_SEQS_NON_GAP_LENGTH_OPTION + "=<n> remove sequences with less than n non-gap characters (after column operations)" )
      puts( "           -" + REMOVE_MATCHING_SEQUENCES_OPTION + "=<s> remove all sequences with names containing s" )
      puts( "           -" + KEEP_MATCHING_SEQUENCES_OPTION + "=<s> keep only sequences with names containing s" )
      puts( "           -" + SPLIT + "=<n> split a fasta file into n files of equal number of sequences (expect for " )
      puts( "            last one), cannot be used with other options" )
      puts( "           -" + REM_RED_OPTION + ": remove redundant sequences" )
      puts()
    end





  end # class MsaProcessor


end # module Evoruby

#
# = lib/evo/msa/msa.rb - Msa class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: msa.rb,v 1.11 2009/01/03 00:42:08 cmzmasek Exp $
#


require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/sequence/sequence'

module Evoruby

  class Msa

    def initialize()
      @sequences = Array.new
      @identical_seqs_detected = Array.new
    end


    def add_sequence( sequence )
      @sequences.push( sequence )
    end

    def add( name, molecular_sequence_str )
      add_sequence( Sequence.new( name, molecular_sequence_str ) )
    end

    def get_sequence( index )
      if ( index < 0 || index > get_number_of_seqs() - 1 )
        error_msg = "attempt to get sequence " <<
         index.to_s << " in alignment of " << get_number_of_seqs().to_s <<
         " sequences"
        raise ArgumentError, error_msg
      end
      return @sequences[ index ]
    end

    def remove_sequence!( index )
      if ( index < 0 || index > get_number_of_seqs() - 1 )
        error_msg = "attempt to remove sequence " <<
         index.to_s << " in alignment of " << get_number_of_seqs().to_s <<
         " sequences"
        raise ArgumentError, error_msg
      end
      @sequences.delete_at( index )
    end

    def get_identical_seqs_detected
      @identical_seqs_detected
    end


    def is_aligned()
      if ( get_number_of_seqs < 1 )
        return false
      else
        l = @sequences[ 0 ].get_length()
        for i in 0 ... get_number_of_seqs()
          if ( get_sequence( i ).get_length() != l )
            return false
          end
        end
      end
      return true
    end

    def find_by_name( name, case_sensitive, partial_match )
      indices = Array.new()
      for i in 0 ... get_number_of_seqs()
        current_name = get_sequence( i ).get_name()
        if !case_sensitive
          current_name = current_name.downcase
          name = name.downcase
        end
        if current_name == name ||
           ( partial_match && current_name.include?( name ) )
          indices.push( i )
        end
      end
      indices
    end

    def find_by_name_pattern( name_re, avoid_similar_to = true )
      indices = []
      for i in 0 ... get_number_of_seqs()
        if avoid_similar_to
          m = name_re.match( get_sequence( i ).get_name() )
          if m && !m.pre_match.downcase.include?( "similar to " )
            indices.push( i )
          end
        else
          if name_re.match( get_sequence( i ).get_name() )
            indices.push( i )
          end
        end
      end
      indices
    end

    # throws ArgumentError
    def get_by_name_pattern( name_re , avoid_similar_to = true )
      indices = find_by_name_pattern( name_re, avoid_similar_to  )
      if ( indices.length > 1 )
        error_msg = "pattern  " + name_re.to_s + " not unique"
        raise ArgumentError, error_msg
      elsif ( indices.length < 1 )
        error_msg = "pattern " + name_re.to_s + " not found"
        raise ArgumentError, error_msg
      end
      get_sequence( indices[ 0 ] )
    end

    def find_by_name_start( name, case_sensitive )
      indices = []
      for i in 0 ... get_number_of_seqs()
        get_sequence( i ).get_name() =~ /^\s*(\S+)/
        current_name = $1
        if !case_sensitive
          current_name = current_name.downcase
          name = name.downcase
        end
        if  ( current_name == name )
          indices.push( i )
        end
      end
      indices
    end

    def has?( name, case_sensitive = true, partial_match = false )
      for i in 0 ... get_number_of_seqs()
        current_name = get_sequence( i ).get_name()
        if !case_sensitive
          current_name = current_name.downcase
          name = name.downcase
        end
        if current_name == name ||
           ( partial_match && current_name.include?( name ) )
          return true
        end
      end
      false
    end

    # throws ArgumentError
    def get_by_name( name, case_sensitive = true, partial_match = false )
      indices = find_by_name( name, case_sensitive, partial_match )
      if ( indices.length > 1 )
        error_msg = "\"" + name + "\" not unique"
        raise ArgumentError, error_msg
      elsif ( indices.length < 1 )
        error_msg = "\"" + name + "\" not found"
        raise ArgumentError, error_msg
      end
      get_sequence( indices[ 0 ] )
    end

    # throws ArgumentError
    def get_by_name_start( name, case_sensitive = true )
      indices = find_by_name_start( name, case_sensitive )
      if ( indices.length > 1 )
        error_msg = "\"" + name + "\" not unique"
        raise ArgumentError, error_msg
      elsif ( indices.length < 1 )
        error_msg = "\"" + name + "\" not found"
        raise ArgumentError, error_msg
      end
      get_sequence( indices[ 0 ] )
    end


    def get_sub_alignment( seq_numbers )
      msa = Msa.new()
      for i in 0 ... seq_numbers.length()
        msa.add_sequence( get_sequence( seq_numbers[ i ] ).copy() )
      end
      return msa
    end

    def get_number_of_seqs()
      @sequences.length
    end

    def get_length
      if !is_aligned()
        error_msg = "attempt to get length of unaligned msa"
        raise StandardError, error_msg, caller
      end
      if get_number_of_seqs() < 1
        -1
      else
        @sequences[ 0 ].get_length()
      end
    end

    def to_str
      to_fasta
    end

    def to_fasta
      s = String.new
      for i in 0...get_number_of_seqs
        s += @sequences[ i ].to_fasta + Constants::LINE_DELIMITER
      end
      s
    end


    def print_overlap_diagram( min_overlap = 1, io = STDOUT, max_name_length = 10 )
      if ( !is_aligned() )
        error_msg = "attempt to get overlap diagram of unaligned msa"
        raise StandardError, error_msg, caller
      end
      for i in 0 ... get_number_of_seqs()
        io.print( Util.normalize_seq_name( get_sequence( i ).get_name(), max_name_length ) )
        for j in 0 ... get_number_of_seqs()
          if i == j
            io.print( " " )
          else
            if overlap?( i, j, min_overlap )
              io.print( "+" )
            else
              io.print( "-" )
            end
          end
        end
        io.print( Evoruby::Constants::LINE_DELIMITER )
      end
    end

    #returns array of Msa with an overlap of min_overlap
    def split_into_overlapping_msa( min_overlap = 1 )
      if ( !is_aligned() )
        error_msg = "attempt to split into overlapping msas of unaligned msa"
        raise StandardError, error_msg, caller
      end
      msas = Array.new()
      bins = get_overlaps( min_overlap )
      for i in 0 ... bins.length
        msas.push( get_sub_alignment( bins[ i ] ) )
      end
      msas
    end

    def overlap?( index_1, index_2, min_overlap = 1 )
      seq_1 = get_sequence( index_1 )
      seq_2 = get_sequence( index_2 )
      overlap_count = 0
      for i in 0...seq_1.get_length()
        if !Util.is_aa_gap_character?( seq_1.get_character_code( i ) ) &&
           !Util.is_aa_gap_character?( seq_2.get_character_code( i ) )
          overlap_count += 1
          if overlap_count >= min_overlap
            return true
          end
        end
      end
      return false
    end

    def calculate_overlap( index_1, index_2 )
      seq_1 = get_sequence( index_1 )
      seq_2 = get_sequence( index_2 )
      overlap_count = 0
      for i in 0...seq_1.get_length
        if !Util.is_aa_gap_character?( seq_1.get_character_code( i ) ) &&
           !Util.is_aa_gap_character?( seq_2.get_character_code( i ) )
          overlap_count += 1
        end
      end
      overlap_count
    end

    def calculate_identities( index_1, index_2 )
      seq_1 = get_sequence( index_1 )
      seq_2 = get_sequence( index_2 )
      identities_count = 0
      for i in 0...seq_1.get_length
        if !Util.is_aa_gap_character?( seq_1.get_character_code( i ) ) &&
           !Util.is_aa_gap_character?( seq_2.get_character_code( i ) ) &&
           seq_1.get_character_code( i ) != 63 &&
           ( seq_1.get_residue( i ).downcase() ==
             seq_2.get_residue( i ).downcase() )
          identities_count += 1
        end
      end
      identities_count
    end

    def remove_gap_only_columns!()
      remove_columns!( get_gap_only_columns() )
    end

    def remove_gap_columns!()
      remove_columns!( get_gap_columns() )
    end

    # removes columns for which seqs with gap / number of sequences > gap_ratio
    def remove_gap_columns_w_gap_ratio!( gap_ratio )
      remove_columns!( get_gap_columns_w_gap_ratio( gap_ratio ) )
    end


    def remove_sequences_by_gap_ratio!( gap_ratio )
      if ( !is_aligned() )
        error_msg = "attempt to remove sequences by gap ratio on unaligned msa"
        raise StandardError, error_msg, caller
      end
      n = get_number_of_seqs
      removed = Array.new
      for s in 0 ... n
        if ( get_sequence( ( n - 1 ) - s  ).get_gap_ratio() > gap_ratio )
          if ( Evoruby::Constants::VERBOSE )
            puts( "removed: " + get_sequence( ( n - 1 ) - s  ).get_name )
          end
          removed << get_sequence( ( n - 1 ) - s  ).get_name
          remove_sequence!( ( n - 1 ) - s  )
        end
      end
      removed
    end


    def remove_redundant_sequences!( consider_taxonomy = false, verbose = false )
      n = get_number_of_seqs
      removed = Array.new
      to_be_removed = Set.new
      @identical_seqs_detected = Array.new
      for i in 0 ... ( n - 1 )
        for j in ( i + 1 ) ... n
          if !to_be_removed.include?( i ) && !to_be_removed.include?( j )
            if  !consider_taxonomy ||
               ( ( get_sequence( i ).get_taxonomy == nil && get_sequence( j ).get_taxonomy == nil ) ||
                 ( get_sequence( i ).get_taxonomy == get_sequence( j ).get_taxonomy ) )
              if Util.clean_seq_str( get_sequence( i ).get_sequence_as_string ) ==
                 Util.clean_seq_str( get_sequence( j ).get_sequence_as_string )
                to_be_removed.add( j )

                tax_i = ""
                tax_j = ""
                if get_sequence( i ).get_taxonomy != nil
                  tax_i = get_sequence( i ).get_taxonomy.get_name
                end
                if get_sequence( j ).get_taxonomy != nil
                  tax_j = get_sequence( j ).get_taxonomy.get_name
                end
                identical_pair = get_sequence( i ).get_name + " [#{tax_i}] == " + get_sequence( j ).get_name + " [#{tax_j}]"
                @identical_seqs_detected.push( identical_pair )
                if verbose
                  puts identical_pair
                end
              end
            end
          end
        end
      end

      to_be_removed_ary = to_be_removed.to_a.sort.reverse

      to_be_removed_ary.each { | index |
        removed.push( get_sequence( index ).get_name )
        remove_sequence!( index )
      }
      removed
    end


    def remove_sequences_by_non_gap_length!( min_non_gap_length )
      if ( !is_aligned() )
        error_msg = "attempt to remove sequences by non gap length on unaligned msa"
        raise StandardError, error_msg, caller
      end
      n = get_number_of_seqs
      l = get_length
      removed = Array.new
      for s in 0 ... n
        if ( ( l - get_sequence( ( n - 1 ) - s ).get_gap_length ) < min_non_gap_length )
          if ( Evoruby::Constants::VERBOSE )
            puts( "removed: " + get_sequence( ( n - 1 ) - s  ).get_name )
          end
          removed << get_sequence( ( n - 1 ) - s  ).get_name
          remove_sequence!( ( n - 1 ) - s )
        end
      end
      removed
    end

    def trim!( first, last )
      cols = Array.new()
      for i in 0 ... get_length()
        if ( i < first || i > last )
          cols.push( i )
        end
      end
      remove_columns!( cols )
    end

    def get_gap_only_columns()
      if ( !is_aligned() )
        error_msg = "attempt to get gap only columns of unaligned msa"
        raise StandardError, error_msg, caller
      end
      cols = Array.new()
      for c in 0 ... get_length
        nogap_char_found = false
        for s in 0 ... get_number_of_seqs
          unless Util.is_aa_gap_character?( get_sequence( s ).get_character_code( c ) )
            nogap_char_found = true
            break
          end
        end
        unless nogap_char_found
          cols.push( c )
        end
      end
      return cols
    end

    def calculate_gap_proportion()
      if ( !is_aligned() )
        error_msg = "attempt to get gap only columns of unaligned msa"
        raise StandardError, error_msg, caller
      end
      total_sum = 0.0
      gap_sum = 0.0
      for c in 0 ... get_length
        for s in 0 ... get_number_of_seqs
          total_sum = total_sum + 1
          if Util.is_aa_gap_character?( get_sequence( s ).get_character_code( c ) )
            gap_sum = gap_sum  + 1
          end
        end

      end
      return gap_sum / total_sum
    end

    def get_gap_columns()
      if ( !is_aligned() )
        error_msg = "attempt to get gap columns of unaligned msa"
        raise StandardError, error_msg, caller
      end
      cols = Array.new()
      for c in 0 ... get_length
        gap_char_found = false
        for s in 0 ... get_number_of_seqs
          if Util.is_aa_gap_character?( get_sequence( s ).get_character_code( c ) )
            gap_char_found = true
            break
          end
        end
        if gap_char_found
          cols.push( c )
        end
      end
      return cols
    end

    # gap_ratio = seqs with gap / number of sequences
    # returns column indices for which seqs with gap / number of sequences > gap_ratio
    def get_gap_columns_w_gap_ratio( gap_ratio )
      if ( !is_aligned() )
        error_msg = "attempt to get gap columns with gap_ratio of unaligned msa"
        raise StandardError, error_msg, caller
      end
      if ( gap_ratio < 0 || gap_ratio > 1 )
        error_msg = "gap ratio must be between 0 and 1 inclusive"
        raise ArgumentError, error_msg, caller
      end
      cols = Array.new()
      for c in 0 ... get_length
        gap_chars_found = 0
        for s in 0 ... get_number_of_seqs
          if Util.is_aa_gap_character?( get_sequence( s ).get_character_code( c ) )
            gap_chars_found += 1
          end
        end
        if ( ( gap_chars_found.to_f / get_number_of_seqs ) > gap_ratio )
          cols.push( c )
        end
      end
      return cols
    end


    # Split an alignment into n alignemnts of equal size, except last one
    def split( n, verbose = false )
      if ( n < 2 || n > get_number_of_seqs )
        error_msg = "attempt to split into less than two or more than the number of sequences"
        raise StandardError, error_msg, caller
      end
      msas = Array.new()
      r = get_number_of_seqs % n
      x = get_number_of_seqs / n
      for i in 0 ... n
        msa = Msa.new()
        s = 0

        if ( ( r > 0 ) && ( i == ( n - 1 ) ) )
          y = x + r
          if ( verbose )
            puts( i.to_s + ": " + y.to_s )
          end
          for j in 0 ... y
            msa.add_sequence( get_sequence( ( i * x ) + j ) )
          end
        else
          if ( verbose )
            puts( i.to_s + ": " + x.to_s )
          end
          for j in 0 ... x
            msa.add_sequence( get_sequence( ( i * x ) + j ) )
          end
        end
        msas.push( msa )
      end
      msas
    end


    private

    def get_overlaps( min_overlap = 1 )
      if ( !is_aligned() )
        error_msg = "attempt to get overlaps of unaligned msa"
        raise StandardError, error_msg, caller
      end
      bins = Array.new()
      for i in 0 ... get_number_of_seqs()
        found_bin = false
        for j in 0 ... bins.length
          current_bin = bins[ j ]
          # does seq i overlap with all seqs in current_bin?
          all_overlap = true
          for z in 0 ... current_bin.length
            unless overlap?( i, current_bin[ z ], min_overlap )
              all_overlap = false
              break
            end
          end
          if all_overlap
            current_bin.push( i )
            found_bin = true
          end
        end
        unless found_bin
          new_bin = Array.new()
          new_bin.push( i )
          bins.push( new_bin )
        end
      end
      return bins
    end

    def remove_columns!( cols )
      if ( !is_aligned() )
        error_msg = "attempt to remove columns of unaligned msa"
        raise StandardError, error_msg, caller
      end
      cols.reverse!()
      for c in 0 ... cols.length()
        col = cols[ c ]
        for s in 0 ... get_number_of_seqs()
          get_sequence( s ).delete_residue!( col )
        end
      end
      return self
    end


  end # class Msa

end # module Evoruby


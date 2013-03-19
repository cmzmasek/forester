#
# = lib/evo/sequence/sequence.rb - Sequence class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: sequence.rb,v 1.10 2009/01/07 02:48:20 cmzmasek Exp $

require 'set'

module Evoruby

    class Sequence

        def initialize( name,
                molecular_sequence_str,
                accession = nil,
                accession_source = nil,
                taxonomy = nil,
                symbol = nil,
                secondary_accession = nil,
                secondary_accession_source = nil )
            @name               = String.new( name.strip() )
            @molecular_sequence = String.new( molecular_sequence_str )
            if ( accession == nil )
                @accession = String.new()
            else
                @accession = String.new( accession.strip() )
            end
            if ( accession_source == nil )
                @accession_source = String.new()
            else
                @accession_source = String.new( accession_source.strip() )
            end
            @taxonomy = taxonomy
            if ( symbol == nil )
                @symbol = String.new()
            else
                @symbol = String.new( symbol.strip() )
            end
            if ( secondary_accession == nil )
                @secondary_accession = String.new()
            else
                @secondary_accession = String.new( secondary_accession.strip() )
            end
            if ( secondary_accession_source == nil )
                @secondary_accession_source = String.new()
            else
                @secondary_accession_source = String.new( secondary_accession_source.strip() )
            end
        end

        def copy
            if get_taxonomy == nil
                Sequence.new( get_name, get_sequence_as_string, get_accession, get_accession_source, nil, get_symbol, get_secondary_accession, get_secondary_accession_source )
            else
                Sequence.new( get_name, get_sequence_as_string, get_accession, get_accession_source, get_taxonomy.copy, get_symbol, get_secondary_accession, get_secondary_accession_source )
            end
        end

        def get_name()
            @name
        end

        def set_name( name )
            @name = name
        end

        def get_sequence_as_string()
            @molecular_sequence
        end

        def get_accession()
            @accession
        end

        def get_accession_source()
            @accession_source
        end

        def get_secondary_accession()
            @secondary_accession
        end

        def get_secondary_accession_source()
            @secondary_accession_source
        end

        def get_symbol()
            @symbol
        end

        def get_taxonomy()
            @taxonomy
        end

        def get_length()
            @molecular_sequence.length
        end

        def get_residue( position )
            get_slice( position, 1 )
        end

        def get_character_code( position )
            @molecular_sequence.getbyte( position )
        end

        def get_gap_ratio()
            return get_gap_length().to_f / get_length()
        end

        def get_gap_length()
            counter = 0
            for i in 0 ... get_length()
                if ( Util.is_aa_gap_character?( get_character_code( i ) ) )
                    counter += 1
                end
            end
            return counter;
        end

        def delete_residue!( position )
            if ( position < 0 || position >= get_length() )
                error_msg = "attempt to delete residue at postion out of range"
                raise ArgumentError, error_msg
            end
            @molecular_sequence.slice!( position )
        end

        def get_slice( start, length )
            if ( start < 0 || start + length > get_length() )
                error_msg = "attempt to get sequence residue(s) at postion out of range"
                raise ArgumentError, error_msg
            end
            @molecular_sequence.slice( start, length )
        end

        def get_slice!( start, length )
            if ( start < 0 || start + length > get_length() )
                error_msg = "attempt to get sequence residue(s) at postion out of range"
                raise ArgumentError, error_msg
            end
            @molecular_sequence.slice!( start, length )
        end

        def get_subsequence( first, last )
            if ( last < first )
                error_msg = "attempt to get subsequence from " + first + " to " + last
                raise ArgumentError, error_msg
            end
            return Sequence.new( get_name, @molecular_sequence.slice( first, last - first + 1 ) )
        end

        def append!( molecular_sequence_str )
            @molecular_sequence.concat( molecular_sequence_str )
        end

        def to_str
            return "[" + @name + "] " + @molecular_sequence
        end

        def to_fasta
            return ">" + @name + Constants::LINE_DELIMITER  + @molecular_sequence
        end


    end # class Sequence

end # module Evoruby

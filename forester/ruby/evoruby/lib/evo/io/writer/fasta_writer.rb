#
# = lib/evo/io/writer/fasta_writer.rb - FastaWriter class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: fasta_writer.rb,v 1.6 2008/09/12 23:52:11 cmzmasek Exp $
#
# last modified: 05/16/2007

require 'lib/evo/io/writer/msa_writer'

module Evoruby

    class FastaWriter < MsaWriter

        LINE_WIDTH_DEFAULT      = 60
        MAX_NAME_LENGTH_DEFAULT = 0

        def initialize()
            @line_width          = LINE_WIDTH_DEFAULT
            @max_name_length     = MAX_NAME_LENGTH_DEFAULT
            @remove_gap_chars    = false
            @clean               = false
            @ex_if_name_too_long = false
        end


        def set_line_width( line_width = LINE_WIDTH_DEFAULT )
            if ( line_width < 1 )
                line_width = LINE_WIDTH_DEFAULT
            end
            @line_width = line_width
        end

        def set_max_name_length( length = MAX_NAME_LENGTH_DEFAULT )
            if ( length < 1 )
                length = MAX_NAME_LENGTH_DEFAULT
            end
            @max_name_length = length
        end

        def remove_gap_chars( remove_gap_chars = true )
            @remove_gap_chars = remove_gap_chars
        end

        def clean( clean = true )
            @clean = clean
        end

        def set_exception_if_name_too_long( exception_if_name_too_long )
          @ex_if_name_too_long = exception_if_name_too_long
        end

        def write( msa, path )
            Util.check_file_for_writability( path )
            f = File.open( path, "a" )
            for i in 0 ... msa.get_number_of_seqs()
                seq_obj = msa.get_sequence( i )
                name = seq_obj.get_name()
                f.print( ">" )
                if ( @max_name_length != MAX_NAME_LENGTH_DEFAULT )
                    name = Util.normalize_seq_name( name, @max_name_length, @ex_if_name_too_long )
                end
                f.print( name )
                counter = 0
                for j in 0 ... seq_obj.get_length()
                    unless @remove_gap_chars && Util.is_aa_gap_character?( seq_obj.get_character_code( j ) )
                        char = seq_obj.get_residue( j )
                        if ( @clean )
                            char = Util.clean_seq_str( char )
                            if ( char.length < 1 )
                                next
                            end
                        end
                        if counter % @line_width == 0
                            f.print( Evoruby::Constants::LINE_DELIMITER )
                        end
                        f.print( char )
                        counter += 1
                    end
                end
                f.print( Evoruby::Constants::LINE_DELIMITER )
            end
            f.close()
        end

    end # class FastaWriter

end # module Evoruby


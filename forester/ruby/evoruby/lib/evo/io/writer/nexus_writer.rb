#
# = lib/evo/io/writer/nexus_writer.rb - NexusWriter class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: nexus_writer.rb,v 1.4 2009/11/04 01:50:59 cmzmasek Exp $
#
# last modified: 05/16/2007

require 'lib/evo/io/writer/msa_writer'

module Evoruby

    class NexusWriter < MsaWriter

        MAX_NAME_LENGTH_DEFAULT = 10

        def initialize()
            @max_name_length     = MAX_NAME_LENGTH_DEFAULT
            @clean               = false
            @ex_if_name_too_long = false
        end

        def set_max_name_length( length = MAX_NAME_LENGTH_DEFAULT )
            if length < 1
                length = MAX_NAME_LENGTH_DEFAULT
            end
            @max_name_length = length
        end

        def clean( clean = true )
            @clean = clean
        end

        def set_exception_if_name_too_long( exception_if_name_too_long )
          @ex_if_name_too_long = exception_if_name_too_long
        end

        def write( msa, path )
            if ( !msa.is_aligned() )
                error_msg = "attempt to write unaligned msa in nexus format"
                raise StandardError, error_msg, caller
            end

            Util.check_file_for_writability( path )

            f = File.open( path, "a" )

            f.print( "Begin Data;" )
            f.print( Evoruby::Constants::LINE_DELIMITER )
            f.print( "   Dimensions NTax=" )
            f.print( msa.get_number_of_seqs().to_s() )
            f.print( " NChar=" )
            f.print( msa.get_length().to_s() )
            f.print( ";" )
            f.print( Evoruby::Constants::LINE_DELIMITER )
            f.print( "   Format DataType=Protein Interleave=No gap=-;" )
            f.print( Evoruby::Constants::LINE_DELIMITER )
            f.print( "   Matrix" )
            f.print( Evoruby::Constants::LINE_DELIMITER )
            for i in 0 ... msa.get_number_of_seqs()
                seq_obj = msa.get_sequence( i )
                name = seq_obj.get_name()
                seq  = seq_obj.get_sequence_as_string()
                name = name.gsub( /\s+$/, '')
                name = name.gsub( /\s+/, '_')
                name = Util.normalize_seq_name( name, @max_name_length, @ex_if_name_too_long )
                f.print( "      " )
                f.print( name )
                f.print( " " )
                if ( @clean )
                    seq = Util.clean_seq_str( seq )
                end
                f.print( seq )
                f.print( Evoruby::Constants::LINE_DELIMITER )
            end
            f.print( "   ;" )
            f.print( Evoruby::Constants::LINE_DELIMITER )
            f.print( "End;" )
            f.print( Evoruby::Constants::LINE_DELIMITER )
            f.close()
        end

    end # class NexusWriter

end # module Evoruby

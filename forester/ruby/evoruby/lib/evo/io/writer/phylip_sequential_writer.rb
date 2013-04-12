#
# = lib/evo/io/writer/phylip_sequential_writer.rb - PhylipSequentialWriter class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: phylip_sequential_writer.rb,v 1.4 2008/09/03 00:31:38 cmzmasek Exp $
#
# last modified: 05/16/2007

require 'lib/evo/io/writer/msa_writer'

module Evoruby

  class PhylipSequentialWriter < MsaWriter

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
        error_msg = "attempt to write unaligned msa in phylip sequential format"
        raise StandardError, error_msg, caller
      end


      Util.check_file_for_writability( path )

      f = File.open( path, "a" )

      f.print( msa.get_number_of_seqs().to_s() )
      f.print( " " )
      f.print( msa.get_length().to_s() )
      f.print( Evoruby::Constants::LINE_DELIMITER )
      for i in 0 ... msa.get_number_of_seqs()
        seq_obj = msa.get_sequence( i )
        name = seq_obj.get_name()
        seq  = seq_obj.get_sequence_as_string()
        name = name.gsub( /\s+$/, '')
        name = name.gsub( /\s+/, '_')
        name = Util.normalize_seq_name( name, @max_name_length, @ex_if_name_too_long )
        f.print( name )
        f.print( " " )
        if ( @clean )
          seq = Util.clean_seq_str( seq )
        end
        f.print( seq )
        f.print( Evoruby::Constants::LINE_DELIMITER )
      end
      f.close()
    end

  end # class PhylipSequentialWriter

end # module Evoruby

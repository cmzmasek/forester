# = lib/evo/table/basic_table.rb - BasicTable class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)'s
#
# $Id: basic_table.rb,v 1.3 2007/09/28 03:12:10 cmzmasek Exp $
#
# last modified: 05/16/2007

#require 'lib/evo/util/constants'

module Evoruby

  class BasicTable

        def initialize()
            @rows = Hash.new
            @max_row = 0
            @max_col = 0
        end

        # raises ArgumentError
        def set_value( row, col, value )
            if ( ( row < 0 ) || ( col < 0 ) )
                raise( ArgumentError, "attempt to use negative values for row or column" )
            end
            if ( row > get_max_row() )
                set_max_row( row )
            end
            if ( col > get_max_col() )
                set_max_col( col )
            end
            row_map = nil
            if ( @rows.has_key?( row ) )
                row_map = @rows[ row ]
            else
                row_map = Hash.new
                @rows[ row ] = row_map
            end
            row_map[ col ] = value
        end

        # raises ArgumentError
        def get_value_as_string( row, col )
            return ( get_value( row, col ) ).to_s
        end

        # raises ArgumentError
        def get_value( row, col )
            if ( ( row > get_max_row() ) || ( row < 0 ) )
                raise( ArgumentError, "value for row (" + row.to_s +
                         ") is out of range [max row: " + get_max_row().to_s + "]" )
            elsif ( ( col > get_max_col() ) || ( row < 0 ) )
                raise( ArgumentError, "value for column (" + col.to_s +
                         ") is out of range [max column: " + get_max_col().to_s + "]" )
            end
            row_map = @rows[ row ]
            if ( ( row_map == nil ) || ( row_map.length < 1 ) )
                return nil
            end
            return row_map[ col ]
        end

        def get_max_col()
            return @max_col
        end

        def get_max_row()
            return @max_row
        end

        # raises ArgumentError
        def get_columns_as_map( key_col, value_col )
            map = Hash.new
            for row in 0 .. get_max_row
                key = get_value( row, key_col )
                value = get_value( row, value_col )
                if ( ( key != nil ) && ( value != nil ) )
                    if ( map.has_key?( key ) )
                        raise( ArgumentError, "attempt to use non-unique table value as key [" +
                                        + key + "]" )
                    end
                    map[ key ] = value
                end
            end
            return map
        end

        def to_s
            str = String.new
            for row in 0 .. get_max_row
               for col in 0 .. get_max_col
                   str << col.to_s << " "
               end
               str << LEvoruby::Constants::LINE_DELIMITER
               for col in 0 .. get_max_col
                   str << row.to_s << ": "
                   str << get_value( row, col ) << " "
               end
               str << Evoruby::Constants::LINE_DELIMITER
            end
            return str
        end


        private

        def get_row( row )
            return @rows[ row ]
        end

        def set_max_col( max_col )
            @max_col = max_col
        end

        def set_max_row( max_row )
            @max_row = max_row
        end

    end # class BasicTable

end # module Evoruby

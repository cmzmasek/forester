#
# = lib/evo/io/parser/basic_table_parser - BasicTableParser class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: basic_table_parser.rb,v 1.3 2007/09/28 03:12:10 cmzmasek Exp $
#
# last modified: 05/16/2007

module Evoruby

    class BasicTableParser

        START_OF_COMMENT_LINE_CHAR = "#"

        # raises ArgumentError
        def BasicTableParser.parse( path, column_delimiter )
            Util.check_file_for_readability( path )
            table = BasicTable.new
            row = 0
            File.open( path ) do | file |
                while line = file.gets
                    if ( !Util.is_string_empty?( line ) &&
                         !line.slice( 0, 1 ).eql?( START_OF_COMMENT_LINE_CHAR ) )
                        values = line.split( column_delimiter )
                        col = 0
                        values.each { | value | 
                            table.set_value( row, col, value.strip! )
                            col += 1
                        }
                        row += 1
                    end
                end
            end
            return table
        end

    end # class BasicTableParser

end # module Evoruby

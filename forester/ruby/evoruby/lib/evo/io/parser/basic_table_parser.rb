# forester -- software libraries and applications
# for evolutionary biology and genomics.
# Copyright (C) 2026 Christian M. Zmasek
# All rights reserved
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact: czmasek at jcvi dot org

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

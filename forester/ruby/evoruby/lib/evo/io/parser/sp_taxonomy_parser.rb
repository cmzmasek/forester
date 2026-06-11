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

    require 'lib/evo/taxonomy/sp_taxonomy'

    class SpTaxonomyParser

        START_OF_COMMENT_LINE_CHAR = "#"

        # raises ArgumentError
        def SpTaxonomyParser.parse( path )
            Util.check_file_for_readability( path )
            row = 0
            sp_taxonomies = Array.new
            File.open( path ) do | file |
                while line = file.gets
                    row += 1
                    if !Util.is_string_empty?( line )
                        if line =~ /([A-Z0-9]{3,5})\s+[A-Z]\s+(\d+):\s+N=(.+)/
                            code = $1
                            id = $2
                            sci_name = $3
                            tax = SpTaxonomy.new(code, id, sci_name )
                            #puts tax.to_str
                            sp_taxonomies.push( tax )
                        end
                    end
                end
            end
            sp_taxonomies
        end
    end # class SpTaxonomyParser

end # module Evoruby
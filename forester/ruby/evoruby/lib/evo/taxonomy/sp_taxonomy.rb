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

    class SpTaxonomy

        attr_accessor :code, :id, :scientific_name, :common_name
        
        def initialize( code, id, scientific_name, common_name = nil )
            @code = String.new( code.strip() )
            @id = String.new( id.strip() )
            @scientific_name = String.new( scientific_name.strip() )
            if ( common_name == nil )
                @common_name = String.new()
            else
                @common_name = String.new( common_name.strip() )
            end
        end

        def copy
            return Taxonomy.new( code, id, scientific_name, common_name  )
        end

        def to_str()
            code + " " + id + ": N=" + scientific_name
        end

    end # class SpTaxonomy

end # module Evoruby
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

    class Taxonomy

        def initialize( name, id = nil, id_source = nil )
            @name = String.new( name.strip() )
            if ( id == nil )
                @id = String.new()
            else
                @id = String.new( id.strip() )
            end
            if ( id_source == nil )
                @id_source = String.new()
            else
                @id_source = String.new( id_source.strip() )
            end
        end

        def == ( taxonomy )
            if taxonomy == nil
                return false
            else
                return ( ( get_name == taxonomy.get_name ) &&
                     ( get_id == taxonomy.get_id ) &&
                     ( get_id_source == taxonomy.get_id_source ) )
            end
        end

        def copy
            return Taxonomy.new( get_name, get_id, get_id_source )
        end

        def get_name()
            @name
        end

        def get_id()
            @id
        end

        def get_id_source()
            @id_source
        end

        def to_str()
            if Util.is_string_empty?( get_id )
                @name
            else
                "[" + get_id + "] " + @name
            end
        end

    end # class Taxonomy

end # module Evoruby
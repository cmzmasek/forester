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

    class ProteinDomain

        def initialize( name, from, to, id, confidence )
            @name       = String.new( name )
            @from       = from
            @to         = to
            @id         = String.new( id )
            @confidence = confidence
        end

        def get_name()
            return @name
        end

        def get_from()
            return @from
        end

        def get_to()
            return @to
        end

        def get_id()
            return @id
        end

        def get_confidence()
            return @confidence
        end

    end # class ProteinDomain

end # module Evoruby

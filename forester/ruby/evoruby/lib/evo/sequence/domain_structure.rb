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

require 'lib/evo/util/constants'

module Evoruby

    class DomainStructure

        def initialize( total_length )
            @domains = Hash.new
            @total_length = total_length
        end

        def add_domain( domain, overwrite_if_same_from_to )
            key = domain.get_from
            if ( @domains.has_key?( key ) )
                prev_domain = @domains[ key ]
                if ( prev_domain.get_to == domain.get_to )
                    puts( "WARNING: more than one domain at the same location [" +
                        key.to_s + "-" + domain.get_to.to_s + "]: " + prev_domain.get_name + " and " + domain.get_name)
                    if ( overwrite_if_same_from_to )
                        puts( "         ignored the one with higher E-value [" +
                        prev_domain.get_confidence().to_s + " vs " + domain.get_confidence().to_s + "]" )
                        if prev_domain.get_confidence() < domain.get_confidence()
                            return # keep previous one
                        else
                            @domains[ key ] = domain
                            return
                        end
                    end
                end

                while ( @domains.has_key?( key ) )
                    key = key + 0.0001
                end

            end
            @domains[ key ] = domain
        end

        def to_NHX
            str = String.new
            str << get_total_length.to_s
            a = @domains.sort
            for d in a
                domain = d[ 1 ]
                str << Evoruby::Constants::DOMAIN_STRUCTURE_NHX_SEPARATOR
                str << domain.get_from.to_s
                str << Evoruby::Constants::DOMAIN_STRUCTURE_NHX_SEPARATOR
                str << domain.get_to.to_s
                str << Evoruby::Constants::DOMAIN_STRUCTURE_NHX_SEPARATOR
                str << domain.get_confidence.to_s
                str << Evoruby::Constants::DOMAIN_STRUCTURE_NHX_SEPARATOR
                str << domain.get_name
            end
            return str
        end

        def get_total_length
            return @total_length
        end

    end # class DomainStructure

end # module Evoruby

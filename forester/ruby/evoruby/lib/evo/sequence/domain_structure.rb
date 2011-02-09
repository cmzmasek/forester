#
# = lib/evo/sequence/domain_structure.rb - DomainStructure class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: domain_structure.rb,v 1.2 2007/06/12 04:51:33 cmzmasek Exp $
#
# last modified: 05/16/2007

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

#
# = lib/evo/sequence/protein_domain.rb - ProteinDomain class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: protein_domain.rb,v 1.2 2007/06/12 04:51:33 cmzmasek Exp $
#
# last modified: 05/16/2007

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

#
# = lib/evo/taxonomy/sp_taxonomy.rb - SpTaxonomy class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: sp_taxonomy.rb,v 1.1 2008/12/30 05:28:00 cmzmasek Exp $



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
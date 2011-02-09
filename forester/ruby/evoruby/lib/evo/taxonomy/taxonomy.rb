#
# = lib/evo/taxonomy/taxonomy.rb - Taxonomy class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: taxonomy.rb,v 1.2 2009/01/03 00:19:08 cmzmasek Exp $



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
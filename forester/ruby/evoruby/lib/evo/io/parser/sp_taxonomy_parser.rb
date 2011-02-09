#
# = lib/evo/io/parser/sp_taxonomy_parser - SpTaxonomyParser class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: sp_taxonomy_parser.rb,v 1.2 2008/12/31 03:21:45 cmzmasek Exp $


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
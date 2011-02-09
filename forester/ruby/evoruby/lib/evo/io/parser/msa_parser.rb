#
# = lib/evo/io/parser/msa_parser - MsaParser class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: msa_parser.rb,v 1.2 2007/06/12 04:51:34 cmzmasek Exp $
#
# last modified: 05/16/2007

module Evoruby

    class MsaParser
        def initialize()
            raise TypeError, "Cannot instanciate abstract class MsaParser"
        end

        def parse( path )
        end
    end

end # module Evoruby

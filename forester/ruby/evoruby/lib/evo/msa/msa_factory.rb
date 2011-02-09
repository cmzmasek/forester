#
# = lib/evo/msa/msa_factory.rb - MsaFactory class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: msa_factory.rb,v 1.2 2007/06/12 04:51:34 cmzmasek Exp $
#
# last modified: 05/16/2007

module Evoruby

    class MsaFactory

        def initialize
        end

        def create_msa_from_file( path, msa_parser )
            msa_parser.parse( path )
        end

    end # class MsaFactory

end # module Evoruby

#
# = lib/evo/io/msa_io.rb - MsaIO class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: msa_io.rb,v 1.2 2007/06/12 04:51:35 cmzmasek Exp $
#
# last modified: 05/16/2007

module Evoruby

  class MsaIO

        def initialize()
        end

        def write_to_file( msa, path, msa_writer )
            msa_writer.write( msa, path )
        end

  end # module Evoruby

end # class MsaIO

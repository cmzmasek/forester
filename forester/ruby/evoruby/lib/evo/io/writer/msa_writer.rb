#
# = lib/evo/io/writer/msa_writer.rb - MsaWriter class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: msa_writer.rb,v 1.2 2007/06/12 04:51:35 cmzmasek Exp $
#
# last modified: 05/16/2007

require 'lib/evo/util/constants'
require 'lib/evo/util/util'

module Evoruby

  class MsaWriter

    def initialize()
      raise TypeError, "Cannot instanciate abstract class MsaWriter"
    end

    def set_max_name_length( length )
    end

    def set_exception_if_name_too_long( exception_if_name_too_long )
    end

    def write( msa, path )
    end

  end # class MsaWriter

end # module Evoruby

#
# = lib/evo/util/constants.rb - Constants class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: constants.rb,v 1.3 2007/12/21 04:13:33 cmzmasek Exp $
#
# last modified: 05/11/2007

module Evoruby

    class Constants

        VERBOSE = true

        EVORUBY_VERSION = '1.0'

        FORESTER_HOME_ENV_VARIABLE = 'FORESTER_HOME'
        JAVA_HOME_ENV_VARIABLE     = 'JAVA_HOME'

        EVORUBY         = 'evoruby'

        LINE_DELIMITER  = "\n"

        FILE_SEPARATOR  = File::SEPARATOR

        DOMAIN_STRUCTURE_NHX_SEPARATOR = '>'


    end # class Constants

end # module Evoruby

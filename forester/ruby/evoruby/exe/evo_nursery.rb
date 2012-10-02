#!/usr/local/bin/ruby -W0
#
# = exe/evo_nursery
#
# Copyright::  Copyright (C) 2009-2010 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: evo_nursery.rb,v 1.1 2009/10/07 21:59:41 cmzmasek Exp $
#
# last modified: 2009/10/07

require 'lib/evo/tool/evo_nursery'

module Evoruby

    en = EvoNursery.new()

    en.run()

end  # module Evoruby
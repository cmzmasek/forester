#!/usr/local/bin/ruby -w
#
# = exe/msa_pro
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: msa_pro.rb,v 1.4 2008/08/28 17:09:06 cmzmasek Exp $
#


require 'lib/evo/tool/msa_processor'

module Evoruby

    mp = MsaProcessor.new()

    mp.run()

end  # module Evoruby

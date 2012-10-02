#!/usr/local/bin/ruby -w
#
# = exe/d2f
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: mse.rb,v 1.2 2008/08/28 17:09:06 cmzmasek Exp $


require 'lib/evo/tool/multi_sequence_extractor'

module Evoruby
    
    mse = MultiSequenceExtractor.new()
    
    mse.run()
    
end  # module Evoruby
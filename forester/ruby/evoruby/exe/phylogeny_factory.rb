#!/usr/local/bin/ruby -w
# 
# = exe/phylogeny_factory
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: phylogeny_factory.rb,v 1.5 2008/08/28 17:09:06 cmzmasek Exp $
#
 

require 'lib/evo/tools/phylogeny_factory'

module Evoruby

    pf = PhylogenyFactory.new

    pf.run

end  # module Evoruby

#!/usr/local/bin/ruby -w
#
# = exe/phylogeny_factory
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: phylogenies_decorator.rb,v 1.3 2008/08/28 17:09:06 cmzmasek Exp $
#
 

require 'lib/evo/tools/phylogenies_decorator'

module Evoruby

    fd = PhylogeniesDecorator.new

    fd.run

end  # module Evoruby

#!/usr/local/bin/ruby -w
#
# = exe/tap
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: tap.rb,v 1.3 2008/08/28 17:09:06 cmzmasek Exp $
#
# last modified: 05/18/2007

require 'lib/evo/tool/taxonomy_processor'

module Evoruby

    tap = TaxonomyProcessor.new()

    tap.run()

end  # module Evoruby

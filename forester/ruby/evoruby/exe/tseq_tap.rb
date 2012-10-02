#!/usr/local/bin/ruby -w
#
# = exe/tseq_tap
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: tseq_tap.rb,v 1.1 2008/12/31 06:00:08 cmzmasek Exp $


require 'lib/evo/tool/tseq_taxonomy_processor'

module Evoruby

    tap = TseqTaxonomyProcessor.new()

    tap.run()

end  # module Evoruby

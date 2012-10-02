#!/usr/local/bin/ruby -w
#
# = exe/fasta_tap
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: fasta_tap.rb,v 1.1 2009/01/20 20:44:54 cmzmasek Exp $


require 'lib/evo/tool/fasta_taxonomy_processor'

module Evoruby

    tap = FastaTaxonomyProcessor.new()

    tap.run()

end  # module Evoruby

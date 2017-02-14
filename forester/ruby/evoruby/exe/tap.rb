#!/usr/local/bin/ruby -w
#
# = exe/tap
#
# Copyright::    Copyright (C) 2017 Christian M Zmasek
# License::      GNU Lesser General Public License (LGPL)

require 'lib/evo/tool/taxonomy_processor'

module Evoruby

    tap = TaxonomyProcessor.new()

    tap.run()

end  # module Evoruby

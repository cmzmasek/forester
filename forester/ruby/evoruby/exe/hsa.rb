#!/usr/local/bin/ruby -w
#
# = exe/hsp
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)

require 'lib/evo/tool/hmmscan_analysis'

module Evoruby

    hsa = HmmscanAnalysis.new()

    hsa.run()

end  # module Evoruby

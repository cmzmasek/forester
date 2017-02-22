#!/usr/local/bin/ruby -w
#
# = exe/dsx
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)


require 'lib/evo/tool/multi_domain_seq_extractor'

module Evoruby

    mdsx = MultiDomainSeqExtractor.new()

    mdsx.run()

end  # module Evoruby

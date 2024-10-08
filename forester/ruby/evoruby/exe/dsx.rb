#!/usr/local/bin/ruby -w
#
# = exe/dsx
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)


require 'lib/evo/tool/domain_sequence_extractor'

module Evoruby

    dsx = DomainSequenceExtractor.new()

    dsx.run()

end  # module Evoruby

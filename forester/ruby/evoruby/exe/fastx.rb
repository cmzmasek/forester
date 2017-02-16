#!/usr/local/bin/ruby -w
#
# = exe/fastx
#
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)

require 'lib/evo/tool/fasta_extractor'

module Evoruby
    
    mse = FastaExtractor.new()
    
    mse.run()
    
end  # module Evoruby
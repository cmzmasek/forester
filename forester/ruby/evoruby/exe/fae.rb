#!/usr/local/bin/ruby -w
#
# = exe/fae
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: fae.rb,v 1.1 2008/09/10 02:16:34 cmzmasek Exp $


require 'lib/evo/tools/fasta_extractor'

module Evoruby
    
    mse = FastaExtractor.new()
    
    mse.run()
    
end  # module Evoruby
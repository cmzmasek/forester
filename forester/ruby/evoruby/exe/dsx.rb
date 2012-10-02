#!/usr/local/bin/ruby -w
#
# = exe/dsx
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: dsx.rb,v 1.3 2008/08/28 17:09:06 cmzmasek Exp $
#
# last modified: 06/11/2007

require 'lib/evo/tool/domain_sequence_extractor'

module Evoruby

    dsx = DomainSequenceExtractor.new()

    dsx.run()

end  # module Evoruby

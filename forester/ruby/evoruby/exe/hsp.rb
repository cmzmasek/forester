#!/usr/local/bin/ruby -w
#
# = exe/hsp
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: hsp.rb,v 1.1 2009/11/25 05:42:04 cmzmasek Exp $
#
# last modified: 11/24/2009

require 'lib/evo/tool/hmmscan_summary'

module Evoruby

    hsp = HmmscanSummary.new()

    hsp.run()

end  # module Evoruby

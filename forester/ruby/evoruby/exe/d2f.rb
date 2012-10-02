#!/usr/local/bin/ruby -w
#
# = exe/d2f
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: d2f.rb,v 1.3 2008/08/28 17:09:06 cmzmasek Exp $
#
# last modified: 06/11/2007

require 'lib/evo/tool/domains_to_forester'

module Evoruby

    dtf = DomainsToForester.new()

    dtf.run()

end  # module Evoruby

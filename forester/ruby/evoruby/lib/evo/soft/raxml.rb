#
# = lib/soft/raxml - Raxml class
#
# Copyright::  Copyright (C) 2009 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: raxml.rb,v 1.1 2009/10/07 00:08:35 cmzmasek Exp $
#
# last modified: 2009/10/06

require 'lib/evo/soft/resource_locations'
require 'lib/evo/util/util'

module Evoruby
    
    class Raxml 
        
        VERBOSE = true
        
        def initialize
            @fast_me_home = Util.get_env_variable_value( ResourceLocations::FASTME_HOME_ENV_VARIABLE )
            Util.check_file_for_readability( @fast_me_home )
        end
        
        def run( pwd_file, bootstrap_number, initial_tree ) 
            Util.check_file_for_readability( pwd_file )
            if bootstrap_number == nil || bootstrap_number < 0
                error_msg = "illegal bootstrap number: " + bootstrap_number
                raise ArgumentError, error_msg
            end
            if initial_tree == nil || (!initial_tree.eql?( "BME" ) && !initial_tree.eql?( "GME" ) && !initial_tree.eql?( "NJ" ) )
                error_msg = "illegal initial tree: " + initial_tree
                raise ArgumentError, error_msg
            end
            input = String.new()
            if bootstrap_number > 1 
                input = '-b #{initial_tree} -i #{pwd_file} -n #{bootstrap_number} -s b'
            else 
                input = '-b #{initial_tree} -i #{pwd_file} -s b'
            end
            if VERBOSE
                puts @fast_me_home + " " + input
            end
            IO.popen( @fast_me_home, 'r+' ) do |io|
                io.puts input
                io.close_write
                return io.read
            end
        end
    end # class Raxml 
    
end # module Evoruby
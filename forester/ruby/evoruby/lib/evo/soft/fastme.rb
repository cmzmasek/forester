#
# = lib/soft/fastme - FastMe class
#
# Copyright::  Copyright (C) 2009 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: fastme.rb,v 1.3 2009/10/08 22:44:54 cmzmasek Exp $
#
# last modified: 2009/10/06

require 'lib/evo/soft/resource_locations'
require 'lib/evo/util/util'

module Evoruby

    class FastMe

        VERBOSE = false
       
        OUTTREE      = 'output.tre'
        OUTPUT_D     = 'output.d'
        VERSION      = '2.0'
        
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
            init_tree_option = determine_initial_tree( initial_tree )
            input = String.new()
            if bootstrap_number > 1
                input = "-b #{init_tree_option} -i #{pwd_file} -n #{bootstrap_number} -s b"
            else
                input = "-b #{init_tree_option} -i #{pwd_file} -s b"
            end
            if VERBOSE
                puts @fast_me_home + " " + input
            end
            IO.popen( @fast_me_home + " " + input, 'r+' ) do |io|
                io.close_write
                return io.read
            end
        end

        private

        def determine_initial_tree( initial_tree )
            opt = nil
            if ( initial_tree == :BME )
                opt = "BME"
            elsif ( initial_tree == :GME )
                opt = "GME"
            elsif ( initial_tree == :NJ )
                opt = "NJ"
            else
                error_msg = "unknown initial tree"
                raise ArgumentError, error_msg
            end
            return opt
        end

    end # class FastMe

end # module Evoruby
#
# = lib/evo/util/command_line_arguments.rb - CommandLineArguments class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: command_line_arguments.rb,v 1.2 2007/06/12 04:51:34 cmzmasek Exp $
#
# last modified: 05/16/2007

module Evoruby

    class CommandLineArguments

        OPTIONS_PREFIX          = "-"
        EXTENDED_OPTIONS_PREFIX = "--"
        OPTIONS_SEPARATOR       = "="

        # raises ArgumentError
        def initialize( args )
            @options  = Hash.new
            @extended_options = Hash.new
            @file_names = Array.new
            parse_arguments( args )
        end

        def get_file_names
            return @file_names
        end

        def get_file_name( i )
            return @file_names[ i ]
        end

        def get_number_of_files()
            return @file_names.length
        end

        def is_option_set?( option_name )
            o = get_all_options
            return ( o.has_key?( option_name ) )
        end

        # raises ArgumentError
        def get_option_value( option_name )
            o = get_all_options
            if ( o.has_key?( option_name ) )
                value = o[ option_name ]
                if ( !Util.is_string_empty?( value ) )
                    return value
                else
                    raise( ArgumentError, "value for option \"" +
                         option_name + "\" is not set", caller )
                end
            else
                raise( ArgumentError, "option \"" + option_name +
                     "\" is not set", caller )
            end
        end

        def get_option_value_as_int( option_name )
            return get_option_value( option_name ).to_i
        end

        def get_option_value_as_float( option_name )
            return get_option_value( option_name ).to_f
        end

        # mandatory_options (Array)
        #
        def validate_mandatory_options( mandatory_options )
            o = get_all_options
            missing = Array.new
            for ma in mandatory_options
                if ( !o.has_key?( ma ) )
                    missing.push( ma )
                end
            end
            return missing
        end

        # mandatory_options (Array)
        #
        def validate_mandatory_options_as_str( mandatory_options )
            missing = validate_mandatory_options( mandatory_options )
            return missing.join( ", " )
        end

        # allowed_options (Array)
        #
        def validate_allowed_options( allowed_options )
            o = get_all_options
            disallowed = Array.new
            o.each_key { |op|
                if ( !allowed_options.include?( op ) )
                    disallowed.push( op )
                end
            }
            return disallowed
        end

        # allowed_options (Array)
        #
        def validate_allowed_options_as_str( allowed_options )
            disallowed = validate_allowed_options( allowed_options )
            return disallowed.join( ", " )
        end

        private

        def get_all_options
            o = Hash.new
            o.merge!( get_options_list )
            o.merge!( get_extended_options_list )
            return o
        end

        def parse_arguments( args )
            for arg in args
                if ( arg.index( EXTENDED_OPTIONS_PREFIX ) == 0 )
                    parse_option( arg.slice( EXTENDED_OPTIONS_PREFIX.length, arg.length() - 1 ),
                                  get_extended_options_list )

                elsif ( arg.index( OPTIONS_PREFIX ) == 0 )
                    parse_option( arg.slice( OPTIONS_PREFIX.length, arg.length() - 1 ),
                                  get_options_list )

                else
                    get_file_names.push( arg )
                end
            end
        end

        # raises ArgumentError
        def parse_option( option, options_map )
            sep_index = option.index( OPTIONS_SEPARATOR )
            if ( sep_index == nil )
                if ( Util.is_string_empty?( option ) )
                    raise( ArgumentError, "attempt to set option with an empty name" )
                end
                if ( get_all_options.has_key?( option ) )
                     raise( ArgumentError, "attempt to set option \"" +
                            option + "\" mutiple times" )
                end
                options_map[ option ] = ""
            else
                key = option.slice( 0, sep_index )
                value = option.slice( sep_index + 1, option.length() - 1 )
                if ( Util.is_string_empty?( key ) )
                    raise( ArgumentError, "attempt to set option with an empty name" )
                end
                if ( Util.is_string_empty?( value ) )
                    raise( ArgumentError, "attempt to set option with an empty value" )
                end
                if ( get_all_options.has_key?( key ) )
                    raise( ArgumentError, "attempt to set option \"" +
                            key + "\" mutiple times [" + option + "]" )
                end
                options_map[ key ] = value
            end
        end

        def get_file_names_list
            return @file_names
        end

        def get_options_list
            return @options
        end

        def get_extended_options_list
            return @extended_options
        end

    end # class CommandLineArguments

end # module Evoruby

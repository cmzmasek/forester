#
# = lib/evo/apps/fasta_extractor.rb - FastaExtractor class
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: fasta_extractor.rb,v 1.2 2010/12/13 19:00:11 cmzmasek Exp $


require 'lib/evo/util/util'
require 'lib/evo/util/constants'
require 'lib/evo/util/command_line_arguments'


module Evoruby

    class FastaExtractor

        PRG_NAME                           = "fae"
        PRG_VERSION                        = "1.0.0"
        PRG_DESC                           = "extraction of nucleotide sequences from a fasta file by names from wublast search"
        PRG_DATE                           = "2008.08.09"
        COPYRIGHT                          = "2008-2009 Christian M Zmasek"
        CONTACT                            = "phylosoft@gmail.com"
        WWW                                = "www.phylosoft.org"
        HELP_OPTION_1                      = 'help'
        HELP_OPTION_2                      = 'h'


        def run()

            Util.print_program_information( PRG_NAME,
                PRG_VERSION,
                PRG_DESC ,
                PRG_DATE,
                COPYRIGHT,
                CONTACT,
                WWW,
                STDOUT )

            ld = Constants::LINE_DELIMITER

            begin
                cla = CommandLineArguments.new( ARGV )
            rescue ArgumentError => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s )
            end

            if ( cla.is_option_set?( HELP_OPTION_1 ) ||
                     cla.is_option_set?( HELP_OPTION_2 ) )
                print_help
                exit( 0 )
            end

            if ( cla.get_number_of_files != 3 )
                print_help
                exit( -1 )
            end

            allowed_opts = Array.new

            disallowed = cla.validate_allowed_options_as_str( allowed_opts )
            if ( disallowed.length > 0 )
                Util.fatal_error( PRG_NAME,
                    "unknown option(s): " + disallowed,
                    STDOUT )
            end

            input_file  = cla.get_file_name( 0 )
            names_file  = cla.get_file_name( 1 )
            output_file = cla.get_file_name( 2 )

            if  !File.exist?( input_file )
                Util.fatal_error( PRG_NAME, "error: input file [#{input_file}] does not exist" )
            end
            if  !File.exist?( names_file )
                Util.fatal_error( PRG_NAME, "error: names file [#{names_file}] does not exist" )
            end
            if File.exist?( output_file   )
                Util.fatal_error( PRG_NAME, "error: [#{output_file }] already exists" )
            end

            names = extract_names_with_frames( names_file )

            extract_sequences( names, input_file, output_file )

            puts
            Util.print_message( PRG_NAME, "OK" )
            puts

        end


        def extract_names_with_frames( names_file )
            names = Hash.new()
            File.open( names_file ) do | file |
                while line = file.gets
                    if ( !Util.is_string_empty?( line ) && !(line =~ /\s*#/ ) )
                        if ( line =~ /(\S+)\s+([+|-]\d)\s+\d+\s+(\S+)/ )
                            name  = $1
                            frame = $2
                            e     = $3
                            names[ name ] =  "[" + frame + "] [" + e + "]"
                        end
                    end
                end
            end
            names
        end

        def extract_sequences( names, fasta_file, output_file )
            output = File.open( output_file, "a" )
            matching_state = false
            counter = 0
            File.open( fasta_file ) do | file |
                while line = file.gets
                    if !Util.is_string_empty?( line )
                        if ( line =~ /\s*>\s*(.+)/ )
                            name = $1
                            if names.has_key?( name )
                                matching_state = true
                                counter += 1
                                puts counter.to_s + ". " +name + " " + names[ name ]
                                output.print( ">" + name + " " + names[ name ] )
                                output.print( Evoruby::Constants::LINE_DELIMITER )
                            else
                                matching_state = false
                            end
                        elsif matching_state
                            output.print( line )
                        end
                    end
                end
            end
            output.close()
        end

        def print_help()
            puts( "Usage:" )
            puts()
            puts( "  " + PRG_NAME + ".rb <input fasta file> <names file based on blast output> <output file>" )
            puts()
        end

    end # class FastaExtractor
end
#
# = lib/evo/apps/ -  class
#
# Copyright::  Copyright (C) 2009 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: new_tap.rb,v 1.4 2010/12/13 19:00:11 cmzmasek Exp $


require 'lib/evo/util/util'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/msa/msa'
require 'lib/evo/io/msa_io'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/io/parser/general_msa_parser'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/writer/phylip_sequential_writer'
require 'lib/evo/util/command_line_arguments'

module Evoruby

    class TaxonomyProcessor

        PRG_NAME       = ""
        PRG_DATE       = "2009.10.09"
        PRG_DESC       = "replacement of labels in multiple sequence files"
        PRG_VERSION    = "1.00"
        COPYRIGHT      = "2009 Christian M Zmasek"
        CONTACT        = "phylosoft@gmail.com"
        WWW            = "www.phylosoft.org"

        REMOVE_REDUNDANT_SEQS_OPTION = "rr"
        
        def initialize()
            @taxonomies = Hash.new()
        end

        def run()

            Util.print_program_information( PRG_NAME,
                PRG_VERSION,
                PRG_DESC,
                PRG_DATE,
                COPYRIGHT,
                CONTACT,
                WWW,
                STDOUT )

            if ( ARGV == nil || ( ARGV.length != 3 && ARGV.length != 4 ) )
                puts( "Usage: #{PRG_NAME}.rb <input sequences> <output sequences> <output map>" )
                puts()
                puts( "  options: -" + REMOVE_REDUNDANT_SEQS_OPTION + ": to remove redundant sequences" )
                puts()
                exit( -1 )
            end

            begin
                cla = CommandLineArguments.new( ARGV )
            rescue ArgumentError => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s )
            end
            
            input     = cla.get_file_name( 0 )
            output    = cla.get_file_name( 1 )
            map_file = cla.get_file_name( 2 )

            allowed_opts = Array.new
            allowed_opts.push( REMOVE_REDUNDANT_SEQS_OPTION ) 
            
            disallowed = cla.validate_allowed_options_as_str( allowed_opts )
            if ( disallowed.length > 0 )
                Util.fatal_error( PRG_NAME, "unknown option(s): " + disallowed )
            end

            
            remove_redudant = false
            if ( cla.is_option_set?( REMOVE_REDUNDANT_SEQS_OPTION ) )
                remove_redudant = true
            end

            if ( File.exists?( output ) )
                Util.fatal_error( PRG_NAME, "outfile [" + output + "] already exists" )
            end
            if ( File.exists?( map_file ) )
                Util.fatal_error( PRG_NAME, "map file [" + map_file + "] already exists" )
            end
            if ( !File.exists?( input) )
                Util.fatal_error( PRG_NAME, "infile [" + input + "] does not exist" )
            end
           
            fasta_like = Util.looks_like_fasta?( input )

            puts()
            puts( "Input alignment : " + input )
            puts( "Output alignment: " + output )
            puts( "Output map      : " + map_file )
            if ( fasta_like )
                puts( "Format          : Fasta"  )
            else
                puts( "Format          : Phylip like" )
            end
            puts()

            species_map = Hash.new
           
            f = MsaFactory.new()
            begin
                if ( fasta_like )
                    msa = f.create_msa_from_file( input, FastaParser.new() )
                else
                    msa = f.create_msa_from_file( input, GeneralMsaParser.new() )
                end
            rescue Exception => e
                Util.fatal_error( PRG_NAME, "failed to read file: " + e.to_s )
            end

            if ( msa == nil || msa.get_number_of_seqs() < 1 )
                Util.fatal_error( PRG_NAME, "failed to read MSA" )
            end
            begin
                Util.check_file_for_writability( map_file )
            rescue Exception => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_, STDOUT )
            end

            if ( remove_redudant ) 
                removed = msa.remove_redundant_sequences!( true )
                if removed.size > 0
                    Util.print_message( PRG_NAME, "going to ignore the following " + removed.size.to_s + " redundant sequences:" )
                    removed.each { | seq_name |
                        puts seq_name
                    }
                    Util.print_message( PRG_NAME, "will process " + msa.get_number_of_seqs.to_s + " non redundant sequences" )
                end
            end

            lf = File.open( map_file, "a" )
            for i in 0 ... msa.get_number_of_seqs
                seq  = msa.get_sequence( i )
            end

            io = MsaIO.new()
            w = nil
            if ( fasta_like )
                w = FastaWriter.new()
            else
                w = PhylipSequentialWriter.new()
            end
            w.set_max_name_length( 10 )
            w.clean( true )
            begin
                io.write_to_file( msa, output, w )
            rescue Exception => e
                Util.fatal_error( PRG_NAME, "failed to write file: " + e.to_s )
            end
            lf.close()
            if ( @taxonomies.length > 0 )
                Util.print_message( PRG_NAME, "number of unique taxonomies: " + @taxonomies.length.to_s )
            end
            Util.print_message( PRG_NAME, "wrote: " + map_file )
            Util.print_message( PRG_NAME, "wrote: " + output )
            Util.print_message( PRG_NAME, "OK" )
        end

    end # class 

end # module Evoruby
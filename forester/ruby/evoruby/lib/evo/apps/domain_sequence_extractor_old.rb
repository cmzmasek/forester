#
# = lib/evo/apps/domain_sequence_extractor.rb - DomainSequenceExtractor class
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: domain_sequence_extractor.rb,v 1.19 2010/12/13 19:00:11 cmzmasek Exp $


require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/io/parser/hmmsearch_domain_extractor'

module Evoruby

    class DomainSequenceExtractor

        PRG_NAME       = "dsx"
        PRG_VERSION    = "1.1.0"
        PRG_DESC       = "extraction of domain sequences from hmmsearch output"
        PRG_DATE       = "2008.01.03"
        COPYRIGHT      = "2008-2009 Christian M Zmasek"
        CONTACT        = "phylosoft@gmail.com"
        WWW            = "www.phylosoft.org"

        E_VALUE_THRESHOLD_OPTION           = 'e'
        LENGTH_THRESHOLD_OPTION            = 'l'
        ADD_POSITION_OPTION                = 'p'
        ADD_DOMAIN_NUMBER_OPTION           = 'd'
        ADD_DOMAIN_NUMBER_OPTION_AS_DIGIT  = 'dd'
        ADD_DOMAIN_NUMBER_OPTION_AS_LETTER = 'dl'
        TRIM_OPTION                        = 't'
        LOG_FILE_SUFFIX                    = '_domain_seq_extr.log'
        PASSED_SEQS_SUFFIX                 = '_domain_seq_extr_passed'
        FAILED_SEQS_SUFFIX                 = '_domain_seq_extr_failed'
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
            rescue ArgumentError
                Util.fatal_error( PRG_NAME, "error: " + $!, STDOUT )
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
            allowed_opts.push( E_VALUE_THRESHOLD_OPTION )
            allowed_opts.push( ADD_POSITION_OPTION )
            allowed_opts.push( ADD_DOMAIN_NUMBER_OPTION )
            allowed_opts.push( LENGTH_THRESHOLD_OPTION )
            allowed_opts.push( ADD_DOMAIN_NUMBER_OPTION_AS_DIGIT )
            allowed_opts.push( ADD_DOMAIN_NUMBER_OPTION_AS_LETTER )
            allowed_opts.push( TRIM_OPTION )

            disallowed = cla.validate_allowed_options_as_str( allowed_opts )
            if ( disallowed.length > 0 )
                Util.fatal_error( PRG_NAME,
                    "unknown option(s): " + disallowed,
                    STDOUT )
            end

            hmmsearch_output    = cla.get_file_name( 0 )
            fasta_sequence_file = cla.get_file_name( 1 )
            outfile             = cla.get_file_name( 2 )

            add_position = false
            if ( cla.is_option_set?( ADD_POSITION_OPTION ) )
                add_position = true
            end

            trim = false
            if ( cla.is_option_set?( TRIM_OPTION ) )
                trim = true
            end

            add_domain_number           = false
            add_domain_number_as_letter = false
            add_domain_number_as_digit  = false

            if ( cla.is_option_set?( ADD_DOMAIN_NUMBER_OPTION ) )
                add_domain_number = true
            end
            if ( cla.is_option_set?( ADD_DOMAIN_NUMBER_OPTION_AS_LETTER ) )
                add_domain_number_as_letter = true
            end
            if ( cla.is_option_set?( ADD_DOMAIN_NUMBER_OPTION_AS_DIGIT ) )
                add_domain_number_as_digit = true
            end

            if ( add_domain_number_as_letter && add_domain_number_as_digit )
                puts( "attempt to add domain number as letter and digit at the same time" )
                print_help
                exit( -1 )
            end

            e_value_threshold = -1.0
            if ( cla.is_option_set?( E_VALUE_THRESHOLD_OPTION ) )
                begin
                    e_value_threshold = cla.get_option_value_as_float( E_VALUE_THRESHOLD_OPTION )
                rescue ArgumentError => e
                    Forester::Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
                end
                if ( e_value_threshold < 0.0 )
                    Forester::Util.fatal_error( PRG_NAME, "attempt to use a negative E-value threshold", STDOUT )
                end
            end

            length_threshold = -1
            if ( cla.is_option_set?( LENGTH_THRESHOLD_OPTION ) )
                begin
                    length_threshold = cla.get_option_value_as_int( LENGTH_THRESHOLD_OPTION )
                rescue ArgumentError => e
                    Forester::Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
                end
                if ( length_threshold < 0)
                    Forester::Util.fatal_error( PRG_NAME, "attempt to use a negative length threshold", STDOUT )
                end
            end

            log = String.new

            puts()
            puts( "Hmmsearch outputfile                   : " + hmmsearch_output )
            log << "Hmmsearch outputfile                   : " + hmmsearch_output + ld
            puts( "Fasta sequencefile (complete sequences): " + fasta_sequence_file )
            log << "Fasta sequencefile (complete sequences): " + fasta_sequence_file + ld
            puts( "Outputfile                             : " + outfile )
            log << "Outputfile                             : " + outfile + ld
            puts( "Passed sequences outfile (fasta)       : " + outfile + PASSED_SEQS_SUFFIX )
            log << "Passed sequences outfile (fasta)       : " + outfile + PASSED_SEQS_SUFFIX + ld
            puts( "Failed sequences outfile (fasta)       : " + outfile + FAILED_SEQS_SUFFIX )
            log << "Failed sequences outfile (fasta)       : " + outfile + FAILED_SEQS_SUFFIX + ld
            puts( "Logfile                                : " + outfile + LOG_FILE_SUFFIX )
            log <<  "Logfile                                : " + outfile + LOG_FILE_SUFFIX + ld
            if ( e_value_threshold >= 0.0 )
                puts( "E-value threshold : " + e_value_threshold.to_s )
                log << "E-value threshold : " + e_value_threshold.to_s + ld
            else
                puts( "E-value threshold : no threshold" )
                log << "E-value threshold : no threshold" + ld
            end
            if ( length_threshold > 0 )
                puts( "Length threshold  : " + length_threshold.to_s )
                log << "Length threshold  : " + length_threshold.to_s + ld
            else
                puts( "Length threshold  : no threshold" )
                log << "Length threshold  : no threshold" + ld
            end

            if ( trim )
                puts( "Trim last 2 chars : true" )
                log << "Trim last 2 chars : true" + ld
            else
                puts( "Trim names        : false" )
                log << "Trim names        : false" + ld
            end


            if ( add_position )
                puts( "Add positions (rel to complete seq) to extracted domains: true" )
                log << "Add positions (rel to complete seq) to extracted domains: true" + ld
            else
                puts( "Add positions (rel to complete seq) to extracted domains: false" )
                log << "Add positions (rel to complete seq) to extracted domains: false" + ld
            end

            if ( add_domain_number || add_domain_number_as_digit || add_domain_number_as_letter )
                puts( "Add numbers to extracted domains (in case of more than one domain per complete seq): true" )
                log << "Add numbers to extracted domains (in case of more than one domain per complete seq): true" + ld
            else
                puts( "Add numbers to extracted domains (in case of more than one domain per complete seq): false" )
                log << "Add numbers to extracted domains (in case of more than one domain per complete seq): false" + ld
            end

            puts

            domain_count = 0
            begin
                parser = HmmsearchDomainExtractor.new()
                domain_count = parser.parse( hmmsearch_output,
                    fasta_sequence_file,
                    outfile,
                    outfile + PASSED_SEQS_SUFFIX,
                    outfile + FAILED_SEQS_SUFFIX,
                    e_value_threshold,
                    length_threshold,
                    add_position,
                    add_domain_number,
                    add_domain_number_as_digit,
                    add_domain_number_as_letter,
                    trim,
                    log )
            rescue ArgumentError, IOError, StandardError => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
            rescue Exception => e
                Util.fatal_error( PRG_NAME, "unexpected exception!: " + e.to_s, STDOUT )
            end

            puts
            Util.print_message( PRG_NAME, "extracted a total of " + domain_count.to_s + " domains" )
            Util.print_message( PRG_NAME, "wrote;               " + outfile )
            Util.print_message( PRG_NAME, "wrote:               " + outfile + LOG_FILE_SUFFIX )
            Util.print_message( PRG_NAME, "(wrote:              " + outfile + PASSED_SEQS_SUFFIX + ")" )
            Util.print_message( PRG_NAME, "(wrote:              " + outfile + FAILED_SEQS_SUFFIX + ")" )

            begin
                f = File.open( outfile + LOG_FILE_SUFFIX, 'a' )
                f.print( log )
                f.close
            rescue Exception => e
                Util.fatal_error( PRG_NAME, "error: " + e.to_s )
            end

            puts
            Util.print_message( PRG_NAME, "OK" )
            puts

        end

        def print_help()
            puts()
            puts( "Usage:" )
            puts()
            puts( "  " + PRG_NAME + ".rb [options] <hmmsearch outputfile> <file containing complete sequences in fasta format> <outputfile>" )
            puts()
            puts( "  options: -" + E_VALUE_THRESHOLD_OPTION  + "=<f>: E-value threshold, default is no threshold" )
            puts( "           -" + LENGTH_THRESHOLD_OPTION   + "=<i>: length threshold, default is no threshold" )
            puts( "           -" + ADD_POSITION_OPTION  + ": to add positions (rel to complete seq) to extracted domains" )
            puts( "           -" + ADD_DOMAIN_NUMBER_OPTION  + ": to add numbers to extracted domains (in case of more than one domain per complete seq) (example \"domain~2-3\")" )
            puts( "           -" + ADD_DOMAIN_NUMBER_OPTION_AS_DIGIT  + ": to add numbers to extracted domains as digit (example \"domain2\")" )
            puts( "           -" + ADD_DOMAIN_NUMBER_OPTION_AS_LETTER  + ": to add numbers to extracted domains as letter (example \"domaina\")" )
            puts( "           -" + TRIM_OPTION  + ": to remove the last 2 characters from sequence names" )
            puts()
        end

    end # class DomainSequenceExtractor

end # module Evoruby

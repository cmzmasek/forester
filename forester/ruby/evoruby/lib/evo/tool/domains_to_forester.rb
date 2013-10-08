#
# = lib/evo/apps/domains_to_forester - DomainsToForester class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: Exp $
#
# last modified: 06/11/2007

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/parser/fasta_parser'
require 'lib/evo/sequence/protein_domain'
require 'lib/evo/sequence/domain_structure'

module Evoruby

  class DomainsToForester

    PRG_NAME       = "d2f"
    PRG_DESC       = "parsed hmmpfam output to forester format"
    PRG_VERSION    = "1.001"
    PRG_DATE       = "20120807"
    COPYRIGHT      = "2012 Christian M Zmasek"
    CONTACT        = "phylosoft@gmail.com"
    WWW            = "www.phylosoft.org"

    E_VALUE_THRESHOLD_OPTION         = "e"
    OVERWRITE_IF_SAME_FROM_TO_OPTION = "o"
    HELP_OPTION_1                    = "help"
    HELP_OPTION_2                    = "h"

    def parse( domains_list_file,
        original_seqs_file,
        outfile,
        column_delimiter,
        e_value_threshold,
        overwrite_if_same_from_to )
      Util.check_file_for_readability( domains_list_file )
      Util.check_file_for_readability( original_seqs_file )
      Util.check_file_for_writability( outfile )

      domain_structures = Hash.new() # protein name is key, domain structure is value

      f = MsaFactory.new

      original_seqs = f.create_msa_from_file( original_seqs_file, FastaParser.new )
      if ( original_seqs.get_number_of_seqs < 1 )
        error_msg = "\"" + original_seqs_file + "\" appears devoid of sequences in fasta-format"
        raise ArgumentError, error_msg
      end

      File.open( domains_list_file ) do | file |
        while line = file.gets
          if !is_ignorable?( line )
            
            a = line.split( column_delimiter )
            l = a.length
            if ( ( l < 4 ) || ( e_value_threshold >= 0.0 && l < 5 ) )
              error_msg = "unexpected format at line: " + line
              raise IOError, error_msg
            end
            protein_name = a[ 0 ]
            domain_name  = a[ 1 ]
            seq_from     = -1
            seq_to       = -1
            ##########################################
            if domain_name =~ /RRM_\d/
              puts "ignoring " + line 
              next
            end
            ##########################################
            
            
            begin
              seq_from = a[ 2 ].to_i
            rescue Exception
              error_msg = "failed to parse seq from from \"" + a[ 2 ] + "\" [line: " + line + "]"
              raise IOError, error_msg
            end
            begin
              seq_to = a[ 3 ].to_i
            rescue Exception
              error_msg = "failed to parse seq to from \"" + a[ 3 ] + "\" [line: " + line + "]"
              raise IOError, error_msg
            end

            e_value = -1
            if l > 4
              begin
                e_value = a[ 4 ].to_f
              rescue Exception
                error_msg = "failed to parse E-value from \"" + a[ 4 ] + "\" [line: " + line + "]"
                raise IOError, error_msg
              end
            end

            seq = original_seqs.get_by_name_start( protein_name )

            total_length = seq.get_length

            if ( ( ( e_value_threshold < 0.0 ) || ( e_value <= e_value_threshold ) )  )
              pd = ProteinDomain.new( domain_name, seq_from, seq_to, "", e_value )
              ds = nil
              if ( domain_structures.has_key?( protein_name ) )
                ds = domain_structures[ protein_name ]
              else
                ds = DomainStructure.new( total_length )
                domain_structures[ protein_name ] = ds
              end
              ds.add_domain( pd, overwrite_if_same_from_to )
            end

          end
        end
      end

      out = File.open( outfile, "a" )
      ds = domain_structures.sort
      for d in ds
        protein_name     = d[ 0 ]
        domain_structure = d[ 1 ]
        out.print( protein_name.to_s )
        out.print( "\t" )
        out.print( domain_structure.to_NHX )
        out.print( Constants::LINE_DELIMITER  )
      end

      out.flush()
      out.close()

    end # parse




    def run()

      Util.print_program_information( PRG_NAME,
        PRG_VERSION,
        PRG_DESC,
        PRG_DATE,
        COPYRIGHT,
        CONTACT,
        WWW,
        STDOUT )

      begin
        cla = CommandLineArguments.new( ARGV )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      end

      if ( cla.is_option_set?( HELP_OPTION_1 ) ||
           cla.is_option_set?( HELP_OPTION_2 ) )
        print_help
        exit( 0 )
      end

      if cla.get_number_of_files != 3
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push( E_VALUE_THRESHOLD_OPTION )
      allowed_opts.push( OVERWRITE_IF_SAME_FROM_TO_OPTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
          "unknown option(s): " + disallowed,
          STDOUT )
      end

      domains_list_file       = cla.get_file_name( 0 )
      original_sequences_file = cla.get_file_name( 1 )
      outfile                 = cla.get_file_name( 2 )


      e_value_threshold = -1.0
      if cla.is_option_set?( E_VALUE_THRESHOLD_OPTION )
        begin
          e_value_threshold = cla.get_option_value_as_float( E_VALUE_THRESHOLD_OPTION )
        rescue ArgumentError => e
          Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
        end
        if ( e_value_threshold < 0.0 )
          Util.fatal_error( PRG_NAME, "attempt to use a negative E-value threshold", STDOUT )
        end
      end
      overwrite_if_same_from_to = false
      if ( cla.is_option_set?( OVERWRITE_IF_SAME_FROM_TO_OPTION ) )
        overwrite_if_same_from_to = true
      end

      puts
      puts( "Domains list file                      : " + domains_list_file )
      puts( "Fasta sequencefile (complete sequences): " + original_sequences_file )
      puts( "Outputfile                             : " + outfile )
      if ( e_value_threshold >= 0.0 )
        puts( "E-value threshold                      : " + e_value_threshold.to_s )
      else
        puts( "E-value threshold                      : no threshold" )
      end
      if ( overwrite_if_same_from_to )
        puts( "Overwrite if same from and to          : true" )
      else
        puts( "Overwrite if same from and to          : false" )
      end

      puts

      begin
        parse( domains_list_file,
          original_sequences_file,
          outfile,
          " ",
          e_value_threshold,
          overwrite_if_same_from_to )

      rescue ArgumentError, IOError, StandardError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s, STDOUT )
      rescue Exception => e
        Util.fatal_error( PRG_NAME, "unexpected exception: " + e.to_s, STDOUT )
      end


      puts
      Util.print_message( PRG_NAME, 'OK' )
      puts

    end

    private

    def print_help()
      puts
      puts( "Usage:" )
      puts
      puts( "  " + PRG_NAME + ".rb [options] <domains list file (parsed hmmpfam output)> <file containing complete sequences in fasta format> <outputfile>" )
      puts()
      puts( "  options: -" + E_VALUE_THRESHOLD_OPTION  + "=<f> : E-value threshold, default is no threshold" )
      puts( "               -" + OVERWRITE_IF_SAME_FROM_TO_OPTION  + " : overwrite domain with same start and end with domain with better E-value" )
      puts
    end



    def is_ignorable?( line )
      return ( line !~ /[A-Za-z0-9-]/ || line =~ /^\s*#/)
    end


  end # class DomainsToForester


end # module Evoruby

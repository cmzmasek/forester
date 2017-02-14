#
# = lib/evo/apps/domains_to_forester - DomainsToForester class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)

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
    PRG_DESC       = "converting of parsed hmmpfam output to forester format"
    PRG_VERSION    = "1.002"
    PRG_DATE       = "20170213"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

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
      WWW,
      STDOUT )

      if ( ARGV == nil || ( ARGV.length < 1 )  )
        print_help
        exit( -1 )
      end

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

      unless ( cla.get_number_of_files == 1 || cla.get_number_of_files == 2 || cla.get_number_of_files == 3 )
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

      domains_list_file = cla.get_file_name( 0 )
      original_sequences_file = ""
      outfile = ""
      if (cla.get_number_of_files == 3)
        original_sequences_file = cla.get_file_name( 1 )
        outfile = cla.get_file_name( 2 )
      elsif (cla.get_number_of_files == 1 || cla.get_number_of_files == 2 )
        if ( cla.get_number_of_files == 2 )
          original_sequences_file = cla.get_file_name( 1 )
        else
          hmmscan_index = domains_list_file.index("hmmscan")
          if ( hmmscan_index != nil )
            prefix = domains_list_file[0 .. hmmscan_index-1 ]
            suffix = Constants::ID_NORMALIZED_FASTA_FILE_SUFFIX
            files = Dir.entries( "." )
            matching_files = Util.get_matching_files( files, prefix, suffix)
            if matching_files.length < 1
              Util.fatal_error( PRG_NAME, 'no file matching [' + prefix +
              '...' + suffix + '] present in current directory: need to indicate <file containing complete sequences in fasta format> as second argument' )
            end
            if matching_files.length > 1
              Util.fatal_error( PRG_NAME, 'more than one file matching [' +
              prefix  + '...' + suffix + '] present in current directory: need to indicate <file containing complete sequences in fasta format> as second argument' )
            end
            original_sequences_file = matching_files[ 0 ]
          end
        end
        outfile = domains_list_file
        if (outfile.end_with?(Constants::DOMAIN_TABLE_SUFFIX) )
          outfile = outfile.chomp(Constants::DOMAIN_TABLE_SUFFIX)
        end
        if ( e_value_threshold >= 0.0 )
          outfile = outfile + Constants::DOMAINS_TO_FORESTER_EVALUE_CUTOFF_SUFFIX + e_value_threshold.to_s
        end
        outfile = outfile + Constants::DOMAINS_TO_FORESTER_OUTFILE_SUFFIX
      end

      overwrite_if_same_from_to = false
      if ( cla.is_option_set?( OVERWRITE_IF_SAME_FROM_TO_OPTION ) )
        overwrite_if_same_from_to = true
      end

      puts
      puts( "Domain table                            : " + domains_list_file )
      puts( "Fasta sequence file (complete sequences): " + original_sequences_file )
      puts( "Outputfile                              : " + outfile )
      if ( e_value_threshold >= 0.0 )
        puts( "E-value threshold                       : " + e_value_threshold.to_s )
      else
        puts( "E-value threshold                       : no threshold" )
      end
      if ( overwrite_if_same_from_to )
        puts( "Overwrite if same from and to           : true" )
      else
        puts( "Overwrite if same from and to           : false" )
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
      Util.print_message( PRG_NAME, "wrote: " + outfile )
      Util.print_message( PRG_NAME, "next steps in standard analysis pipeline: hmmsearch followed by dsx.rb")
      Util.print_message( PRG_NAME, 'OK' )
      puts

    end

    private

    def print_help()
      puts
      puts( "Usage:" )
      puts
      puts( "  " + PRG_NAME + ".rb [options] <domain table (parsed hmmpfam output)> [file containing complete sequences in fasta format] [outputfile]" )
      puts()
      puts( "  options: -" + E_VALUE_THRESHOLD_OPTION  + "=<f> : E-value threshold, default is no threshold" )
      puts( "               -" + OVERWRITE_IF_SAME_FROM_TO_OPTION  + " : overwrite domain with same start and end with domain with better E-value" )
      puts
      puts( "Examples:" )
      puts
      puts( "  " + PRG_NAME + ".rb P53_hmmscan_#{Constants::PFAM_V_FOR_EX}_10_domain_table P53_ni.fasta P53_hmmscan_300_10.dff" )
      puts
      puts( "  " + PRG_NAME + ".rb P53_hmmscan_#{Constants::PFAM_V_FOR_EX}_10_domain_table P53_ni.fasta" )
      puts
      puts( "  " + PRG_NAME + ".rb P53_hmmscan_#{Constants::PFAM_V_FOR_EX}_10_domain_table" )
      puts()
    end

    def is_ignorable?( line )
      return ( line !~ /[A-Za-z0-9-]/ || line =~ /^\s*#/)
    end

  end # class DomainsToForester

end # module Evoruby

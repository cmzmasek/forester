#!/usr/local/bin/ruby -w
#
# = lib/evo/apps/phylogenies_decorator
#
# Copyright (C) 2018 Christian M. Zmasek
# Copyright (C) 2018 J. Craig Venter Institute
# License::    GNU Lesser General Public License (LGPL)
#
# decoration of phylogenies with sequence/species names and domain architectures
#
# Environment variable FORESTER_HOME needs to point to the appropriate
# directory (e.g. setenv FORESTER_HOME $HOME/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/)

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'date'
require 'fileutils'

module Evoruby
  class PhylogeniesDecorator

    DECORATOR_OPTIONS_SEQ_NAMES = '-p -t -mp -or'
    DECORATOR_OPTIONS_DOMAINS   = '-p -t'
    SLEEP_TIME                  = 0.01
    REMOVE_NI                   = false
    TMP_FILE_1                  = '___PD1___'
    TMP_FILE_2                  = '___PD2___'
    LOG_FILE                    = '00_phylogenies_decorator.log'
    FORESTER_HOME               = ENV[Constants::FORESTER_HOME_ENV_VARIABLE]
    JAVA_HOME                   = ENV[Constants::JAVA_HOME_ENV_VARIABLE]

    PRG_NAME       = "phylogenies_decorator"
    PRG_DATE       = "170428"
    PRG_DESC       = "decoration of phylogenies with sequence/species names and domain architectures"
    PRG_VERSION    = "1.05"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    HELP_OPTION_1                           = "help"
    HELP_OPTION_2                           = "h"
    NO_DOMAINS_OPTION                       = 'nd'
    NO_SEQS_OPTION                          = 'ns'
    VERBOSE_OPTION                          = 'v'
    EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION = 'tc'

    NL = Constants::LINE_DELIMITER
    def run

      Util.print_program_information( PRG_NAME,
      PRG_VERSION,
      PRG_DESC,
      PRG_DATE,
      WWW,
      STDOUT )

      if ( ARGV == nil || ARGV.length < 2  )
        print_help
        exit( -1 )
      end

      if FORESTER_HOME == nil || FORESTER_HOME.length < 1
        Util.fatal_error( PRG_NAME, "apparently environment variable #{Constants::FORESTER_HOME_ENV_VARIABLE} has not been set" )
      end
      if JAVA_HOME == nil ||  JAVA_HOME.length < 1
        Util.fatal_error( PRG_NAME, "apparently environment variable #{Constants::JAVA_HOME_ENV_VARIABLE} has not been set" )
      end

      if !File.exist?( FORESTER_HOME )
        Util.fatal_error( PRG_NAME, '[' + FORESTER_HOME + '] does not exist' )
      end
      if !File.exist?( JAVA_HOME )
        Util.fatal_error( PRG_NAME, '[' + JAVA_HOME + '] does not exist' )
      end

      decorator = JAVA_HOME + '/bin/java -cp ' + FORESTER_HOME + '/java/forester.jar org.forester.application.decorator'

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

      if ( cla.get_number_of_files != 2 && cla.get_number_of_files != 3 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push(NO_DOMAINS_OPTION)
      allowed_opts.push(NO_SEQS_OPTION)
      allowed_opts.push(EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION)
      allowed_opts.push(VERBOSE_OPTION)

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
        "unknown option(s): " + disallowed,
        STDOUT )
      end

      no_domains = false
      if cla.is_option_set?(NO_DOMAINS_OPTION)
        no_domains = true
      end

      no_seqs_files = false
      if cla.is_option_set?(NO_SEQS_OPTION)
        no_seqs_files = true
      end

      extr_bracketed_tc = false
      if cla.is_option_set?(EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION)
        extr_bracketed_tc = true
      end

      verbose = false
      if cla.is_option_set?(VERBOSE_OPTION)
        verbose = true
      end

      if File.exist? LOG_FILE
        Util.fatal_error( PRG_NAME, 'logfile [' + LOG_FILE + '] already exists' )
      end

      in_suffix = cla.get_file_name( 0 )
      out_suffix = cla.get_file_name( 1 )

      mapping_files_dir = nil

      if cla.get_number_of_files == 3
        mapping_files_dir = cla.get_file_name( 2 )
      else
        mapping_files_dir = Dir.getwd
      end
      unless File.exist? mapping_files_dir
        Util.fatal_error( PRG_NAME, 'mapping files directory [' + mapping_files_dir + '] does not exist' )
      end
      unless File.directory? mapping_files_dir
        Util.fatal_error( PRG_NAME, '[' + mapping_files_dir + '] is not a directory' )
      end
      if Dir.entries(mapping_files_dir).length <= 2
        Util.fatal_error( PRG_NAME, 'mapping files directory [' + mapping_files_dir + '] is empty' )
      end

      mapping_files_dir = Util.canonical_path( mapping_files_dir.to_s )

      log = String.new

      now = DateTime.now
      log << "Program              : " + PRG_NAME + NL
      log << "Version              : " + PRG_VERSION + NL
      log << "Program date         : " + PRG_DATE + NL
      log << "Input/Output dir     : " + Dir.getwd + NL
      log << "Mappings file dir    : " + mapping_files_dir + NL
      log << "Input suffix         : " + in_suffix + NL
      log << "Output suffix        : " + out_suffix + NL
      log << "No domains data      : " + no_domains.to_s + NL
      log << "No mol seq data      : " + no_seqs_files.to_s + NL
      log << "Extract tax codes    : " + extr_bracketed_tc.to_s + NL
      log << "Date/time: " + now.to_s + NL + NL

      Util.print_message( PRG_NAME, 'Input/Output dir : ' + Dir.getwd )
      Util.print_message( PRG_NAME, 'Mappings file dir: ' + mapping_files_dir )
      Util.print_message( PRG_NAME, 'Input suffix     : ' + in_suffix )
      Util.print_message( PRG_NAME, 'Output suffix    : ' + out_suffix )
      Util.print_message( PRG_NAME, 'No domains data  : ' + no_domains.to_s )
      Util.print_message( PRG_NAME, 'No mol seq data  : ' + no_seqs_files.to_s )
      Util.print_message( PRG_NAME, 'Extract tax codes: ' + extr_bracketed_tc.to_s )

      if ( File.exist?( TMP_FILE_1 ) )
        File.delete( TMP_FILE_1 )
      end
      if ( File.exist?( TMP_FILE_2 ) )
        File.delete( TMP_FILE_2 )
      end

      files = Dir.entries( "." )

      counter = 0

      files.each { | phylogeny_file |
        if ( !File.directory?( phylogeny_file ) &&
        phylogeny_file !~ /^\./ &&
        phylogeny_file !~ /^00/ &&
        phylogeny_file !~ /#{out_suffix}$/ &&
        phylogeny_file =~ /#{in_suffix}$/ )
          begin
            Util.check_file_for_readability( phylogeny_file )
          rescue ArgumentError
            Util.fatal_error( PRG_NAME, 'can not read from: ' + phylogeny_file + ': '+ $!.to_s )
          end

          counter += 1

          outfile = phylogeny_file.sub( /#{in_suffix}$/, out_suffix )

          if REMOVE_NI
            outfile = outfile.sub( /_ni_/, '_' )
          end

          if File.exist?( outfile )
            msg = counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile +
            ' : already exists, skipping'
            Util.print_message( PRG_NAME, msg  )
            log << msg + NL
            next
          end

          if verbose
            puts
          end
          Util.print_message( PRG_NAME, counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile )
          log << counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile + NL

          phylogeny_id = phylogeny_file
          if  phylogeny_id == nil || phylogeny_id.size < 1
            Util.fatal_error( PRG_NAME, 'could not get id from ' + phylogeny_file.to_s )
          end
          if verbose
            Util.print_message( PRG_NAME, "Id: " + phylogeny_id )
          end
          log << "Id: " + phylogeny_id + NL

          ids_mapfile_path = nil
          domains_mapfile_path = nil
          seqs_file_path = nil

          ids_mapfile_name = get_file( mapping_files_dir, phylogeny_id, Constants::ID_MAP_FILE_SUFFIX )
          ids_mapfile_path = Util.canonical_path(mapping_files_dir, ids_mapfile_name)

          begin
            Util.check_file_for_readability( ids_mapfile_path)
          rescue IOError
            Util.fatal_error( PRG_NAME, "failed to read from [#{ids_mapfile_path}]: " + $!.to_s )
          end
          if verbose
            Util.print_message( PRG_NAME, "Ids mapfile: " + ids_mapfile_path )
          end
          log << "Ids mapfile: " + ids_mapfile_path + NL

          unless no_seqs_files
            seqs_file_name = get_file( mapping_files_dir, phylogeny_id, Constants::ID_NORMALIZED_FASTA_FILE_SUFFIX )
            seqs_file_path = Util.canonical_path(mapping_files_dir, seqs_file_name)
            begin
              Util.check_file_for_readability( seqs_file_path  )
            rescue IOError
              Util.fatal_error( PRG_NAME, "failed to read from [#{seqs_file_path}]: " + $!.to_s )
            end
            if verbose
              Util.print_message( PRG_NAME, "Seq file: " + seqs_file_path )
            end
            log << "Seq file: " + seqs_file_path + NL
          end

          unless no_domains
            domains_mapfile_name = get_file( mapping_files_dir , phylogeny_id, Constants::DOMAINS_TO_FORESTER_OUTFILE_SUFFIX  )
            domains_mapfile_path = Util.canonical_path(mapping_files_dir, domains_mapfile_name)
            begin
              Util.check_file_for_readability( domains_mapfile_path )
            rescue IOError
              Util.fatal_error( PRG_NAME, "failed to read from [#{domains_mapfile_path}]: " + $!.to_s )
            end
            if verbose
              Util.print_message( PRG_NAME, "Domains file: " + domains_mapfile_path )
            end
            log << "Domains file: " + domains_mapfile_path + NL
          end

          log <<  NL + NL

          if no_seqs_files
            FileUtils.cp(phylogeny_file, TMP_FILE_1)
          else
            cmd = decorator +
            ' -t -p -f=m ' + phylogeny_file + ' ' +
            seqs_file_path  + ' ' + TMP_FILE_1
            if verbose
              puts cmd
            end
            begin
              execute_cmd( cmd, log )
            rescue Exception
              Util.fatal_error( PRG_NAME, 'error: ' + $!.to_s )
            end
          end

          unless no_domains
            cmd = decorator + ' ' + DECORATOR_OPTIONS_DOMAINS + ' ' +
            '-f=d ' + TMP_FILE_1 + ' ' +
            domains_mapfile_path + ' ' + TMP_FILE_2
            if verbose
              puts cmd
            end
            begin
              execute_cmd( cmd, log )
            rescue Exception
              Util.fatal_error( PRG_NAME, 'error: ' + $!.to_s )
            end
          end

          opts = DECORATOR_OPTIONS_SEQ_NAMES
          if extr_bracketed_tc
            opts += ' -tc'
          end

          if no_domains
            cmd = decorator + ' ' + opts + ' -f=n ' + TMP_FILE_1 + ' ' +
            ids_mapfile_path + ' ' + outfile
            if verbose
              puts cmd
            end
            begin
              execute_cmd( cmd, log )
            rescue Exception
              Util.fatal_error( PRG_NAME, 'error: ' + $!.to_s )
            end
            File.delete( TMP_FILE_1 )
          else
            cmd = decorator + ' ' + opts + ' -f=n ' + TMP_FILE_2 + ' ' +
            ids_mapfile_path + ' ' + outfile
            if verbose
              puts cmd
            end
            begin
              execute_cmd( cmd, log )
            rescue Exception
              Util.fatal_error( PRG_NAME, 'error: ' + $!.to_s )
            end
            File.delete( TMP_FILE_1 )
            File.delete( TMP_FILE_2 )
          end
        end
      }
      open( LOG_FILE, 'w' ) do | f |
        f.write( log )
      end
      if verbose
        puts
      end
      Util.print_message( PRG_NAME, 'OK' )
      puts
    end # def run

    def execute_cmd( cmd, log )
      log << 'executing ' + cmd + NL
      IO.popen( cmd , 'r+' ) do | pipe |
        pipe.close_write
        log << pipe.read + NL + NL
      end
      if $?.to_i != 0
        raise StandardError, "failed to execute " + cmd
      end
      sleep( SLEEP_TIME )
    end

    def get_file( files_in_dir, phylogeny_id, suffix_pattern )
      begin
        Util.get_matching_file( files_in_dir, phylogeny_id, suffix_pattern )
      rescue Exception
        Util.fatal_error( PRG_NAME, 'error: ' + $!.to_s )
      end
    end

    def print_help()
      puts "Usage:"
      puts
      puts "   " + PRG_NAME + ".rb [options] <suffix of in-trees to be decorated> <suffix for decorated out-trees> [mapping files directory, default: current dir]"
      puts
      puts "   " + PRG_NAME + ".rb [options] <input directory> <output directory> <mapping files directory>"
      puts
      puts "   required file  (in mapping files directory): " + "name mappings       : #{Constants::ID_MAP_FILE_SUFFIX}"
      puts "   optional files (in mapping files directory): " + "sequences           : #{Constants::ID_NORMALIZED_FASTA_FILE_SUFFIX}"
      puts "                                                " + "domain architectures: #{Constants::DOMAINS_TO_FORESTER_OUTFILE_SUFFIX}"
      puts
      puts "   options: -" + NO_DOMAINS_OPTION  + ": to not add domain architecture information (#{Constants::DOMAINS_TO_FORESTER_OUTFILE_SUFFIX} file)"
      puts "            -" + NO_SEQS_OPTION   + ": to not add molecular sequence information (#{Constants::ID_NORMALIZED_FASTA_FILE_SUFFIX} file)"
      puts "            -" + EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION  + ": to extract bracketed taxonomic codes, e.g. [NEMVE]"
      puts "            -" + VERBOSE_OPTION  + " : verbose"
      puts
      puts "Examples: " + PRG_NAME + ".rb .xml _d.xml"
      puts "          " + PRG_NAME + ".rb -#{NO_DOMAINS_OPTION} -#{NO_SEQS_OPTION} .xml _d.xml"
      puts "          " + PRG_NAME + ".rb -#{NO_DOMAINS_OPTION} -#{NO_SEQS_OPTION} .xml _d.xml mappings_dir"
      # puts
      # puts "          " + PRG_NAME + ".rb in_trees_dir out_dir mappings_dir"
      puts
    end
  end # class PhylogenyiesDecorator

end # module Evoruby

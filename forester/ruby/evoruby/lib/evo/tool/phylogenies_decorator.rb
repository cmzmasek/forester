#!/usr/local/bin/ruby -w
#
# = lib/evo/apps/phylogenies_decorator
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# decoration of phylogenies with sequence/species names and domain architectures
#
# $Id: phylogenies_decorator.rb,v 1.34 2010/12/13 19:00:11 cmzmasek Exp $
#
# Environment variable FORESTER_HOME needs to point to the appropriate
# directory (e.g. setenv FORESTER_HOME $HOME/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/)

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'

require 'date'

module Evoruby

  class PhylogeniesDecorator

    #DECORATOR_OPTIONS_SEQ_NAMES = '-r=1 -mdn'
    #DECORATOR_OPTIONS_SEQ_NAMES = '-p -t -sn'
    DECORATOR_OPTIONS_SEQ_NAMES = '-p -t -c -tc -mp -or'
    # -mdn is a hidden expert option to rename e.g. "6_ORYLA3" to "6_[3]_ORYLA"
    #DECORATOR_OPTIONS_SEQ_NAMES = '-sn -r=1'
    #DECORATOR_OPTIONS_DOMAINS = '-r=1'
    DECORATOR_OPTIONS_DOMAINS = '-p -t'
    IDS_MAPFILE_SUFFIX        = '.nim'
    DOMAINS_MAPFILE_SUFFIX    = '.dff'
    SLEEP_TIME                = 0.1
    REMOVE_NI                 = true
    TMP_FILE                  = '___PD___'
    LOG_FILE                  = '00_phylogenies_decorator.log'
    FORESTER_HOME             = ENV[Constants::FORESTER_HOME_ENV_VARIABLE]
    JAVA_HOME                 = ENV[Constants::JAVA_HOME_ENV_VARIABLE]

    PRG_NAME       = "phylogenies_decorator"
    PRG_DATE       = "2012.10.11"
    PRG_DESC       = "decoration of phylogenies with sequence/species names and domain architectures"
    PRG_VERSION    = "1.02"
    COPYRIGHT      = "2012 Christian M Zmasek"
    CONTACT        = "phylosoft@gmail.com"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    IDS_ONLY_OPTION     = "n"
    DOMAINS_ONLY_OPTION = "d"
    HELP_OPTION_1       = "help"
    HELP_OPTION_2       = "h"

    NL = Constants::LINE_DELIMITER

    def run

      Util.print_program_information( PRG_NAME,
        PRG_VERSION,
        PRG_DESC,
        PRG_DATE,
        COPYRIGHT,
        CONTACT,
        WWW,
        STDOUT )

      if ( ARGV == nil || ARGV.length > 3 || ARGV.length < 2  )
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

      if File.exist?( LOG_FILE )
        Util.fatal_error( PRG_NAME, 'logfile [' + LOG_FILE + '] already exists' )
      end

      allowed_opts = Array.new
      allowed_opts.push( IDS_ONLY_OPTION )
      allowed_opts.push( DOMAINS_ONLY_OPTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME, "unknown option(s): " + disallowed )
      end

      ids_only = false
      domains_only = false

      in_suffix = cla.get_file_name( 0 )
      out_suffix = cla.get_file_name( 1 )

      if cla.is_option_set?( IDS_ONLY_OPTION )
        ids_only = true
      end
      if cla.is_option_set?( DOMAINS_ONLY_OPTION )
        domains_only = true
      end

      if ( ids_only && domains_only )
        Util.fatal_error( PRG_NAME, 'attempt to use ids only and domains only at the same time' )
      end

      log = String.new

      now = DateTime.now
      log << "Program              : " + PRG_NAME + NL
      log << "Version              : " + PRG_VERSION + NL
      log << "Program date         : " + PRG_DATE + NL
      log << "Options for seq names: " + DECORATOR_OPTIONS_SEQ_NAMES + NL
      log << "Options for domains  : " + DECORATOR_OPTIONS_DOMAINS + NL
      log << "FORESTER_HOME        : " + FORESTER_HOME + NL
      log << "JAVA_HOME            : " + JAVA_HOME + NL + NL
      log << "Date/time: " + now.to_s + NL
      log << "Directory: " + Dir.getwd  + NL + NL

      Util.print_message( PRG_NAME, 'input suffix     : ' + in_suffix )
      Util.print_message( PRG_NAME, 'output suffix    : ' + out_suffix )

      log << 'input suffix     : ' + in_suffix + NL
      log << 'output suffix    : ' + out_suffix + NL

      if ( File.exists?( TMP_FILE ) )
        File.delete( TMP_FILE )
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
            Util.fatal_error( PRG_NAME, 'can not read from: ' + phylogeny_file + ': '+ $! )
          end

          counter += 1

          outfile = phylogeny_file.sub( /#{in_suffix}$/, out_suffix )

          if REMOVE_NI
            outfile = outfile.sub( /_ni_/, '_' )
          end

          if File.exists?( outfile )
            msg = counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile +
             ' : already exists, skipping'
            Util.print_message( PRG_NAME, msg  )
            log << msg + NL
            next
          end

          Util.print_message( PRG_NAME, counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile )
          log << counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile + NL

          phylogeny_id = get_id( phylogeny_file )

          ids_mapfile_name = nil
          domains_mapfile_name = nil

          if ids_only
            ids_mapfile_name = get_file( files, phylogeny_id, IDS_MAPFILE_SUFFIX )
          elsif domains_only
            domains_mapfile_name = get_file( files, phylogeny_id, DOMAINS_MAPFILE_SUFFIX )
          else
            ids_mapfile_name = get_file( files, phylogeny_id, IDS_MAPFILE_SUFFIX )
            domains_mapfile_name = get_file( files, phylogeny_id, DOMAINS_MAPFILE_SUFFIX )
          end

          if domains_mapfile_name != nil
            begin
              Util.check_file_for_readability( domains_mapfile_name )
            rescue ArgumentError
              Util.fatal_error( PRG_NAME, 'failed to read from [#{domains_mapfile_name}]: ' + $! )
            end
          end

          if ids_mapfile_name != nil
            begin
              Util.check_file_for_readability( ids_mapfile_name )
            rescue ArgumentError
              Util.fatal_error( PRG_NAME, 'failed to read from [#{ids_mapfile_name}]: ' + $! )
            end
          end

          if domains_mapfile_name != nil
            if ids_mapfile_name != nil
              my_outfile = TMP_FILE
            else
              my_outfile = outfile
            end
            cmd = decorator + ' ' + DECORATOR_OPTIONS_DOMAINS + ' ' +
             '-f=d ' + phylogeny_file + ' ' +
             domains_mapfile_name + ' ' + my_outfile
            puts cmd
            execute_cmd( cmd, log )
          end

          if ids_mapfile_name != nil
            if domains_mapfile_name != nil
              my_infile = TMP_FILE
            else
              my_infile = phylogeny_file
            end
            cmd = decorator + ' ' +  DECORATOR_OPTIONS_SEQ_NAMES + ' ' +
             '-f=n ' + my_infile + ' ' +
             ids_mapfile_name + ' ' + outfile
            puts cmd
            execute_cmd( cmd, log )
          end

          if ( File.exists?( TMP_FILE ) )
            File.delete( TMP_FILE )
          end
        end
      }
      open( LOG_FILE, 'w' ) do | f |
        f.write( log )
      end
      puts
      Util.print_message( PRG_NAME, 'OK' )
      puts
    end # def run

    def execute_cmd( cmd, log )
      log << 'excuting ' + cmd + NL
      IO.popen( cmd , 'r+' ) do | pipe |
        pipe.close_write
        log << pipe.read + NL + NL
      end
      sleep( SLEEP_TIME )
    end


    def get_id( phylogeny_file_name )
      phylogeny_file_name =~ /^([^_]+)/
      $1
    end

    def get_file( files_in_dir, phylogeny_id, suffix_pattern )
      matching_files = Array.new
      matching_suffix_files = Array.new
      files_in_dir.each { | file |

        if ( !File.directory?( file ) &&
             file !~ /^\./ &&
             file !~ /^00/ &&
             file =~ /^#{phylogeny_id}.*#{suffix_pattern}$/ )
          matching_files << file
        end
        if ( !File.directory?( file ) &&
             file !~ /^\./ &&
             file !~ /^00/ &&
             file =~ /#{suffix_pattern}$/ )
          matching_suffix_files << file
        end
      }
      if matching_files.length < 1 && matching_suffix_files.length == 1
        return matching_suffix_files[ 0 ]
      end

      if matching_files.length < 1 && matching_suffix_files.length < 1
        Util.fatal_error( PRG_NAME, 'no file matching [' + phylogeny_id +
           '_] [' + suffix_pattern + '] present in current directory' )
      end
      if matching_files.length > 1
        Util.fatal_error( PRG_NAME, 'more than one file matching [' + phylogeny_id +
           '_] [' + suffix_pattern + '] present in current directory' )
      end
      matching_files[ 0 ]
    end

    def print_help()
      puts( "Usage:" )
      puts()
      puts( "  " + PRG_NAME + ".rb [options] <suffix of intrees to be decorated> <suffix for decorated outtrees> " )
      puts()
      puts( "  options: -" + IDS_ONLY_OPTION + ": decorate with sequence/species names only" )
      puts( "           -" + DOMAINS_ONLY_OPTION + ": decorate with domain structures only" )
      puts()
    end
  end # class PhylogenyiesDecorator

end # module Evoruby

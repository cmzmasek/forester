#!/usr/local/bin/ruby -w
#
# = lib/evo/apps/phylogenies_decorator
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)
#
# decoration of phylogenies with sequence/species names and domain architectures
#
# Environment variable FORESTER_HOME needs to point to the appropriate
# directory (e.g. setenv FORESTER_HOME $HOME/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/)

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'
require 'date'

module Evoruby
  class PhylogeniesDecorator

    DECORATOR_OPTIONS_SEQ_NAMES = '-p -t -mp -or'
    DECORATOR_OPTIONS_DOMAINS = '-p -t'
    IDS_MAPFILE_SUFFIX        = '.nim'
    DOMAINS_MAPFILE_SUFFIX    = '.dff'
    SLEEP_TIME                = 0.01
    REMOVE_NI                 = true
    TMP_FILE_1                  = '___PD1___'
    TMP_FILE_2                  = '___PD2___'
    LOG_FILE                  = '00_phylogenies_decorator.log'
    FORESTER_HOME             = ENV[Constants::FORESTER_HOME_ENV_VARIABLE]
    JAVA_HOME                 = ENV[Constants::JAVA_HOME_ENV_VARIABLE]

    PRG_NAME       = "phylogenies_decorator"
    PRG_DATE       = "170329"
    PRG_DESC       = "decoration of phylogenies with sequence/species names and domain architectures"
    PRG_VERSION    = "1.02"
    WWW            = "https://sites.google.com/site/cmzmasek/home/software/forester"

    HELP_OPTION_1                           = "help"
    HELP_OPTION_2                           = "h"
    NO_DOMAINS_OPTION                       = 'nd'
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

      if ( cla.get_number_of_files != 2 )
        print_help
        exit( -1 )
      end

      allowed_opts = Array.new
      allowed_opts.push(NO_DOMAINS_OPTION)
      allowed_opts.push(EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION)

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

      extr_bracketed_tc = false
      if cla.is_option_set?(EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION)
        extr_bracketed_tc = true
      end

      if File.exist?( LOG_FILE )
        Util.fatal_error( PRG_NAME, 'logfile [' + LOG_FILE + '] already exists' )
      end

      in_suffix = cla.get_file_name( 0 )
      out_suffix = cla.get_file_name( 1 )

      log = String.new

      now = DateTime.now
      log << "Program              : " + PRG_NAME + NL
      log << "Version              : " + PRG_VERSION + NL
      log << "Program date         : " + PRG_DATE + NL
      log << "No domains           : " + no_domains.to_s + NL
      log << "Extract taxo codes   : " + extr_bracketed_tc.to_s + NL
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
            Util.fatal_error( PRG_NAME, 'can not read from: ' + phylogeny_file + ': '+ $! )
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

          Util.print_message( PRG_NAME, counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile )
          log << counter.to_s + ': ' + phylogeny_file + ' -> ' +  outfile + NL

          phylogeny_id = get_id( phylogeny_file )
          if  phylogeny_id == nil || phylogeny_id.size < 1
            Util.fatal_error( PRG_NAME, 'could not get id from ' + phylogeny_file.to_s )
          end
          puts
          Util.print_message( PRG_NAME, "Id: " + phylogeny_id )
          log << "Id: " + phylogeny_id + NL

          ids_mapfile_name = nil
          domains_mapfile_name = nil
          seqs_file_name = nil

          ids_mapfile_name = get_file( files, phylogeny_id, IDS_MAPFILE_SUFFIX )

          begin
            Util.check_file_for_readability( ids_mapfile_name )
          rescue IOError
            Util.fatal_error( PRG_NAME, 'failed to read from [#{ids_mapfile_name}]: ' + $! )
          end
          Util.print_message( PRG_NAME, "Ids mapfile: " + ids_mapfile_name )
          log << "Ids mapfile: " + ids_mapfile_name + NL

          seqs_file_name = get_seq_file( files, phylogeny_id )
          begin
            Util.check_file_for_readability( seqs_file_name  )
          rescue IOError
            Util.fatal_error( PRG_NAME, 'failed to read from [#{seqs_file_name }]: ' + $! )
          end
          Util.print_message( PRG_NAME, "Seq file: " + seqs_file_name )
          log << "Seq file: " + seqs_file_name + NL

          unless no_domains
            domains_mapfile_name = get_file( files, phylogeny_id, DOMAINS_MAPFILE_SUFFIX )
            begin
              Util.check_file_for_readability( domains_mapfile_name )
            rescue IOError
              Util.fatal_error( PRG_NAME, 'failed to read from [#{domains_mapfile_name}]: ' + $! )
            end
            Util.print_message( PRG_NAME, "Domains file: " + domains_mapfile_name )
            log << "Domains file: " + domains_mapfile_name + NL
          end

          cmd = decorator +
          ' -t -p -f=m ' + phylogeny_file + ' ' +
          seqs_file_name  + ' ' + TMP_FILE_1
          puts cmd
          begin
            execute_cmd( cmd, log )
          rescue Error
            Util.fatal_error( PRG_NAME, 'error: ' + $! )
          end

          unless no_domains
            cmd = decorator + ' ' + DECORATOR_OPTIONS_DOMAINS + ' ' +
            '-f=d ' + TMP_FILE_1 + ' ' +
            domains_mapfile_name + ' ' + TMP_FILE_2
            puts cmd
            begin
              execute_cmd( cmd, log )
            rescue Error
              Util.fatal_error( PRG_NAME, 'error: ' + $! )
            end
          end

          opts = DECORATOR_OPTIONS_SEQ_NAMES
          if extr_bracketed_tc
            opts += ' -tc'
          end

          if no_domains
            cmd = decorator + ' ' + opts + ' -f=n ' + TMP_FILE_1 + ' ' +
            ids_mapfile_name + ' ' + outfile
            puts cmd
            begin
              execute_cmd( cmd, log )
            rescue Error
              Util.fatal_error( PRG_NAME, 'error: ' + $! )
            end
            File.delete( TMP_FILE_1 )
          else
            cmd = decorator + ' ' + opts + ' -f=n ' + TMP_FILE_2 + ' ' +
            ids_mapfile_name + ' ' + outfile
            puts cmd
            begin
              execute_cmd( cmd, log )
            rescue Error
              Util.fatal_error( PRG_NAME, 'error: ' + $! )
            end
            File.delete( TMP_FILE_1 )
            File.delete( TMP_FILE_2 )
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
      log << 'executing ' + cmd + NL
      IO.popen( cmd , 'r+' ) do | pipe |
        pipe.close_write
        log << pipe.read + NL + NL
      end
      sleep( SLEEP_TIME )
    end

    def get_id( phylogeny_file_name )
      if phylogeny_file_name =~ /^(.+?_DA)_/
        return $1
      elsif phylogeny_file_name =~ /^(.+?)_/
        return $1
      end
      nil
    end

    def get_file( files_in_dir, phylogeny_id, suffix_pattern )
      matching_files = Util.get_matching_files( files_in_dir, phylogeny_id, suffix_pattern )
      if matching_files.length < 1
        Util.fatal_error( PRG_NAME, 'no file matching [' + phylogeny_id +
        '...' + suffix_pattern + '] present in current directory' )
      end
      if matching_files.length > 1
        Util.fatal_error( PRG_NAME, 'more than one file matching [' +
        phylogeny_id  + '...' + suffix_pattern + '] present in current directory' )
      end
      matching_files[ 0 ]
    end

    def get_seq_file( files_in_dir, phylogeny_id )
      matching_files = Array.new

      files_in_dir.each { | file |

        if ( !File.directory?( file ) &&
        file !~ /^\./ &&
        file !~ /^00/ &&
        ( file =~ /^#{phylogeny_id}__.+\d$/ || file =~ /^#{phylogeny_id}_.*\.fasta$/ ) )
          matching_files << file
        end
      }

      if matching_files.length < 1
        Util.fatal_error( PRG_NAME, 'no seq file matching [' +
        phylogeny_id + '_] present in current directory' )
      end
      if matching_files.length > 1
        Util.fatal_error( PRG_NAME, 'more than one seq file matching [' +
        phylogeny_id + '_] present in current directory' )
      end
      matching_files[ 0 ]
    end

    def print_help()
      puts 'Usage:'
      puts
      puts "   " + PRG_NAME + ".rb <suffix of in-trees to be decorated> <suffix for decorated out-trees> "
      puts
      puts "   required files (in this dir): " + "name mappings       : .nim"
      puts "                                 " + "sequences           : _ni.fasta"
      puts "                                 " + "domain architectures: .dff"
      puts
      puts "   options: -" + NO_DOMAINS_OPTION  + ": to not add domain architecture information (.dff file)"
      puts "            -" + EXTRACT_BRACKETED_TAXONOMIC_CODE_OPTION  + ": to extract bracketed taxonomic codes, e.g. [NEMVE]"
      puts
      puts "Example: " + PRG_NAME + ".rb .xml _d.xml"
      puts
    end
  end # class PhylogenyiesDecorator

end # module Evoruby

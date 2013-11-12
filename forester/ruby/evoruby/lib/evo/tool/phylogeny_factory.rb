#
# = lib/evo/apps/phylogeny_factory - PhylogenyFactory class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: phylogeny_factory.rb,v 1.32 2010/12/13 19:00:11 cmzmasek Exp $

require 'lib/evo/util/constants'
require 'lib/evo/util/util'
require 'lib/evo/util/command_line_arguments'

require 'set'
require 'date'

module Evoruby

  class PhylogenyFactory

    PRG_NAME       = "phylogeny_factory"
    PRG_DATE       = "1301111"
    PRG_DESC       = "automated phylogeny reconstruction using queing system"
    PRG_VERSION    = "1.100"
    COPYRIGHT      = "2013 Christian M Zmasek"
    CONTACT        = "phylosoft@gmail.com"
    WWW            = "www.phylosoft.org"

    USE_JOB_SUBMISSION_SYSTEM_OPTION  = 's'
    BS_OPTION                         = 'b'
    LOG_FILE                          = '00_phylogeny_factory.log'
    TEMPLATE_FILE                     = '00_phylogeny_factory.template'
    PBS_O_WORKDIR                     = '$PBS_O_WORKDIR/'
    MIN_LENGTH_DEFAULT                = 40
    PFAM_HHMS                         = "/home/czmasek/DATA/PFAM/PFAM270X/PFAM_A_HMMs/"
    WALLTIME                          = '100:00:00'
    QUEUE                             = 'default'

    TMP_CMD_FILE_SUFFIX = '_QSUB'

    RSL                 = 'RSL'
    HMM                 = 'HMM'
    PHYLO_PL            = 'PHYLO_PL'

    OPTION_OPEN          = '%['
    OPTION_CLOSE          = ']%'

    WAIT                 = 1.0

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

      begin
        cla = CommandLineArguments.new( ARGV )
      rescue ArgumentError => e
        Util.fatal_error( PRG_NAME, "error: " + e.to_s )
      end

      allowed_opts = Array.new
      allowed_opts.push( USE_JOB_SUBMISSION_SYSTEM_OPTION )
      allowed_opts.push( BS_OPTION )

      disallowed = cla.validate_allowed_options_as_str( allowed_opts )
      if ( disallowed.length > 0 )
        Util.fatal_error( PRG_NAME,
          "unknown option(s): " + disallowed,
          STDOUT )
      end

      if File.exists?( LOG_FILE )
        puts( '[' + PRG_NAME + '] > log file [' + LOG_FILE + '] already exists' )
        exit( -1 )
      end

      if !File.exists?( TEMPLATE_FILE )
        puts( '[' + PRG_NAME + '] > template file [' + TEMPLATE_FILE + '] not found' )
        exit( -1 )
      end

      use_job_submission_system = false
      if cla.is_option_set?( USE_JOB_SUBMISSION_SYSTEM_OPTION )
        use_job_submission_system = true
      end

      bootstraps = 1
      if cla.is_option_set?( BS_OPTION )
        bootstraps = cla.get_option_value_as_int( BS_OPTION )
      end
      if bootstraps < 0
        puts( '[' + PRG_NAME + '] > negative bootstrap value' )
        exit( -1 )
      end
      if bootstraps == 0
        bootstraps = 1
      end

      log = String.new

      now = DateTime.now
      log << "Program     : " + PRG_NAME + NL
      log << "Version     : " + PRG_VERSION + NL
      log << "Program date: " + PRG_DATE + NL + NL
      log << "Bootstraps  : " + bootstraps.to_s + NL
      log << "Date/time   : " + now.to_s + NL
      log << "Directory   : " + Dir.getwd  + NL + NL

      puts( '[' + PRG_NAME + '] > reading ' + TEMPLATE_FILE )

      paths       = Hash.new  # path placeholder -> full path
      min_lengths = Hash.new  # alignment id -> minimal length
      options     = Hash.new  # option placeholder -> option

      commands    = Array.new

      log <<  "////////////////////////////////////////////////////////////////// #{NL}"
      log << "Template file [" + TEMPLATE_FILE + "]:#{NL}"

      command = String.new

      open( TEMPLATE_FILE ).each { | line |
        log << line
        if ( line =~ /^#/ )
        elsif ( line =~ /^\$\s*(\S+)\s*=\s*(\S+)/ )
          paths[ $1 ] = $2
          puts( '[' + PRG_NAME + '] > paths      : ' + $1 + ' => ' + $2 )

        elsif ( line =~ /^%\s*#{RSL}\s*(\S+)\s*=\s*(\S+)/ )
          min_lengths[ $1 ] = $2
          puts( '[' + PRG_NAME + '] > min lengths: ' + $1 + ' => ' + $2 )

        elsif ( line =~ /^%\s*(\S+)\s*=\s*(\S+)/ )
          key = $1
          value = $2
          if key == PHYLO_PL
            value = update_phylo_pl_options( value, bootstraps )
          end
          options[ key ] = value
          puts( '[' + PRG_NAME + '] > options    : ' + key + ' => ' + value )

        elsif ( line =~ /^>\s*(.+)/ )
          command = command + $1 + ";#{NL}"

        elsif ( line =~ /^-/  )
          commands << prepare( command, paths )
          command = String.new
        end
      }
      log << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ #{NL}#{NL}"

      files = Dir.entries( "." )

      files.each { | file |
        if ( !File.directory?( file ) &&
             file !~ /^\./ &&
             file !~ /#{TEMPLATE_FILE}/ &&
             file !~ /.bck$/ &&
             file !~ /.log$/ &&
             file !~ /nohup/ &&
             file !~ /^00/ )
          aln_name = file.to_str
          id = get_id( aln_name )
          if id != nil && id.length > 0
            puts
            puts( '[' + PRG_NAME + '] > file [id]  : ' + aln_name + ' [' + id + ']' )
          else
            puts
            puts( '[' + PRG_NAME + '] > file [id]  : ' + aln_name + ' [WARNING: could not get id!]' )
          end
          commands.each do | cmd |
            cmd = subst_hmm( cmd, id )
            cmd = subst_min_length( cmd, id, min_lengths )
            cmd = subst_options( cmd, options )
            if use_job_submission_system
              cmd = subst_aln_name( cmd, PBS_O_WORKDIR + aln_name )
            else
              cmd = subst_aln_name( cmd, aln_name )
            end

            if ( cmd =~ /%/ )
              cmd =~ /(%.*?%)/
              problem = $1
              puts( '[' + PRG_NAME + '] > WARNING    : command still contains placeholder: ' + problem )
              log << "WARNING: command still contains placeholder: " + cmd + NL
            else
              tmp_cmd_file = file.to_str[ 0..4 ] + TMP_CMD_FILE_SUFFIX
              if ( File.exists?( tmp_cmd_file ) )
                File.delete( tmp_cmd_file )
              end
              if use_job_submission_system
                open( tmp_cmd_file, 'w' ) do |f|
                  f.write( cmd )
                end
              end

              log << cmd + NL

              if use_job_submission_system
                IO.popen( 'qsub -q ' + QUEUE  + ' -l walltime=' + WALLTIME + ' ' + tmp_cmd_file , 'r+' ) do | pipe |
                  pipe.close_write
                end
              else
                spawn( 'nohup ' + cmd + ' &', STDERR => "/dev/null" )
              end

              sleep( WAIT )
              if ( File.exists?( tmp_cmd_file ) )
                File.delete( tmp_cmd_file )
              end
            end
          end
        end
      }

      open( LOG_FILE, 'w' ) do | f |
        f.write( log )
      end

      puts()
      puts( '[' + PRG_NAME + '] > OK' )
      puts()

    end # def run

    def prepare( command, paths )
      paths.each_pair{ | name, full |
        command = command.gsub( name, full )
      }
      command
    end

    def subst_options( command, options )
      opt_placeholders = command.scan( /%\[\S+\]%/ )
      opt_placeholders.each { | opt_placeholder |
        opt_placeholder = opt_placeholder.gsub( OPTION_OPEN , '' )
        opt_placeholder = opt_placeholder.gsub( OPTION_CLOSE, '' )
        opt_value = options[ opt_placeholder ]
        if ( opt_value != nil && opt_value.size > 0 )
          command = command.gsub( OPTION_OPEN + opt_placeholder + OPTION_CLOSE, opt_value )
        end
      }
      command
    end

    def subst_aln_name( command, aln_name )
      command = command.gsub( '$', aln_name )
      command
    end

    def subst_hmm( command, id )
      if id != nil && id.length > 0
        hmm = PFAM_HHMS + id + ".hmm"
        command = command.gsub( OPTION_OPEN + HMM + OPTION_CLOSE, hmm )
      end
      command
    end

    def update_phylo_pl_options( phylo_pl_options, bootstraps )
      opts = phylo_pl_options
      unless opts  =~ /B\d/
        opts = 'B' + bootstraps.to_s + opts
      end
      opts = '-' + opts
      opts
    end

    def subst_min_length( command, id, min_lengths )
      min_length = nil
      if id != nil && id.length > 0
        min_length = min_lengths[ id ]
      end
      if  min_length != nil && min_length.size > 0
        command = command.gsub( OPTION_OPEN + RSL + OPTION_CLOSE, min_length )
      else
        command = command.gsub( OPTION_OPEN + RSL + OPTION_CLOSE, MIN_LENGTH_DEFAULT.to_s )
      end
      command
    end

    def get_id( aln_name )
      if aln_name =~ /_{2}(.+)_{2}/
        return $1
      end
      nil
    end

  end # class PhylogenyFactory

end # module Evoruby


=begin

# Name convention if alignment specific parameters
# are to be used:
#  the substring between the first two double underscores is a
#  unique identifier and needs to match the identifiers
#  in '% <parameter-type> <unique-id>=<value>' statements
#  Example:
#  alignment name     : 'x__bcl2__e1'
#  parameter statments: '% RSL bcl2=60'
$ PROBCONS=/home/czmasek/SOFTWARE/PROBCONS/probcons_v1_12/probcons
$ DIALIGN_TX=/home/czmasek/SOFTWARE/DIALIGNTX/DIALIGN-TX_1.0.2/source/dialign-tx
$ DIALIGN_CONF=/home/czmasek/SOFTWARE/DIALIGNTX/DIALIGN-TX_1.0.2/conf
$ MAFFT=/home/czmasek/SOFTWARE/MAFFT/mafft-7.017-without-extensions/scripts/mafft
$ MUSCLE=/home/czmasek/SOFTWARE/MUSCLE/muscle3.8.31/src/muscle
$ CLUSTALO=/home/czmasek/SOFTWARE/CLUSTAL_OMEGA/clustal-omega-1.1.0/src/clustalo
$ KALIGN=/home/czmasek/SOFTWARE/KALIGN/kalign203/kalign
$ HMMALIGN=/home/czmasek/SOFTWARE/HMMER/hmmer-3.0/src/hmmalign
$ MSA_PRO=/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/msa_pro.rb
$ PHYLO_PL=/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/archive/perl/phylo_pl.pl


% RSL Hormone_recep=60
%
% RSL Y_phosphatase=100
% RSL Y_phosphatase2=75
% RSL Y_phosphatase3=50
% RSL Y_phosphatase3C=40

% PHYLO_OPT=B100q@1r4j2IGS21X

% TMP_DIR  = /home/czmasek/tmp/


> KALIGN $ > $_kalign
> MSA_PRO -o=p -n=10 -d -rr=0.5 -c -rsl=%[RSL]% $_kalign $_kalign_05_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_kalign_05_%[RSL]%.aln $_kalign_05_%[RSL]% %[TMP_DIR]%
-

> KALIGN $ > $_kalign_
> MSA_PRO -o=p -n=10 -d -rr=0.9 -c -rsl=%[RSL]% $_kalign_ $_kalign_09_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_kalign_09_%[RSL]%.aln $_kalign_09_%[RSL]% %[TMP_DIR]%
-


> HMMALIGN --amino --trim --outformat Pfam -o $_hmmalign %[HMM]% $ > /dev/null
> MSA_PRO -o=p -n=10 -d -rr=0.5 -c -rsl=%[RSL]% $_hmmalign $_hmmalign_05_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_hmmalign_05_%[RSL]%.aln $_hmmalign_05_%[RSL]% %[TMP_DIR]%
-

> HMMALIGN --amino --trim --outformat Pfam -o $_hmmalign_ %[HMM]% $ > /dev/null
> MSA_PRO -o=p -n=10 -d -rr=0.9 -c -rsl=%[RSL]% $_hmmalign_ $_hmmalign_09_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_hmmalign_09_%[RSL]%.aln $_hmmalign_09_%[RSL]% %[TMP_DIR]%
-


> MAFFT --maxiterate 1000 --localpair $ > $_mafft
> MSA_PRO -o=p -n=10 -d -rr=0.5 -c -rsl=%[RSL]% $_mafft $_mafft_05_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_mafft_05_%[RSL]%.aln $_mafft_05_%[RSL]% %[TMP_DIR]%
-

> MAFFT --maxiterate 1000 --localpair $ > $_mafft_
> MSA_PRO -o=p -n=10 -d -rr=0.9 -c -rsl=%[RSL]% $_mafft_ $_mafft_09_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_mafft_09_%[RSL]%.aln $_mafft_09_%[RSL]% %[TMP_DIR]%
-


> MUSCLE  -maxiters 1000 -maxtrees 100 -in $ -out $_muscle
> MSA_PRO -o=p -n=10 -d -rr=0.5 -c -rsl=%[RSL]% $_muscle $_muscle_05_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_muscle_05_%[RSL]%.aln $_muscle_05_%[RSL]% %[TMP_DIR]%
-

> MUSCLE  -maxiters 1000 -maxtrees 100 -in $ -out $_muscle_
> MSA_PRO -o=p -n=10 -d -rr=0.9 -c -rsl=%[RSL]% $_muscle_ $_muscle_09_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_muscle_09_%[RSL]%.aln $_muscle_09_%[RSL]% %[TMP_DIR]%
-


> CLUSTALO --full --full-iter --iter=5 -i $ -o $_clustalo
> MSA_PRO -o=p -n=10 -d -rr=0.5 -c -rsl=%[RSL]% $_clustalo $_clustalo_05_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_clustalo_05_%[RSL]%.aln $_clustalo_05_%[RSL]% %[TMP_DIR]%
-

> CLUSTALO --full --full-iter --iter=5 -i $ -o $_clustalo_
> MSA_PRO -o=p -n=10 -d -rr=0.9 -c -rsl=%[RSL]% $_clustalo_ $_clustalo_09_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_clustalo_09_%[RSL]%.aln $_clustalo_09_%[RSL]% %[TMP_DIR]%
-


> PROBCONS $ > $_probcons
> MSA_PRO -o=p -n=10 -d -rem_red -rr=0.5 -c -rsl=%[RSL]% $_probcons $_probcons_05_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_probcons_05_%[RSL]%.aln $_probcons_05_%[RSL]% %[TMP_DIR]%
-

> PROBCONS $ > $_probcons_
> MSA_PRO -o=p -n=10 -d -rem_red -rr=0.9 -c -rsl=%[RSL]% $_probcons_ $_probcons_09_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_probcons_09_%[RSL]%.aln $_probcons_09_%[RSL]% %[TMP_DIR]%
-


> DIALIGN_TX DIALIGN_CONF $ $_dialigntx
> MSA_PRO -o=p -n=10 -d -rem_red -rr=0.5 -c -rsl=%[RSL]% $_dialigntx $_dialigntx_05_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_dialigntx_05_%[RSL]%.aln $_dialigntx_05_%[RSL]% %[TMP_DIR]%
-

> DIALIGN_TX DIALIGN_CONF $ $_dialigntx_
> MSA_PRO -o=p -n=10 -d -rem_red -rr=0.9 -c -rsl=%[RSL]% $_dialigntx_ $_dialigntx_09_%[RSL]%.aln
> PHYLO_PL %[PHYLO_OPT]% $_dialigntx_09_%[RSL]%.aln $_dialigntx_09_%[RSL]% %[TMP_DIR]%
-

=end


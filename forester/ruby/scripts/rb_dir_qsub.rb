#!/usr/local/bin/ruby -w
#
# = rb_dir_qsub
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: rb_dir_qsub.rb,v 1.15 2009/11/07 02:06:59 cmzmasek Exp $
#
# To execute qsub commands.
# Submits PRG for every file in the current directory.
#
# Examples for PARAMETER_FILE:
#
# PRG:     /home/user/SOFTWARE/HMMER/hmmer-2.3.2/src/hmmpfam
# OPT:     -E 20 -A 0 /home/user/DATA/PFAM/Pfam_ls
# SUFFIX:  _hmmpfam_22_20_ls
#
# PRG:     /home/user/SOFTWARE/WUBLAST/tblastn
# OPT:
# VOPT:    AMPQU
# VOPT:    HYDMA
# SUFFIX:  _blast
#
# PRG:     /home/czmasek/SOFTWARE/HMMER/hmmer-3.0b2/binaries/intel-linux-x86_64/hmmscan
# OPT:     -E 2 --notextw --qformat fasta /home/czmasek/DATA/PFAM/PFAM240/Pfam-A.hmm
# SUFFIX:  .hmmscan30b2_240
# OUTPUT:  --domtblout


module ForesterScripts

  if RUBY_VERSION !~ /1.9/
    puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
    exit( -1 )
  end

  PARAMETER_FILE = 'parameters.rb_dir_qsub'
  SLEEP          = 1.0
  REMOVE_SUFFIX  = true
  LIMIT_TO_TAX_CODE  = true

  PRG         = 'PRG:'
  OPT         = 'OPT:'
  VOPT        = 'VOPT:'
  OUTPUT_OPT  = 'OUTPUT:'
  SUFFIX      = 'SUFFIX:'
  INPUT_PART  = 'INPUT_PART:'


  PBS_O_WORKDIR       = '$PBS_O_WORKDIR/'
  TMP_CMD_FILE_SUFFIX = '__QSUB'
  NAME                = 'rb_dir_qsub'

  if ( !File.exists?( PARAMETER_FILE ) )
    puts( '[' + NAME + '] > parameters file "' + PARAMETER_FILE + '" not found' )
    Process.exit!
  end
  puts( '[' + NAME + '] > reading ' + PARAMETER_FILE )

  prg = ''
  opt = ''
  vopts = Array.new
  suffix = ''
  input_part = ''
  output_opt = ''
  open( PARAMETER_FILE ).each { |line|
    if ( line.length > 1 && line =~ /^[^#]\S+/ )
      if line =~ /^#{PRG}\s+(\S+)/
        prg = $1
      end
      if line =~ /^\s*#{OPT}\s+(\S+.+)/
        opt = $1
      end
      if line =~ /^\s*#{VOPT}\s+(\S+.+)/
        vopts.push( $1 )
      end
      if line =~ /^\s*#{SUFFIX}\s+(\S+)/
        suffix = $1
      end
      if line =~ /^\s*#{INPUT_PART}\s+(\S+)/
        input_part = $1
      end
      if line =~ /^\s*#{OUTPUT_OPT}\s+(\S+.+)/
        output_opt = $1
      end
    end
  }
  if ( prg.length < 1 )
    puts( '[' + NAME + '] > no program name found in parameters file "' + PARAMETER_FILE + '"' )
    Process.exit!
  end
  puts( '[' + NAME + '] > program: ' + prg )
  puts( '[' + NAME + '] > option :  ' + opt )
  vopts.each { |vopt|
    puts( '[' + NAME + '] > voption:  ' + vopt )
  }
  puts( '[' + NAME + '] > suffix :  ' + suffix )
  if ( input_part.length > 0 )
    puts( '[' + NAME + '] > input:  ' + input_part )
  end
  if ( output_opt.length > 0 )
    puts( '[' + NAME + '] > output opt :  ' + output_opt )
  end
  if vopts.empty?
    vopts.push( "" )
  end

  files = Dir.entries( "." )

  files.each { |file|
    if ( !File.directory?( file ) && file !~ /^\./ && file !~ /#{PARAMETER_FILE}/ )

      if ( input_part.length > 0 && file !~ /#{input_part}/ )
        next
      end
      vopts.each { |vopt|
        cmd = ""
        outputfile = file.to_str
        if LIMIT_TO_TAX_CODE
          if outputfile =~ /^([A-Z0-9]{3,5})[_\.]/
            outputfile = $1
          end
        end
        if REMOVE_SUFFIX
          if outputfile =~ /(.+)\..{1,5}/
            outputfile = $1
          end
        end
        if output_opt.length > 0
          cmd = prg + ' ' +
           output_opt + ' ' + PBS_O_WORKDIR + outputfile + suffix + ' ' +
           opt + ' ' +
           PBS_O_WORKDIR + file.to_str +
           ' > /dev/null'
        elsif vopt.length > 0
          cmd = prg + ' ' + opt + ' ' + vopt + ' ' + PBS_O_WORKDIR + file.to_str +
           ' > ' + PBS_O_WORKDIR + vopt + "_" + outputfile + suffix
        else
          cmd = prg + ' ' + opt + ' ' + PBS_O_WORKDIR + file.to_str +
           ' > ' + PBS_O_WORKDIR + outputfile + suffix
        end
        tmp_cmd_file = file.to_str + TMP_CMD_FILE_SUFFIX
        if File.exists?( tmp_cmd_file )
          File.delete( tmp_cmd_file )
        end
        open( tmp_cmd_file, 'w' ) do |f|
          f.write( cmd )
        end
        puts( '[' + NAME + '] > excuting ' + cmd )
        IO.popen( 'qsub ' + tmp_cmd_file , 'r+' ) do |pipe|
          pipe.close_write
          puts pipe.read
        end
        sleep( SLEEP )
        if File.exists?( tmp_cmd_file )
          File.delete( tmp_cmd_file )
        end
      }
    end
  }
  puts( '[' + NAME + '] > OK.' )
  puts

end

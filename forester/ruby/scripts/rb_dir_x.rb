#!/usr/local/bin/ruby -w
#
# = rb_x_qsub
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: rb_dir_x.rb,v 1.8 2008/09/16 23:31:39 cmzmasek Exp $
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


module ForesterScripts

    if RUBY_VERSION !~ /1.9/
        puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
        exit( -1 )
    end

    PARAMETER_FILE    = 'parameters.rb_dir_x'
    SLEEP = 1.0
    SPAWN = true

    PRG         = 'PRG:'
    OPT         = 'OPT:'
    VOPT        = 'VOPT:'
    OUTPUT_OPT  = 'OUTPUT:' # TODO e.g. > or -o
    SUFFIX      = 'SUFFIX:'
    INPUT_PART  = 'INPUT_PART:'

    NAME        = 'rb_dir_x'

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
        puts( '[' + NAME + '] > input  :  ' + input_part )
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
                if vopt.length > 0
                    cmd = 'nohup ' + prg + ' ' + opt + ' ' + vopt + ' ' + file.to_str +
                     ' > ' + vopt + "_" + file.to_str + suffix + ' &'
                else
                    cmd = 'nohup ' + prg + ' ' + opt + ' ' + file.to_str +
                     ' > ' + file.to_str + suffix + ' &'
                end

                puts( '[' + NAME + '] > excuting ' + cmd )
                if SPAWN
                    spawn( cmd, STDERR => "/dev/null" )
                else
                    IO.popen( cmd , 'r+' ) do |pipe|
                        pipe.close_write
                        puts pipe.read
                    end
                end
                sleep( SLEEP )

            }
        end
    }
    puts( '[' + NAME + '] > OK.' )
    puts

end

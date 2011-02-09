#!/usr/local/bin/ruby -w
#
# = rb_qsub 
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: rb_qsub.rb,v 1.6 2008/08/30 19:57:59 cmzmasek Exp $
#
# last modified: 11/13/2007
#
#
# To execute qsub commands.
# Each line l (unless precded by a # or space) in file
# 'commands.qsub' is executed as 'qsub l'


module ForesterScripts

    if RUBY_VERSION !~ /1.9/
        puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
        exit( -1 )
    end     
    
    CMDS_FILE    = 'commands.qsub'
    TMP_CMD_FILE = '__QSUB_RB_CMD__'
    PRG_NAME     = 'rb_qsub'

    if ( !File.exists?( CMDS_FILE ) ) 
        puts( '[' +PRG_NAME + '] > commands file "' + CMDS_FILE + '" not found' )
        Process.exit!          
    end    
    
    puts( '[' +PRG_NAME + '] > reading ' + CMDS_FILE )

    open( CMDS_FILE ).each { |line| 
        if ( line.length > 1 && line =~ /^[^#]\S+/ )
            if ( File.exists?( TMP_CMD_FILE ) ) 
                File.delete( TMP_CMD_FILE ) 
            end
            open( TMP_CMD_FILE, 'w' ) do |f|
                f.write( line )
            end 
            puts( '[' +PRG_NAME + '] > excuting ' + line )
            IO.popen( 'qsub ' + TMP_CMD_FILE , 'r+' ) do |pipe|
                pipe.close_write
                puts pipe.read
            end
            if ( File.exists?( TMP_CMD_FILE ) ) 
                File.delete( TMP_CMD_FILE ) 
            end
            sleep( 10.0 )
        end
    }
    puts( '[' +PRG_NAME + '] > OK.' )
    puts

end


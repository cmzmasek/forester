#!/usr/local/bin/ruby -w

# forester -- software libraries and applications
# for evolutionary biology and genomics.
# Copyright (C) 2026 Christian M. Zmasek
# All rights reserved
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Contact: czmasek at jcvi dot org

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


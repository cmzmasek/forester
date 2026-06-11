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

    NUMBER_OF_LETTERS = 3

    if ( ARGV == nil || ARGV.length != 3 )
        puts( "usage: replace_id.rb <map-file> <infile> <outfile>" )         
        exit( -1 )
    end    
    mapfile = ARGV[ 0 ]
    infile  = ARGV[ 1 ]
    outfile = ARGV[ 2 ]
    
    
    if ( File.exists?( outfile ) ) 
        puts( "outfile [" + outfile + "] already exists" )
        exit( -1 )  
    end
    if ( !File.exists?( infile) )
        puts( "infile [" + infile + "] does not exist" )
        exit( -1 ) 
    end 
    if ( !File.exists?( mapfile ) ) 
        puts( "mapfile [" + mapfile + "] does not exist" )
        exit( -1 ) 
    end                
    
    number_to_complete_id_map = Hash.new
    
    File.open( mapfile ) do | file |
        while line = file.gets
            if ( line =~ /(\d+_\S+)\s*:/ )
                complete_id = $1
                complete_id =~ /(\d+)_\S+/
                number_to_complete_id_map[ $1 ] = complete_id     
                puts( $1 + ' => ' + complete_id )
            end
        end
    end 
    
    if ( number_to_complete_id_map.size < 1 ) 
        puts( "mapping file was empty" )         
        exit( -1 )    
    end   
    
    data_str = String.new
    
    File.open( infile ) do | file |
        while line = file.gets
            data_str = data_str + line.chomp
        end 
    end     
    
    replacements = 0
    number_to_complete_id_map.each_pair{ |number, id|
        data_str.gsub!( /\b#{number}_[A-Z]{#{NUMBER_OF_LETTERS}}/, id )         
    }
    
    open( outfile, 'w' ) do |file|
        file.write( data_str )
    end      
    
    puts( "wrote " + outfile )
    puts( "OK" )
    
end


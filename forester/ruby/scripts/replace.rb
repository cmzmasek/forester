#!/usr/local/bin/ruby -w
#
# = replace 
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: replace.rb,v 1.5 2008/08/28 17:09:07 cmzmasek Exp $
#
# To replace multiple strings in file.
# Map file contains intructions for replacement (one on each line)
# in the following format (by example): old#new
#


module ForesterScripts
    
    if RUBY_VERSION !~ /1.9/
                      puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
                      exit( -1 )
                end     
    
    if ( ARGV == nil || ARGV.length != 3 )
        puts( "usage: replace.rb <map-file> <infile> <outfile>" )         
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
    
    old_new_map = Hash.new
    
    File.open( mapfile ) do | file |
        while line = file.gets
            if ( line =~/(\S+)\s*#\s*(\S+)/ )
                old_new_map[ $1 ] = $2      
                puts( $1 + ' => ' + $2 )     
            end
        end
    end 
    
    if ( old_new_map.size < 1 ) 
        puts( "mapping file was empty" )         
        exit( -1 )    
    end   
    
    data_str = String.new
    
    File.open( infile ) do | file |
        while line = file.gets
            data_str =  data_str + line.chomp
        end 
    end     
    
    old_new_map.each_pair{ |old, new|
        data_str = data_str.gsub( old, new )
    }
    
    open( outfile, 'w' ) do |file|
        file.write( data_str )
    end      
    
    puts( "wrote " + outfile )
    
end

    
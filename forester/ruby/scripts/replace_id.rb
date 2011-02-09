#!/usr/local/bin/ruby -w
#
# = replace_id 
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: replace_id.rb,v 1.8 2008/08/28 17:09:07 cmzmasek Exp $
#
# To replace ()by way of example '123_CHI5' with '123_CHICK5'
# given a mapping file containing '123_CHICKEN'
# (in the form '123_CHICKEN: some description which is ignored').
#
# Note. This will break if the species id ends with a number (as is 
# in the case for many bacteria).


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


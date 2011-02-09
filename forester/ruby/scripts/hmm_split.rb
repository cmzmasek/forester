#!/usr/local/bin/ruby -w
#
# = hmm_split 
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: hmm_split.rb,v 1.5 2008/11/17 22:32:43 cmzmasek Exp $
#
# To split a Pfam HMM file into one file for each HMM.
#


module ForesterScripts
    
    if RUBY_VERSION !~ /1.9/
                      puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
                      exit( -1 )
                end     
    
    
    if ( ARGV == nil || ARGV.length != 3 )
        puts( "usage: hmm_split.rb <Pfam HMM file> <outfile suffix> <outdir>" )         
        exit( -1 )
    end    
       
       hmmfile = ARGV[ 0 ]
       suffix  = ARGV[ 1 ]
       outdir  = ARGV[ 2 ]
     
       if ( !File.exists?( outdir ) )
           puts( "outdir [" + outdir + "] does not exist" )
           exit( -1 ) 
       end 
       if ( !File.exists?( hmmfile ) ) 
           puts( "Pfam HMM file [" + hmmfile + "] does not exist" )
           exit( -1 ) 
       end                
       
       data = String.new
       name = String.new
       line_count = 0
       count = 0
       
       File.open( hmmfile ) do | file |
           while line = file.gets
               data = data + line
               line_count += 1
               if ( line =~ /NAME\s+(.+)/ )
                   if name.length > 0
                       puts( "Pfam HMM file [" + hmmfile + "] format error [line: " + line + "]" )
                       exit( -1 )                        
                   end    
                   name = $1    
               elsif ( line =~ /\/\// )
                   if name.length < 1
                       puts( "Pfam HMM file [" + hmmfile + "] format error [line: " + line + "]" )
                       exit( -1 )                        
                   end
                                       
                   outfile = outdir + '/' + name + suffix
                   if ( File.exists?( outfile ) ) 
                       puts( "file [" + outfile + "] already exists" )
                       exit( -1 )  
                   end                   
                   open( outfile, 'w' ) do | out |
                       out.write( data )
                   end
                   count += 1
                   puts( count.to_s + ": " + name )
                   data = String.new
                   name = String.new                                      
               end
           end
       end 
         
       puts()
       puts( "wrote " + count.to_s + " individual HMM files to " + outdir )    
    
end    
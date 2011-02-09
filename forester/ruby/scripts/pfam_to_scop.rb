#!/usr/local/bin/ruby -w
#
# = pfam_to_scop
#
# Copyright::  Copyright (C) 2008-2009 Christian M. Zmasek. All rights reserved.
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: pfam_to_scop.rb,v 1.2 2008/08/28 17:09:07 cmzmasek Exp $
#
# This extracts ID and SCOP fa (or fa and sf) from Pfam data files.
#
# Created 2008-06-25 in San Diego, CA, USA by CMZ
#
# Usage: pfam_to_scop.rb <infile: Pfam data file such as Pfam-A.full> <outfile>

require 'iconv'

module ForesterScripts

    if RUBY_VERSION !~ /1.9/
                      puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
                      exit( -1 )
                end      
    
    SF = true

    SEP = "\t"
    LINE_DELIMITER  = "\n"

    if ( ARGV == nil || ARGV.length != 2 )
        puts( "usage: pfam_to_scop.rb <infile: Pfam data file such as Pfam-A.full> <outfile>" )
        exit( -1 )
    end

    pfamfile = ARGV[ 0 ]
    outfile  = ARGV[ 1 ]

    if ( !File.exists?( pfamfile ) )
        puts( "Pfam data file [" + pfamfile + "] does not exist" )
        exit( -1 )
    end
    if ( File.exists?( outfile ) )
        puts( "outfile [" + outfile + "] already exists" )
        exit( -1 )
    end

    ic = Iconv.new( 'UTF-8//IGNORE', 'UTF-8' )

    id = nil
    scops = Array.new()
    line_count = 0
    count = 0
    scop_count = 0

    out = File.open( outfile, 'w' )

    File.open( pfamfile ) do | file |
        while line = file.gets
            line_count += 1

            line = ic.iconv( line )

            if ( line =~ /#=GF ID\s+(.+)/ )
                if ( id != nil )
                    puts( "Pfam data file [" + pfamfile + "] format error [line: " + line + "]" )
                    exit( -1 )
                end
                id = $1
            elsif ( line =~ /#=GF\s+DR\s+SCOP;\s+(\w+);\s+fa/ )
                scops.push( $1 )
            elsif ( SF && line =~ /#=GF\s+DR\s+SCOP;\s+(\w+);\s+sf/ )
                scops.push( $1 )
            elsif ( line =~ /^\/\// )
                if ( id == nil )
                    puts( "Pfam data file [" + pfamfile + "] format error [line: " + line + "]" )
                    exit( -1 )
                end
                scops.each { |s|
                    out.write( id )
                    out.write( SEP )
                    out.write( s )
                    out.write( LINE_DELIMITER )
                    scop_count += 1
                }
                id = nil
                scops = Array.new()
                count += 1
            end
        end
    end

    out.close

    puts()
    if ( SF )
        puts( "Extracted #{scop_count} scop fa and sf identifiers for #{count.to_s} individual Pfams to " + outfile )
    else
        puts( "Extracted #{scop_count} scop fa identifiers for #{count.to_s} individual Pfams to " + outfile )
    end
    puts( "OK" )
    puts()

end # module ForesterScripts
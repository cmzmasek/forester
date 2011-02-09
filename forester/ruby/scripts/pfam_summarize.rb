#!/usr/local/bin/ruby -w
#
# = pfam_summarize
#
# Copyright::  Copyright (C) 2008-2009 Christian M. Zmasek. All rights reserved.
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: pfam_summarize.rb,v 1.2 2008/08/28 17:09:07 cmzmasek Exp $
#
# This extracts ID, AC, DE, TP, and DR values from Pfam data files.
#
# Created 2008-06-25 in San Diego, CA, USA by CMZ
#
# Usage: pfam_summarize.rb <infile: Pfam data file such as Pfam-A.full> <outfile>

require 'iconv'

module ForesterScripts
    if RUBY_VERSION !~ /1.9/
                      puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
                      exit( -1 )
                end     
    SEP = "\t"
    LINE_DELIMITER  = "\n"

    if ( ARGV == nil || ARGV.length != 2 )
        puts( "usage: pfam_summarize.rb <infile: Pfam data file such as Pfam-A.full> <outfile>" )
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
    ac = nil
    de = nil
    tp = nil
    dr = Array.new()
    line_count = 0
    count = 0

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
            elsif ( line =~ /#=GF AC\s+(.+)/ )
                ac = $1
            elsif ( line =~ /#=GF DE\s+(.+)/ )
                de = $1
            elsif ( line =~ /#=GF TP\s+(.+)/ )
                tp = $1
            elsif ( line =~ /#=GF DR\s+(.+)/ )
                dr.push( $1 )
            elsif ( line =~ /^\/\// )
                if ( id == nil || ac == nil )
                    puts( "Pfam data file [" + pfamfile + "] format error [line: " + line + "]" )
                    exit( -1 )
                end
                out.write( id )
                out.write( SEP )
                out.write( ac )
                out.write( SEP )
                out.write( tp )
                out.write( SEP )
                out.write( '[' )
                out.write( de )
                out.write( ']' )
                out.write( SEP )
                out.write( '[' )
                dr.each { |d|
                    out.write( d )
                    out.write( ' ' )
                }
                out.write( ']' )
                out.write( LINE_DELIMITER )

                id = nil
                ac = nil
                de = nil
                tp = nil
                dr = Array.new()
                count += 1
            end
        end
    end

    out.close

    puts()
    puts( "Summarized data for " + count.to_s + " individual Pfams to " + outfile )
    puts( "OK" )
    puts()

end # module ForesterScripts


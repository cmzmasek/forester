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
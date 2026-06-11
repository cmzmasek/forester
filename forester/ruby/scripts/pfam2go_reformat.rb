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

    require 'set'

    if RUBY_VERSION !~ /1.9/
        puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
        exit( -1 )
    end

    if ( ARGV == nil || ARGV.length != 2 )
        puts( "usage: pfam2go_reformat.rb <pfam2go file> <outfiles base>" )
        exit( -1 )
    end

    infile  = ARGV[ 0 ]
    outfilebase = ARGV[ 1 ]
    outfile_sgd_style = outfilebase + "_sgd_style_associations"
    outfile_simple_map = outfilebase + "_basic_associations"
    outfile_all_pfams = outfilebase + "_all_associated_pfams"

    pfams = SortedSet.new

    if ( File.exists?( outfile_sgd_style ) )
        puts( "outfile [" +  outfile_sgd_style + "] already exists" )
        exit( -1 )
    end
    if ( File.exists?( outfile_simple_map ) )
        puts( "outfile [" +  outfile_simple_map + "] already exists" )
        exit( -1 )
    end
    if ( File.exists?( outfile_all_pfams ) )
        puts( "outfile [" + outfile_all_pfams + "] already exists" )
        exit( -1 )
    end
    if ( !File.exists?( infile) )
        puts( "infile [" + infile + "] does not exist" )
        exit( -1 )
    end

    out_str_sgd = String.new
    out_str_basic = String.new

    File.open( infile ) do | file |
        while line = file.gets
            if line =~ /^\s*Pfam:PF(\d+)\s+(\S+)\s.+(GO:\d+)\s*$/
                pfam_id = $1
                pfam_name = $2
                go_id = $3
                new_line = "PFAM" + "\t" + pfam_name + "\t" + pfam_name + "\t\t" + go_id + "\t" + "PF:" + pfam_id + "\t\t\t\t\t\t\t\t\t"
                out_str_sgd = out_str_sgd + new_line + "\n"
                out_str_basic = out_str_basic + pfam_name + "\t" + go_id + "\n"
                pfams.add( pfam_name )
            end
        end
    end

    open(  outfile_sgd_style, 'w' ) do |file|
        file.write( out_str_sgd )
    end
    open( outfile_simple_map, 'w' ) do |file|
        file.write( out_str_basic )
    end
    open( outfile_all_pfams, 'w' ) do |file|
        pfams.each { |pfam|
            file.write( pfam )
            file.write( "\n" )
        }
    end
    puts( "number of associated pfams         : " + pfams.size.to_s )
    puts( "wrote assocations in sgd style to  : " + outfile_sgd_style )
    puts( "wrote assocations in basic style to: " + outfile_simple_map )
    puts( "wrote all associated pfams to      : " + outfile_all_pfams )
    puts( "OK")

end


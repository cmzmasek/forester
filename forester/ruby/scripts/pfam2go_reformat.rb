#!/usr/local/bin/ruby -w
#
# = pfam2go_reformat
#
# Copyright::  Copyright (C) 2006-2008 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: pfam2go_reformat.rb,v 1.4 2008/11/27 01:41:36 cmzmasek Exp $
#
# Reformat pfam2go to a "association" file suitable as input
# for microarray GO enrichment/overrepresentation-type analyses,
# and create a file listing all mapped Pfams as well.


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


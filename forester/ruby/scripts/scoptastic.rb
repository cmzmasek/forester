#!/usr/local/bin/ruby -w
#
# = scoptastic
#
# Copyright::  Copyright (C) 2008-2009 Christian M. Zmasek.
#              All rights reserved.
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: scoptastic.rb,v 1.3 2008/08/28 17:09:07 cmzmasek Exp $
#
# To create Pfam id to SCOP mappings, one for each of four levels of SCOP
# classification.
#
# Created 2008-06-25 in San Diego, CA, USA by CMZ
#
# Usage: scoptastic.rb scoptastic.rb <Pfam id to ac map file, e.g.
# pfam_summarize.rb output> <Pfam ac to SCOP classification map file> <Pfam id
# to SCOP outfile root>


require 'iconv'

module ForesterScripts

    if RUBY_VERSION !~ /1.9/
                      puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
                      exit( -1 )
                end     
    
    CLASS_LEVEL_SUFFIX       = "_SCOP_2_CLASS"
    FOLD_LEVEL_SUFFIX        = "_SCOP_3_FOLD"
    SUPERFAMILY_LEVEL_SUFFIX = "_SCOP_4_SUPERFAMILY"
    FAMILY_LEVEL_SUFFIX      = "_SCOP_5_FAMILY"

    SEP = "\t"
    LINE_DELIMITER  = "\n"

    if ( ARGV == nil || ARGV.length != 3 )
        puts( "usage: scoptastic.rb <Pfam id to ac map file, e.g. pfam_summarize.rb output> <Pfam ac to SCOP classification map file> <Pfam id to SCOP outfile root>" )
        exit( -1 )
    end

    pfam_id_to_ac   = ARGV[ 0 ]
    pfam_ac_to_scop = ARGV[ 1 ]
    outfile         = ARGV[ 2 ]

    if ( !File.exists?( pfam_id_to_ac ) )
        puts( "Pfam id to ac map file [" + pfam_id_to_ac + "] does not exist" )
        exit( -1 )
    end
    if ( !File.exists?( pfam_ac_to_scop ) )
        puts( "Pfam ac to SCOP classification map file [" + pfam_ac_to_scop + "] does not exist" )
        exit( -1 )
    end
    if ( File.exists?( outfile + CLASS_LEVEL_SUFFIX ) )
        puts( "Outfile [" + outfile + CLASS_LEVEL_SUFFIX + "] already exists" )
        exit( -1 )
    end
    if ( File.exists?( outfile +  FOLD_LEVEL_SUFFIX ) )
        puts( "Outfile [" + outfile +  FOLD_LEVEL_SUFFIX + "] already exists" )
        exit( -1 )
    end
    if ( File.exists?( outfile + SUPERFAMILY_LEVEL_SUFFIX ) )
        puts( "Outfile [" + outfile + SUPERFAMILY_LEVEL_SUFFIX + "] already exists" )
        exit( -1 )
    end
    if ( File.exists?( outfile + FAMILY_LEVEL_SUFFIX ) )
        puts( "Outfile [" + outfile + FAMILY_LEVEL_SUFFIX + "] already exists" )
        exit( -1 )
    end

    ic = Iconv.new( 'UTF-8//IGNORE', 'UTF-8' )

    pfam_ac_to_id_map = Hash.new

    pfam_ac_to_scop_map = Hash.new

    count = 0

    File.open( pfam_id_to_ac  ) do | file |
        while line = file.gets
            line = ic.iconv( line )
            if ( line !~ /^#/ && line =~ /\S/ )
                if ( line =~ /^(\S+)\s+(PF\d+)/ )
                    pfam_ac_to_id_map[ $2 ] = $1
                    count += 1
                else
                    puts( "Pfam id to ac map file [" + pfam_id_to_ac + "] format error [line: " + line + "]" )
                    exit( -1 )
                end
            end
        end
    end
    puts()
    puts( "Extracted #{count} Pfam id to ac mappings from file [#{pfam_id_to_ac}]" )

    count = 0
    File.open( pfam_ac_to_scop ) do | file |
        while line = file.gets
            line = ic.iconv( line )
            if ( line !~ /^#/ && line =~ /\S/ )
                if ( line =~ /^(PF\d+)\.?\d*\s+([a-z]\.\d+\.\d+\.\d+)/ )
                    pfam_ac_to_scop_map[ $1 ] = $2
                    count += 1
                else
                    puts( "Pfam ac to SCOP classification map file [" + pfam_ac_to_scop + "] format error [line: " + line + "]" )
                    exit( -1 )
                end
            end
        end
    end

    puts( "Extracted #{count} Pfam ac to SCOP classification mappings from file [#{pfam_ac_to_scop}]" )

    out_class_level = File.open( outfile + CLASS_LEVEL_SUFFIX, 'w' )
    out_fold_level = File.open( outfile + FOLD_LEVEL_SUFFIX  , 'w' )
    out_superfamily_level = File.open( outfile + SUPERFAMILY_LEVEL_SUFFIX, 'w' )
    out_family_level = File.open( outfile + FAMILY_LEVEL_SUFFIX, 'w' )

    count = 0
    pfam_ac_to_scop_map.each { | pfam_ac,scop |
        if ( pfam_ac_to_id_map.has_key?( pfam_ac ) )
            pfam_id = pfam_ac_to_id_map[ pfam_ac ]
            scop_split = scop.split( "\." )

            out_class_level.write( pfam_id )
            out_fold_level.write( pfam_id )
            out_superfamily_level.write( pfam_id )
            out_family_level.write( pfam_id )

            out_class_level.write( SEP )
            out_fold_level.write( SEP )
            out_superfamily_level.write( SEP )
            out_family_level.write( SEP )

            out_class_level.write( scop_split[ 0 ] )
            out_fold_level.write( scop_split[ 0 ] + "." + scop_split[ 1 ] )
            out_superfamily_level.write( scop_split[ 0 ] + "." + scop_split[ 1 ] + "." + scop_split[ 2 ] )
            out_family_level.write( scop )

            out_class_level.write( LINE_DELIMITER )
            out_fold_level.write( LINE_DELIMITER )
            out_superfamily_level.write( LINE_DELIMITER )
            out_family_level.write( LINE_DELIMITER )
            count += 1
        else
            puts( "Pfam ac #{pfam_ac} not found in Pfam id to ac map file [" + pfam_id_to_ac + "]" )
            exit( -1 )
        end
    }

    out_class_level.close
    out_fold_level.close
    out_superfamily_level.close
    out_family_level.close

    puts()
    puts( "Wrote #{count} Pfam id to SCOP mappings to files '#{outfile + CLASS_LEVEL_SUFFIX}', '#{outfile + FOLD_LEVEL_SUFFIX}', '#{outfile + SUPERFAMILY_LEVEL_SUFFIX}', and '#{ outfile + FAMILY_LEVEL_SUFFIX}'" )
    puts( "OK" )
    puts()

end # module ForesterScripts


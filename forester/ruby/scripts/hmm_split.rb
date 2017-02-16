#!/usr/local/bin/ruby -w
#
# = hmm_split
#
# Copyright::    Copyright (C) 2017 Christian M Zmasek
# License::      GNU Lesser General Public License (LGPL)

module ForesterScripts

  if ( ARGV == nil || ARGV.length != 3 )
    puts( "usage: hmm_split.rb <Pfam HMM file> <outfile suffix> <outdir>" )
    exit( -1 )
  end

  hmmfile = ARGV[ 0 ]
  suffix  = ARGV[ 1 ]
  outdir  = ARGV[ 2 ]

  if ( !File.exist?( outdir ) )
    puts( "outdir [" + outdir + "] does not exist" )
    exit( -1 )
  end
  if ( !File.exist?( hmmfile ) )
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
        if ( File.exist?( outfile ) )
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
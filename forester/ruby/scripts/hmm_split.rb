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
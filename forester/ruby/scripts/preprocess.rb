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

if ( ARGV == nil || ARGV.length != 2 )
  puts( "usage: preprocess.rb <input> <path to Pfam A HMMs>" )         
  exit( -1 )
end 

input = ARGV[ 0 ]

pfam = ARGV[ 1 ]

system( "hmmscan --nobias --domtblout " + input + "_hmmscan_260_10 -E 10 " + pfam + " "  + input + ".ni.fasta" )

system( "hsp " + input + "_hmmscan_260_10 " + input + "_hmmscan_260_10_domain_table" )

system( "d2f " + input + "_hmmscan_260_10_domain_table " + input + ".ni.fasta " + input + "_hmmscan_260_10.dff" )
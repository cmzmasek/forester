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

require 'lib/evo/sequence/sequence'
require 'lib/evo/msa/msa'
require 'lib/evo/msa/msa_factory'
require 'lib/evo/io/writer/fasta_writer'
require 'lib/evo/io/parser/fasta_parser'

module Evoruby

  input1 = ARGV[ 0 ]
  input2 = ARGV[ 1 ]

  f = MsaFactory.new()

  msa1 = nil
  msa2 = nil

  begin
    msa1 = f.create_msa_from_file( input1, FastaParser.new() )
    msa2 = f.create_msa_from_file( input2, FastaParser.new() )
  rescue Exception => e
    puts "error: " + e.to_s
    exit
  end

  only_in_1 = Msa.new
  only_in_2 = Msa.new
  in_both = Msa.new

  for i in 0 ... msa1.get_number_of_seqs()
    seq = msa1.get_sequence( i )
    name = seq.get_name
    if msa2.has?( name )
      in_both.add_sequence(seq)
    else
      only_in_1.add_sequence(seq)
    end
  end

  for j in 0 ... msa2.get_number_of_seqs()
    seq = msa2.get_sequence( j )
    name = seq.get_name
    unless msa1.has?( name )
      only_in_2.add_sequence(seq)
    end
  end

  puts "only in 1: " + only_in_1.get_number_of_seqs.to_s
  puts "only in 2: " + only_in_2.get_number_of_seqs.to_s
  puts "in both  : " + in_both.get_number_of_seqs.to_s

  w = FastaWriter.new
  w.write(only_in_1, "only_in_1.fasta")
  w.write(only_in_2, "only_in_2.fasta")
  w.write(in_both, "in_both.fasta")

end

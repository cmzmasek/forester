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

require 'lib/evo/io/parser/msa_parser'
require 'lib/evo/msa/msa'

module Evoruby

  class FastaParser < MsaParser

    def initialize
    end

    def parse( path )
      Util.check_file_for_readability( path )
      msa = Msa.new
      current_seq = String.new()
      name        = String.new()
      saw_first_seq = false
      File.open( path ) do | file |
        while line = file.gets
          
          line.encode!("UTF-8", :invalid => :replace, :undef => :replace, :replace => "?")          
          if can_ignore?( line, saw_first_seq )

          elsif line =~ /^\s*>\s*(.+)/
            saw_first_seq = true
            add_seq( name, current_seq, msa )
            name = $1
            current_seq = String.new()
          elsif line =~ /^\s*(.+)/
            if name.length < 1
              error_msg = "format error at: " + line
              raise IOError, error_msg
            end
            # was: seq = $1.rstrip
            seq =  $1.gsub(/\s+/, '')
            current_seq << seq
          else
            error_msg = "Unexpected line: " + line
            raise IOError, error_msg
          end
        end
      end
      add_seq( name, current_seq, msa )
      return msa
    end

    private

    def add_seq( name, seq, msa )
      if name.length > 0 && seq.length > 0
        msa.add( name, seq )
      end
    end

    def can_ignore?( line, saw_first_seq )
      return ( line !~ /\S/  ||
         line =~ /^\s*#/ ||
         ( !saw_first_seq && line =~/^\s*[^>]/ ) )
    end

  end # class FastaParser

end # module Evoruby

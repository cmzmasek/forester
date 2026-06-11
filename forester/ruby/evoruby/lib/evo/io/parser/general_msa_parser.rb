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
  class GeneralMsaParser < MsaParser
    def initialize
    end

    def parse( path )
      Util.check_file_for_readability( path )
      block                       = -1
      current_seq_index_per_block = -1
      current_name                = nil
      saw_ignorable = true
      is_first      = true
      msa = Msa.new

      File.open( path ) do | file |
        while line = file.gets
          line.encode!("UTF-8", :invalid => :replace, :undef => :replace, :replace => "?")
          if can_ignore?( line )
            saw_ignorable = true
          elsif ( is_first && is_program_name_line?( line ) )
          elsif( line =~ /^\S+\s+.+\s*$/ || line =~ /^\s+.+\s*$/ || line =~ /^\S+\s*$/ )
            if ( saw_ignorable )
              block += 1
              current_seq_index_per_block = -1
              saw_ignorable = false
            end
            current_seq_index_per_block += 1
            if ( line =~ /^(\S+)\s+(.+?)\s*$/ )
              name = $1
              seq  = $2.gsub( /\s/, '.' )
              a = msa.find_by_name( name, false, false )
              if ( a.length < 1 )
                msa.add( name, seq )
              elsif ( a.length == 1 )
                msa.get_sequence( a[ 0 ] ).append!( seq )
              else
                error_msg = "Unexpected error at line: " + line
                raise IOError, error_msg
              end
              current_name = name
            elsif ( line =~ /^\s+(.+?)\s*$/ )
              seq = $1.gsub( /\s/, '.' )
              a = msa.find_by_name( current_name, false, false )
              if ( a.length != 1  )
                error_msg = "Unexpected error at line: " + line
                raise IOError, error_msg
              else
                msa.get_sequence( a[ 0 ] ).append!( seq )
              end

            elsif ( line =~ /^(\S+)\s*$/ )
              seq = $1
              if block == 0
                error_msg = "First block cannot contain unnamed sequences: "  + line
                raise IOError, error_msg
              else
                msa.get_sequence( current_seq_index_per_block ).append!( seq )
              end
              current_name = nil
            end
          else
            error_msg = "Unexpected line: " + line
            raise IOError, error_msg
          end
          if ( is_first )
            is_first = false
          end
        end
      end
      return msa
    end # def parse( path )

    private

    def can_ignore?( line )
      return ( line !~ /[A-Za-z\-?\*_\.]/ ||
      line =~ /^\s+[*\.:]/ ||
      line =~ /^\s*#/ ||
      line =~ /^\s*%/ ||
      line =~ /^\s*\/\// ||
      line =~ /^\s*!!/  )
    end

    def is_program_name_line?( line )
      return ( line =~ /^CLUSTAL\s/ ||
      line =~ /^MUSCLE\s\(/ ||
      line =~ /^PROBCONS\s/ )
    end
  end # class GeneralMsaParser

end # module Evoruby

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

require 'lib/evo/util/util'

module Evoruby

  class UniprotParser

    ID = "ID"
    DE = "DE"
    DR = "DR"
    LAST = '//'

    def initialize

    end


    def parse( lines )
      de = []
      dr = []
      id = nil
      lines.each do | line |

        if line.include?( ID ) && line.index( ID ) == 0
          id = line.split[ 1 ]
        elsif id != nil
          if line.include?( LAST ) && line.index( LAST ) == 0
            e = UniprotEntry.new
            e.id = id
            e.de = de
            e.dr = dr
            return e
          else
            if line.include?( DE ) && line.index( DE ) == 0
              add( line, de )
            elsif line.include?( DR ) && line.index( DR ) == 0
              add( line, dr )
            end
          end
        end
      end
      return nil
    end

    private

    def add( line, ary )
      line =~/[A-Z]{2}\s+(.+)/
      ary << $1
    end


  end # class UniprotParser

  class UniprotEntry

    attr_accessor :id
    attr_accessor :ac
    attr_accessor :de
    attr_accessor :gn
    attr_accessor :os
    attr_accessor :ox
    attr_accessor :dr
    attr_accessor :pe
    attr_accessor :kw

    def get_pdb_ids
      ids = []
      if dr != nil
        dr.each do | dr |
          if dr != nil
            if dr =~ /PDB;\s+([A-Z0-9]{4});/
              ids << $1
            end
          end
        end
      end
      ids
    end

    def get_go_descriptions
      gos = []
      if dr != nil
        dr.each do | dr |
          if dr != nil
            if dr =~ /GO;\s+GO:\d+(.+);\s+([^;]+)/
              gos << $1
            end
          end
        end
      end
      gos
    end

    def get_full_name
     # DE   RecName: Full=Apoptosis regulator Bcl-2;
    end


     def get_reactome_descriptions
      s = []
      if dr != nil
        dr.each do | dr |
          if dr != nil
            if dr =~ /Reactome;\s+REACT_\d+;\s+([^.]+)/
              s << $1
            end
          end
        end
      end
      s
    end


  end








end # module Evoruby

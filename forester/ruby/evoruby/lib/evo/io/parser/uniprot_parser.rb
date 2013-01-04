#
# = lib/evo/io/parser/uniprot_parser - UniprotParser class
#
# Copyright::  Copyright (C) 2012 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id:  Exp $
#
# last modified: 121003


#require 'iconv'


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

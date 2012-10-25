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

    def initialize file
      Util.check_file_for_readability file
      @file = file
    end


    def parse( ids )
      entries = Hash.new
      de = []
      dr = []
      id = nil
      File.open( @file ).each do | line |
        if line.index( ID ) == 0
          #   puts line
          ids.each do | i |
            #puts " " + i
            if line.include?( i ) && line.split[ 1 ] == i
              id = i
              break
            end
          end
        end
        if id != nil
          if line.include?( LAST ) && line.index( LAST ) == 0
            e = UniprotEntry.new
            e.de = de
            e.dr = dr
            entries[ id ] = e
            puts id
            id = nil
            de = []
            dr = []
          else
            if line.include?( DE ) && line.index( DE ) == 0
              add( line, de )
            elsif line.include?( DR ) && line.index( DR ) == 0
              add( line, dr )
            end
          end
        end
      end
      entries
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

  end


end # module Evoruby

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
      #ic = Iconv.new( 'UTF-8//IGNORE', 'UTF-8' )
      entries = []
      de = []
      dr = []
      read = false
      File.open( @file ).each do | line |
        if line.index( ID ) == 0
          puts line 
          ids.each do | id |
            puts " " + id
            if line.index( id ) == 0
              read = true
              break
            end
          end
        end
        if read
          if line.index LAST == 0
            read = false
            e = UniprotEntry.new
            e.de = de
            e.dr = dr
            entries[ id ] = e
            de = []
            dr = []
          else
            if line.index DE == 0
              add( line, de )
            elsif line.index DR == 0
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

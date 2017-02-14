#
# = lib/evo/io/parser/fasta_parser - FastaParser class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)
#
# Last modified: 2017/02/07

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

#
# = lib/evo/io/parser/general_msa_parser - GeneralMsaParser class
#
# Copyright::  Copyright (C) 2006-2007 Christian M. Zmasek
# License::    GNU Lesser General Public License (LGPL)
#
# $Id: general_msa_parser.rb,v 1.8 2009/10/08 22:44:54 cmzmasek Exp $
#
# last modified: 2009/10/08

require 'lib/evo/io/parser/msa_parser'
require 'lib/evo/msa/msa'

require 'iconv'

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
            ic = Iconv.new( 'UTF-8//IGNORE', 'UTF-8' )
            File.open( path ) do | file |
                while line = file.gets
                    line = ic.iconv( line )
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
                                error_msg = "First block cannot contain unnamed sequences"
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

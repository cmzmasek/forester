#
# = lib/evo/util/util.rb - Util class
#
# Copyright::    Copyright (C) 2023 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)
#
# Last modified: 2023/10/10

require 'pathname'
require 'lib/evo/util/constants'

module Evoruby
  class Util
    def Util.canonical_path( parent, child = nil )
      if child == nil
        return  File.expand_path(Pathname.new(parent).cleanpath.to_s).to_s
      end

      s = nil
      if parent.end_with?('/')
        s = parent + child
      else
        s = parent + '/' + child
      end
      File.expand_path(Pathname.new(s).cleanpath.to_s).to_s
    end

    def Util.get_matching_files( files, prefix_pattern, suffix_pattern )
      matching_files = Array.new
      files.each { | file |
        if ( !File.directory?( file ) &&
        file !~ /^\./ &&
        file =~ /^#{prefix_pattern}.*#{suffix_pattern}$/ )
          matching_files << file
        end
      }
      matching_files
    end

    def Util.get_matching_file( dir_name, prefix, suffix )
      unless Dir.exist? dir_name
        raise IOError, "directory [#{dir_name}] does not exist"
      end

      all_files = Dir.entries(dir_name)

      if ( all_files.size <= 2 )
        raise IOError, "directory [#{dir_name}] is empty"
      end
      matching_files = Array.new
      all_files.each { | file |
        if  (( !File.directory?( file )) && (file.end_with? suffix))
          matching_files << file
        end
      }

      if ( matching_files.size == 0 )
        raise IOError, "no files ending with \"" + suffix + "\" found in [" + dir_name + "]"
      end
      my_prefix = prefix
      done = false
      more_than_one = false
      the_one = nil;

      loop do
        puts my_prefix
        matches = 0
        matching_files.each { | file |
          if file.start_with?( my_prefix + "." ) # was  if file.start_with?( my_prefix )
            matches += 1
            if matches > 1
              the_one = nil
              break
            end
            the_one = file
          end
        }
        if matches > 1
          more_than_one = true
          done = true
        end
        if matches == 1
          done = true
        else
          if my_prefix.length <= 1
            raise IOError, "no file matching \"" + prefix + "\" and ending with \"" + suffix + "\" found in [" + dir_name + "]"
          end
          my_prefix = my_prefix[ 0 ... ( my_prefix.length - 1 ) ]

        end
        break if done
      end

      if more_than_one
        raise IOError, "multiple files matching \"" +  prefix  +
        "\" and ending with \"" + suffix + "\" found in [" + dir_name + "]"
      elsif the_one != nil
      else
        raise IOError, "no file matching \"" + prefix  + "\" and ending with \"" +
        suffix + "\" found in [" + dir_name + "]"
      end
      the_one
    end

    def Util.normalize_seq_name( name, length, exception_if_too_long = false )
      if name.length > length
        if exception_if_too_long
          error_msg = "sequence name \"#{name}\" is too long (>#{length})"
          raise StandardError, error_msg
        end
        name = name[ 0, length ]
      elsif name.length < length
        t = length - name.length
        t.times do
          name = name + " "
        end
      end
      name
    end

    # Returns true if char_code corresponds to: space * - . _
    def Util.is_aa_gap_character?( char_code )
      return ( char_code <= 32  || char_code == 42 || char_code == 45 || char_code == 46 ||char_code == 95  )
    end

    # Deletes *, digits, and whitespace, replaces BJOUZ? with X, and replaces non-(letters, -) with -
    def Util.clean_seq_str( seq_str )
      seq_str = seq_str.upcase
      seq_str = seq_str.gsub( /\s+/, '' )
      seq_str = seq_str.gsub( /\d+/, '' )
      seq_str = seq_str.gsub( '*', '' )
      seq_str = seq_str.gsub( /[BJOUZ?]/, 'X' )
      seq_str = seq_str.gsub( /[^A-Z\-]/, '-' )
      seq_str
    end

    # raises ArgumentError
    def Util.check_file_for_readability( path )
      unless ( File.exist?( path ) )
        error_msg = "file [#{path}] does not exist"
        raise IOError, error_msg
      end
      unless ( File.file?( path ) )
        error_msg = "file [#{path}] is not a regular file"
        raise IOError, error_msg
      end
      unless ( File.readable?( path ) )
        error_msg = "file [#{path}] is not a readable file"
        raise IOError, error_msg
      end
      if ( File.zero?( path ) )
        error_msg = "file [#{path}] is empty"
        raise IOError, error_msg
      end
    end

    # raises ArgumentError
    def Util.check_file_for_writability( path )
      if File.directory?( path )
        error_msg = "file [#{path}] is an existing directory"
        raise IOError, error_msg
      elsif File.exist?( path )
        error_msg = "file [#{path}] already exists"
        raise IOError, error_msg
      elsif File.writable?( path )
        error_msg = "file [#{path}] is not writeable"
        raise IOError, error_msg
      end
    end

    def Util.fatal_error_if_not_writable( prg_name, path )
      begin
        Util.check_file_for_writability( path )
      rescue IOError => e
        Util.fatal_error( prg_name, e.to_s )
      end
    end

    def Util.fatal_error_if_not_readable( prg_name, path )
      begin
        Util.check_file_for_readability( path )
      rescue IOError => e
        Util.fatal_error( prg_name, e.to_s )
      end
    end

    def Util.get_env_variable_value( env_variable )
      value = ENV[env_variable]
      if value == nil || value.empty?
        error_msg = "apparently environment variable #{env_variable} has not been set"
        raise StandardError, error_msg
      end
      value
    end

    # raises ArgumentError
    def Util.file2array( path, split_by_semicolon )
      Util.check_file_for_readability( path )
      a = Array.new()
      c = 0
      File.open( path ) do | file |
        while line = file.gets
          if ( line =~ /^\s*(\S.*?)\s*$/ )
            s = $1
            if ( split_by_semicolon && s =~/;/ )
              sa = s.split( /;/ )
              for i in 0 ... sa.length()
                a[ c ] = sa[ i ].strip!
              end
            else
              a[ c ] = s
            end
            c += 1
          end
        end
      end
      return a
    end

    def Util.print_program_information( prg_name,
      prg_version,
      prg_desc,
      date,
      www,
      io = STDOUT )

      ruby_version = RUBY_VERSION
      l = prg_name.length + prg_version.length + date.length + ruby_version.length + 12
      io.print( Evoruby::Constants::LINE_DELIMITER )
      io.print( prg_name + " " + prg_version + " [" + date + "] [ruby " + ruby_version + "]")
      io.print( Evoruby::Constants::LINE_DELIMITER )
      l.times {
        io.print( "_" )
      }
      io.print( Constants::LINE_DELIMITER )
      io.print( Constants::LINE_DELIMITER )
      io.print( prg_desc )
      io.print( Constants::LINE_DELIMITER )
      io.print( Constants::LINE_DELIMITER )
      io.print( "Website: " + www )
      io.print( Constants::LINE_DELIMITER )
      io.print( Constants::LINE_DELIMITER )
    end

    def Util.fatal_error( prg_name, message, io = STDOUT )
      io.print( Constants::LINE_DELIMITER )
      if ( !Util.is_string_empty?( prg_name ) )
        io.print( "[" + prg_name + "] > " + message )
      else
        io.print( " > " + message )
      end
      io.print( Constants::LINE_DELIMITER )
      io.print( Constants::LINE_DELIMITER )
      exit( -1 )
    end

    def Util.print_message( prg_name, message, io = STDOUT )
      if ( !Util.is_string_empty?( prg_name ) )
        io.print( "[" + prg_name + "] > " + message )
      else
        io.print( " > " + message )
      end
      io.print( Constants::LINE_DELIMITER )
    end

    def Util.print_warning_message( prg_name, message, io = STDOUT )
      if ( !Util.is_string_empty?( prg_name ) )
        io.print( "[" + prg_name + "] > WARNING: " + message )
      else
        io.print( " > " + message )
      end
      io.print( Constants::LINE_DELIMITER )
    end

    def Util.is_string_empty?( s )
      return ( s == nil || s.length < 1 )
    end

    # From "Ruby Cookbook"
    # counts_hash: key is a "name", value is the count (integer)
    def Util.draw_histogram( counts_hash, char = "#" )
      pairs = counts_hash.keys.collect { |x| [ x.to_s, counts_hash[ x ] ] }.sort
      largest_key_size = pairs.max { |x, y| x[ 0 ].size <=> y[ 0 ].size }[ 0 ].size
      pairs.inject( "" ) do | s, kv |
        s << "#{ kv[ 0 ].ljust( largest_key_size ) }  | #{ char*kv[ 1 ] }" + Constants::LINE_DELIMITER
      end
    end

    def Util.looks_like_fasta?( path )
      Util.check_file_for_readability( path )
      File.open( path ) do | file |
        while line = file.gets
          if ( line !~ /\S/ || line =~ /^\s*#/ )
          elsif line =~ /^\s*>\s*(.+)/
            return true
          else
            return false
          end
        end
      end
      error_msg = "unexpected format"
      raise IOError, error_msg
    end

  end # class Util

end # module Evoruby

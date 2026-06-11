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

module ForesterScripts

  if RUBY_VERSION !~ /1.9/
    puts( "Your ruby version is #{RUBY_VERSION}, expected 1.9.x " )
    exit( -1 )
  end

  PARAMETER_FILE = 'parameters.rb_dir_qsub'
  SLEEP          = 1.0
  REMOVE_SUFFIX  = true
  LIMIT_TO_TAX_CODE  = true

  PRG         = 'PRG:'
  OPT         = 'OPT:'
  VOPT        = 'VOPT:'
  OUTPUT_OPT  = 'OUTPUT:'
  SUFFIX      = 'SUFFIX:'
  INPUT_PART  = 'INPUT_PART:'


  PBS_O_WORKDIR       = '$PBS_O_WORKDIR/'
  TMP_CMD_FILE_SUFFIX = '__QSUB'
  NAME                = 'rb_dir_qsub'

  if ( !File.exists?( PARAMETER_FILE ) )
    puts( '[' + NAME + '] > parameters file "' + PARAMETER_FILE + '" not found' )
    Process.exit!
  end
  puts( '[' + NAME + '] > reading ' + PARAMETER_FILE )

  prg = ''
  opt = ''
  vopts = Array.new
  suffix = ''
  input_part = ''
  output_opt = ''
  open( PARAMETER_FILE ).each { |line|
    if ( line.length > 1 && line =~ /^[^#]\S+/ )
      if line =~ /^#{PRG}\s+(\S+)/
        prg = $1
      end
      if line =~ /^\s*#{OPT}\s+(\S+.+)/
        opt = $1
      end
      if line =~ /^\s*#{VOPT}\s+(\S+.+)/
        vopts.push( $1 )
      end
      if line =~ /^\s*#{SUFFIX}\s+(\S+)/
        suffix = $1
      end
      if line =~ /^\s*#{INPUT_PART}\s+(\S+)/
        input_part = $1
      end
      if line =~ /^\s*#{OUTPUT_OPT}\s+(\S+.+)/
        output_opt = $1
      end
    end
  }
  if ( prg.length < 1 )
    puts( '[' + NAME + '] > no program name found in parameters file "' + PARAMETER_FILE + '"' )
    Process.exit!
  end
  puts( '[' + NAME + '] > program: ' + prg )
  puts( '[' + NAME + '] > option :  ' + opt )
  vopts.each { |vopt|
    puts( '[' + NAME + '] > voption:  ' + vopt )
  }
  puts( '[' + NAME + '] > suffix :  ' + suffix )
  if ( input_part.length > 0 )
    puts( '[' + NAME + '] > input:  ' + input_part )
  end
  if ( output_opt.length > 0 )
    puts( '[' + NAME + '] > output opt :  ' + output_opt )
  end
  if vopts.empty?
    vopts.push( "" )
  end

  files = Dir.entries( "." )

  files.each { |file|
    if ( !File.directory?( file ) && file !~ /^\./ && file !~ /#{PARAMETER_FILE}/ )

      if ( input_part.length > 0 && file !~ /#{input_part}/ )
        next
      end
      vopts.each { |vopt|
        cmd = ""
        outputfile = file.to_str
        if LIMIT_TO_TAX_CODE
          if outputfile =~ /^([A-Z0-9]{3,5})[_\.]/
            outputfile = $1
          end
        end
        if REMOVE_SUFFIX
          if outputfile =~ /(.+)\..{1,5}/
            outputfile = $1
          end
        end
        if output_opt.length > 0
          cmd = prg + ' ' +
           output_opt + ' ' + PBS_O_WORKDIR + outputfile + suffix + ' ' +
           opt + ' ' +
           PBS_O_WORKDIR + file.to_str +
           ' > /dev/null'
        elsif vopt.length > 0
          cmd = prg + ' ' + opt + ' ' + vopt + ' ' + PBS_O_WORKDIR + file.to_str +
           ' > ' + PBS_O_WORKDIR + vopt + "_" + outputfile + suffix
        else
          cmd = prg + ' ' + opt + ' ' + PBS_O_WORKDIR + file.to_str +
           ' > ' + PBS_O_WORKDIR + outputfile + suffix
        end
        tmp_cmd_file = file.to_str + TMP_CMD_FILE_SUFFIX
        if File.exists?( tmp_cmd_file )
          File.delete( tmp_cmd_file )
        end
        open( tmp_cmd_file, 'w' ) do |f|
          f.write( cmd )
        end
        puts( '[' + NAME + '] > excuting ' + cmd )
        IO.popen( 'qsub ' + tmp_cmd_file , 'r+' ) do |pipe|
          pipe.close_write
          puts pipe.read
        end
        sleep( SLEEP )
        if File.exists?( tmp_cmd_file )
          File.delete( tmp_cmd_file )
        end
      }
    end
  }
  puts( '[' + NAME + '] > OK.' )
  puts

end

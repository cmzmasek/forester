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

NEW_DIR = '/home/czmasek/WORK/GENOME_HMMPFAM/HMMSCAN3_no_bias_260/'
OLD_DIR = '/home/czmasek/WORK/GENOME_HMMPFAM/HMMSCAN3_no_bias_250/'

# this need to be run in the dir where the links to be moved are
# need to create a dir named 'newdir' first

Dir.foreach('.') { |f| 
  if File.symlink?( f )
     link = File.readlink( f )
     puts f + ' -> ' + link 
   
     link =~ /\/([^\/]+)\/([^\/]+)\.hmmscan_250/
     group = $1
     
     core_name = $2
     puts '  => ' + group + ' / ' + core_name
     new_link = NEW_DIR + group + '/' + core_name + '.hmmscan_260'
     puts '  => ' + new_link
     puts
     if ( !File.exists?( new_link ) )
        puts 'ERROR ' + new_link.to_s + ' does not exist'  
        exit
     end
     if ( !File.exists?( link ) )
        puts 'ERROR ' + link.to_s + ' does not exist' 
        exit
     end
     #what is called new_link will by linked to by 'newdir/' + f.to_s:
     File.symlink( new_link, 'newdir/' + f.to_s )
  end
}


 

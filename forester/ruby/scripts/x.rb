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

PRG = "/home/czmasek/SOFTWARE/FORESTER/DEV/forester/forester/ruby/evoruby/exe/hsp.rb"
Dir.foreach(".") do |f|
  next if f == '.' or f == '..'
  if f =~ /^([A-Z0-9]{3,6})\.hmmscan_260/
    taxon = $1
    cmd = "#{PRG} -ie=10 -a -pe=1 -m=RRMa -s=#{taxon} #{taxon}.hmmscan_260 #{taxon}_ie10_domain_table > #{taxon}_ie10_summary_table.txt"
    puts cmd
    output = IO.popen( cmd )
  elsif f =~ /^([A-Z0-9]{3,6})_cdhit090\.hmmscan_260/
    taxon = $1
    cmd = "#{PRG} -ie=10 -a -pe=1 -m=RRMa -s=#{taxon} #{taxon}_cdhit090.hmmscan_260 #{taxon}_cdhit090_ie10_domain_table > #{taxon}_cdhit090_ie10_summary_table.txt"
    puts cmd
    output = IO.popen( cmd )
  end
end
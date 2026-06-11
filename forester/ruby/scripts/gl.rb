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

File.open( "g.txt" ) do | file |
  while line = file.gets
    if line =~ /^([A-Z0-9]{3,5})/
      puts $1 + "\t" + "/home/czmasek/DATA/SEQUENCES/GENOMES/PROTEIN_PREDICTIONS_CDHIT_095/ALL/" + $1 + ".fasta"
    else
      puts line
    end
  end
end
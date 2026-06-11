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

require 'net/ftp'

EMAIL           = '?'
PUB_RELEASE_DIR = '/pub/release-65/fasta'
PEP_DIR         = '/pep'

ftp = Net::FTP.new('ftp.ensembl.org', 'anonymous', EMAIL)
ftp.passive = true # To avoid "No route to host" error.
ftp.chdir( PUB_RELEASE_DIR )
files = ftp.list('*_*') # To only look at files with an underscore.
count = 0
files.each do | file |
  species = file.split().last
  begin
    ftp.chdir(species + PEP_DIR)
    pepfiles = ftp.list()
    pepfiles.each do | pepfile |
      pepfile = pepfile.split().last
      if pepfile =~ /all.fa.gz/ # Only want the "all.fa.gz" files (and not the
                                # "abinitio" files).
        ftp.getbinaryfile(pepfile)
        puts 'downloaded "' + pepfile + '"'
        count += 1
      end
    end
  rescue Exception
    puts 'ignoring "' + species + '"'
  end
  ftp.chdir(PUB_RELEASE_DIR) # To go back to the starting directory.
end
ftp.close
puts 'done (downloaded ' + count.to_s + ' files)'
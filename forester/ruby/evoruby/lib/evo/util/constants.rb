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

module Evoruby
  class Constants

    VERBOSE = true

    EVORUBY_VERSION = '1.1'

    ID_NORMALIZED_FASTA_FILE_SUFFIX          = "_ni.fasta"
    ID_MAP_FILE_SUFFIX                       = ".nim"
    DOMAIN_TABLE_SUFFIX                      = "_domain_table"
    HMMSCAN                                  = "_hmmscan"
    DOMAINS_TO_FORESTER_OUTFILE_SUFFIX       = ".dff"
    DOMAINS_TO_FORESTER_EVALUE_CUTOFF_SUFFIX = "_dtfE"
    
    PFAM_V_FOR_EX                             = "300" # Pfam version for examples

    FORESTER_HOME_ENV_VARIABLE = 'FORESTER_HOME'

    EVORUBY         = 'evoruby'

    LINE_DELIMITER  = "\n"

    FILE_SEPARATOR  = File::SEPARATOR

    DOMAIN_STRUCTURE_NHX_SEPARATOR = '>'

  end # class Constants

end # module Evoruby

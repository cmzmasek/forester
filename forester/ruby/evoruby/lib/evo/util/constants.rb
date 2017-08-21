#
# = lib/evo/util/constants.rb - Constants class
#
# Copyright::    Copyright (C) 2017 Christian M. Zmasek
# License::      GNU Lesser General Public License (LGPL)

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
    JAVA_HOME_ENV_VARIABLE     = 'JAVA_HOME'

    EVORUBY         = 'evoruby'

    LINE_DELIMITER  = "\n"

    FILE_SEPARATOR  = File::SEPARATOR

    DOMAIN_STRUCTURE_NHX_SEPARATOR = '>'

  end # class Constants

end # module Evoruby

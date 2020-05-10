
package org.forester.surfacing;

final class GlobalOptions {

    private static int     _verbosity                                                  = 0;
    private static boolean _uniprot_priority_for_accessor_parsing                      = false;
    private static int     _obtain_names_for_das_from_db_max_ids_to_search_per_species = 0;
    static int getVerbosity() {
        return _verbosity;
    }

    static void setVerbosity( final int verbosity ) {
        if ( verbosity < 0 ) {
            _verbosity = 0;
        }
        else if ( verbosity > 3 ) {
            _verbosity = 3;
        }
        else {
            _verbosity = verbosity;
        }
    }

    static boolean isUniprotPriorityForAccessorParsing() {
        return _uniprot_priority_for_accessor_parsing;
    }

    static void setUniprotPriorityForAccessorParsing( final boolean uniprot_priority_for_accessor_parsing ) {
        _uniprot_priority_for_accessor_parsing = uniprot_priority_for_accessor_parsing;
    }

    static int getObtainNamesForDasFromDbMaxIdsToSearchPerSpecies() {
        return _obtain_names_for_das_from_db_max_ids_to_search_per_species;
    }

    static void setObtainNamesForDasFromDbMaxIdsToSearchPerSpecies( final int obtain_names_for_das_from_db_max_ids_to_search_per_species ) {
        if ( obtain_names_for_das_from_db_max_ids_to_search_per_species < 1 ) {
            _obtain_names_for_das_from_db_max_ids_to_search_per_species = 1;
        }
        else if ( obtain_names_for_das_from_db_max_ids_to_search_per_species > 1000 ) {
            _obtain_names_for_das_from_db_max_ids_to_search_per_species = 1000;
        }
        else {
            _obtain_names_for_das_from_db_max_ids_to_search_per_species = obtain_names_for_das_from_db_max_ids_to_search_per_species;
        }
    }
}

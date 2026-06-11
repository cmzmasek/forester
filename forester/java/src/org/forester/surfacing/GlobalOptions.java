// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

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

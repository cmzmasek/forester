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

package org.forester.ws.seqdb;

/**
 * Just the organism a sequence accession belongs to: its NCBI tax-id and scientific name (either may be
 * {@code null}). The deliberately tiny return of {@link OrganismSource} -- so resolving "which organism
 * is this protein?" never has to retrieve or retain the full sequence record. Immutable.
 */
public final class Organism {

    /** Nothing resolved. */
    public static final Organism EMPTY = new Organism( null, null );

    private final String _tax_id;
    private final String _scientific_name;

    public Organism( final String tax_id, final String scientific_name ) {
        _tax_id = tax_id;
        _scientific_name = scientific_name;
    }

    /** The organism's NCBI tax-id (numeric string), or {@code null}. */
    public String getTaxId() {
        return _tax_id;
    }

    public String getScientificName() {
        return _scientific_name;
    }

    public boolean isEmpty() {
        return isBlank( _tax_id ) && isBlank( _scientific_name );
    }

    private static boolean isBlank( final String s ) {
        return ( s == null ) || s.trim().isEmpty();
    }

    @Override
    public String toString() {
        return "Organism[" + _scientific_name + " (" + _tax_id + ")]";
    }
}

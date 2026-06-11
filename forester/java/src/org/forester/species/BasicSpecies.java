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

package org.forester.species;

import org.forester.util.ForesterUtil;

public class BasicSpecies implements Species {

    final private String _species_id;

    public BasicSpecies( final String species_id ) {
        if ( ForesterUtil.isEmpty( species_id ) ) {
            throw new IllegalArgumentException( "attempt to create new species from empty or null string" );
        }
        _species_id = species_id.trim();
    }

    @Override
    public int compareTo( final Species species ) {
        if ( this == species ) {
            return 0;
        }
        return getSpeciesId().toLowerCase().compareTo( species.getSpeciesId().toLowerCase() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return getSpeciesId().equals( ( ( Species ) o ).getSpeciesId() );
        }
    }

    /* (non-Javadoc)
     * @see org.forester.surfacing.Species#getSpeciesId()
     */
    @Override
    public String getSpeciesId() {
        return _species_id;
    }

    @Override
    public int hashCode() {
        return getSpeciesId().hashCode();
    }

    @Override
    public String toString() {
        return getSpeciesId();
    }
}

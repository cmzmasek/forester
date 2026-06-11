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

import java.util.ArrayList;
import java.util.List;

import org.forester.protein.BinaryDomainCombination;
import org.forester.species.Species;

final class DirectedCombinableDomains extends BasicCombinableDomains {

    public DirectedCombinableDomains( final String n_terminal_key_domain, final Species species ) {
        super( n_terminal_key_domain, species );
    }

    @Override
    public List<BinaryDomainCombination> toBinaryDomainCombinations() {
        final List<BinaryDomainCombination> binary_combinations = new ArrayList<BinaryDomainCombination>( getNumberOfCombinableDomains() );
        for( final String domain : getCombiningDomains().keySet() ) {
            // Precondition (!): key domain is most upstream domain.
            //TODO ensure this is true.
            binary_combinations.add( DirectedBinaryDomainCombination.obtainInstance( getKeyDomain(), domain ) );
        }
        return binary_combinations;
    }
}

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

import java.util.List;

final class CombinationsBasedPairwiseDomainSimilarityCalculator implements PairwiseDomainSimilarityCalculator {

    @Override
    public PairwiseDomainSimilarity calculateSimilarity( final CombinableDomains domains_1,
                                                         final CombinableDomains domains_2 ) {
        if ( !domains_1.getKeyDomain().equals( domains_2.getKeyDomain() ) ) {
            throw new IllegalArgumentException( "attempt to calculate similarity between domain collection with different keys" );
        }
        final List<String> d1 = domains_1.getCombinableDomains();
        final List<String> d2 = domains_2.getCombinableDomains();
        int same = 0;
        int different = 0;
        for( final String domain : d1 ) {
            if ( d2.contains( domain ) ) {
                same++;
            }
            else {
                different++;
            }
        }
        for( final String domain : d2 ) {
            if ( !( d1.contains( domain ) ) ) {
                different++;
            }
        }
        final int difference = domains_1.getNumberOfCombinableDomains() - domains_2.getNumberOfCombinableDomains();
        return new CombinationsBasedPairwiseDomainSimilarity( same, different, difference );
    }
}

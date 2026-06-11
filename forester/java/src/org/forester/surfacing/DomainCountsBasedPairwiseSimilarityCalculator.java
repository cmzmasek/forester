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

final class DomainCountsBasedPairwiseSimilarityCalculator implements PairwiseDomainSimilarityCalculator {

    @Override
    public PairwiseDomainSimilarity calculateSimilarity( final CombinableDomains domains_1,
                                                         final CombinableDomains domains_2 ) {
        if ( !domains_1.getKeyDomain().equals( domains_2.getKeyDomain() ) ) {
            throw new IllegalArgumentException( "attempt to calculate similarity between domain collection with different keys" );
        }
        if ( ( domains_1.getKeyDomainCount() > Short.MAX_VALUE ) || ( domains_2.getKeyDomainCount() > Short.MAX_VALUE )
                || ( ( domains_1.getKeyDomainCount() + domains_2.getKeyDomainCount() ) > Short.MAX_VALUE ) ) {
            throw new IllegalArgumentException( "too large for short!" );
        }
        final short dc1 = ( short ) domains_1.getKeyDomainCount();
        final short dc2 = ( short ) domains_2.getKeyDomainCount();
        return new CountsBasedPairwiseDomainSimilarity( ( short ) ( dc1 - dc2 ), ( short ) ( dc1 + dc2 ) );
    }
}

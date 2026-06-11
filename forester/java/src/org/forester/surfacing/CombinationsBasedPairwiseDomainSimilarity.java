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

final class CombinationsBasedPairwiseDomainSimilarity implements PairwiseDomainSimilarity {

    private final int    _difference_in_counts;
    private final int    _different_domains;
    private final int    _same_domains;
    private final double _score;

    public CombinationsBasedPairwiseDomainSimilarity( final int same_domains,
                                                      final int different_domains,
                                                      final int difference_in_counts ) {
        if ( ( same_domains < 0 ) || ( different_domains < 0 ) ) {
            throw new IllegalArgumentException( "attempt to use domain counts less than 0" );
        }
        _difference_in_counts = difference_in_counts;
        _same_domains = same_domains;
        _different_domains = different_domains;
        if ( _different_domains == 0 ) {
            _score = 1.0;
        }
        else {
            _score = ( double ) _same_domains / ( _different_domains + _same_domains );
        }
    }

    @Override
    public int getDifferenceInCounts() {
        return _difference_in_counts;
    }

    public int getNumberOfDifferentDomains() {
        return _different_domains;
    }

    public int getNumberOfSameDomains() {
        return _same_domains;
    }

    @Override
    public double getSimilarityScore() {
        return _score;
    }
}

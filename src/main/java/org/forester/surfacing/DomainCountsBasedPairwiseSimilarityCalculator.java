// $Id:
// 04:20:19 cmzmasek Exp $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.surfacing;

public class DomainCountsBasedPairwiseSimilarityCalculator implements PairwiseDomainSimilarityCalculator {

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

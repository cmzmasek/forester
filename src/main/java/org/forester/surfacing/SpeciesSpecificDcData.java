// $Id:
// cmzmasek Exp $
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

import java.util.SortedMap;
import java.util.SortedSet;

/*
 * A helper class for PrintableDomainSimilarity.
 */
interface SpeciesSpecificDcData {

    public void addProteinsExhibitingCombinationCount( final String domain_id, final int count );

    /**
     * This should return a sorted map mapping domain ids to their corresponding
     * counts
     *
     * @return a sorted map mapping domain ids to their corresponding counts
     */
    public SortedMap<String, Integer> getCombinableDomainIdToCountsMap();

    public SortedSet<String> getKeyDomainProteins();

    public int getNumberOfProteinsExhibitingCombinationWith( final String domain_id );

    public StringBuffer toStringBuffer( final DomainSimilarityCalculator.Detailedness detailedness, boolean html );

    void addKeyDomainProtein( String protein );
}

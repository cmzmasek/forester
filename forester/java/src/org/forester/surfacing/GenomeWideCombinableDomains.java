// $Id:
// $
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

import org.forester.protein.BinaryDomainCombination;
import org.forester.protein.BinaryDomainCombination.DomainCombinationType;
import org.forester.species.Species;
import org.forester.util.DescriptiveStatistics;

public interface GenomeWideCombinableDomains {

    public boolean contains( String key_id );

    public CombinableDomains get( String key_id );

    public SortedMap<String, CombinableDomains> getAllCombinableDomainsIds();

    /**
     * This should return all domains ids present in the genome.
     * 
     * @return a sorted set of domains ids
     */
    public SortedSet<String> getAllDomainIds();

    public DomainCombinationType getDomainCombinationType();

    /**
     * This should return a statistic for per domain 
     * promiscuity in a genome.
     * 
     * @return descriptive statistics for per domain promiscuity in a genome
     */
    public DescriptiveStatistics getPerGenomeDomainPromiscuityStatistics();

    public int getSize();

    public Species getSpecies();

    /**
     * This should return all binary domain combinations present in the genome.
     * 
     * @return a sorted set of binary domain combinations
     */
    public SortedSet<BinaryDomainCombination> toBinaryDomainCombinations();

    public StringBuilder toStringBuilder( GenomeWideCombinableDomainsSortOrder order );

    SortedSet<String> getMostPromiscuosDomain();

    public static enum GenomeWideCombinableDomainsSortOrder {
        ALPHABETICAL_KEY_ID, COMBINATIONS_COUNT, KEY_DOMAIN_COUNT, KEY_DOMAIN_PROTEINS_COUNT
    }
}

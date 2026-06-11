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

import java.util.SortedMap;
import java.util.SortedSet;

import org.forester.protein.BinaryDomainCombination;
import org.forester.protein.BinaryDomainCombination.DomainCombinationType;
import org.forester.species.Species;
import org.forester.util.DescriptiveStatistics;

interface GenomeWideCombinableDomains {

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

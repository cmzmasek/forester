// $Id:
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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.protein.BinaryDomainCombination;
import org.forester.species.Species;
import org.forester.util.ForesterUtil;

public class BasicCombinableDomains implements CombinableDomains {

    final private String                   _key_domain;
    private int                            _key_domain_count;
    final private Species                  _species;
    final private TreeMap<String, Integer> _combining_domains;
    final private Set<String>              _proteins_with_key_domain;

    public BasicCombinableDomains( final String key_domain, final Species species ) {
        _key_domain = key_domain;
        _species = species;
        _combining_domains = new TreeMap<String, Integer>();
        _proteins_with_key_domain = new HashSet<String>();
        _key_domain_count = 0;
    }

    @Override
    public void addCombinableDomain( final String protein_domain ) {
        if ( getCombiningDomains().containsKey( protein_domain ) ) {
            getCombiningDomains().put( protein_domain, getCombiningDomains().get( protein_domain ) + 1 );
        }
        else {
            getCombiningDomains().put( protein_domain, 1 );
        }
    }

    @Override
    public void addKeyDomainProtein( final String protein ) {
        if ( ForesterUtil.isEmpty( protein ) ) {
            throw new IllegalArgumentException( "attempt to add null or empty protein" );
        }
        getKeyDomainProteins().add( protein );
    }

    @Override
    public List<String> getAllDomains() {
        final List<String> domains = getCombinableDomains();
        if ( !domains.contains( getKeyDomain() ) ) {
            domains.add( getKeyDomain() );
        }
        return domains;
    }

    @Override
    public List<String> getCombinableDomains() {
        final List<String> domains = new ArrayList<String>( getNumberOfCombinableDomains() );
        for( final String domain : getCombiningDomains().keySet() ) {
            domains.add( domain );
        }
        return domains;
    }

    @Override
    public SortedMap<String, Integer> getCombinableDomainsIds() {
        final SortedMap<String, Integer> ids = new TreeMap<String, Integer>();
        for( final String domain : getCombiningDomains().keySet() ) {
            final String pd = domain;
            ids.put( pd, getCombiningDomains().get( pd ) );
        }
        return ids;
    }

    @Override
    public StringBuilder getCombiningDomainIdsAsStringBuilder() {
        final StringBuilder sb = new StringBuilder();
        for( final Iterator<String> iter = getCombiningDomains().keySet().iterator(); iter.hasNext(); ) {
            final String key = iter.next();
            sb.append( key.toString() );
            sb.append( " [" );
            final int count = getCombiningDomains().get( key );
            sb.append( count );
            sb.append( "]" );
            if ( iter.hasNext() ) {
                sb.append( ", " );
            }
        }
        return sb;
    }

    protected TreeMap<String, Integer> getCombiningDomains() {
        return _combining_domains;
    }

    @Override
    public String getKeyDomain() {
        return _key_domain;
    }

    @Override
    public int getKeyDomainCount() {
        return _key_domain_count;
    }

    @Override
    public int getKeyDomainProteinsCount() {
        return getKeyDomainProteins().size();
    }

    @Override
    public int getNumberOfCombinableDomains() {
        return _combining_domains.size();
    }

    @Override
    public int getNumberOfProteinsExhibitingCombination( final String protein_domain ) {
        if ( getCombiningDomains().containsKey( protein_domain ) ) {
            return getCombiningDomains().get( protein_domain );
        }
        else {
            return 0;
        }
    }

    @Override
    public Species getSpecies() {
        return _species;
    }

    @Override
    public boolean isCombinable( final String protein_domain ) {
        return getCombiningDomains().containsKey( protein_domain );
    }

    @Override
    public void setKeyDomainCount( final int key_domain_count ) {
        _key_domain_count = key_domain_count;
    }

    @Override
    public List<BinaryDomainCombination> toBinaryDomainCombinations() {
        final List<BinaryDomainCombination> binary_combinations = new ArrayList<BinaryDomainCombination>( getNumberOfCombinableDomains() );
        for( final String domain : getCombiningDomains().keySet() ) {
            binary_combinations.add( new BasicBinaryDomainCombination( getKeyDomain(), domain ) );
        }
        return binary_combinations;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( getKeyDomain() );
        sb.append( " [" );
        sb.append( getKeyDomainCount() );
        sb.append( ", " );
        sb.append( getKeyDomainProteinsCount() );
        sb.append( ", " );
        sb.append( getNumberOfCombinableDomains() );
        sb.append( "]: " );
        sb.append( getCombiningDomainIdsAsStringBuilder() );
        return sb.toString();
    }

    @Override
    public Set<String> getKeyDomainProteins() {
        return _proteins_with_key_domain;
    }
}

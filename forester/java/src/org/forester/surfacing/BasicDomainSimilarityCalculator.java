// $Id:
// Exp $
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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.species.Species;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class BasicDomainSimilarityCalculator implements DomainSimilarityCalculator {

    final DomainSimilarity.DomainSimilaritySortField _sort;
    private final boolean                            _calc_similarity_score;
    private final boolean                            _sort_by_species_count_first;
    private final boolean                            _treat_as_binary_comparison;

    public BasicDomainSimilarityCalculator( final DomainSimilarity.DomainSimilaritySortField sort,
                                            final boolean sort_by_species_count_first,
                                            final boolean treat_as_binary_comparison,
                                            final boolean calc_similarity_score ) {
        _sort = sort;
        _sort_by_species_count_first = sort_by_species_count_first;
        _treat_as_binary_comparison = treat_as_binary_comparison;
        _calc_similarity_score = calc_similarity_score;
    }

    @Override
    public SortedSet<DomainSimilarity> calculateSimilarities( final PairwiseDomainSimilarityCalculator pairwise_calculator,
                                                              final List<GenomeWideCombinableDomains> cdc_list,
                                                              final boolean ignore_domains_without_combinations_in_any_genome,
                                                              final boolean ignore_domains_specific_to_one_genome ) {
        if ( cdc_list.size() < 2 ) {
            throw new IllegalArgumentException( "attempt to calculate multiple combinable domains similarity for less than two combinale domains collections" );
        }
        final SortedSet<DomainSimilarity> similarities = new TreeSet<DomainSimilarity>();
        final SortedSet<String> keys = new TreeSet<String>();
        for( final GenomeWideCombinableDomains cdc : cdc_list ) {
            keys.addAll( ( cdc ).getAllCombinableDomainsIds().keySet() );
        }
        final DecimalFormat pf = new java.text.DecimalFormat( "000000" );
        int counter = 1;
        System.out.println( keys.size() );
        for( final String key : keys ) {
            ForesterUtil.updateProgress( counter, pf );
            counter++;
            final List<CombinableDomains> same_id_cd_list = new ArrayList<CombinableDomains>( cdc_list.size() );
            final List<Species> species_with_key_id_domain = new ArrayList<Species>();
            for( final GenomeWideCombinableDomains cdc : cdc_list ) {
                if ( cdc.contains( key ) ) {
                    same_id_cd_list.add( cdc.get( key ) );
                    species_with_key_id_domain.add( cdc.getSpecies() );
                }
            }
            if ( ignore_domains_without_combinations_in_any_genome ) { //TODO: test me..........................................<<<<<<<<<<<<<
                boolean without_combinations = true;
                for( final CombinableDomains cd : same_id_cd_list ) {
                    if ( cd.getNumberOfCombinableDomains() > 0 ) {
                        without_combinations = false;
                        break;
                    }
                }
                if ( without_combinations ) {
                    continue;
                }
            }
            if ( same_id_cd_list.size() > 0 ) {
                if ( !ignore_domains_specific_to_one_genome || ( same_id_cd_list.size() > 1 ) ) {
                    final DomainSimilarity s = calculateSimilarity( pairwise_calculator, same_id_cd_list );
                    if ( s != null ) {
                        similarities.add( s );
                    }
                    else {
                        throw new RuntimeException( "similarity is null: this should not have happened" );
                    }
                }
            }
            else {
                throw new RuntimeException( "this should not have happened" );
            }
        }
        System.out.println();
        return similarities;
    }

    public boolean isCalcSimilarityScore() {
        return _calc_similarity_score;
    }

    private DomainSimilarity calculateSimilarity( final PairwiseDomainSimilarityCalculator pairwise_calculator,
                                                  final List<CombinableDomains> domains_list ) {
        if ( domains_list.size() == 1 ) {
            final SortedMap<Species, SpeciesSpecificDcData> species_data = new TreeMap<Species, SpeciesSpecificDcData>();
            species_data.put( domains_list.get( 0 ).getSpecies(),
                              createSpeciesSpecificDomainSimilariyData( domains_list.get( 0 ) ) );
            if ( !isCalcSimilarityScore() ) {
                return new DomainSimilarity( domains_list.get( 0 ),
                                             0,
                                             0,
                                             species_data,
                                             isSortBySpeciesCountFirst(),
                                             isTreatAsBinaryComparison() );
            }
            else {
                return new DomainSimilarity( domains_list.get( 0 ),
                                             1.0,
                                             1.0,
                                             1.0,
                                             1.0,
                                             0.0,
                                             0,
                                             0,
                                             0,
                                             species_data,
                                             isSortBySpeciesCountFirst(),
                                             isTreatAsBinaryComparison() );
            }
        }
        DescriptiveStatistics stat = null;
        if ( isCalcSimilarityScore() ) {
            stat = new BasicDescriptiveStatistics();
        }
        final SortedMap<Species, SpeciesSpecificDcData> species_data = new TreeMap<Species, SpeciesSpecificDcData>();
        species_data.put( domains_list.get( 0 ).getSpecies(),
                          createSpeciesSpecificDomainSimilariyData( domains_list.get( 0 ) ) );
        int max_difference_in_counts = 0;
        int max_difference = 0;
        final boolean is_domain_combination_based = pairwise_calculator instanceof CombinationsBasedPairwiseDomainSimilarityCalculator;
        for( int i = 1; i < domains_list.size(); ++i ) {
            species_data.put( domains_list.get( i ).getSpecies(),
                              createSpeciesSpecificDomainSimilariyData( domains_list.get( i ) ) );
            final CombinableDomains domains_i = domains_list.get( i );
            for( int j = 0; j < i; ++j ) {
                final PairwiseDomainSimilarity pairwise_similarity = pairwise_calculator
                        .calculateSimilarity( domains_i, domains_list.get( j ) );
                final int difference_in_counts = pairwise_similarity.getDifferenceInCounts();
                int difference = 0;
                if ( is_domain_combination_based ) {
                    difference = ( ( CombinationsBasedPairwiseDomainSimilarity ) pairwise_similarity )
                            .getNumberOfDifferentDomains();
                }
                else {
                    difference = difference_in_counts;
                }
                if ( Math.abs( difference_in_counts ) > Math.abs( max_difference_in_counts ) ) {
                    max_difference_in_counts = difference_in_counts;
                }
                if ( Math.abs( difference ) > Math.abs( max_difference ) ) {
                    max_difference = difference;
                }
                if ( isCalcSimilarityScore() ) {
                    stat.addValue( pairwise_similarity.getSimilarityScore() );
                }
            }
        }
        if ( isCalcSimilarityScore() ) {
            if ( stat.getN() < 1 ) {
                throw new RuntimeException( "empty descriptive statistics: this should not have happened" );
            }
            if ( ( stat.getN() != 1 ) && isTreatAsBinaryComparison() ) {
                throw new IllegalArgumentException( "attmpt to treat similarity with N not equal to one as binary comparison" );
            }
        }
        if ( !isTreatAsBinaryComparison() && ( max_difference_in_counts < 0 ) ) {
            max_difference_in_counts = Math.abs( max_difference_in_counts );
            if ( !is_domain_combination_based ) {
                max_difference = Math.abs( max_difference );
            }
        }
        DomainSimilarity similarity = null;
        if ( !isCalcSimilarityScore() ) {
            similarity = new DomainSimilarity( domains_list.get( 0 ),
                                               max_difference_in_counts,
                                               max_difference,
                                               species_data,
                                               isSortBySpeciesCountFirst(),
                                               isTreatAsBinaryComparison() );
        }
        else {
            if ( stat.getN() == 1 ) {
                similarity = new DomainSimilarity( domains_list.get( 0 ),
                                                   stat.getMin(),
                                                   stat.getMax(),
                                                   stat.arithmeticMean(),
                                                   stat.median(),
                                                   0.0,
                                                   stat.getN(),
                                                   max_difference_in_counts,
                                                   max_difference,
                                                   species_data,
                                                   isSortBySpeciesCountFirst(),
                                                   isTreatAsBinaryComparison() );
            }
            else {
                similarity = new DomainSimilarity( domains_list.get( 0 ),
                                                   stat.getMin(),
                                                   stat.getMax(),
                                                   stat.arithmeticMean(),
                                                   stat.median(),
                                                   stat.sampleStandardDeviation(),
                                                   stat.getN(),
                                                   max_difference_in_counts,
                                                   max_difference,
                                                   species_data,
                                                   isSortBySpeciesCountFirst(),
                                                   isTreatAsBinaryComparison() );
            }
        }
        return similarity;
    }

    private boolean isSortBySpeciesCountFirst() {
        return _sort_by_species_count_first;
    }

    private boolean isTreatAsBinaryComparison() {
        return _treat_as_binary_comparison;
    }

    private static SpeciesSpecificDcData createSpeciesSpecificDomainSimilariyData( final CombinableDomains cd ) {
        final SpeciesSpecificDcData sd = new PrintableSpeciesSpecificDcData( cd.getKeyDomainCount(),
                                                                             cd.getNumberOfCombinableDomains() );
        for( final String prot : cd.getKeyDomainProteins() ) {
            sd.addKeyDomainProtein( prot );
        }
        for( final String domain : cd.getCombinableDomains() ) {
            sd.addProteinsExhibitingCombinationCount( domain, cd.getNumberOfProteinsExhibitingCombination( domain ) );
        }
        return sd;
    }
}

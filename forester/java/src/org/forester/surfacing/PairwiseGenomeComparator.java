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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.evoinference.matrix.distance.DistanceMatrix;
import org.forester.go.GoId;
import org.forester.go.GoNameSpace;
import org.forester.go.GoTerm;
import org.forester.species.Species;
import org.forester.surfacing.DomainSimilarityCalculator.Detailedness;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public class PairwiseGenomeComparator {

    private List<DistanceMatrix> _domain_distance_scores_means;
    private List<DistanceMatrix> _shared_domains_based_distances;
    private List<DistanceMatrix> _shared_binary_combinations_based_distances;

    //private List<HistogramData>  _histogram_datas;
    public PairwiseGenomeComparator() {
        init();
    }

    public List<DistanceMatrix> getDomainDistanceScoresMeans() {
        return _domain_distance_scores_means;
    }

    //public List<HistogramData> getHistogramDatas() {
    //    return _histogram_datas;
    //}
    public List<DistanceMatrix> getSharedBinaryCombinationsBasedDistances() {
        return _shared_binary_combinations_based_distances;
    }

    public List<DistanceMatrix> getSharedDomainsBasedDistances() {
        return _shared_domains_based_distances;
    }

    private void init() {
        //_histogram_datas = new ArrayList<HistogramData>();
        _domain_distance_scores_means = new ArrayList<DistanceMatrix>();
        _shared_domains_based_distances = new ArrayList<DistanceMatrix>();
        _shared_binary_combinations_based_distances = new ArrayList<DistanceMatrix>();
    }

    public void performPairwiseComparisons( final StringBuilder html_desc,
                                            final boolean sort_by_species_count_first,
                                            final Detailedness detailedness,
                                            final boolean ignore_domains_without_combs_in_all_spec,
                                            final boolean ignore_domains_specific_to_one_species,
                                            final DomainSimilarity.DomainSimilaritySortField domain_similarity_sort_field,
                                            final PrintableDomainSimilarity.PRINT_OPTION domain_similarity_print_option,
                                            final DomainSimilarity.DomainSimilarityScoring scoring,
                                            final Map<String, List<GoId>> domain_id_to_go_ids_map,
                                            final Map<GoId, GoTerm> go_id_to_term_map,
                                            final GoNameSpace go_namespace_limit,
                                            final Species[] species,
                                            final int number_of_genomes,
                                            final List<GenomeWideCombinableDomains> list_of_genome_wide_combinable_domains,
                                            final PairwiseDomainSimilarityCalculator pw_calc,
                                            final String automated_pairwise_comparison_suffix,
                                            final boolean verbose,
                                            final String automated_pairwise_comparison_prefix,
                                            final String command_line_prg_name,
                                            final File out_dir,
                                            final boolean write_pairwise_comparisons,
                                            final Map<String, Integer> tax_code_to_id_map ) {
        init();
        final BasicSymmetricalDistanceMatrix domain_distance_scores_means = new BasicSymmetricalDistanceMatrix( number_of_genomes );
        final BasicSymmetricalDistanceMatrix shared_domains_based_distances = new BasicSymmetricalDistanceMatrix( number_of_genomes );
        final BasicSymmetricalDistanceMatrix shared_binary_combinations_based_distances = new BasicSymmetricalDistanceMatrix( number_of_genomes );
        if ( verbose ) {
            System.out.println();
            System.out.println( "Pairwise genome distances:" );
            System.out.print( "[species-i - species-j:" );
            System.out.print( " mean-score-based" );
            System.out.print( " (sd)" );
            System.out.print( " [N]" );
            System.out.print( " | shared-domains-based" );
            System.out.println( " | shared-binary-combinations-based]" );
            System.out.println();
        }
        for( int i = 0; i < number_of_genomes; ++i ) {
            final String species_i = species[ i ].getSpeciesId();
            domain_distance_scores_means.setIdentifier( i, species_i );
            shared_domains_based_distances.setIdentifier( i, species_i );
            shared_binary_combinations_based_distances.setIdentifier( i, species_i );
            if ( verbose ) {
                System.out.println( ( i + 1 ) + "/" + number_of_genomes );
            }
            for( int j = 0; j < i; ++j ) {
                if ( ( list_of_genome_wide_combinable_domains.get( i ).getSize() < 1 )
                        || ( list_of_genome_wide_combinable_domains.get( j ).getSize() < 1 ) ) {
                    domain_distance_scores_means
                            .setValue( i, j, DomainArchitectureBasedGenomeSimilarityCalculator.MAX_SIMILARITY_SCORE );
                    shared_domains_based_distances
                            .setValue( i, j, DomainArchitectureBasedGenomeSimilarityCalculator.MAX_SIMILARITY_SCORE );
                    shared_binary_combinations_based_distances
                            .setValue( i, j, DomainArchitectureBasedGenomeSimilarityCalculator.MAX_SIMILARITY_SCORE );
                    continue;
                }
                final List<GenomeWideCombinableDomains> genome_pair = new ArrayList<GenomeWideCombinableDomains>( 2 );
                genome_pair.add( list_of_genome_wide_combinable_domains.get( i ) );
                genome_pair.add( list_of_genome_wide_combinable_domains.get( j ) );
                DomainSimilarityCalculator.GoAnnotationOutput go_annotation_output = DomainSimilarityCalculator.GoAnnotationOutput.NONE;
                if ( domain_id_to_go_ids_map != null ) {
                    go_annotation_output = DomainSimilarityCalculator.GoAnnotationOutput.ALL;
                }
                final DomainSimilarityCalculator calc = new BasicDomainSimilarityCalculator( domain_similarity_sort_field,
                                                                                             sort_by_species_count_first,
                                                                                             true );
                final SortedSet<DomainSimilarity> similarities = calc
                        .calculateSimilarities( pw_calc,
                                                genome_pair,
                                                ignore_domains_without_combs_in_all_spec,
                                                ignore_domains_specific_to_one_species );
                SurfacingUtil.decoratePrintableDomainSimilarities( similarities, detailedness );
                final DescriptiveStatistics stats = SurfacingUtil
                        .calculateDescriptiveStatisticsForMeanValues( similarities );
                final String species_j = species[ j ].getSpeciesId();
                final DomainArchitectureBasedGenomeSimilarityCalculator genome_similarity_calculator = new DomainArchitectureBasedGenomeSimilarityCalculator( list_of_genome_wide_combinable_domains
                                                                                                                                                                      .get( i ),
                                                                                                                                                              list_of_genome_wide_combinable_domains
                                                                                                                                                                      .get( j ) );
                genome_similarity_calculator.setAllowDomainsToBeIgnored( false );
                double dissimilarity_score_mean;
                if ( stats.getN() < 1 ) {
                    // No domains in common
                    dissimilarity_score_mean = 1.0;
                }
                else {
                    dissimilarity_score_mean = 1.0 - stats.arithmeticMean();
                }
                final double shared_domains_based_genome_distance = 1.0 - genome_similarity_calculator
                        .calculateSharedDomainsBasedGenomeSimilarityScore();
                final double shared_binary_combinations_based_genome_distance = 1.0 - genome_similarity_calculator
                        .calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore();
                domain_distance_scores_means.setValue( i, j, dissimilarity_score_mean );
                shared_domains_based_distances.setValue( i, j, shared_domains_based_genome_distance );
                shared_binary_combinations_based_distances.setValue( i,
                                                                     j,
                                                                     shared_binary_combinations_based_genome_distance );
                if ( verbose ) {
                    System.out.print( species_i + "-" );
                    System.out.print( species_j + ": " );
                    System.out.print( ForesterUtil.round( dissimilarity_score_mean, 2 ) );
                    if ( stats.getN() > 1 ) {
                        System.out.print( " (" + ForesterUtil.round( stats.sampleStandardDeviation(), 2 ) + ")" );
                    }
                    else {
                        System.out.print( " (n/a)" );
                    }
                    System.out.print( " [" + stats.getN() + "]" );
                    System.out.print( " | " );
                    System.out.print( ForesterUtil.round( shared_domains_based_genome_distance, 2 ) );
                    System.out.print( " | " );
                    System.out.println( ForesterUtil.round( shared_binary_combinations_based_genome_distance, 2 ) );
                }
                String pairwise_similarities_output_file_str = automated_pairwise_comparison_prefix + species_i + "_"
                        + species_j + automated_pairwise_comparison_suffix;
                switch ( domain_similarity_print_option ) {
                    case HTML:
                        if ( !pairwise_similarities_output_file_str.endsWith( ".html" ) ) {
                            pairwise_similarities_output_file_str += ".html";
                        }
                        break;
                }
                DescriptiveStatistics pw_stats = null;
                if ( write_pairwise_comparisons ) {
                    try {
                        final Writer writer = new BufferedWriter( new FileWriter( out_dir == null ? pairwise_similarities_output_file_str
                                : out_dir + ForesterUtil.FILE_SEPARATOR + pairwise_similarities_output_file_str ) );
                        pw_stats = SurfacingUtil.writeDomainSimilaritiesToFile( html_desc,
                                                                                new StringBuilder( species_i + "-"
                                                                                        + species_j ),
                                                                                writer,
                                                                                null,
                                                                                similarities,
                                                                                true,
                                                                                null,
                                                                                domain_similarity_print_option,
                                                                                domain_similarity_sort_field,
                                                                                scoring,
                                                                                false,
                                                                                tax_code_to_id_map );
                    }
                    catch ( final IOException e ) {
                        ForesterUtil.fatalError( command_line_prg_name, "Failed to write similarites to: \""
                                + pairwise_similarities_output_file_str + "\" [" + e.getMessage() + "]" );
                    }
                }
                if ( pw_stats != null ) {
                    if ( pw_stats.getMin() >= pw_stats.getMax() ) {
                        ForesterUtil
                                .printWarningMessage( command_line_prg_name, "for [" + species_i + "-" + species_j
                                        + "] score minimum is [" + pw_stats.getMin() + "] while score maximum is ["
                                        + pw_stats.getMax()
                                        + "], possibly indicating that a genome is compared to itself" );
                    }
                }
            }
        }
        getDomainDistanceScoresMeans().add( domain_distance_scores_means );
        getSharedDomainsBasedDistances().add( shared_domains_based_distances );
        getSharedBinaryCombinationsBasedDistances().add( shared_binary_combinations_based_distances );
        if ( verbose ) {
            System.out.println();
        }
    }

    public void performPairwiseComparisonsJacknifed( final Species[] species,
                                                     final int number_of_genomes,
                                                     final List<GenomeWideCombinableDomains> list_of_genome_wide_combinable_domains,
                                                     final boolean verbose,
                                                     final int number_of_resamplings,
                                                     final double jacknife_ratio,
                                                     final long random_seed ) {
        init();
        if ( number_of_resamplings < 2 ) {
            throw new IllegalArgumentException( "attempt to perform jacknife resampling with less than 2 resamplings" );
        }
        if ( jacknife_ratio <= 0.0 ) {
            throw new IllegalArgumentException( "attempt to perform jacknife resampling with jacknife ratio of 0.0 or less" );
        }
        else if ( jacknife_ratio >= 1.0 ) {
            throw new IllegalArgumentException( "attempt to perform jacknife resampling with jacknife ratio 1.0 or more" );
        }
        final String[] all_unique_domain_ids = getAllUniqueDomainIdAsArray( list_of_genome_wide_combinable_domains );
        if ( verbose ) {
            System.out.println();
            System.out.println( "Jacknife: total of domains: " + all_unique_domain_ids.length );
        }
        if ( verbose ) {
            System.out.print( "resampling " );
        }
        final Random generator = new Random( random_seed );
        for( int r = 0; r < number_of_resamplings; ++r ) {
            if ( verbose ) {
                System.out.print( " " + r );
            }
            final SortedSet<String> domain_ids_to_ignore = randomlyPickDomainIds( all_unique_domain_ids,
                                                                                  jacknife_ratio,
                                                                                  generator );
            final BasicSymmetricalDistanceMatrix shared_domains_based_distances = new BasicSymmetricalDistanceMatrix( number_of_genomes );
            final BasicSymmetricalDistanceMatrix shared_binary_combinations_based_distances = new BasicSymmetricalDistanceMatrix( number_of_genomes );
            for( int i = 0; i < number_of_genomes; ++i ) {
                final String species_i = species[ i ].getSpeciesId();
                shared_domains_based_distances.setIdentifier( i, species_i );
                shared_binary_combinations_based_distances.setIdentifier( i, species_i );
                for( int j = 0; j < i; ++j ) {
                    final List<GenomeWideCombinableDomains> genome_pair = new ArrayList<GenomeWideCombinableDomains>( 2 );
                    genome_pair.add( list_of_genome_wide_combinable_domains.get( i ) );
                    genome_pair.add( list_of_genome_wide_combinable_domains.get( j ) );
                    final DomainArchitectureBasedGenomeSimilarityCalculator genome_simiarity_calculator = new DomainArchitectureBasedGenomeSimilarityCalculator( list_of_genome_wide_combinable_domains
                                                                                                                                                                         .get( i ),
                                                                                                                                                                 list_of_genome_wide_combinable_domains
                                                                                                                                                                         .get( j ) );
                    genome_simiarity_calculator.setAllowDomainsToBeIgnored( true );
                    genome_simiarity_calculator.setDomainIdsToIgnore( domain_ids_to_ignore );
                    shared_domains_based_distances.setValue( i, j, 1.0 - genome_simiarity_calculator
                            .calculateSharedDomainsBasedGenomeSimilarityScore() );
                    shared_binary_combinations_based_distances.setValue( i, j, 1.0 - genome_simiarity_calculator
                            .calculateSharedBinaryDomainCombinationBasedGenomeSimilarityScore() );
                }
            }
            getSharedDomainsBasedDistances().add( shared_domains_based_distances );
            getSharedBinaryCombinationsBasedDistances().add( shared_binary_combinations_based_distances );
        }
        if ( verbose ) {
            System.out.println();
        }
    }

    static private String[] getAllUniqueDomainIdAsArray( final List<GenomeWideCombinableDomains> list_of_genome_wide_combinable_domains ) {
        String[] all_domain_ids_array;
        final SortedSet<String> all_domain_ids = new TreeSet<String>();
        for( final GenomeWideCombinableDomains genome_wide_combinable_domains : list_of_genome_wide_combinable_domains ) {
            final SortedSet<String> all_domains = genome_wide_combinable_domains.getAllDomainIds();
            for( final String domain : all_domains ) {
                all_domain_ids.add( domain );
            }
        }
        all_domain_ids_array = new String[ all_domain_ids.size() ];
        int n = 0;
        for( final String domain_id : all_domain_ids ) {
            all_domain_ids_array[ n++ ] = domain_id;
        }
        return all_domain_ids_array;
    }

    static private SortedSet<String> randomlyPickDomainIds( final String[] all_domain_ids_array,
                                                            final double jacknife_ratio,
                                                            final Random generator ) {
        final int size = all_domain_ids_array.length;
        final SortedSet<String> random_domain_ids = new TreeSet<String>();
        final int number_of_ids_pick = ForesterUtil.roundToInt( jacknife_ratio * size );
        while ( random_domain_ids.size() < number_of_ids_pick ) {
            final int r = generator.nextInt( size );
            random_domain_ids.add( all_domain_ids_array[ r ] );
        }
        return random_domain_ids;
    }
}

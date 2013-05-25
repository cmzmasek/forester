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
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.application.surfacing;
import org.forester.evoinference.distance.NeighborJoining;
import org.forester.evoinference.matrix.character.BasicCharacterStateMatrix;
import org.forester.evoinference.matrix.character.CharacterStateMatrix;
import org.forester.evoinference.matrix.character.CharacterStateMatrix.BinaryStates;
import org.forester.evoinference.matrix.character.CharacterStateMatrix.Format;
import org.forester.evoinference.matrix.character.CharacterStateMatrix.GainLossStates;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.evoinference.matrix.distance.DistanceMatrix;
import org.forester.go.GoId;
import org.forester.go.GoNameSpace;
import org.forester.go.GoTerm;
import org.forester.go.PfamToGoMapping;
import org.forester.io.parsers.nexus.NexusConstants;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.PhylogenyNode.NH_CONVERSION_SUPPORT_VALUE_STYLE;
import org.forester.phylogeny.data.BinaryCharacters;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.protein.BasicDomain;
import org.forester.protein.BasicProtein;
import org.forester.protein.BinaryDomainCombination;
import org.forester.protein.Domain;
import org.forester.protein.DomainId;
import org.forester.protein.Protein;
import org.forester.species.Species;
import org.forester.surfacing.DomainSimilarityCalculator.Detailedness;
import org.forester.surfacing.DomainSimilarityCalculator.GoAnnotationOutput;
import org.forester.surfacing.GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder;
import org.forester.util.AsciiHistogram;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.DescriptiveStatistics;
import org.forester.util.ForesterUtil;

public final class SurfacingUtil {

    private final static NumberFormat       FORMATTER                        = new DecimalFormat( "0.0E0" );
    private final static NumberFormat       FORMATTER_3                      = new DecimalFormat( "0.000" );
    private static final Comparator<Domain> ASCENDING_CONFIDENCE_VALUE_ORDER = new Comparator<Domain>() {

                                                                                 @Override
                                                                                 public int compare( final Domain d1,
                                                                                                     final Domain d2 ) {
                                                                                     if ( d1.getPerSequenceEvalue() < d2
                                                                                             .getPerSequenceEvalue() ) {
                                                                                         return -1;
                                                                                     }
                                                                                     else if ( d1
                                                                                             .getPerSequenceEvalue() > d2
                                                                                             .getPerSequenceEvalue() ) {
                                                                                         return 1;
                                                                                     }
                                                                                     else {
                                                                                         return d1.compareTo( d2 );
                                                                                     }
                                                                                 }
                                                                             };
    public final static Pattern             PATTERN_SP_STYLE_TAXONOMY        = Pattern.compile( "^[A-Z0-9]{3,5}$" );
    private static final boolean            USE_LAST                         = true;

    private SurfacingUtil() {
        // Hidden constructor.
    }

    public static void performDomainArchitectureAnalysis( final SortedMap<String, Set<String>> domain_architecutures,
                                                          final SortedMap<String, Integer> domain_architecuture_counts,
                                                          final int min_count,
                                                          final File da_counts_outfile,
                                                          final File unique_da_outfile ) {
        checkForOutputFileWriteability( da_counts_outfile );
        checkForOutputFileWriteability( unique_da_outfile );
        try {
            final BufferedWriter da_counts_out = new BufferedWriter( new FileWriter( da_counts_outfile ) );
            final BufferedWriter unique_da_out = new BufferedWriter( new FileWriter( unique_da_outfile ) );
            final Iterator<Entry<String, Integer>> it = domain_architecuture_counts.entrySet().iterator();
            while ( it.hasNext() ) {
                final Map.Entry<String, Integer> e = it.next();
                final String da = e.getKey();
                final int count = e.getValue();
                if ( count >= min_count ) {
                    da_counts_out.write( da );
                    da_counts_out.write( "\t" );
                    da_counts_out.write( String.valueOf( count ) );
                    da_counts_out.write( ForesterUtil.LINE_SEPARATOR );
                }
                if ( count == 1 ) {
                    final Iterator<Entry<String, Set<String>>> it2 = domain_architecutures.entrySet().iterator();
                    while ( it2.hasNext() ) {
                        final Map.Entry<String, Set<String>> e2 = it2.next();
                        final String genome = e2.getKey();
                        final Set<String> das = e2.getValue();
                        if ( das.contains( da ) ) {
                            unique_da_out.write( genome );
                            unique_da_out.write( "\t" );
                            unique_da_out.write( da );
                            unique_da_out.write( ForesterUtil.LINE_SEPARATOR );
                        }
                    }
                }
            }
            unique_da_out.close();
            da_counts_out.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote distance matrices to \"" + da_counts_outfile + "\"" );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote distance matrices to \"" + unique_da_outfile + "\"" );
        //
    }

    public static int storeDomainArchitectures( final String genome,
                                                final SortedMap<String, Set<String>> domain_architecutures,
                                                final List<Protein> protein_list,
                                                final Map<String, Integer> distinct_domain_architecuture_counts ) {
        final Set<String> da = new HashSet<String>();
        domain_architecutures.put( genome, da );
        for( final Protein protein : protein_list ) {
            final String da_str = ( ( BasicProtein ) protein ).toDomainArchitectureString( "~" );
            if ( !da.contains( da_str ) ) {
                if ( !distinct_domain_architecuture_counts.containsKey( da_str ) ) {
                    distinct_domain_architecuture_counts.put( da_str, 1 );
                }
                else {
                    distinct_domain_architecuture_counts.put( da_str,
                                                              distinct_domain_architecuture_counts.get( da_str ) + 1 );
                }
                da.add( da_str );
            }
        }
        return da.size();
    }

    public static void addAllBinaryDomainCombinationToSet( final GenomeWideCombinableDomains genome,
                                                           final SortedSet<BinaryDomainCombination> binary_domain_combinations ) {
        final SortedMap<DomainId, CombinableDomains> all_cd = genome.getAllCombinableDomainsIds();
        for( final DomainId domain_id : all_cd.keySet() ) {
            binary_domain_combinations.addAll( all_cd.get( domain_id ).toBinaryDomainCombinations() );
        }
    }

    public static void addAllDomainIdsToSet( final GenomeWideCombinableDomains genome,
                                             final SortedSet<DomainId> domain_ids ) {
        final SortedSet<DomainId> domains = genome.getAllDomainIds();
        for( final DomainId domain : domains ) {
            domain_ids.add( domain );
        }
    }

    public static void addHtmlHead( final Writer w, final String title ) throws IOException {
        w.write( SurfacingConstants.NL );
        w.write( "<head>" );
        w.write( "<title>" );
        w.write( title );
        w.write( "</title>" );
        w.write( SurfacingConstants.NL );
        w.write( "<style>" );
        w.write( SurfacingConstants.NL );
        w.write( "a:visited { color : #6633FF; text-decoration : none; }" );
        w.write( SurfacingConstants.NL );
        w.write( "a:link { color : #6633FF; text-decoration : none; }" );
        w.write( SurfacingConstants.NL );
        w.write( "a:active { color : #99FF00; text-decoration : none; }" );
        w.write( SurfacingConstants.NL );
        w.write( "a:hover { color : #FFFFFF; background-color : #99FF00; text-decoration : none; }" );
        w.write( SurfacingConstants.NL );
        w.write( "td { text-align: left; vertical-align: top; font-family: Verdana, Arial, Helvetica; font-size: 8pt}" );
        w.write( SurfacingConstants.NL );
        w.write( "h1 { color : #0000FF; font-family: Verdana, Arial, Helvetica; font-size: 18pt; font-weight: bold }" );
        w.write( SurfacingConstants.NL );
        w.write( "h2 { color : #0000FF; font-family: Verdana, Arial, Helvetica; font-size: 16pt; font-weight: bold }" );
        w.write( SurfacingConstants.NL );
        w.write( "</style>" );
        w.write( SurfacingConstants.NL );
        w.write( "</head>" );
        w.write( SurfacingConstants.NL );
    }

    public static DescriptiveStatistics calculateDescriptiveStatisticsForMeanValues( final Set<DomainSimilarity> similarities ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final DomainSimilarity similarity : similarities ) {
            stats.addValue( similarity.getMeanSimilarityScore() );
        }
        return stats;
    }

    private static void calculateIndependentDomainCombinationGains( final Phylogeny local_phylogeny_l,
                                                                    final String outfilename_for_counts,
                                                                    final String outfilename_for_dc,
                                                                    final String outfilename_for_dc_for_go_mapping,
                                                                    final String outfilename_for_dc_for_go_mapping_unique,
                                                                    final String outfilename_for_rank_counts,
                                                                    final String outfilename_for_ancestor_species_counts,
                                                                    final String outfilename_for_protein_stats,
                                                                    final Map<String, DescriptiveStatistics> protein_length_stats_by_dc,
                                                                    final Map<String, DescriptiveStatistics> domain_number_stats_by_dc,
                                                                    final Map<String, DescriptiveStatistics> domain_length_stats_by_domain ) {
        try {
            //
            //            if ( protein_length_stats_by_dc != null ) {
            //                for( final Entry<?, DescriptiveStatistics> entry : protein_length_stats_by_dc.entrySet() ) {
            //                    System.out.print( entry.getKey().toString() );
            //                    System.out.print( ": " );
            //                    double[] a = entry.getValue().getDataAsDoubleArray();
            //                    for( int i = 0; i < a.length; i++ ) {
            //                        System.out.print( a[ i ] + " " );
            //                    }
            //                    System.out.println();
            //                }
            //            }
            //            if ( domain_number_stats_by_dc != null ) {
            //                for( final Entry<?, DescriptiveStatistics> entry : domain_number_stats_by_dc.entrySet() ) {
            //                    System.out.print( entry.getKey().toString() );
            //                    System.out.print( ": " );
            //                    double[] a = entry.getValue().getDataAsDoubleArray();
            //                    for( int i = 0; i < a.length; i++ ) {
            //                        System.out.print( a[ i ] + " " );
            //                    }
            //                    System.out.println();
            //                }
            //            }
            //
            final BufferedWriter out_counts = new BufferedWriter( new FileWriter( outfilename_for_counts ) );
            final BufferedWriter out_dc = new BufferedWriter( new FileWriter( outfilename_for_dc ) );
            final BufferedWriter out_dc_for_go_mapping = new BufferedWriter( new FileWriter( outfilename_for_dc_for_go_mapping ) );
            final BufferedWriter out_dc_for_go_mapping_unique = new BufferedWriter( new FileWriter( outfilename_for_dc_for_go_mapping_unique ) );
            final SortedMap<String, Integer> dc_gain_counts = new TreeMap<String, Integer>();
            for( final PhylogenyNodeIterator it = local_phylogeny_l.iteratorPostorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                final Set<String> gained_dc = n.getNodeData().getBinaryCharacters().getGainedCharacters();
                for( final String dc : gained_dc ) {
                    if ( dc_gain_counts.containsKey( dc ) ) {
                        dc_gain_counts.put( dc, dc_gain_counts.get( dc ) + 1 );
                    }
                    else {
                        dc_gain_counts.put( dc, 1 );
                    }
                }
            }
            final SortedMap<Integer, Integer> histogram = new TreeMap<Integer, Integer>();
            final SortedMap<Integer, StringBuilder> domain_lists = new TreeMap<Integer, StringBuilder>();
            final SortedMap<Integer, DescriptiveStatistics> dc_reapp_counts_to_protein_length_stats = new TreeMap<Integer, DescriptiveStatistics>();
            final SortedMap<Integer, DescriptiveStatistics> dc_reapp_counts_to_domain_number_stats = new TreeMap<Integer, DescriptiveStatistics>();
            final SortedMap<Integer, DescriptiveStatistics> dc_reapp_counts_to_domain_lengths_stats = new TreeMap<Integer, DescriptiveStatistics>();
            final SortedMap<Integer, PriorityQueue<String>> domain_lists_go = new TreeMap<Integer, PriorityQueue<String>>();
            final SortedMap<Integer, SortedSet<String>> domain_lists_go_unique = new TreeMap<Integer, SortedSet<String>>();
            final Set<String> dcs = dc_gain_counts.keySet();
            final SortedSet<String> more_than_once = new TreeSet<String>();
            final DescriptiveStatistics gained_once_lengths_stats = new BasicDescriptiveStatistics();
            final DescriptiveStatistics gained_once_domain_count_stats = new BasicDescriptiveStatistics();
            final DescriptiveStatistics gained_multiple_times_lengths_stats = new BasicDescriptiveStatistics();
            final DescriptiveStatistics gained_multiple_times_domain_count_stats = new BasicDescriptiveStatistics();
            long gained_multiple_times_domain_length_sum = 0;
            long gained_once_domain_length_sum = 0;
            long gained_multiple_times_domain_length_count = 0;
            long gained_once_domain_length_count = 0;
            for( final String dc : dcs ) {
                final int count = dc_gain_counts.get( dc );
                if ( histogram.containsKey( count ) ) {
                    histogram.put( count, histogram.get( count ) + 1 );
                    domain_lists.get( count ).append( ", " + dc );
                    domain_lists_go.get( count ).addAll( splitDomainCombination( dc ) );
                    domain_lists_go_unique.get( count ).addAll( splitDomainCombination( dc ) );
                }
                else {
                    histogram.put( count, 1 );
                    domain_lists.put( count, new StringBuilder( dc ) );
                    final PriorityQueue<String> q = new PriorityQueue<String>();
                    q.addAll( splitDomainCombination( dc ) );
                    domain_lists_go.put( count, q );
                    final SortedSet<String> set = new TreeSet<String>();
                    set.addAll( splitDomainCombination( dc ) );
                    domain_lists_go_unique.put( count, set );
                }
                if ( protein_length_stats_by_dc != null ) {
                    if ( !dc_reapp_counts_to_protein_length_stats.containsKey( count ) ) {
                        dc_reapp_counts_to_protein_length_stats.put( count, new BasicDescriptiveStatistics() );
                    }
                    dc_reapp_counts_to_protein_length_stats.get( count ).addValue( protein_length_stats_by_dc.get( dc )
                            .arithmeticMean() );
                }
                if ( domain_number_stats_by_dc != null ) {
                    if ( !dc_reapp_counts_to_domain_number_stats.containsKey( count ) ) {
                        dc_reapp_counts_to_domain_number_stats.put( count, new BasicDescriptiveStatistics() );
                    }
                    dc_reapp_counts_to_domain_number_stats.get( count ).addValue( domain_number_stats_by_dc.get( dc )
                            .arithmeticMean() );
                }
                if ( domain_length_stats_by_domain != null ) {
                    if ( !dc_reapp_counts_to_domain_lengths_stats.containsKey( count ) ) {
                        dc_reapp_counts_to_domain_lengths_stats.put( count, new BasicDescriptiveStatistics() );
                    }
                    final String[] ds = dc.split( "=" );
                    dc_reapp_counts_to_domain_lengths_stats.get( count ).addValue( domain_length_stats_by_domain
                            .get( ds[ 0 ] ).arithmeticMean() );
                    dc_reapp_counts_to_domain_lengths_stats.get( count ).addValue( domain_length_stats_by_domain
                            .get( ds[ 1 ] ).arithmeticMean() );
                }
                if ( count > 1 ) {
                    more_than_once.add( dc );
                    if ( protein_length_stats_by_dc != null ) {
                        final DescriptiveStatistics s = protein_length_stats_by_dc.get( dc );
                        for( final double element : s.getData() ) {
                            gained_multiple_times_lengths_stats.addValue( element );
                        }
                    }
                    if ( domain_number_stats_by_dc != null ) {
                        final DescriptiveStatistics s = domain_number_stats_by_dc.get( dc );
                        for( final double element : s.getData() ) {
                            gained_multiple_times_domain_count_stats.addValue( element );
                        }
                    }
                    if ( domain_length_stats_by_domain != null ) {
                        final String[] ds = dc.split( "=" );
                        final DescriptiveStatistics s0 = domain_length_stats_by_domain.get( ds[ 0 ] );
                        final DescriptiveStatistics s1 = domain_length_stats_by_domain.get( ds[ 1 ] );
                        for( final double element : s0.getData() ) {
                            gained_multiple_times_domain_length_sum += element;
                            ++gained_multiple_times_domain_length_count;
                        }
                        for( final double element : s1.getData() ) {
                            gained_multiple_times_domain_length_sum += element;
                            ++gained_multiple_times_domain_length_count;
                        }
                    }
                }
                else {
                    if ( protein_length_stats_by_dc != null ) {
                        final DescriptiveStatistics s = protein_length_stats_by_dc.get( dc );
                        for( final double element : s.getData() ) {
                            gained_once_lengths_stats.addValue( element );
                        }
                    }
                    if ( domain_number_stats_by_dc != null ) {
                        final DescriptiveStatistics s = domain_number_stats_by_dc.get( dc );
                        for( final double element : s.getData() ) {
                            gained_once_domain_count_stats.addValue( element );
                        }
                    }
                    if ( domain_length_stats_by_domain != null ) {
                        final String[] ds = dc.split( "=" );
                        final DescriptiveStatistics s0 = domain_length_stats_by_domain.get( ds[ 0 ] );
                        final DescriptiveStatistics s1 = domain_length_stats_by_domain.get( ds[ 1 ] );
                        for( final double element : s0.getData() ) {
                            gained_once_domain_length_sum += element;
                            ++gained_once_domain_length_count;
                        }
                        for( final double element : s1.getData() ) {
                            gained_once_domain_length_sum += element;
                            ++gained_once_domain_length_count;
                        }
                    }
                }
            }
            final Set<Integer> histogram_keys = histogram.keySet();
            for( final Integer histogram_key : histogram_keys ) {
                final int count = histogram.get( histogram_key );
                final StringBuilder dc = domain_lists.get( histogram_key );
                out_counts.write( histogram_key + "\t" + count + ForesterUtil.LINE_SEPARATOR );
                out_dc.write( histogram_key + "\t" + dc + ForesterUtil.LINE_SEPARATOR );
                out_dc_for_go_mapping.write( "#" + histogram_key + ForesterUtil.LINE_SEPARATOR );
                final Object[] sorted = domain_lists_go.get( histogram_key ).toArray();
                Arrays.sort( sorted );
                for( final Object domain : sorted ) {
                    out_dc_for_go_mapping.write( domain + ForesterUtil.LINE_SEPARATOR );
                }
                out_dc_for_go_mapping_unique.write( "#" + histogram_key + ForesterUtil.LINE_SEPARATOR );
                for( final String domain : domain_lists_go_unique.get( histogram_key ) ) {
                    out_dc_for_go_mapping_unique.write( domain + ForesterUtil.LINE_SEPARATOR );
                }
            }
            out_counts.close();
            out_dc.close();
            out_dc_for_go_mapping.close();
            out_dc_for_go_mapping_unique.close();
            final SortedMap<String, Integer> lca_rank_counts = new TreeMap<String, Integer>();
            final SortedMap<String, Integer> lca_ancestor_species_counts = new TreeMap<String, Integer>();
            for( final String dc : more_than_once ) {
                final List<PhylogenyNode> nodes = new ArrayList<PhylogenyNode>();
                for( final PhylogenyNodeIterator it = local_phylogeny_l.iteratorExternalForward(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    if ( n.getNodeData().getBinaryCharacters().getGainedCharacters().contains( dc ) ) {
                        nodes.add( n );
                    }
                }
                for( int i = 0; i < ( nodes.size() - 1 ); ++i ) {
                    for( int j = i + 1; j < nodes.size(); ++j ) {
                        final PhylogenyNode lca = PhylogenyMethods.calculateLCA( nodes.get( i ), nodes.get( j ) );
                        String rank = "unknown";
                        if ( lca.getNodeData().isHasTaxonomy()
                                && !ForesterUtil.isEmpty( lca.getNodeData().getTaxonomy().getRank() ) ) {
                            rank = lca.getNodeData().getTaxonomy().getRank();
                        }
                        addToCountMap( lca_rank_counts, rank );
                        String lca_species;
                        if ( lca.getNodeData().isHasTaxonomy()
                                && !ForesterUtil.isEmpty( lca.getNodeData().getTaxonomy().getScientificName() ) ) {
                            lca_species = lca.getNodeData().getTaxonomy().getScientificName();
                        }
                        else if ( lca.getNodeData().isHasTaxonomy()
                                && !ForesterUtil.isEmpty( lca.getNodeData().getTaxonomy().getCommonName() ) ) {
                            lca_species = lca.getNodeData().getTaxonomy().getCommonName();
                        }
                        else {
                            lca_species = lca.getName();
                        }
                        addToCountMap( lca_ancestor_species_counts, lca_species );
                    }
                }
            }
            final BufferedWriter out_for_rank_counts = new BufferedWriter( new FileWriter( outfilename_for_rank_counts ) );
            final BufferedWriter out_for_ancestor_species_counts = new BufferedWriter( new FileWriter( outfilename_for_ancestor_species_counts ) );
            ForesterUtil.map2writer( out_for_rank_counts, lca_rank_counts, "\t", ForesterUtil.LINE_SEPARATOR );
            ForesterUtil.map2writer( out_for_ancestor_species_counts,
                                     lca_ancestor_species_counts,
                                     "\t",
                                     ForesterUtil.LINE_SEPARATOR );
            out_for_rank_counts.close();
            out_for_ancestor_species_counts.close();
            if ( !ForesterUtil.isEmpty( outfilename_for_protein_stats )
                    && ( ( domain_length_stats_by_domain != null ) || ( protein_length_stats_by_dc != null ) || ( domain_number_stats_by_dc != null ) ) ) {
                final BufferedWriter w = new BufferedWriter( new FileWriter( outfilename_for_protein_stats ) );
                w.write( "Domain Lengths: " );
                w.write( "\n" );
                if ( domain_length_stats_by_domain != null ) {
                    for( final Entry<Integer, DescriptiveStatistics> entry : dc_reapp_counts_to_domain_lengths_stats
                            .entrySet() ) {
                        w.write( entry.getKey().toString() );
                        w.write( "\t" + entry.getValue().arithmeticMean() );
                        w.write( "\t" + entry.getValue().median() );
                        w.write( "\n" );
                    }
                }
                w.flush();
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Protein Lengths: " );
                w.write( "\n" );
                if ( protein_length_stats_by_dc != null ) {
                    for( final Entry<Integer, DescriptiveStatistics> entry : dc_reapp_counts_to_protein_length_stats
                            .entrySet() ) {
                        w.write( entry.getKey().toString() );
                        w.write( "\t" + entry.getValue().arithmeticMean() );
                        w.write( "\t" + entry.getValue().median() );
                        w.write( "\n" );
                    }
                }
                w.flush();
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Number of domains: " );
                w.write( "\n" );
                if ( domain_number_stats_by_dc != null ) {
                    for( final Entry<Integer, DescriptiveStatistics> entry : dc_reapp_counts_to_domain_number_stats
                            .entrySet() ) {
                        w.write( entry.getKey().toString() );
                        w.write( "\t" + entry.getValue().arithmeticMean() );
                        w.write( "\t" + entry.getValue().median() );
                        w.write( "\n" );
                    }
                }
                w.flush();
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Gained once, domain lengths:" );
                w.write( "\n" );
                w.write( "N: " + gained_once_domain_length_count );
                w.write( "\n" );
                w.write( "Avg: " + ( ( double ) gained_once_domain_length_sum / gained_once_domain_length_count ) );
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Gained multiple times, domain lengths:" );
                w.write( "\n" );
                w.write( "N: " + gained_multiple_times_domain_length_count );
                w.write( "\n" );
                w.write( "Avg: "
                        + ( ( double ) gained_multiple_times_domain_length_sum / gained_multiple_times_domain_length_count ) );
                w.write( "\n" );
                w.write( "\n" );
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Gained once, protein lengths:" );
                w.write( "\n" );
                w.write( gained_once_lengths_stats.toString() );
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Gained once, domain counts:" );
                w.write( "\n" );
                w.write( gained_once_domain_count_stats.toString() );
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Gained multiple times, protein lengths:" );
                w.write( "\n" );
                w.write( gained_multiple_times_lengths_stats.toString() );
                w.write( "\n" );
                w.write( "\n" );
                w.write( "Gained multiple times, domain counts:" );
                w.write( "\n" );
                w.write( gained_multiple_times_domain_count_stats.toString() );
                w.flush();
                w.close();
            }
        }
        catch ( final IOException e ) {
            ForesterUtil.printWarningMessage( surfacing.PRG_NAME, "Failure to write: " + e );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote independent domain combination gains fitch counts to ["
                + outfilename_for_counts + "]" );
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote independent domain combination gains fitch lists to ["
                + outfilename_for_dc + "]" );
        ForesterUtil.programMessage( surfacing.PRG_NAME,
                                     "Wrote independent domain combination gains fitch lists to (for GO mapping) ["
                                             + outfilename_for_dc_for_go_mapping + "]" );
        ForesterUtil.programMessage( surfacing.PRG_NAME,
                                     "Wrote independent domain combination gains fitch lists to (for GO mapping, unique) ["
                                             + outfilename_for_dc_for_go_mapping_unique + "]" );
    }

    private final static void addToCountMap( final Map<String, Integer> map, final String s ) {
        if ( map.containsKey( s ) ) {
            map.put( s, map.get( s ) + 1 );
        }
        else {
            map.put( s, 1 );
        }
    }

    public static int calculateOverlap( final Domain domain, final List<Boolean> covered_positions ) {
        int overlap_count = 0;
        for( int i = domain.getFrom(); i <= domain.getTo(); ++i ) {
            if ( ( i < covered_positions.size() ) && ( covered_positions.get( i ) == true ) ) {
                ++overlap_count;
            }
        }
        return overlap_count;
    }

    public static void checkForOutputFileWriteability( final File outfile ) {
        final String error = ForesterUtil.isWritableFile( outfile );
        if ( !ForesterUtil.isEmpty( error ) ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, error );
        }
    }

    private static SortedSet<String> collectAllDomainsChangedOnSubtree( final PhylogenyNode subtree_root,
                                                                        final boolean get_gains ) {
        final SortedSet<String> domains = new TreeSet<String>();
        for( final PhylogenyNode descendant : PhylogenyMethods.getAllDescendants( subtree_root ) ) {
            final BinaryCharacters chars = descendant.getNodeData().getBinaryCharacters();
            if ( get_gains ) {
                domains.addAll( chars.getGainedCharacters() );
            }
            else {
                domains.addAll( chars.getLostCharacters() );
            }
        }
        return domains;
    }

    public static void collectChangedDomainCombinationsFromBinaryStatesMatrixAsListToFile( final CharacterStateMatrix<CharacterStateMatrix.GainLossStates> matrix,
                                                                                           final BinaryDomainCombination.DomainCombinationType dc_type,
                                                                                           final List<BinaryDomainCombination> all_binary_domains_combination_gained,
                                                                                           final boolean get_gains ) {
        final SortedSet<String> sorted_ids = new TreeSet<String>();
        for( int i = 0; i < matrix.getNumberOfIdentifiers(); ++i ) {
            sorted_ids.add( matrix.getIdentifier( i ) );
        }
        for( final String id : sorted_ids ) {
            for( int c = 0; c < matrix.getNumberOfCharacters(); ++c ) {
                if ( ( get_gains && ( matrix.getState( id, c ) == CharacterStateMatrix.GainLossStates.GAIN ) )
                        || ( !get_gains && ( matrix.getState( id, c ) == CharacterStateMatrix.GainLossStates.LOSS ) ) ) {
                    if ( dc_type == BinaryDomainCombination.DomainCombinationType.DIRECTED_ADJACTANT ) {
                        all_binary_domains_combination_gained.add( AdjactantDirectedBinaryDomainCombination
                                .createInstance( matrix.getCharacter( c ) ) );
                    }
                    else if ( dc_type == BinaryDomainCombination.DomainCombinationType.DIRECTED ) {
                        all_binary_domains_combination_gained.add( DirectedBinaryDomainCombination
                                .createInstance( matrix.getCharacter( c ) ) );
                    }
                    else {
                        all_binary_domains_combination_gained.add( BasicBinaryDomainCombination.createInstance( matrix
                                .getCharacter( c ) ) );
                    }
                }
            }
        }
    }

    private static File createBaseDirForPerNodeDomainFiles( final String base_dir,
                                                            final boolean domain_combinations,
                                                            final CharacterStateMatrix.GainLossStates state,
                                                            final String outfile ) {
        File per_node_go_mapped_domain_gain_loss_files_base_dir = new File( new File( outfile ).getParent()
                + ForesterUtil.FILE_SEPARATOR + base_dir );
        if ( !per_node_go_mapped_domain_gain_loss_files_base_dir.exists() ) {
            per_node_go_mapped_domain_gain_loss_files_base_dir.mkdir();
        }
        if ( domain_combinations ) {
            per_node_go_mapped_domain_gain_loss_files_base_dir = new File( per_node_go_mapped_domain_gain_loss_files_base_dir
                    + ForesterUtil.FILE_SEPARATOR + "DC" );
        }
        else {
            per_node_go_mapped_domain_gain_loss_files_base_dir = new File( per_node_go_mapped_domain_gain_loss_files_base_dir
                    + ForesterUtil.FILE_SEPARATOR + "DOMAINS" );
        }
        if ( !per_node_go_mapped_domain_gain_loss_files_base_dir.exists() ) {
            per_node_go_mapped_domain_gain_loss_files_base_dir.mkdir();
        }
        if ( state == GainLossStates.GAIN ) {
            per_node_go_mapped_domain_gain_loss_files_base_dir = new File( per_node_go_mapped_domain_gain_loss_files_base_dir
                    + ForesterUtil.FILE_SEPARATOR + "GAINS" );
        }
        else if ( state == GainLossStates.LOSS ) {
            per_node_go_mapped_domain_gain_loss_files_base_dir = new File( per_node_go_mapped_domain_gain_loss_files_base_dir
                    + ForesterUtil.FILE_SEPARATOR + "LOSSES" );
        }
        else {
            per_node_go_mapped_domain_gain_loss_files_base_dir = new File( per_node_go_mapped_domain_gain_loss_files_base_dir
                    + ForesterUtil.FILE_SEPARATOR + "PRESENT" );
        }
        if ( !per_node_go_mapped_domain_gain_loss_files_base_dir.exists() ) {
            per_node_go_mapped_domain_gain_loss_files_base_dir.mkdir();
        }
        return per_node_go_mapped_domain_gain_loss_files_base_dir;
    }

    public static Map<DomainId, List<GoId>> createDomainIdToGoIdMap( final List<PfamToGoMapping> pfam_to_go_mappings ) {
        final Map<DomainId, List<GoId>> domain_id_to_go_ids_map = new HashMap<DomainId, List<GoId>>( pfam_to_go_mappings
                .size() );
        for( final PfamToGoMapping pfam_to_go : pfam_to_go_mappings ) {
            if ( !domain_id_to_go_ids_map.containsKey( pfam_to_go.getKey() ) ) {
                domain_id_to_go_ids_map.put( pfam_to_go.getKey(), new ArrayList<GoId>() );
            }
            domain_id_to_go_ids_map.get( pfam_to_go.getKey() ).add( pfam_to_go.getValue() );
        }
        return domain_id_to_go_ids_map;
    }

    public static Map<DomainId, Set<String>> createDomainIdToSecondaryFeaturesMap( final File secondary_features_map_file )
            throws IOException {
        final BasicTable<String> primary_table = BasicTableParser.parse( secondary_features_map_file, '\t' );
        final Map<DomainId, Set<String>> map = new TreeMap<DomainId, Set<String>>();
        for( int r = 0; r < primary_table.getNumberOfRows(); ++r ) {
            final DomainId domain_id = new DomainId( primary_table.getValue( 0, r ) );
            if ( !map.containsKey( domain_id ) ) {
                map.put( domain_id, new HashSet<String>() );
            }
            map.get( domain_id ).add( primary_table.getValue( 1, r ) );
        }
        return map;
    }

    public static Phylogeny createNjTreeBasedOnMatrixToFile( final File nj_tree_outfile, final DistanceMatrix distance ) {
        checkForOutputFileWriteability( nj_tree_outfile );
        final NeighborJoining nj = NeighborJoining.createInstance();
        final Phylogeny phylogeny = nj.execute( ( BasicSymmetricalDistanceMatrix ) distance );
        phylogeny.setName( nj_tree_outfile.getName() );
        writePhylogenyToFile( phylogeny, nj_tree_outfile.toString() );
        return phylogeny;
    }

    private static SortedSet<BinaryDomainCombination> createSetOfAllBinaryDomainCombinationsPerGenome( final GenomeWideCombinableDomains gwcd ) {
        final SortedMap<DomainId, CombinableDomains> cds = gwcd.getAllCombinableDomainsIds();
        final SortedSet<BinaryDomainCombination> binary_combinations = new TreeSet<BinaryDomainCombination>();
        for( final DomainId domain_id : cds.keySet() ) {
            final CombinableDomains cd = cds.get( domain_id );
            binary_combinations.addAll( cd.toBinaryDomainCombinations() );
        }
        return binary_combinations;
    }

    public static void decoratePrintableDomainSimilarities( final SortedSet<DomainSimilarity> domain_similarities,
                                                            final Detailedness detailedness,
                                                            final GoAnnotationOutput go_annotation_output,
                                                            final Map<GoId, GoTerm> go_id_to_term_map,
                                                            final GoNameSpace go_namespace_limit ) {
        if ( ( go_namespace_limit != null ) && ( ( go_id_to_term_map == null ) || go_id_to_term_map.isEmpty() ) ) {
            throw new IllegalArgumentException( "attempt to use a GO namespace limit without a GO id to term map" );
        }
        for( final DomainSimilarity domain_similarity : domain_similarities ) {
            if ( domain_similarity instanceof PrintableDomainSimilarity ) {
                final PrintableDomainSimilarity printable_domain_similarity = ( PrintableDomainSimilarity ) domain_similarity;
                printable_domain_similarity.setDetailedness( detailedness );
                printable_domain_similarity.setGoAnnotationOutput( go_annotation_output );
                printable_domain_similarity.setGoIdToTermMap( go_id_to_term_map );
                printable_domain_similarity.setGoNamespaceLimit( go_namespace_limit );
            }
        }
    }

    public static void executeDomainLengthAnalysis( final String[][] input_file_properties,
                                                    final int number_of_genomes,
                                                    final DomainLengthsTable domain_lengths_table,
                                                    final File outfile ) throws IOException {
        final DecimalFormat df = new DecimalFormat( "#.00" );
        checkForOutputFileWriteability( outfile );
        final BufferedWriter out = new BufferedWriter( new FileWriter( outfile ) );
        out.write( "MEAN BASED STATISTICS PER SPECIES" );
        out.write( ForesterUtil.LINE_SEPARATOR );
        out.write( domain_lengths_table.createMeanBasedStatisticsPerSpeciesTable().toString() );
        out.write( ForesterUtil.LINE_SEPARATOR );
        out.write( ForesterUtil.LINE_SEPARATOR );
        final List<DomainLengths> domain_lengths_list = domain_lengths_table.getDomainLengthsList();
        out.write( "OUTLIER SPECIES PER DOMAIN (Z>=1.5)" );
        out.write( ForesterUtil.LINE_SEPARATOR );
        for( final DomainLengths domain_lengths : domain_lengths_list ) {
            final List<Species> species_list = domain_lengths.getMeanBasedOutlierSpecies( 1.5 );
            if ( species_list.size() > 0 ) {
                out.write( domain_lengths.getDomainId() + "\t" );
                for( final Species species : species_list ) {
                    out.write( species + "\t" );
                }
                out.write( ForesterUtil.LINE_SEPARATOR );
                // DescriptiveStatistics stats_for_domain = domain_lengths
                //         .calculateMeanBasedStatistics();
                //AsciiHistogram histo = new AsciiHistogram( stats_for_domain );
                //System.out.println( histo.toStringBuffer( 40, '=', 60, 4 ).toString() );
            }
        }
        out.write( ForesterUtil.LINE_SEPARATOR );
        out.write( ForesterUtil.LINE_SEPARATOR );
        out.write( "OUTLIER SPECIES (Z 1.0)" );
        out.write( ForesterUtil.LINE_SEPARATOR );
        final DescriptiveStatistics stats_for_all_species = domain_lengths_table
                .calculateMeanBasedStatisticsForAllSpecies();
        out.write( stats_for_all_species.asSummary() );
        out.write( ForesterUtil.LINE_SEPARATOR );
        final AsciiHistogram histo = new AsciiHistogram( stats_for_all_species );
        out.write( histo.toStringBuffer( 40, '=', 60, 4 ).toString() );
        out.write( ForesterUtil.LINE_SEPARATOR );
        final double population_sd = stats_for_all_species.sampleStandardDeviation();
        final double population_mean = stats_for_all_species.arithmeticMean();
        for( final Species species : domain_lengths_table.getSpecies() ) {
            final double x = domain_lengths_table.calculateMeanBasedStatisticsForSpecies( species ).arithmeticMean();
            final double z = ( x - population_mean ) / population_sd;
            out.write( species + "\t" + z );
            out.write( ForesterUtil.LINE_SEPARATOR );
        }
        out.write( ForesterUtil.LINE_SEPARATOR );
        for( final Species species : domain_lengths_table.getSpecies() ) {
            final DescriptiveStatistics stats_for_species = domain_lengths_table
                    .calculateMeanBasedStatisticsForSpecies( species );
            final double x = stats_for_species.arithmeticMean();
            final double z = ( x - population_mean ) / population_sd;
            if ( ( z <= -1.0 ) || ( z >= 1.0 ) ) {
                out.write( species + "\t" + df.format( z ) + "\t" + stats_for_species.asSummary() );
                out.write( ForesterUtil.LINE_SEPARATOR );
            }
        }
        out.close();
        //        final List<HistogramData> histogram_datas = new ArrayList<HistogramData>();
        //        for( int i = 0; i < number_of_genomes; ++i ) {
        //            final Species species = new BasicSpecies( input_file_properties[ i ][ 0 ] );
        //            histogram_datas
        //                    .add( new HistogramData( species.toString(), domain_lengths_table
        //                            .calculateMeanBasedStatisticsForSpecies( species )
        //                            .getDataAsDoubleArray(), 5, 600, null, 60 ) );
        //        }
        //        final HistogramsFrame hf = new HistogramsFrame( histogram_datas );
        //        hf.setVisible( true );
        System.gc();
    }

    /**
     * 
     * @param all_binary_domains_combination_lost_fitch 
     * @param consider_directedness_and_adjacency_for_bin_combinations 
     * @param all_binary_domains_combination_gained if null ignored, otherwise this is to list all binary domain combinations
     * which were gained under unweighted (Fitch) parsimony.
     */
    public static void executeParsimonyAnalysis( final long random_number_seed_for_fitch_parsimony,
                                                 final boolean radomize_fitch_parsimony,
                                                 final String outfile_name,
                                                 final DomainParsimonyCalculator domain_parsimony,
                                                 final Phylogeny phylogeny,
                                                 final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                 final Map<GoId, GoTerm> go_id_to_term_map,
                                                 final GoNameSpace go_namespace_limit,
                                                 final String parameters_str,
                                                 final Map<DomainId, Set<String>>[] domain_id_to_secondary_features_maps,
                                                 final SortedSet<DomainId> positive_filter,
                                                 final boolean output_binary_domain_combinations_for_graphs,
                                                 final List<BinaryDomainCombination> all_binary_domains_combination_gained_fitch,
                                                 final List<BinaryDomainCombination> all_binary_domains_combination_lost_fitch,
                                                 final BinaryDomainCombination.DomainCombinationType dc_type,
                                                 final Map<String, DescriptiveStatistics> protein_length_stats_by_dc,
                                                 final Map<String, DescriptiveStatistics> domain_number_stats_by_dc,
                                                 final Map<String, DescriptiveStatistics> domain_length_stats_by_domain ) {
        final String sep = ForesterUtil.LINE_SEPARATOR + "###################" + ForesterUtil.LINE_SEPARATOR;
        final String date_time = ForesterUtil.getCurrentDateTime();
        final SortedSet<String> all_pfams_encountered = new TreeSet<String>();
        final SortedSet<String> all_pfams_gained_as_domains = new TreeSet<String>();
        final SortedSet<String> all_pfams_lost_as_domains = new TreeSet<String>();
        final SortedSet<String> all_pfams_gained_as_dom_combinations = new TreeSet<String>();
        final SortedSet<String> all_pfams_lost_as_dom_combinations = new TreeSet<String>();
        writeToNexus( outfile_name, domain_parsimony, phylogeny );
        // DOLLO DOMAINS
        // -------------
        Phylogeny local_phylogeny_l = phylogeny.copy();
        if ( ( positive_filter != null ) && ( positive_filter.size() > 0 ) ) {
            domain_parsimony.executeDolloParsimonyOnDomainPresence( positive_filter );
        }
        else {
            domain_parsimony.executeDolloParsimonyOnDomainPresence();
        }
        SurfacingUtil.writeMatrixToFile( domain_parsimony.getGainLossMatrix(), outfile_name
                + surfacing.PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_DOMAINS, Format.FORESTER );
        SurfacingUtil.writeMatrixToFile( domain_parsimony.getGainLossCountsMatrix(), outfile_name
                + surfacing.PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_DOMAINS, Format.FORESTER );
        SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                           CharacterStateMatrix.GainLossStates.GAIN,
                                                           outfile_name + surfacing.PARSIMONY_OUTPUT_DOLLO_GAINS_D,
                                                           sep,
                                                           ForesterUtil.LINE_SEPARATOR,
                                                           null );
        SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                           CharacterStateMatrix.GainLossStates.LOSS,
                                                           outfile_name + surfacing.PARSIMONY_OUTPUT_DOLLO_LOSSES_D,
                                                           sep,
                                                           ForesterUtil.LINE_SEPARATOR,
                                                           null );
        SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(), null, outfile_name
                + surfacing.PARSIMONY_OUTPUT_DOLLO_PRESENT_D, sep, ForesterUtil.LINE_SEPARATOR, null );
        //HTML:
        writeBinaryStatesMatrixToList( domain_id_to_go_ids_map,
                                       go_id_to_term_map,
                                       go_namespace_limit,
                                       false,
                                       domain_parsimony.getGainLossMatrix(),
                                       CharacterStateMatrix.GainLossStates.GAIN,
                                       outfile_name + surfacing.PARSIMONY_OUTPUT_DOLLO_GAINS_HTML_D,
                                       sep,
                                       ForesterUtil.LINE_SEPARATOR,
                                       "Dollo Parsimony | Gains | Domains",
                                       "+",
                                       domain_id_to_secondary_features_maps,
                                       all_pfams_encountered,
                                       all_pfams_gained_as_domains,
                                       "_dollo_gains_d" );
        writeBinaryStatesMatrixToList( domain_id_to_go_ids_map,
                                       go_id_to_term_map,
                                       go_namespace_limit,
                                       false,
                                       domain_parsimony.getGainLossMatrix(),
                                       CharacterStateMatrix.GainLossStates.LOSS,
                                       outfile_name + surfacing.PARSIMONY_OUTPUT_DOLLO_LOSSES_HTML_D,
                                       sep,
                                       ForesterUtil.LINE_SEPARATOR,
                                       "Dollo Parsimony | Losses | Domains",
                                       "-",
                                       domain_id_to_secondary_features_maps,
                                       all_pfams_encountered,
                                       all_pfams_lost_as_domains,
                                       "_dollo_losses_d" );
        writeBinaryStatesMatrixToList( domain_id_to_go_ids_map,
                                       go_id_to_term_map,
                                       go_namespace_limit,
                                       false,
                                       domain_parsimony.getGainLossMatrix(),
                                       null,
                                       outfile_name + surfacing.PARSIMONY_OUTPUT_DOLLO_PRESENT_HTML_D,
                                       sep,
                                       ForesterUtil.LINE_SEPARATOR,
                                       "Dollo Parsimony | Present | Domains",
                                       "",
                                       domain_id_to_secondary_features_maps,
                                       all_pfams_encountered,
                                       null,
                                       "_dollo_present_d" );
        preparePhylogeny( local_phylogeny_l,
                          domain_parsimony,
                          date_time,
                          "Dollo parsimony on domain presence/absence",
                          "dollo_on_domains_" + outfile_name,
                          parameters_str );
        SurfacingUtil.writePhylogenyToFile( local_phylogeny_l, outfile_name
                + surfacing.DOMAINS_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO );
        try {
            writeAllDomainsChangedOnAllSubtrees( local_phylogeny_l, true, outfile_name, "_dollo_all_gains_d" );
            writeAllDomainsChangedOnAllSubtrees( local_phylogeny_l, false, outfile_name, "_dollo_all_losses_d" );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getLocalizedMessage() );
        }
        if ( domain_parsimony.calculateNumberOfBinaryDomainCombination() > 0 ) {
            // FITCH DOMAIN COMBINATIONS
            // -------------------------
            local_phylogeny_l = phylogeny.copy();
            String randomization = "no";
            if ( radomize_fitch_parsimony ) {
                domain_parsimony.executeFitchParsimonyOnBinaryDomainCombintion( random_number_seed_for_fitch_parsimony );
                randomization = "yes, seed = " + random_number_seed_for_fitch_parsimony;
            }
            else {
                domain_parsimony.executeFitchParsimonyOnBinaryDomainCombintion( USE_LAST );
            }
            SurfacingUtil.writeMatrixToFile( domain_parsimony.getGainLossMatrix(), outfile_name
                    + surfacing.PARSIMONY_OUTPUT_GL_SUFFIX_FITCH_BINARY_COMBINATIONS, Format.FORESTER );
            SurfacingUtil.writeMatrixToFile( domain_parsimony.getGainLossCountsMatrix(), outfile_name
                    + surfacing.PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_FITCH_BINARY_COMBINATIONS, Format.FORESTER );
            SurfacingUtil
                    .writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                          CharacterStateMatrix.GainLossStates.GAIN,
                                                          outfile_name + surfacing.PARSIMONY_OUTPUT_FITCH_GAINS_BC,
                                                          sep,
                                                          ForesterUtil.LINE_SEPARATOR,
                                                          null );
            SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                               CharacterStateMatrix.GainLossStates.LOSS,
                                                               outfile_name
                                                                       + surfacing.PARSIMONY_OUTPUT_FITCH_LOSSES_BC,
                                                               sep,
                                                               ForesterUtil.LINE_SEPARATOR,
                                                               null );
            SurfacingUtil.writeBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(), null, outfile_name
                    + surfacing.PARSIMONY_OUTPUT_FITCH_PRESENT_BC, sep, ForesterUtil.LINE_SEPARATOR, null );
            if ( all_binary_domains_combination_gained_fitch != null ) {
                collectChangedDomainCombinationsFromBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                                                    dc_type,
                                                                                    all_binary_domains_combination_gained_fitch,
                                                                                    true );
            }
            if ( all_binary_domains_combination_lost_fitch != null ) {
                collectChangedDomainCombinationsFromBinaryStatesMatrixAsListToFile( domain_parsimony.getGainLossMatrix(),
                                                                                    dc_type,
                                                                                    all_binary_domains_combination_lost_fitch,
                                                                                    false );
            }
            if ( output_binary_domain_combinations_for_graphs ) {
                SurfacingUtil
                        .writeBinaryStatesMatrixAsListToFileForBinaryCombinationsForGraphAnalysis( domain_parsimony
                                                                                                           .getGainLossMatrix(),
                                                                                                   null,
                                                                                                   outfile_name
                                                                                                           + surfacing.PARSIMONY_OUTPUT_FITCH_PRESENT_BC_OUTPUTFILE_SUFFIX_FOR_GRAPH_ANALYSIS,
                                                                                                   sep,
                                                                                                   ForesterUtil.LINE_SEPARATOR,
                                                                                                   BinaryDomainCombination.OutputFormat.DOT );
            }
            // HTML:
            writeBinaryStatesMatrixToList( domain_id_to_go_ids_map,
                                           go_id_to_term_map,
                                           go_namespace_limit,
                                           true,
                                           domain_parsimony.getGainLossMatrix(),
                                           CharacterStateMatrix.GainLossStates.GAIN,
                                           outfile_name + surfacing.PARSIMONY_OUTPUT_FITCH_GAINS_HTML_BC,
                                           sep,
                                           ForesterUtil.LINE_SEPARATOR,
                                           "Fitch Parsimony | Gains | Domain Combinations",
                                           "+",
                                           null,
                                           all_pfams_encountered,
                                           all_pfams_gained_as_dom_combinations,
                                           "_fitch_gains_dc" );
            writeBinaryStatesMatrixToList( domain_id_to_go_ids_map,
                                           go_id_to_term_map,
                                           go_namespace_limit,
                                           true,
                                           domain_parsimony.getGainLossMatrix(),
                                           CharacterStateMatrix.GainLossStates.LOSS,
                                           outfile_name + surfacing.PARSIMONY_OUTPUT_FITCH_LOSSES_HTML_BC,
                                           sep,
                                           ForesterUtil.LINE_SEPARATOR,
                                           "Fitch Parsimony | Losses | Domain Combinations",
                                           "-",
                                           null,
                                           all_pfams_encountered,
                                           all_pfams_lost_as_dom_combinations,
                                           "_fitch_losses_dc" );
            writeBinaryStatesMatrixToList( domain_id_to_go_ids_map,
                                           go_id_to_term_map,
                                           go_namespace_limit,
                                           true,
                                           domain_parsimony.getGainLossMatrix(),
                                           null,
                                           outfile_name + surfacing.PARSIMONY_OUTPUT_FITCH_PRESENT_HTML_BC,
                                           sep,
                                           ForesterUtil.LINE_SEPARATOR,
                                           "Fitch Parsimony | Present | Domain Combinations",
                                           "",
                                           null,
                                           all_pfams_encountered,
                                           null,
                                           "_fitch_present_dc" );
            writeAllEncounteredPfamsToFile( domain_id_to_go_ids_map,
                                            go_id_to_term_map,
                                            outfile_name,
                                            all_pfams_encountered );
            writePfamsToFile( outfile_name + surfacing.ALL_PFAMS_GAINED_AS_DOMAINS_SUFFIX, all_pfams_gained_as_domains );
            writePfamsToFile( outfile_name + surfacing.ALL_PFAMS_LOST_AS_DOMAINS_SUFFIX, all_pfams_lost_as_domains );
            writePfamsToFile( outfile_name + surfacing.ALL_PFAMS_GAINED_AS_DC_SUFFIX,
                              all_pfams_gained_as_dom_combinations );
            writePfamsToFile( outfile_name + surfacing.ALL_PFAMS_LOST_AS_DC_SUFFIX, all_pfams_lost_as_dom_combinations );
            preparePhylogeny( local_phylogeny_l,
                              domain_parsimony,
                              date_time,
                              "Fitch parsimony on binary domain combination presence/absence randomization: "
                                      + randomization,
                              "fitch_on_binary_domain_combinations_" + outfile_name,
                              parameters_str );
            SurfacingUtil.writePhylogenyToFile( local_phylogeny_l, outfile_name
                    + surfacing.BINARY_DOMAIN_COMBINATIONS_PARSIMONY_TREE_OUTPUT_SUFFIX_FITCH );
            calculateIndependentDomainCombinationGains( local_phylogeny_l,
                                                        outfile_name
                                                                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_COUNTS_OUTPUT_SUFFIX,
                                                        outfile_name
                                                                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_DC_OUTPUT_SUFFIX,
                                                        outfile_name
                                                                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_OUTPUT_SUFFIX,
                                                        outfile_name
                                                                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_OUTPUT_UNIQUE_SUFFIX,
                                                        outfile_name + "_indep_dc_gains_fitch_lca_ranks.txt",
                                                        outfile_name + "_indep_dc_gains_fitch_lca_taxonomies.txt",
                                                        outfile_name + "_indep_dc_gains_fitch_protein_statistics.txt",
                                                        protein_length_stats_by_dc,
                                                        domain_number_stats_by_dc,
                                                        domain_length_stats_by_domain );
        }
    }

    public static void executeParsimonyAnalysisForSecondaryFeatures( final String outfile_name,
                                                                     final DomainParsimonyCalculator secondary_features_parsimony,
                                                                     final Phylogeny phylogeny,
                                                                     final String parameters_str,
                                                                     final Map<Species, MappingResults> mapping_results_map ) {
        final String sep = ForesterUtil.LINE_SEPARATOR + "###################" + ForesterUtil.LINE_SEPARATOR;
        final String date_time = ForesterUtil.getCurrentDateTime();
        System.out.println();
        writeToNexus( outfile_name + surfacing.NEXUS_SECONDARY_FEATURES,
                      secondary_features_parsimony.createMatrixOfSecondaryFeaturePresenceOrAbsence( null ),
                      phylogeny );
        Phylogeny local_phylogeny_copy = phylogeny.copy();
        secondary_features_parsimony.executeDolloParsimonyOnSecondaryFeatures( mapping_results_map );
        SurfacingUtil.writeMatrixToFile( secondary_features_parsimony.getGainLossMatrix(), outfile_name
                + surfacing.PARSIMONY_OUTPUT_GL_SUFFIX_DOLLO_SECONDARY_FEATURES, Format.FORESTER );
        SurfacingUtil.writeMatrixToFile( secondary_features_parsimony.getGainLossCountsMatrix(), outfile_name
                + surfacing.PARSIMONY_OUTPUT_GL_COUNTS_SUFFIX_DOLLO_SECONDARY_FEATURES, Format.FORESTER );
        SurfacingUtil
                .writeBinaryStatesMatrixAsListToFile( secondary_features_parsimony.getGainLossMatrix(),
                                                      CharacterStateMatrix.GainLossStates.GAIN,
                                                      outfile_name
                                                              + surfacing.PARSIMONY_OUTPUT_DOLLO_GAINS_SECONDARY_FEATURES,
                                                      sep,
                                                      ForesterUtil.LINE_SEPARATOR,
                                                      null );
        SurfacingUtil
                .writeBinaryStatesMatrixAsListToFile( secondary_features_parsimony.getGainLossMatrix(),
                                                      CharacterStateMatrix.GainLossStates.LOSS,
                                                      outfile_name
                                                              + surfacing.PARSIMONY_OUTPUT_DOLLO_LOSSES_SECONDARY_FEATURES,
                                                      sep,
                                                      ForesterUtil.LINE_SEPARATOR,
                                                      null );
        SurfacingUtil
                .writeBinaryStatesMatrixAsListToFile( secondary_features_parsimony.getGainLossMatrix(),
                                                      null,
                                                      outfile_name
                                                              + surfacing.PARSIMONY_OUTPUT_DOLLO_PRESENT_SECONDARY_FEATURES,
                                                      sep,
                                                      ForesterUtil.LINE_SEPARATOR,
                                                      null );
        preparePhylogeny( local_phylogeny_copy,
                          secondary_features_parsimony,
                          date_time,
                          "Dollo parsimony on secondary feature presence/absence",
                          "dollo_on_secondary_features_" + outfile_name,
                          parameters_str );
        SurfacingUtil.writePhylogenyToFile( local_phylogeny_copy, outfile_name
                + surfacing.SECONDARY_FEATURES_PARSIMONY_TREE_OUTPUT_SUFFIX_DOLLO );
        // FITCH DOMAIN COMBINATIONS
        // -------------------------
        local_phylogeny_copy = phylogeny.copy();
        final String randomization = "no";
        secondary_features_parsimony.executeFitchParsimonyOnBinaryDomainCombintionOnSecondaryFeatures( USE_LAST );
        preparePhylogeny( local_phylogeny_copy,
                          secondary_features_parsimony,
                          date_time,
                          "Fitch parsimony on secondary binary domain combination presence/absence randomization: "
                                  + randomization,
                          "fitch_on_binary_domain_combinations_" + outfile_name,
                          parameters_str );
        SurfacingUtil.writePhylogenyToFile( local_phylogeny_copy, outfile_name
                + surfacing.BINARY_DOMAIN_COMBINATIONS_PARSIMONY_TREE_OUTPUT_SUFFIX_FITCH_MAPPED );
        calculateIndependentDomainCombinationGains( local_phylogeny_copy, outfile_name
                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_COUNTS_MAPPED_OUTPUT_SUFFIX, outfile_name
                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_DC_MAPPED_OUTPUT_SUFFIX, outfile_name
                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_MAPPED_OUTPUT_SUFFIX, outfile_name
                + surfacing.INDEPENDENT_DC_GAINS_FITCH_PARS_DC_FOR_GO_MAPPING_MAPPED_OUTPUT_UNIQUE_SUFFIX, outfile_name
                + "_MAPPED_indep_dc_gains_fitch_lca_ranks.txt", outfile_name
                + "_MAPPED_indep_dc_gains_fitch_lca_taxonomies.txt", null, null, null, null );
    }

    public static void doit( final List<Protein> proteins,
                             final List<DomainId> query_domain_ids_nc_order,
                             final Writer out,
                             final String separator,
                             final String limit_to_species,
                             final Map<String, List<Integer>> average_protein_lengths_by_dc ) throws IOException {
        for( final Protein protein : proteins ) {
            if ( ForesterUtil.isEmpty( limit_to_species )
                    || protein.getSpecies().getSpeciesId().equalsIgnoreCase( limit_to_species ) ) {
                if ( protein.contains( query_domain_ids_nc_order, true ) ) {
                    out.write( protein.getSpecies().getSpeciesId() );
                    out.write( separator );
                    out.write( protein.getProteinId().getId() );
                    out.write( separator );
                    out.write( "[" );
                    final Set<DomainId> visited_domain_ids = new HashSet<DomainId>();
                    boolean first = true;
                    for( final Domain domain : protein.getProteinDomains() ) {
                        if ( !visited_domain_ids.contains( domain.getDomainId() ) ) {
                            visited_domain_ids.add( domain.getDomainId() );
                            if ( first ) {
                                first = false;
                            }
                            else {
                                out.write( " " );
                            }
                            out.write( domain.getDomainId().getId() );
                            out.write( " {" );
                            out.write( "" + domain.getTotalCount() );
                            out.write( "}" );
                        }
                    }
                    out.write( "]" );
                    out.write( separator );
                    if ( !( ForesterUtil.isEmpty( protein.getDescription() ) || protein.getDescription()
                            .equals( SurfacingConstants.NONE ) ) ) {
                        out.write( protein.getDescription() );
                    }
                    out.write( separator );
                    if ( !( ForesterUtil.isEmpty( protein.getAccession() ) || protein.getAccession()
                            .equals( SurfacingConstants.NONE ) ) ) {
                        out.write( protein.getAccession() );
                    }
                    out.write( SurfacingConstants.NL );
                }
            }
        }
        out.flush();
    }

    public static void extractProteinNames( final List<Protein> proteins,
                                            final List<DomainId> query_domain_ids_nc_order,
                                            final Writer out,
                                            final String separator,
                                            final String limit_to_species ) throws IOException {
        for( final Protein protein : proteins ) {
            if ( ForesterUtil.isEmpty( limit_to_species )
                    || protein.getSpecies().getSpeciesId().equalsIgnoreCase( limit_to_species ) ) {
                if ( protein.contains( query_domain_ids_nc_order, true ) ) {
                    out.write( protein.getSpecies().getSpeciesId() );
                    out.write( separator );
                    out.write( protein.getProteinId().getId() );
                    out.write( separator );
                    out.write( "[" );
                    final Set<DomainId> visited_domain_ids = new HashSet<DomainId>();
                    boolean first = true;
                    for( final Domain domain : protein.getProteinDomains() ) {
                        if ( !visited_domain_ids.contains( domain.getDomainId() ) ) {
                            visited_domain_ids.add( domain.getDomainId() );
                            if ( first ) {
                                first = false;
                            }
                            else {
                                out.write( " " );
                            }
                            out.write( domain.getDomainId().getId() );
                            out.write( " {" );
                            out.write( "" + domain.getTotalCount() );
                            out.write( "}" );
                        }
                    }
                    out.write( "]" );
                    out.write( separator );
                    if ( !( ForesterUtil.isEmpty( protein.getDescription() ) || protein.getDescription()
                            .equals( SurfacingConstants.NONE ) ) ) {
                        out.write( protein.getDescription() );
                    }
                    out.write( separator );
                    if ( !( ForesterUtil.isEmpty( protein.getAccession() ) || protein.getAccession()
                            .equals( SurfacingConstants.NONE ) ) ) {
                        out.write( protein.getAccession() );
                    }
                    out.write( SurfacingConstants.NL );
                }
            }
        }
        out.flush();
    }

    public static void extractProteinNames( final SortedMap<Species, List<Protein>> protein_lists_per_species,
                                            final DomainId domain_id,
                                            final Writer out,
                                            final String separator,
                                            final String limit_to_species,
                                            final double domain_e_cutoff ) throws IOException {
        System.out.println( "Per domain E-value: " + domain_e_cutoff );
        for( final Species species : protein_lists_per_species.keySet() ) {
            System.out.println( species + ":" );
            for( final Protein protein : protein_lists_per_species.get( species ) ) {
                if ( ForesterUtil.isEmpty( limit_to_species )
                        || protein.getSpecies().getSpeciesId().equalsIgnoreCase( limit_to_species ) ) {
                    final List<Domain> domains = protein.getProteinDomains( domain_id );
                    if ( domains.size() > 0 ) {
                        out.write( protein.getSpecies().getSpeciesId() );
                        out.write( separator );
                        out.write( protein.getProteinId().getId() );
                        out.write( separator );
                        out.write( domain_id.toString() );
                        out.write( separator );
                        int prev_to = -1;
                        for( final Domain domain : domains ) {
                            if ( ( domain_e_cutoff < 0 ) || ( domain.getPerDomainEvalue() <= domain_e_cutoff ) ) {
                                out.write( "/" );
                                out.write( domain.getFrom() + "-" + domain.getTo() );
                                if ( prev_to >= 0 ) {
                                    final int l = domain.getFrom() - prev_to;
                                    System.out.println( l );
                                }
                                prev_to = domain.getTo();
                            }
                        }
                        out.write( "/" );
                        out.write( separator );
                        final List<Domain> domain_list = new ArrayList<Domain>();
                        for( final Domain domain : protein.getProteinDomains() ) {
                            if ( ( domain_e_cutoff < 0 ) || ( domain.getPerDomainEvalue() <= domain_e_cutoff ) ) {
                                domain_list.add( domain );
                            }
                        }
                        final Domain domain_ary[] = new Domain[ domain_list.size() ];
                        for( int i = 0; i < domain_list.size(); ++i ) {
                            domain_ary[ i ] = domain_list.get( i );
                        }
                        Arrays.sort( domain_ary, new DomainComparator( true ) );
                        out.write( "{" );
                        boolean first = true;
                        for( final Domain domain : domain_ary ) {
                            if ( first ) {
                                first = false;
                            }
                            else {
                                out.write( "," );
                            }
                            out.write( domain.getDomainId().toString() );
                            out.write( ":" + domain.getFrom() + "-" + domain.getTo() );
                            out.write( ":" + domain.getPerDomainEvalue() );
                        }
                        out.write( "}" );
                        if ( !( ForesterUtil.isEmpty( protein.getDescription() ) || protein.getDescription()
                                .equals( SurfacingConstants.NONE ) ) ) {
                            out.write( protein.getDescription() );
                        }
                        out.write( separator );
                        if ( !( ForesterUtil.isEmpty( protein.getAccession() ) || protein.getAccession()
                                .equals( SurfacingConstants.NONE ) ) ) {
                            out.write( protein.getAccession() );
                        }
                        out.write( SurfacingConstants.NL );
                    }
                }
            }
        }
        out.flush();
    }

    public static SortedSet<DomainId> getAllDomainIds( final List<GenomeWideCombinableDomains> gwcd_list ) {
        final SortedSet<DomainId> all_domains_ids = new TreeSet<DomainId>();
        for( final GenomeWideCombinableDomains gwcd : gwcd_list ) {
            final Set<DomainId> all_domains = gwcd.getAllDomainIds();
            //    for( final Domain domain : all_domains ) {
            all_domains_ids.addAll( all_domains );
            //    }
        }
        return all_domains_ids;
    }

    public static SortedMap<String, Integer> getDomainCounts( final List<Protein> protein_domain_collections ) {
        final SortedMap<String, Integer> map = new TreeMap<String, Integer>();
        for( final Protein protein_domain_collection : protein_domain_collections ) {
            for( final Object name : protein_domain_collection.getProteinDomains() ) {
                final BasicDomain protein_domain = ( BasicDomain ) name;
                final String id = protein_domain.getDomainId().getId();
                if ( map.containsKey( id ) ) {
                    map.put( id, map.get( id ) + 1 );
                }
                else {
                    map.put( id, 1 );
                }
            }
        }
        return map;
    }

    public static int getNumberOfNodesLackingName( final Phylogeny p, final StringBuilder names ) {
        final PhylogenyNodeIterator it = p.iteratorPostorder();
        int c = 0;
        while ( it.hasNext() ) {
            final PhylogenyNode n = it.next();
            if ( ForesterUtil.isEmpty( n.getName() )
                    && ( !n.getNodeData().isHasTaxonomy() || ForesterUtil.isEmpty( n.getNodeData().getTaxonomy()
                            .getScientificName() ) )
                    && ( !n.getNodeData().isHasTaxonomy() || ForesterUtil.isEmpty( n.getNodeData().getTaxonomy()
                            .getCommonName() ) ) ) {
                if ( n.getParent() != null ) {
                    names.append( " " );
                    names.append( n.getParent().getName() );
                }
                final List l = n.getAllExternalDescendants();
                for( final Object object : l ) {
                    System.out.println( l.toString() );
                }
                ++c;
            }
        }
        return c;
    }

    /**
     * Returns true is Domain domain falls in an uninterrupted stretch of
     * covered positions.
     * 
     * @param domain
     * @param covered_positions
     * @return
     */
    public static boolean isEngulfed( final Domain domain, final List<Boolean> covered_positions ) {
        for( int i = domain.getFrom(); i <= domain.getTo(); ++i ) {
            if ( ( i >= covered_positions.size() ) || ( covered_positions.get( i ) != true ) ) {
                return false;
            }
        }
        return true;
    }

    public static void preparePhylogeny( final Phylogeny p,
                                         final DomainParsimonyCalculator domain_parsimony,
                                         final String date_time,
                                         final String method,
                                         final String name,
                                         final String parameters_str ) {
        domain_parsimony.decoratePhylogenyWithDomains( p );
        final StringBuilder desc = new StringBuilder();
        desc.append( "[Method: " + method + "] [Date: " + date_time + "] " );
        desc.append( "[Cost: " + domain_parsimony.getCost() + "] " );
        desc.append( "[Gains: " + domain_parsimony.getTotalGains() + "] " );
        desc.append( "[Losses: " + domain_parsimony.getTotalLosses() + "] " );
        desc.append( "[Unchanged: " + domain_parsimony.getTotalUnchanged() + "] " );
        desc.append( "[Parameters: " + parameters_str + "]" );
        p.setName( name );
        p.setDescription( desc.toString() );
        p.setConfidence( new Confidence( domain_parsimony.getCost(), "parsimony" ) );
        p.setRerootable( false );
        p.setRooted( true );
    }

    /*
     * species | protein id | n-terminal domain | c-terminal domain | n-terminal domain per domain E-value | c-terminal domain per domain E-value
     * 
     * 
     */
    static public StringBuffer proteinToDomainCombinations( final Protein protein,
                                                            final String protein_id,
                                                            final String separator ) {
        final StringBuffer sb = new StringBuffer();
        if ( protein.getSpecies() == null ) {
            throw new IllegalArgumentException( "species must not be null" );
        }
        if ( ForesterUtil.isEmpty( protein.getSpecies().getSpeciesId() ) ) {
            throw new IllegalArgumentException( "species id must not be empty" );
        }
        final List<Domain> domains = protein.getProteinDomains();
        if ( domains.size() > 1 ) {
            final Map<String, Integer> counts = new HashMap<String, Integer>();
            for( final Domain domain : domains ) {
                final String id = domain.getDomainId().getId();
                if ( counts.containsKey( id ) ) {
                    counts.put( id, counts.get( id ) + 1 );
                }
                else {
                    counts.put( id, 1 );
                }
            }
            final Set<String> dcs = new HashSet<String>();
            for( int i = 1; i < domains.size(); ++i ) {
                for( int j = 0; j < i; ++j ) {
                    Domain domain_n = domains.get( i );
                    Domain domain_c = domains.get( j );
                    if ( domain_n.getFrom() > domain_c.getFrom() ) {
                        domain_n = domains.get( j );
                        domain_c = domains.get( i );
                    }
                    final String dc = domain_n.getDomainId().getId() + domain_c.getDomainId().getId();
                    if ( !dcs.contains( dc ) ) {
                        dcs.add( dc );
                        sb.append( protein.getSpecies() );
                        sb.append( separator );
                        sb.append( protein_id );
                        sb.append( separator );
                        sb.append( domain_n.getDomainId().getId() );
                        sb.append( separator );
                        sb.append( domain_c.getDomainId().getId() );
                        sb.append( separator );
                        sb.append( domain_n.getPerDomainEvalue() );
                        sb.append( separator );
                        sb.append( domain_c.getPerDomainEvalue() );
                        sb.append( separator );
                        sb.append( counts.get( domain_n.getDomainId().getId() ) );
                        sb.append( separator );
                        sb.append( counts.get( domain_c.getDomainId().getId() ) );
                        sb.append( ForesterUtil.LINE_SEPARATOR );
                    }
                }
            }
        }
        else if ( domains.size() == 1 ) {
            sb.append( protein.getSpecies() );
            sb.append( separator );
            sb.append( protein_id );
            sb.append( separator );
            sb.append( domains.get( 0 ).getDomainId().getId() );
            sb.append( separator );
            sb.append( separator );
            sb.append( domains.get( 0 ).getPerDomainEvalue() );
            sb.append( separator );
            sb.append( separator );
            sb.append( 1 );
            sb.append( separator );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        else {
            sb.append( protein.getSpecies() );
            sb.append( separator );
            sb.append( protein_id );
            sb.append( separator );
            sb.append( separator );
            sb.append( separator );
            sb.append( separator );
            sb.append( separator );
            sb.append( separator );
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        return sb;
    }

    /**
     * 
     * Example regarding engulfment: ------------0.1 ----------0.2 --0.3 =>
     * domain with 0.3 is ignored
     * 
     * -----------0.1 ----------0.2 --0.3 => domain with 0.3 is ignored
     * 
     * 
     * ------------0.1 ----------0.3 --0.2 => domains with 0.3 and 0.2 are _not_
     * ignored
     * 
     * @param max_allowed_overlap
     *            maximal allowed overlap (inclusive) to be still considered not
     *            overlapping (zero or negative value to allow any overlap)
     * @param remove_engulfed_domains
     *            to remove domains which are completely engulfed by coverage of
     *            domains with better support
     * @param protein
     * @return
     */
    public static Protein removeOverlappingDomains( final int max_allowed_overlap,
                                                    final boolean remove_engulfed_domains,
                                                    final Protein protein ) {
        final Protein pruned_protein = new BasicProtein( protein.getProteinId().getId(), protein.getSpecies()
                .getSpeciesId(), protein.getLength() );
        final List<Domain> sorted = SurfacingUtil.sortDomainsWithAscendingConfidenceValues( protein );
        final List<Boolean> covered_positions = new ArrayList<Boolean>();
        for( final Domain domain : sorted ) {
            if ( ( ( max_allowed_overlap < 0 ) || ( SurfacingUtil.calculateOverlap( domain, covered_positions ) <= max_allowed_overlap ) )
                    && ( !remove_engulfed_domains || !isEngulfed( domain, covered_positions ) ) ) {
                final int covered_positions_size = covered_positions.size();
                for( int i = covered_positions_size; i < domain.getFrom(); ++i ) {
                    covered_positions.add( false );
                }
                final int new_covered_positions_size = covered_positions.size();
                for( int i = domain.getFrom(); i <= domain.getTo(); ++i ) {
                    if ( i < new_covered_positions_size ) {
                        covered_positions.set( i, true );
                    }
                    else {
                        covered_positions.add( true );
                    }
                }
                pruned_protein.addProteinDomain( domain );
            }
        }
        return pruned_protein;
    }

    public static List<Domain> sortDomainsWithAscendingConfidenceValues( final Protein protein ) {
        final List<Domain> domains = new ArrayList<Domain>();
        for( final Domain d : protein.getProteinDomains() ) {
            domains.add( d );
        }
        Collections.sort( domains, SurfacingUtil.ASCENDING_CONFIDENCE_VALUE_ORDER );
        return domains;
    }

    private static List<String> splitDomainCombination( final String dc ) {
        final String[] s = dc.split( "=" );
        if ( s.length != 2 ) {
            ForesterUtil.printErrorMessage( surfacing.PRG_NAME, "Stringyfied domain combination has illegal format: "
                    + dc );
            System.exit( -1 );
        }
        final List<String> l = new ArrayList<String>( 2 );
        l.add( s[ 0 ] );
        l.add( s[ 1 ] );
        return l;
    }

    public static void writeAllDomainsChangedOnAllSubtrees( final Phylogeny p,
                                                            final boolean get_gains,
                                                            final String outdir,
                                                            final String suffix_for_filename ) throws IOException {
        CharacterStateMatrix.GainLossStates state = CharacterStateMatrix.GainLossStates.GAIN;
        if ( !get_gains ) {
            state = CharacterStateMatrix.GainLossStates.LOSS;
        }
        final File base_dir = createBaseDirForPerNodeDomainFiles( surfacing.BASE_DIRECTORY_PER_SUBTREE_DOMAIN_GAIN_LOSS_FILES,
                                                                  false,
                                                                  state,
                                                                  outdir );
        for( final PhylogenyNodeIterator it = p.iteratorPostorder(); it.hasNext(); ) {
            final PhylogenyNode node = it.next();
            if ( !node.isExternal() ) {
                final SortedSet<String> domains = collectAllDomainsChangedOnSubtree( node, get_gains );
                if ( domains.size() > 0 ) {
                    final Writer writer = ForesterUtil.createBufferedWriter( base_dir + ForesterUtil.FILE_SEPARATOR
                            + node.getName() + suffix_for_filename );
                    for( final String domain : domains ) {
                        writer.write( domain );
                        writer.write( ForesterUtil.LINE_SEPARATOR );
                    }
                    writer.close();
                }
            }
        }
    }

    private static void writeAllEncounteredPfamsToFile( final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                        final Map<GoId, GoTerm> go_id_to_term_map,
                                                        final String outfile_name,
                                                        final SortedSet<String> all_pfams_encountered ) {
        final File all_pfams_encountered_file = new File( outfile_name + surfacing.ALL_PFAMS_ENCOUNTERED_SUFFIX );
        final File all_pfams_encountered_with_go_annotation_file = new File( outfile_name
                + surfacing.ALL_PFAMS_ENCOUNTERED_WITH_GO_ANNOTATION_SUFFIX );
        final File encountered_pfams_summary_file = new File( outfile_name + surfacing.ENCOUNTERED_PFAMS_SUMMARY_SUFFIX );
        int biological_process_counter = 0;
        int cellular_component_counter = 0;
        int molecular_function_counter = 0;
        int pfams_with_mappings_counter = 0;
        int pfams_without_mappings_counter = 0;
        int pfams_without_mappings_to_bp_or_mf_counter = 0;
        int pfams_with_mappings_to_bp_or_mf_counter = 0;
        try {
            final Writer all_pfams_encountered_writer = new BufferedWriter( new FileWriter( all_pfams_encountered_file ) );
            final Writer all_pfams_encountered_with_go_annotation_writer = new BufferedWriter( new FileWriter( all_pfams_encountered_with_go_annotation_file ) );
            final Writer summary_writer = new BufferedWriter( new FileWriter( encountered_pfams_summary_file ) );
            summary_writer.write( "# Pfam to GO mapping summary" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Actual summary is at the end of this file." );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Encountered Pfams without a GO mapping:" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            for( final String pfam : all_pfams_encountered ) {
                all_pfams_encountered_writer.write( pfam );
                all_pfams_encountered_writer.write( ForesterUtil.LINE_SEPARATOR );
                final DomainId domain_id = new DomainId( pfam );
                if ( domain_id_to_go_ids_map.containsKey( domain_id ) ) {
                    ++pfams_with_mappings_counter;
                    all_pfams_encountered_with_go_annotation_writer.write( pfam );
                    all_pfams_encountered_with_go_annotation_writer.write( ForesterUtil.LINE_SEPARATOR );
                    final List<GoId> go_ids = domain_id_to_go_ids_map.get( domain_id );
                    boolean maps_to_bp = false;
                    boolean maps_to_cc = false;
                    boolean maps_to_mf = false;
                    for( final GoId go_id : go_ids ) {
                        final GoTerm go_term = go_id_to_term_map.get( go_id );
                        if ( go_term.getGoNameSpace().isBiologicalProcess() ) {
                            maps_to_bp = true;
                        }
                        else if ( go_term.getGoNameSpace().isCellularComponent() ) {
                            maps_to_cc = true;
                        }
                        else if ( go_term.getGoNameSpace().isMolecularFunction() ) {
                            maps_to_mf = true;
                        }
                    }
                    if ( maps_to_bp ) {
                        ++biological_process_counter;
                    }
                    if ( maps_to_cc ) {
                        ++cellular_component_counter;
                    }
                    if ( maps_to_mf ) {
                        ++molecular_function_counter;
                    }
                    if ( maps_to_bp || maps_to_mf ) {
                        ++pfams_with_mappings_to_bp_or_mf_counter;
                    }
                    else {
                        ++pfams_without_mappings_to_bp_or_mf_counter;
                    }
                }
                else {
                    ++pfams_without_mappings_to_bp_or_mf_counter;
                    ++pfams_without_mappings_counter;
                    summary_writer.write( pfam );
                    summary_writer.write( ForesterUtil.LINE_SEPARATOR );
                }
            }
            all_pfams_encountered_writer.close();
            all_pfams_encountered_with_go_annotation_writer.close();
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote all [" + all_pfams_encountered.size()
                    + "] encountered Pfams to: \"" + all_pfams_encountered_file + "\"" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote all [" + pfams_with_mappings_counter
                    + "] encountered Pfams with GO mappings to: \"" + all_pfams_encountered_with_go_annotation_file
                    + "\"" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote summary (including all ["
                    + pfams_without_mappings_counter + "] encountered Pfams without GO mappings) to: \""
                    + encountered_pfams_summary_file + "\"" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Sum of Pfams encountered                : "
                    + all_pfams_encountered.size() );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Pfams without a mapping                 : "
                    + pfams_without_mappings_counter + " ["
                    + ( ( 100 * pfams_without_mappings_counter ) / all_pfams_encountered.size() ) + "%]" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Pfams without mapping to proc. or func. : "
                    + pfams_without_mappings_to_bp_or_mf_counter + " ["
                    + ( ( 100 * pfams_without_mappings_to_bp_or_mf_counter ) / all_pfams_encountered.size() ) + "%]" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Pfams with a mapping                    : "
                    + pfams_with_mappings_counter + " ["
                    + ( ( 100 * pfams_with_mappings_counter ) / all_pfams_encountered.size() ) + "%]" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Pfams with a mapping to proc. or func.  : "
                    + pfams_with_mappings_to_bp_or_mf_counter + " ["
                    + ( ( 100 * pfams_with_mappings_to_bp_or_mf_counter ) / all_pfams_encountered.size() ) + "%]" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Pfams with mapping to biological process: "
                    + biological_process_counter + " ["
                    + ( ( 100 * biological_process_counter ) / all_pfams_encountered.size() ) + "%]" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Pfams with mapping to molecular function: "
                    + molecular_function_counter + " ["
                    + ( ( 100 * molecular_function_counter ) / all_pfams_encountered.size() ) + "%]" );
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Pfams with mapping to cellular component: "
                    + cellular_component_counter + " ["
                    + ( ( 100 * cellular_component_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Sum of Pfams encountered                : " + all_pfams_encountered.size() );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Pfams without a mapping                 : " + pfams_without_mappings_counter
                    + " [" + ( ( 100 * pfams_without_mappings_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Pfams without mapping to proc. or func. : "
                    + pfams_without_mappings_to_bp_or_mf_counter + " ["
                    + ( ( 100 * pfams_without_mappings_to_bp_or_mf_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Pfams with a mapping                    : " + pfams_with_mappings_counter + " ["
                    + ( ( 100 * pfams_with_mappings_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Pfams with a mapping to proc. or func.  : "
                    + pfams_with_mappings_to_bp_or_mf_counter + " ["
                    + ( ( 100 * pfams_with_mappings_to_bp_or_mf_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Pfams with mapping to biological process: " + biological_process_counter + " ["
                    + ( ( 100 * biological_process_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Pfams with mapping to molecular function: " + molecular_function_counter + " ["
                    + ( ( 100 * molecular_function_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.write( "# Pfams with mapping to cellular component: " + cellular_component_counter + " ["
                    + ( ( 100 * cellular_component_counter ) / all_pfams_encountered.size() ) + "%]" );
            summary_writer.write( ForesterUtil.LINE_SEPARATOR );
            summary_writer.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.printWarningMessage( surfacing.PRG_NAME, "Failure to write: " + e );
        }
    }

    public static void writeBinaryDomainCombinationsFileForGraphAnalysis( final String[][] input_file_properties,
                                                                          final File output_dir,
                                                                          final GenomeWideCombinableDomains gwcd,
                                                                          final int i,
                                                                          final GenomeWideCombinableDomainsSortOrder dc_sort_order ) {
        File dc_outfile_dot = new File( input_file_properties[ i ][ 0 ]
                + surfacing.DOMAIN_COMBINITONS_OUTPUTFILE_SUFFIX_FOR_GRAPH_ANALYSIS );
        if ( output_dir != null ) {
            dc_outfile_dot = new File( output_dir + ForesterUtil.FILE_SEPARATOR + dc_outfile_dot );
        }
        checkForOutputFileWriteability( dc_outfile_dot );
        final SortedSet<BinaryDomainCombination> binary_combinations = createSetOfAllBinaryDomainCombinationsPerGenome( gwcd );
        try {
            final BufferedWriter out_dot = new BufferedWriter( new FileWriter( dc_outfile_dot ) );
            for( final BinaryDomainCombination bdc : binary_combinations ) {
                out_dot.write( bdc.toGraphDescribingLanguage( BinaryDomainCombination.OutputFormat.DOT, null, null )
                        .toString() );
                out_dot.write( SurfacingConstants.NL );
            }
            out_dot.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote binary domain combination for \""
                + input_file_properties[ i ][ 0 ] + "\" (" + input_file_properties[ i ][ 1 ] + ", "
                + input_file_properties[ i ][ 2 ] + ") to: \"" + dc_outfile_dot + "\"" );
    }

    public static void writeBinaryStatesMatrixAsListToFile( final CharacterStateMatrix<CharacterStateMatrix.GainLossStates> matrix,
                                                            final CharacterStateMatrix.GainLossStates state,
                                                            final String filename,
                                                            final String indentifier_characters_separator,
                                                            final String character_separator,
                                                            final Map<String, String> descriptions ) {
        final File outfile = new File( filename );
        checkForOutputFileWriteability( outfile );
        final SortedSet<String> sorted_ids = new TreeSet<String>();
        for( int i = 0; i < matrix.getNumberOfIdentifiers(); ++i ) {
            sorted_ids.add( matrix.getIdentifier( i ) );
        }
        try {
            final BufferedWriter out = new BufferedWriter( new FileWriter( outfile ) );
            for( final String id : sorted_ids ) {
                out.write( indentifier_characters_separator );
                out.write( "#" + id );
                out.write( indentifier_characters_separator );
                for( int c = 0; c < matrix.getNumberOfCharacters(); ++c ) {
                    // Not nice:
                    // using null to indicate either UNCHANGED_PRESENT or GAIN.
                    if ( ( matrix.getState( id, c ) == state )
                            || ( ( state == null ) && ( ( matrix.getState( id, c ) == CharacterStateMatrix.GainLossStates.GAIN ) || ( matrix
                                    .getState( id, c ) == CharacterStateMatrix.GainLossStates.UNCHANGED_PRESENT ) ) ) ) {
                        out.write( matrix.getCharacter( c ) );
                        if ( ( descriptions != null ) && !descriptions.isEmpty()
                                && descriptions.containsKey( matrix.getCharacter( c ) ) ) {
                            out.write( "\t" );
                            out.write( descriptions.get( matrix.getCharacter( c ) ) );
                        }
                        out.write( character_separator );
                    }
                }
            }
            out.flush();
            out.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote characters list: \"" + filename + "\"" );
    }

    public static void writeBinaryStatesMatrixAsListToFileForBinaryCombinationsForGraphAnalysis( final CharacterStateMatrix<CharacterStateMatrix.GainLossStates> matrix,
                                                                                                 final CharacterStateMatrix.GainLossStates state,
                                                                                                 final String filename,
                                                                                                 final String indentifier_characters_separator,
                                                                                                 final String character_separator,
                                                                                                 final BinaryDomainCombination.OutputFormat bc_output_format ) {
        final File outfile = new File( filename );
        checkForOutputFileWriteability( outfile );
        final SortedSet<String> sorted_ids = new TreeSet<String>();
        for( int i = 0; i < matrix.getNumberOfIdentifiers(); ++i ) {
            sorted_ids.add( matrix.getIdentifier( i ) );
        }
        try {
            final BufferedWriter out = new BufferedWriter( new FileWriter( outfile ) );
            for( final String id : sorted_ids ) {
                out.write( indentifier_characters_separator );
                out.write( "#" + id );
                out.write( indentifier_characters_separator );
                for( int c = 0; c < matrix.getNumberOfCharacters(); ++c ) {
                    // Not nice:
                    // using null to indicate either UNCHANGED_PRESENT or GAIN.
                    if ( ( matrix.getState( id, c ) == state )
                            || ( ( state == null ) && ( ( matrix.getState( id, c ) == CharacterStateMatrix.GainLossStates.GAIN ) || ( matrix
                                    .getState( id, c ) == CharacterStateMatrix.GainLossStates.UNCHANGED_PRESENT ) ) ) ) {
                        BinaryDomainCombination bdc = null;
                        try {
                            bdc = BasicBinaryDomainCombination.createInstance( matrix.getCharacter( c ) );
                        }
                        catch ( final Exception e ) {
                            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getLocalizedMessage() );
                        }
                        out.write( bdc.toGraphDescribingLanguage( bc_output_format, null, null ).toString() );
                        out.write( character_separator );
                    }
                }
            }
            out.flush();
            out.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote characters list: \"" + filename + "\"" );
    }

    public static void writeBinaryStatesMatrixToList( final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                                      final Map<GoId, GoTerm> go_id_to_term_map,
                                                      final GoNameSpace go_namespace_limit,
                                                      final boolean domain_combinations,
                                                      final CharacterStateMatrix<CharacterStateMatrix.GainLossStates> matrix,
                                                      final CharacterStateMatrix.GainLossStates state,
                                                      final String filename,
                                                      final String indentifier_characters_separator,
                                                      final String character_separator,
                                                      final String title_for_html,
                                                      final String prefix_for_html,
                                                      final Map<DomainId, Set<String>>[] domain_id_to_secondary_features_maps,
                                                      final SortedSet<String> all_pfams_encountered,
                                                      final SortedSet<String> pfams_gained_or_lost,
                                                      final String suffix_for_per_node_events_file ) {
        if ( ( go_namespace_limit != null ) && ( ( go_id_to_term_map == null ) || ( go_id_to_term_map.size() < 1 ) ) ) {
            throw new IllegalArgumentException( "attempt to use GO namespace limit without a GO-id to term map" );
        }
        else if ( ( ( domain_id_to_go_ids_map == null ) || ( domain_id_to_go_ids_map.size() < 1 ) ) ) {
            throw new IllegalArgumentException( "attempt to output detailed HTML without a Pfam to GO map" );
        }
        else if ( ( ( go_id_to_term_map == null ) || ( go_id_to_term_map.size() < 1 ) ) ) {
            throw new IllegalArgumentException( "attempt to output detailed HTML without a GO-id to term map" );
        }
        final File outfile = new File( filename );
        checkForOutputFileWriteability( outfile );
        final SortedSet<String> sorted_ids = new TreeSet<String>();
        for( int i = 0; i < matrix.getNumberOfIdentifiers(); ++i ) {
            sorted_ids.add( matrix.getIdentifier( i ) );
        }
        try {
            final Writer out = new BufferedWriter( new FileWriter( outfile ) );
            final File per_node_go_mapped_domain_gain_loss_files_base_dir = createBaseDirForPerNodeDomainFiles( surfacing.BASE_DIRECTORY_PER_NODE_DOMAIN_GAIN_LOSS_FILES,
                                                                                                                domain_combinations,
                                                                                                                state,
                                                                                                                filename );
            Writer per_node_go_mapped_domain_gain_loss_outfile_writer = null;
            File per_node_go_mapped_domain_gain_loss_outfile = null;
            int per_node_counter = 0;
            out.write( "<html>" );
            out.write( SurfacingConstants.NL );
            addHtmlHead( out, title_for_html );
            out.write( SurfacingConstants.NL );
            out.write( "<body>" );
            out.write( SurfacingConstants.NL );
            out.write( "<h1>" );
            out.write( SurfacingConstants.NL );
            out.write( title_for_html );
            out.write( SurfacingConstants.NL );
            out.write( "</h1>" );
            out.write( SurfacingConstants.NL );
            out.write( "<table>" );
            out.write( SurfacingConstants.NL );
            for( final String id : sorted_ids ) {
                final Matcher matcher = PATTERN_SP_STYLE_TAXONOMY.matcher( id );
                if ( matcher.matches() ) {
                    continue;
                }
                out.write( "<tr>" );
                out.write( "<td>" );
                out.write( "<a href=\"#" + id + "\">" + id + "</a>" );
                out.write( "</td>" );
                out.write( "</tr>" );
                out.write( SurfacingConstants.NL );
            }
            out.write( "</table>" );
            out.write( SurfacingConstants.NL );
            for( final String id : sorted_ids ) {
                final Matcher matcher = PATTERN_SP_STYLE_TAXONOMY.matcher( id );
                if ( matcher.matches() ) {
                    continue;
                }
                out.write( SurfacingConstants.NL );
                out.write( "<h2>" );
                out.write( "<a name=\"" + id + "\">" + id + "</a>" );
                writeTaxonomyLinks( out, id );
                out.write( "</h2>" );
                out.write( SurfacingConstants.NL );
                out.write( "<table>" );
                out.write( SurfacingConstants.NL );
                out.write( "<tr>" );
                out.write( "<td><b>" );
                out.write( "Pfam domain(s)" );
                out.write( "</b></td><td><b>" );
                out.write( "GO term acc" );
                out.write( "</b></td><td><b>" );
                out.write( "GO term" );
                out.write( "</b></td><td><b>" );
                out.write( "GO namespace" );
                out.write( "</b></td>" );
                out.write( "</tr>" );
                out.write( SurfacingConstants.NL );
                out.write( "</tr>" );
                out.write( SurfacingConstants.NL );
                per_node_counter = 0;
                if ( matrix.getNumberOfCharacters() > 0 ) {
                    per_node_go_mapped_domain_gain_loss_outfile = new File( per_node_go_mapped_domain_gain_loss_files_base_dir
                            + ForesterUtil.FILE_SEPARATOR + id + suffix_for_per_node_events_file );
                    SurfacingUtil.checkForOutputFileWriteability( per_node_go_mapped_domain_gain_loss_outfile );
                    per_node_go_mapped_domain_gain_loss_outfile_writer = ForesterUtil
                            .createBufferedWriter( per_node_go_mapped_domain_gain_loss_outfile );
                }
                else {
                    per_node_go_mapped_domain_gain_loss_outfile = null;
                    per_node_go_mapped_domain_gain_loss_outfile_writer = null;
                }
                for( int c = 0; c < matrix.getNumberOfCharacters(); ++c ) {
                    // Not nice:
                    // using null to indicate either UNCHANGED_PRESENT or GAIN.
                    if ( ( matrix.getState( id, c ) == state )
                            || ( ( state == null ) && ( ( matrix.getState( id, c ) == CharacterStateMatrix.GainLossStates.UNCHANGED_PRESENT ) || ( matrix
                                    .getState( id, c ) == CharacterStateMatrix.GainLossStates.GAIN ) ) ) ) {
                        final String character = matrix.getCharacter( c );
                        String domain_0 = "";
                        String domain_1 = "";
                        if ( character.indexOf( BinaryDomainCombination.SEPARATOR ) > 0 ) {
                            final String[] s = character.split( BinaryDomainCombination.SEPARATOR );
                            if ( s.length != 2 ) {
                                throw new AssertionError( "this should not have happened: unexpected format for domain combination: ["
                                        + character + "]" );
                            }
                            domain_0 = s[ 0 ];
                            domain_1 = s[ 1 ];
                        }
                        else {
                            domain_0 = character;
                        }
                        writeDomainData( domain_id_to_go_ids_map,
                                         go_id_to_term_map,
                                         go_namespace_limit,
                                         out,
                                         domain_0,
                                         domain_1,
                                         prefix_for_html,
                                         character_separator,
                                         domain_id_to_secondary_features_maps,
                                         null );
                        all_pfams_encountered.add( domain_0 );
                        if ( pfams_gained_or_lost != null ) {
                            pfams_gained_or_lost.add( domain_0 );
                        }
                        if ( !ForesterUtil.isEmpty( domain_1 ) ) {
                            all_pfams_encountered.add( domain_1 );
                            if ( pfams_gained_or_lost != null ) {
                                pfams_gained_or_lost.add( domain_1 );
                            }
                        }
                        if ( per_node_go_mapped_domain_gain_loss_outfile_writer != null ) {
                            writeDomainsToIndividualFilePerTreeNode( per_node_go_mapped_domain_gain_loss_outfile_writer,
                                                                     domain_0,
                                                                     domain_1 );
                            per_node_counter++;
                        }
                    }
                }
                if ( per_node_go_mapped_domain_gain_loss_outfile_writer != null ) {
                    per_node_go_mapped_domain_gain_loss_outfile_writer.close();
                    if ( per_node_counter < 1 ) {
                        per_node_go_mapped_domain_gain_loss_outfile.delete();
                    }
                    per_node_counter = 0;
                }
                out.write( "</table>" );
                out.write( SurfacingConstants.NL );
                out.write( "<hr>" );
                out.write( SurfacingConstants.NL );
            } // for( final String id : sorted_ids ) {  
            out.write( "</body>" );
            out.write( SurfacingConstants.NL );
            out.write( "</html>" );
            out.write( SurfacingConstants.NL );
            out.flush();
            out.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote characters detailed HTML list: \"" + filename + "\"" );
    }

    public static void writeDomainCombinationsCountsFile( final String[][] input_file_properties,
                                                          final File output_dir,
                                                          final Writer per_genome_domain_promiscuity_statistics_writer,
                                                          final GenomeWideCombinableDomains gwcd,
                                                          final int i,
                                                          final GenomeWideCombinableDomains.GenomeWideCombinableDomainsSortOrder dc_sort_order ) {
        File dc_outfile = new File( input_file_properties[ i ][ 0 ]
                + surfacing.DOMAIN_COMBINITON_COUNTS_OUTPUTFILE_SUFFIX );
        if ( output_dir != null ) {
            dc_outfile = new File( output_dir + ForesterUtil.FILE_SEPARATOR + dc_outfile );
        }
        checkForOutputFileWriteability( dc_outfile );
        try {
            final BufferedWriter out = new BufferedWriter( new FileWriter( dc_outfile ) );
            out.write( gwcd.toStringBuilder( dc_sort_order ).toString() );
            out.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        final DescriptiveStatistics stats = gwcd.getPerGenomeDomainPromiscuityStatistics();
        try {
            per_genome_domain_promiscuity_statistics_writer.write( input_file_properties[ i ][ 0 ] + "\t" );
            per_genome_domain_promiscuity_statistics_writer.write( FORMATTER_3.format( stats.arithmeticMean() ) + "\t" );
            if ( stats.getN() < 2 ) {
                per_genome_domain_promiscuity_statistics_writer.write( "n/a" + "\t" );
            }
            else {
                per_genome_domain_promiscuity_statistics_writer.write( FORMATTER_3.format( stats
                        .sampleStandardDeviation() ) + "\t" );
            }
            per_genome_domain_promiscuity_statistics_writer.write( FORMATTER_3.format( stats.median() ) + "\t" );
            per_genome_domain_promiscuity_statistics_writer.write( ( int ) stats.getMin() + "\t" );
            per_genome_domain_promiscuity_statistics_writer.write( ( int ) stats.getMax() + "\t" );
            per_genome_domain_promiscuity_statistics_writer.write( stats.getN() + "\t" );
            final SortedSet<DomainId> mpds = gwcd.getMostPromiscuosDomain();
            for( final DomainId mpd : mpds ) {
                per_genome_domain_promiscuity_statistics_writer.write( mpd.getId() + " " );
            }
            per_genome_domain_promiscuity_statistics_writer.write( ForesterUtil.LINE_SEPARATOR );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        if ( input_file_properties[ i ].length == 3 ) {
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote domain combination counts for \""
                    + input_file_properties[ i ][ 0 ] + "\" (" + input_file_properties[ i ][ 1 ] + ", "
                    + input_file_properties[ i ][ 2 ] + ") to: \"" + dc_outfile + "\"" );
        }
        else {
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote domain combination counts for \""
                    + input_file_properties[ i ][ 0 ] + "\" (" + input_file_properties[ i ][ 1 ] + ") to: \""
                    + dc_outfile + "\"" );
        }
    }

    private static void writeDomainData( final Map<DomainId, List<GoId>> domain_id_to_go_ids_map,
                                         final Map<GoId, GoTerm> go_id_to_term_map,
                                         final GoNameSpace go_namespace_limit,
                                         final Writer out,
                                         final String domain_0,
                                         final String domain_1,
                                         final String prefix_for_html,
                                         final String character_separator_for_non_html_output,
                                         final Map<DomainId, Set<String>>[] domain_id_to_secondary_features_maps,
                                         final Set<GoId> all_go_ids ) throws IOException {
        boolean any_go_annotation_present = false;
        boolean first_has_no_go = false;
        int domain_count = 2; // To distinguish between domains and binary domain combinations.
        if ( ForesterUtil.isEmpty( domain_1 ) ) {
            domain_count = 1;
        }
        // The following has a difficult to understand logic.  
        for( int d = 0; d < domain_count; ++d ) {
            List<GoId> go_ids = null;
            boolean go_annotation_present = false;
            if ( d == 0 ) {
                final DomainId domain_id = new DomainId( domain_0 );
                if ( domain_id_to_go_ids_map.containsKey( domain_id ) ) {
                    go_annotation_present = true;
                    any_go_annotation_present = true;
                    go_ids = domain_id_to_go_ids_map.get( domain_id );
                }
                else {
                    first_has_no_go = true;
                }
            }
            else {
                final DomainId domain_id = new DomainId( domain_1 );
                if ( domain_id_to_go_ids_map.containsKey( domain_id ) ) {
                    go_annotation_present = true;
                    any_go_annotation_present = true;
                    go_ids = domain_id_to_go_ids_map.get( domain_id );
                }
            }
            if ( go_annotation_present ) {
                boolean first = ( ( d == 0 ) || ( ( d == 1 ) && first_has_no_go ) );
                for( final GoId go_id : go_ids ) {
                    out.write( "<tr>" );
                    if ( first ) {
                        first = false;
                        writeDomainIdsToHtml( out,
                                              domain_0,
                                              domain_1,
                                              prefix_for_html,
                                              domain_id_to_secondary_features_maps );
                    }
                    else {
                        out.write( "<td></td>" );
                    }
                    if ( !go_id_to_term_map.containsKey( go_id ) ) {
                        throw new IllegalArgumentException( "GO-id [" + go_id + "] not found in GO-id to GO-term map" );
                    }
                    final GoTerm go_term = go_id_to_term_map.get( go_id );
                    if ( ( go_namespace_limit == null ) || go_namespace_limit.equals( go_term.getGoNameSpace() ) ) {
                        // final String top = GoUtils.getPenultimateGoTerm( go_term, go_id_to_term_map ).getName();
                        final String go_id_str = go_id.getId();
                        out.write( "<td>" );
                        out.write( "<a href=\"" + SurfacingConstants.AMIGO_LINK + go_id_str
                                + "\" target=\"amigo_window\">" + go_id_str + "</a>" );
                        out.write( "</td><td>" );
                        out.write( go_term.getName() );
                        if ( domain_count == 2 ) {
                            out.write( " (" + d + ")" );
                        }
                        out.write( "</td><td>" );
                        // out.write( top );
                        // out.write( "</td><td>" );
                        out.write( "[" );
                        out.write( go_term.getGoNameSpace().toShortString() );
                        out.write( "]" );
                        out.write( "</td>" );
                        if ( all_go_ids != null ) {
                            all_go_ids.add( go_id );
                        }
                    }
                    else {
                        out.write( "<td>" );
                        out.write( "</td><td>" );
                        out.write( "</td><td>" );
                        out.write( "</td><td>" );
                        out.write( "</td>" );
                    }
                    out.write( "</tr>" );
                    out.write( SurfacingConstants.NL );
                }
            }
        } //  for( int d = 0; d < domain_count; ++d ) 
        if ( !any_go_annotation_present ) {
            out.write( "<tr>" );
            writeDomainIdsToHtml( out, domain_0, domain_1, prefix_for_html, domain_id_to_secondary_features_maps );
            out.write( "<td>" );
            out.write( "</td><td>" );
            out.write( "</td><td>" );
            out.write( "</td><td>" );
            out.write( "</td>" );
            out.write( "</tr>" );
            out.write( SurfacingConstants.NL );
        }
    }

    private static void writeDomainIdsToHtml( final Writer out,
                                              final String domain_0,
                                              final String domain_1,
                                              final String prefix_for_detailed_html,
                                              final Map<DomainId, Set<String>>[] domain_id_to_secondary_features_maps )
            throws IOException {
        out.write( "<td>" );
        if ( !ForesterUtil.isEmpty( prefix_for_detailed_html ) ) {
            out.write( prefix_for_detailed_html );
            out.write( " " );
        }
        out.write( "<a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + domain_0 + "\">" + domain_0 + "</a>" );
        out.write( "</td>" );
    }

    public static DescriptiveStatistics writeDomainSimilaritiesToFile( final StringBuilder html_desc,
                                                                       final StringBuilder html_title,
                                                                       final Writer single_writer,
                                                                       Map<Character, Writer> split_writers,
                                                                       final SortedSet<DomainSimilarity> similarities,
                                                                       final boolean treat_as_binary,
                                                                       final List<Species> species_order,
                                                                       final PrintableDomainSimilarity.PRINT_OPTION print_option,
                                                                       final DomainSimilarity.DomainSimilaritySortField sort_field,
                                                                       final DomainSimilarity.DomainSimilarityScoring scoring,
                                                                       final boolean verbose ) throws IOException {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        String histogram_title = null;
        switch ( sort_field ) {
            case ABS_MAX_COUNTS_DIFFERENCE:
                if ( treat_as_binary ) {
                    histogram_title = "absolute counts difference:";
                }
                else {
                    histogram_title = "absolute (maximal) counts difference:";
                }
                break;
            case MAX_COUNTS_DIFFERENCE:
                if ( treat_as_binary ) {
                    histogram_title = "counts difference:";
                }
                else {
                    histogram_title = "(maximal) counts difference:";
                }
                break;
            case DOMAIN_ID:
                histogram_title = "score mean:";
                break;
            case MIN:
                histogram_title = "score minimum:";
                break;
            case MAX:
                histogram_title = "score maximum:";
                break;
            case MAX_DIFFERENCE:
                if ( treat_as_binary ) {
                    histogram_title = "difference:";
                }
                else {
                    histogram_title = "(maximal) difference:";
                }
                break;
            case MEAN:
                histogram_title = "score mean:";
                break;
            case SD:
                histogram_title = "score standard deviation:";
                break;
            case SPECIES_COUNT:
                histogram_title = "species number:";
                break;
            default:
                throw new AssertionError( "Unknown sort field: " + sort_field );
        }
        for( final DomainSimilarity similarity : similarities ) {
            switch ( sort_field ) {
                case ABS_MAX_COUNTS_DIFFERENCE:
                    stats.addValue( Math.abs( similarity.getMaximalDifferenceInCounts() ) );
                    break;
                case MAX_COUNTS_DIFFERENCE:
                    stats.addValue( similarity.getMaximalDifferenceInCounts() );
                    break;
                case DOMAIN_ID:
                    stats.addValue( similarity.getMeanSimilarityScore() );
                    break;
                case MIN:
                    stats.addValue( similarity.getMinimalSimilarityScore() );
                    break;
                case MAX:
                    stats.addValue( similarity.getMaximalSimilarityScore() );
                    break;
                case MAX_DIFFERENCE:
                    stats.addValue( similarity.getMaximalDifference() );
                    break;
                case MEAN:
                    stats.addValue( similarity.getMeanSimilarityScore() );
                    break;
                case SD:
                    stats.addValue( similarity.getStandardDeviationOfSimilarityScore() );
                    break;
                case SPECIES_COUNT:
                    stats.addValue( similarity.getSpecies().size() );
                    break;
                default:
                    throw new AssertionError( "Unknown sort field: " + sort_field );
            }
        }
        //
        // final HistogramData[] hists = new HistogramData[ 1 ];
        //      
        //        
        // List<HistogramDataItem> data_items = new
        // ArrayList<HistogramDataItem>();
        // double[] values = stats.getDataAsDoubleArray();
        // for( int i = 0; i < values.length; i++ ) {
        // HistogramDataItem data_item = new BasicHistogramDataItem( "", values[
        // i ] );
        // data_items.add( data_item );
        // }
        //        
        //        
        // HistogramData hd0 = new HistogramData( "name",
        // data_items,
        // null, 20,
        // 40 );
        //        
        //        
        //        
        //        
        // hists[ 0 ] = hd0;
        //       
        // final HistogramsFrame hf = new HistogramsFrame( hists );
        // hf.setVisible( true );
        //
        AsciiHistogram histo = null;
        if ( stats.getMin() < stats.getMin() ) {
            histo = new AsciiHistogram( stats, histogram_title );
        }
        if ( verbose ) {
            if ( histo != null ) {
                System.out.println( histo.toStringBuffer( 20, '|', 40, 5 ) );
            }
            System.out.println();
            System.out.println( "N                   : " + stats.getN() );
            System.out.println( "Min                 : " + stats.getMin() );
            System.out.println( "Max                 : " + stats.getMax() );
            System.out.println( "Mean                : " + stats.arithmeticMean() );
            if ( stats.getN() > 1 ) {
                System.out.println( "SD                  : " + stats.sampleStandardDeviation() );
            }
            else {
                System.out.println( "SD                  : n/a" );
            }
            System.out.println( "Median              : " + stats.median() );
            if ( stats.getN() > 1 ) {
                System.out.println( "Pearsonian skewness : " + stats.pearsonianSkewness() );
            }
            else {
                System.out.println( "Pearsonian skewness : n/a" );
            }
        }
        if ( ( single_writer != null ) && ( ( split_writers == null ) || split_writers.isEmpty() ) ) {
            split_writers = new HashMap<Character, Writer>();
            split_writers.put( '_', single_writer );
        }
        switch ( print_option ) {
            case SIMPLE_TAB_DELIMITED:
                break;
            case HTML:
                for( final Character key : split_writers.keySet() ) {
                    final Writer w = split_writers.get( key );
                    w.write( "<html>" );
                    w.write( SurfacingConstants.NL );
                    if ( key != '_' ) {
                        addHtmlHead( w, "DCs (" + html_title + ") " + key.toString().toUpperCase() );
                    }
                    else {
                        addHtmlHead( w, "DCs (" + html_title + ")" );
                    }
                    w.write( SurfacingConstants.NL );
                    w.write( "<body>" );
                    w.write( SurfacingConstants.NL );
                    w.write( html_desc.toString() );
                    w.write( SurfacingConstants.NL );
                    w.write( "<hr>" );
                    w.write( "<br>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<tt><pre>" );
                    w.write( SurfacingConstants.NL );
                    if ( histo != null ) {
                        w.write( histo.toStringBuffer( 20, '|', 40, 5 ).toString() );
                        w.write( SurfacingConstants.NL );
                    }
                    w.write( "</pre></tt>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<table>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<tr><td>N: </td><td>" + stats.getN() + "</td></tr>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<tr><td>Min: </td><td>" + stats.getMin() + "</td></tr>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<tr><td>Max: </td><td>" + stats.getMax() + "</td></tr>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<tr><td>Mean: </td><td>" + stats.arithmeticMean() + "</td></tr>" );
                    w.write( SurfacingConstants.NL );
                    if ( stats.getN() > 1 ) {
                        w.write( "<tr><td>SD: </td><td>" + stats.sampleStandardDeviation() + "</td></tr>" );
                    }
                    else {
                        w.write( "<tr><td>SD: </td><td>n/a</td></tr>" );
                    }
                    w.write( SurfacingConstants.NL );
                    w.write( "<tr><td>Median: </td><td>" + stats.median() + "</td></tr>" );
                    w.write( SurfacingConstants.NL );
                    if ( stats.getN() > 1 ) {
                        w.write( "<tr><td>Pearsonian skewness: </td><td>" + stats.pearsonianSkewness() + "</td></tr>" );
                    }
                    else {
                        w.write( "<tr><td>Pearsonian skewness: </td><td>n/a</td></tr>" );
                    }
                    w.write( SurfacingConstants.NL );
                    w.write( "</table>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<br>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<hr>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<br>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "<table>" );
                    w.write( SurfacingConstants.NL );
                }
                break;
        }
        for( final Writer w : split_writers.values() ) {
            w.write( SurfacingConstants.NL );
        }
        for( final DomainSimilarity similarity : similarities ) {
            if ( ( species_order != null ) && !species_order.isEmpty() ) {
                ( ( PrintableDomainSimilarity ) similarity ).setSpeciesOrder( species_order );
            }
            if ( single_writer != null ) {
                single_writer.write( similarity.toStringBuffer( print_option ).toString() );
            }
            else {
                Writer local_writer = split_writers.get( ( similarity.getDomainId().getId().charAt( 0 ) + "" )
                        .toLowerCase().charAt( 0 ) );
                if ( local_writer == null ) {
                    local_writer = split_writers.get( '0' );
                }
                local_writer.write( similarity.toStringBuffer( print_option ).toString() );
            }
            for( final Writer w : split_writers.values() ) {
                w.write( SurfacingConstants.NL );
            }
        }
        switch ( print_option ) {
            case HTML:
                for( final Writer w : split_writers.values() ) {
                    w.write( SurfacingConstants.NL );
                    w.write( "</table>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "</font>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "</body>" );
                    w.write( SurfacingConstants.NL );
                    w.write( "</html>" );
                    w.write( SurfacingConstants.NL );
                }
                break;
        }
        for( final Writer w : split_writers.values() ) {
            w.close();
        }
        return stats;
    }

    private static void writeDomainsToIndividualFilePerTreeNode( final Writer individual_files_writer,
                                                                 final String domain_0,
                                                                 final String domain_1 ) throws IOException {
        individual_files_writer.write( domain_0 );
        individual_files_writer.write( ForesterUtil.LINE_SEPARATOR );
        if ( !ForesterUtil.isEmpty( domain_1 ) ) {
            individual_files_writer.write( domain_1 );
            individual_files_writer.write( ForesterUtil.LINE_SEPARATOR );
        }
    }

    public static void writeMatrixToFile( final CharacterStateMatrix<?> matrix,
                                          final String filename,
                                          final Format format ) {
        final File outfile = new File( filename );
        checkForOutputFileWriteability( outfile );
        try {
            final BufferedWriter out = new BufferedWriter( new FileWriter( outfile ) );
            matrix.toWriter( out, format );
            out.flush();
            out.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote matrix: \"" + filename + "\"" );
    }

    public static void writeMatrixToFile( final File matrix_outfile, final List<DistanceMatrix> matrices ) {
        checkForOutputFileWriteability( matrix_outfile );
        try {
            final BufferedWriter out = new BufferedWriter( new FileWriter( matrix_outfile ) );
            for( final DistanceMatrix distance_matrix : matrices ) {
                out.write( distance_matrix.toStringBuffer( DistanceMatrix.Format.PHYLIP ).toString() );
                out.write( ForesterUtil.LINE_SEPARATOR );
                out.flush();
            }
            out.close();
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote distance matrices to \"" + matrix_outfile + "\"" );
    }

    private static void writePfamsToFile( final String outfile_name, final SortedSet<String> pfams ) {
        try {
            final Writer writer = new BufferedWriter( new FileWriter( new File( outfile_name ) ) );
            for( final String pfam : pfams ) {
                writer.write( pfam );
                writer.write( ForesterUtil.LINE_SEPARATOR );
            }
            writer.close();
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote " + pfams.size() + " pfams to [" + outfile_name
                    + "]" );
        }
        catch ( final IOException e ) {
            ForesterUtil.printWarningMessage( surfacing.PRG_NAME, "Failure to write: " + e );
        }
    }

    public static void writePhylogenyToFile( final Phylogeny phylogeny, final String filename ) {
        final PhylogenyWriter writer = new PhylogenyWriter();
        try {
            writer.toPhyloXML( new File( filename ), phylogeny, 1 );
        }
        catch ( final IOException e ) {
            ForesterUtil.printWarningMessage( surfacing.PRG_NAME, "failed to write phylogeny to \"" + filename + "\": "
                    + e );
        }
        ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote phylogeny to \"" + filename + "\"" );
    }

    public static void writeTaxonomyLinks( final Writer writer, final String species ) throws IOException {
        if ( ( species.length() > 1 ) && ( species.indexOf( '_' ) < 1 ) ) {
            final Matcher matcher = PATTERN_SP_STYLE_TAXONOMY.matcher( species );
            writer.write( " [" );
            if ( matcher.matches() ) {
                writer.write( "<a href=\"" + SurfacingConstants.UNIPROT_LINK + species
                        + "\" target=\"taxonomy_window\">uniprot</a>" );
            }
            else {
                writer.write( "<a href=\"" + SurfacingConstants.EOL_LINK + species
                        + "\" target=\"taxonomy_window\">eol</a>" );
                writer.write( "|" );
                writer.write( "<a href=\"" + SurfacingConstants.TOL_LINK + species
                        + "\" target=\"taxonomy_window\">tol</a>" );
            }
            writer.write( "]" );
        }
    }

    private static void writeToNexus( final String outfile_name,
                                      final CharacterStateMatrix<BinaryStates> matrix,
                                      final Phylogeny phylogeny ) {
        if ( !( matrix instanceof BasicCharacterStateMatrix ) ) {
            throw new IllegalArgumentException( "can only write matrices of type [" + BasicCharacterStateMatrix.class
                    + "] to nexus" );
        }
        final BasicCharacterStateMatrix<BinaryStates> my_matrix = ( org.forester.evoinference.matrix.character.BasicCharacterStateMatrix<BinaryStates> ) matrix;
        final List<Phylogeny> phylogenies = new ArrayList<Phylogeny>( 1 );
        phylogenies.add( phylogeny );
        try {
            final BufferedWriter w = new BufferedWriter( new FileWriter( outfile_name ) );
            w.write( NexusConstants.NEXUS );
            w.write( ForesterUtil.LINE_SEPARATOR );
            my_matrix.writeNexusTaxaBlock( w );
            my_matrix.writeNexusBinaryChractersBlock( w );
            PhylogenyWriter.writeNexusTreesBlock( w, phylogenies, NH_CONVERSION_SUPPORT_VALUE_STYLE.NONE );
            w.flush();
            w.close();
            ForesterUtil.programMessage( surfacing.PRG_NAME, "Wrote Nexus file: \"" + outfile_name + "\"" );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( surfacing.PRG_NAME, e.getMessage() );
        }
    }

    private static void writeToNexus( final String outfile_name,
                                      final DomainParsimonyCalculator domain_parsimony,
                                      final Phylogeny phylogeny ) {
        writeToNexus( outfile_name + surfacing.NEXUS_EXTERNAL_DOMAINS,
                      domain_parsimony.createMatrixOfDomainPresenceOrAbsence(),
                      phylogeny );
        writeToNexus( outfile_name + surfacing.NEXUS_EXTERNAL_DOMAIN_COMBINATIONS,
                      domain_parsimony.createMatrixOfBinaryDomainCombinationPresenceOrAbsence(),
                      phylogeny );
    }

    public static void domainsPerProteinsStatistics( final String genome,
                                                     final List<Protein> protein_list,
                                                     final DescriptiveStatistics all_genomes_domains_per_potein_stats,
                                                     final SortedMap<Integer, Integer> all_genomes_domains_per_potein_histo,
                                                     final SortedSet<String> domains_which_are_always_single,
                                                     final SortedSet<String> domains_which_are_sometimes_single_sometimes_not,
                                                     final SortedSet<String> domains_which_never_single,
                                                     final Writer writer ) {
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        for( final Protein protein : protein_list ) {
            final int domains = protein.getNumberOfProteinDomains();
            //System.out.println( domains );
            stats.addValue( domains );
            all_genomes_domains_per_potein_stats.addValue( domains );
            if ( !all_genomes_domains_per_potein_histo.containsKey( domains ) ) {
                all_genomes_domains_per_potein_histo.put( domains, 1 );
            }
            else {
                all_genomes_domains_per_potein_histo.put( domains,
                                                          1 + all_genomes_domains_per_potein_histo.get( domains ) );
            }
            if ( domains == 1 ) {
                final String domain = protein.getProteinDomain( 0 ).getDomainId().getId();
                if ( !domains_which_are_sometimes_single_sometimes_not.contains( domain ) ) {
                    if ( domains_which_never_single.contains( domain ) ) {
                        domains_which_never_single.remove( domain );
                        domains_which_are_sometimes_single_sometimes_not.add( domain );
                    }
                    else {
                        domains_which_are_always_single.add( domain );
                    }
                }
            }
            else if ( domains > 1 ) {
                for( final Domain d : protein.getProteinDomains() ) {
                    final String domain = d.getDomainId().getId();
                    // System.out.println( domain );
                    if ( !domains_which_are_sometimes_single_sometimes_not.contains( domain ) ) {
                        if ( domains_which_are_always_single.contains( domain ) ) {
                            domains_which_are_always_single.remove( domain );
                            domains_which_are_sometimes_single_sometimes_not.add( domain );
                        }
                        else {
                            domains_which_never_single.add( domain );
                        }
                    }
                }
            }
        }
        try {
            writer.write( genome );
            writer.write( "\t" );
            if ( stats.getN() >= 1 ) {
                writer.write( stats.arithmeticMean() + "" );
                writer.write( "\t" );
                if ( stats.getN() >= 2 ) {
                    writer.write( stats.sampleStandardDeviation() + "" );
                }
                else {
                    writer.write( "" );
                }
                writer.write( "\t" );
                writer.write( stats.median() + "" );
                writer.write( "\t" );
                writer.write( stats.getN() + "" );
                writer.write( "\t" );
                writer.write( stats.getMin() + "" );
                writer.write( "\t" );
                writer.write( stats.getMax() + "" );
            }
            else {
                writer.write( "\t" );
                writer.write( "\t" );
                writer.write( "\t" );
                writer.write( "0" );
                writer.write( "\t" );
                writer.write( "\t" );
            }
            writer.write( "\n" );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }

    final static class DomainComparator implements Comparator<Domain> {

        final private boolean _ascending;

        public DomainComparator( final boolean ascending ) {
            _ascending = ascending;
        }

        @Override
        public final int compare( final Domain d0, final Domain d1 ) {
            if ( d0.getFrom() < d1.getFrom() ) {
                return _ascending ? -1 : 1;
            }
            else if ( d0.getFrom() > d1.getFrom() ) {
                return _ascending ? 1 : -1;
            }
            return 0;
        }
    }
}

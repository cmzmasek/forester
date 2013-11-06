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

import java.awt.Color;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.species.Species;
import org.forester.surfacing.DomainSimilarityCalculator.Detailedness;
import org.forester.util.ForesterUtil;

public class PrintableDomainSimilarity implements DomainSimilarity {

    final public static String                              SPECIES_SEPARATOR = "  ";
    final private static int                                EQUAL             = 0;
    final private static String                             NO_SPECIES        = "     ";
    final private CombinableDomains                         _combinable_domains;
    private DomainSimilarityCalculator.Detailedness         _detailedness;
    final private double                                    _max;
    private final int                                       _max_difference;
    private final int                                       _max_difference_in_counts;
    final private double                                    _mean;
    final private double                                    _min;
    final private int                                       _n;
    final private double                                    _sd;
    final private SortedMap<Species, SpeciesSpecificDcData> _species_data;
    private List<Species>                                   _species_order;
    private final boolean                                   _treat_as_binary_comparison;

    public PrintableDomainSimilarity( final CombinableDomains combinable_domains,
                                      final double min,
                                      final double max,
                                      final double mean,
                                      final double median,
                                      final double sd,
                                      final int n,
                                      final int max_difference_in_counts,
                                      final int max_difference,
                                      final SortedMap<Species, SpeciesSpecificDcData> species_data,
                                      final boolean sort_by_species_count_first,
                                      final boolean treat_as_binary_comparison ) {
        if ( combinable_domains == null ) {
            throw new IllegalArgumentException( "attempt to use null combinable domains" );
        }
        if ( species_data == null ) {
            throw new IllegalArgumentException( "attempt to use null species data" );
        }
        if ( species_data.size() < 1 ) {
            throw new IllegalArgumentException( "attempt to use empty species data" );
        }
        if ( n < 0 ) {
            throw new IllegalArgumentException( "attempt to use N less than 0" );
        }
        if ( ( species_data.size() > 1 ) && ( n < 1 ) ) {
            throw new IllegalArgumentException( "attempt to use N less than 1" );
        }
        if ( sd < 0.0 ) {
            throw new IllegalArgumentException( "attempt to use negative SD" );
        }
        if ( max < min ) {
            throw new IllegalArgumentException( "attempt to use max smaller than min" );
        }
        init();
        _combinable_domains = combinable_domains;
        _min = min;
        _max = max;
        _mean = mean;
        _sd = sd;
        _n = n;
        _max_difference_in_counts = max_difference_in_counts;
        _max_difference = max_difference;
        _species_data = species_data;
        _treat_as_binary_comparison = treat_as_binary_comparison;
        final int s = species_data.size();
        if ( ( ( s * s ) - s ) != ( getN() * 2 ) ) {
            throw new IllegalArgumentException( "illegal species count and n: species count:" + s + ", n:" + _n
                    + " for domain " + combinable_domains.getKeyDomain() );
        }
        if ( s > 2 ) {
            if ( getMaximalDifferenceInCounts() < 0 ) {
                throw new IllegalArgumentException( "attempt to use negative max difference in counts with more than two species" );
            }
            if ( getMaximalDifference() < 0 ) {
                throw new IllegalArgumentException( "attempt to use negative max difference with more than two species" );
            }
        }
    }

    public PrintableDomainSimilarity( final CombinableDomains combinable_domains,
                                      final int max_difference_in_counts,
                                      final int max_difference,
                                      final SortedMap<Species, SpeciesSpecificDcData> species_data,
                                      final boolean sort_by_species_count_first,
                                      final boolean treat_as_binary_comparison ) {
        if ( combinable_domains == null ) {
            throw new IllegalArgumentException( "attempt to use null combinable domains" );
        }
        if ( species_data == null ) {
            throw new IllegalArgumentException( "attempt to use null species data" );
        }
        if ( species_data.size() < 1 ) {
            throw new IllegalArgumentException( "attempt to use empty species data" );
        }
        init();
        _combinable_domains = combinable_domains;
        _min = -1;
        _max = -1;
        _mean = -1;
        _sd = -1;
        _n = -1;
        _max_difference_in_counts = max_difference_in_counts;
        _max_difference = max_difference;
        _species_data = species_data;
        _treat_as_binary_comparison = treat_as_binary_comparison;
        final int s = species_data.size();
        if ( s > 2 ) {
            if ( getMaximalDifferenceInCounts() < 0 ) {
                throw new IllegalArgumentException( "attempt to use negative max difference in counts with more than two species" );
            }
            if ( getMaximalDifference() < 0 ) {
                throw new IllegalArgumentException( "attempt to use negative max difference with more than two species" );
            }
        }
    }

    @Override
    public int compareTo( final DomainSimilarity domain_similarity ) {
        if ( this == domain_similarity ) {
            return EQUAL;
        }
        else if ( domain_similarity == null ) {
            throw new IllegalArgumentException( "attempt to compare " + this.getClass() + " to null" );
        }
        else if ( domain_similarity.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to compare " + this.getClass() + " to "
                    + domain_similarity.getClass() );
        }
        return compareByDomainId( domain_similarity );
    }

    @Override
    public SortedSet<String> getCombinableDomainIds( final Species species_of_combinable_domain ) {
        final SortedSet<String> sorted_ids = new TreeSet<String>();
        if ( getSpeciesData().containsKey( species_of_combinable_domain ) ) {
            for( final String id : getSpeciesData().get( species_of_combinable_domain )
                    .getCombinableDomainIdToCountsMap().keySet() ) {
                sorted_ids.add( id );
            }
        }
        return sorted_ids;
    }

    @Override
    public String getDomainId() {
        return getCombinableDomains().getKeyDomain();
    }

    @Override
    public int getMaximalDifference() {
        return _max_difference;
    }

    @Override
    public int getMaximalDifferenceInCounts() {
        return _max_difference_in_counts;
    }

    @Override
    public double getMaximalSimilarityScore() {
        return _max;
    }

    @Override
    public double getMeanSimilarityScore() {
        return _mean;
    }

    @Override
    public double getMinimalSimilarityScore() {
        return _min;
    }

    @Override
    public int getN() {
        return _n;
    }

    @Override
    public SortedSet<Species> getSpecies() {
        final SortedSet<Species> species = new TreeSet<Species>();
        for( final Species s : getSpeciesData().keySet() ) {
            species.add( s );
        }
        return species;
    }

    public List<Species> getSpeciesCustomOrder() {
        return _species_order;
    }

    @Override
    public SortedMap<Species, SpeciesSpecificDcData> getSpeciesData() {
        return _species_data;
    }

    @Override
    public double getStandardDeviationOfSimilarityScore() {
        return _sd;
    }

    public void setDetailedness( final Detailedness detailedness ) {
        _detailedness = detailedness;
    }

    public void setSpeciesOrder( final List<Species> species_order ) {
        if ( !species_order.containsAll( getSpeciesData().keySet() ) ) {
            throw new IllegalArgumentException( "list to order species must contain all species of multiple combinable domains similarity" );
        }
        _species_order = species_order;
    }

    @Override
    public StringBuffer toStringBuffer( final PrintableDomainSimilarity.PRINT_OPTION print_option,
                                        final Map<String, Integer> tax_code_to_id_map,
                                        final Phylogeny phy ) {
        switch ( print_option ) {
            case SIMPLE_TAB_DELIMITED:
                return toStringBufferSimpleTabDelimited();
            case HTML:
                return toStringBufferDetailedHTML( tax_code_to_id_map, phy );
            default:
                throw new AssertionError( "Unknown print option: " + print_option );
        }
    }

    private void addSpeciesSpecificDomainData( final StringBuffer sb,
                                               final Species species,
                                               final boolean html,
                                               final Map<String, Integer> tax_code_to_id_map,
                                               final Phylogeny phy ) {
        if ( html ) {
            addTaxWithLink( sb, species.getSpeciesId(), tax_code_to_id_map, phy );
        }
        else {
            sb.append( species.getSpeciesId() );
        }
        if ( getDetaildness() != DomainSimilarityCalculator.Detailedness.BASIC ) {
            if ( html ) {
                sb.append( ":" );
            }
            else {
                sb.append( "\t" );
            }
            sb.append( getSpeciesData().get( species ).toStringBuffer( getDetaildness(), html ) );
        }
        if ( html ) {
            sb.append( "<br>" );
        }
        else {
            sb.append( "\n\t" );
        }
    }

    private void addTaxWithLink( final StringBuffer sb,
                                 final String tax_code,
                                 final Map<String, Integer> tax_code_to_id_map,
                                 final Phylogeny phy ) {
        String hex = null;
        if ( ( phy != null ) && !phy.isEmpty() ) {
            hex = SurfacingUtil.obtainHexColorStringDependingOnTaxonomyGroup( tax_code, phy );
        }
        sb.append( "<b>" );
        if ( !ForesterUtil.isEmpty( tax_code )
                && ( ( tax_code_to_id_map != null ) && tax_code_to_id_map.containsKey( tax_code ) ) ) {
            if ( !ForesterUtil.isEmpty( hex ) ) {
                sb.append( "<a href=\"" );
                sb.append( SurfacingConstants.UNIPROT_TAXONOMY_ID_LINK );
                sb.append( tax_code_to_id_map.get( tax_code ) );
                sb.append( "\" target=\"tw\"><span style=\"color:" );
                sb.append( hex );
                sb.append( "\">" );
                sb.append( tax_code );
                sb.append( "</span></a>" );
            }
            else {
                sb.append( "<a href=\"" );
                sb.append( SurfacingConstants.UNIPROT_TAXONOMY_ID_LINK );
                sb.append( tax_code_to_id_map.get( tax_code ) );
                sb.append( "\" target=\"tw\">" );
                sb.append( tax_code );
                sb.append( "</a>" );
            }
        }
        else {
            sb.append( tax_code );
        }
        sb.append( "</b>" );
    }

    private int compareByDomainId( final DomainSimilarity other ) {
        return getDomainId().compareToIgnoreCase( other.getDomainId() );
    }

    private CombinableDomains getCombinableDomains() {
        return _combinable_domains;
    }

    private DomainSimilarityCalculator.Detailedness getDetaildness() {
        return _detailedness;
    }

    private StringBuffer getDomainDataInAlphabeticalOrder() {
        final SortedMap<String, SortedSet<String>> m = new TreeMap<String, SortedSet<String>>();
        final StringBuffer sb = new StringBuffer();
        for( final Species species : getSpeciesData().keySet() ) {
            for( final String combable_dom : getCombinableDomainIds( species ) ) {
                if ( !m.containsKey( combable_dom ) ) {
                    m.put( combable_dom, new TreeSet<String>() );
                }
                m.get( combable_dom ).add( species.getSpeciesId() );
            }
        }
        for( final Map.Entry<String, SortedSet<String>> e : m.entrySet() ) {
            sb.append( "<a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + e.getKey() + "\">" + e.getKey() + "</a>" );
            sb.append( ": " );
            sb.append( "<span style=\"font-size:7px\">" );
            for( final String tax : e.getValue() ) {
                final String hex = SurfacingUtil.obtainHexColorStringDependingOnTaxonomyGroup( tax, null );
                if ( !ForesterUtil.isEmpty( hex ) ) {
                    sb.append( "<span style=\"color:" );
                    sb.append( hex );
                    sb.append( "\">" );
                    sb.append( tax );
                    sb.append( "</span>" );
                }
                else {
                    sb.append( tax );
                }
                sb.append( " " );
            }
            sb.append( "</span>" );
            sb.append( "<br>\n" );
        }
        return sb;
    }

    private StringBuffer getTaxonomyGroupDistribution( final Phylogeny tol ) {
        //TODO work on me    
        final SortedMap<String, Set<String>> domain_to_species_set_map = new TreeMap<String, Set<String>>();
        for( final Species species : getSpeciesData().keySet() ) {
            for( final String combable_dom : getCombinableDomainIds( species ) ) {
                if ( !domain_to_species_set_map.containsKey( combable_dom ) ) {
                    domain_to_species_set_map.put( combable_dom, new HashSet<String>() );
                }
                domain_to_species_set_map.get( combable_dom ).add( species.getSpeciesId() );
            }
        }
        final StringBuffer sb = new StringBuffer();
        sb.append( "<table>" );
        for( final Map.Entry<String, Set<String>> domain_to_species_set : domain_to_species_set_map.entrySet() ) {
            final Map<String, Integer> counts = new HashMap<String, Integer>();
            //  final ValueComparator bvc = new ValueComparator( counts );
            //  final SortedMap<String, Integer> sorted_counts = new TreeMap<String, Integer>( bvc );
            for( final String tax_code : domain_to_species_set.getValue() ) {
                final String group = SurfacingUtil.obtainTaxonomyGroup( tax_code, tol );
                if ( !ForesterUtil.isEmpty( group ) ) {
                    if ( !counts.containsKey( group ) ) {
                        counts.put( group, 1 );
                    }
                    else {
                        counts.put( group, counts.get( group ) + 1 );
                    }
                }
                else {
                    return null;
                }
            }
            final SortedMap<Integer, SortedSet<String>> counts_to_groups = new TreeMap<Integer, SortedSet<String>>( new Comparator<Integer>() {

                @Override
                public int compare( final Integer first, final Integer second ) {
                    return second.compareTo( first );
                }
            } );
            for( final Map.Entry<String, Integer> group_to_counts : counts.entrySet() ) {
                final int c = group_to_counts.getValue();
                if ( !counts_to_groups.containsKey( c ) ) {
                    counts_to_groups.put( c, new TreeSet<String>() );
                }
                counts_to_groups.get( c ).add( group_to_counts.getKey() );
            }
            // sorted_counts.putAll( counts );
            sb.append( "<tr>" );
            sb.append( "<td>" );
            sb.append( "<a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + domain_to_species_set.getKey() + "\">"
                    + domain_to_species_set.getKey() + "</a>" );
            sb.append( ": " );
            sb.append( "</td>" );
            // sb.append( "<span style=\"font-size:9px\">" );
            boolean first = true;
            for( final Entry<Integer, SortedSet<String>> count_to_groups : counts_to_groups.entrySet() ) {
                if ( first ) {
                    first = false;
                }
                else {
                    sb.append( "<tr>" );
                    sb.append( "<td>" );
                    sb.append( "</td>" );
                }
                sb.append( "<td>" );
                final SortedSet<String> groups = count_to_groups.getValue();
                sb.append( count_to_groups.getKey() );
                sb.append( ":" );
                for( final String group : groups ) {
                    final Color color = ForesterUtil.obtainColorDependingOnTaxonomyGroup( group );
                    if ( color == null ) {
                        throw new IllegalArgumentException( "no color found for taxonomy group\"" + group + "\"" );
                    }
                    final String hex = String.format( "#%02x%02x%02x",
                                                      color.getRed(),
                                                      color.getGreen(),
                                                      color.getBlue() );
                    sb.append( "<span style=\"color:" );
                    sb.append( hex );
                    sb.append( "\">" );
                    sb.append( " " );
                    sb.append( group );
                    sb.append( "</span>" );
                }
                sb.append( "</td>" );
                sb.append( "</tr>" );
            }
            // sb.append( "</span>" );
            sb.append( ForesterUtil.getLineSeparator() );
        }
        sb.append( "</table>" );
        // i am just a template and need to be modified for "printout" TODO
        //        for( final Map.Entry<String, SortedSet<String>> e : m.entrySet() ) {
        //            sb.append( "<a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + e.getKey() + "\">" + e.getKey() + "</a>" );
        //            sb.append( ": " );
        //            sb.append( "<span style=\"font-size:8px\">" );
        //            for( final String tax : e.getValue() ) {
        //                final String hex = SurfacingUtil.obtainHexColorStringDependingOnTaxonomyGroup( tax, null );
        //                if ( !ForesterUtil.isEmpty( hex ) ) {
        //                    sb.append( "<span style=\"color:" );
        //                    sb.append( hex );
        //                    sb.append( "\">" );
        //                    sb.append( tax );
        //                    sb.append( "</span>" );
        //                }
        //                else {
        //                    sb.append( tax );
        //                }
        //                sb.append( " " );
        //            }
        //            sb.append( "</span>" );
        //            sb.append( "<br>\n" );
        //        }
        return sb;
    }

    /*
     public class Testing {

    public static void main(String[] args) {

        HashMap<String,Double> map = new HashMap<String,Double>();
        ValueComparator bvc =  new ValueComparator(map);
        TreeMap<String,Double> sorted_map = new TreeMap<String,Double>(bvc);

        map.put("A",99.5);
        map.put("B",67.4);
        map.put("C",67.4);
        map.put("D",67.3);

        System.out.println("unsorted map: "+map);

        sorted_map.putAll(map);

        System.out.println("results: "+sorted_map);
    }
    }

       
      
     */
    private StringBuffer getSpeciesDataInAlphabeticalOrder( final boolean html,
                                                            final Map<String, Integer> tax_code_to_id_map,
                                                            final Phylogeny phy ) {
        final StringBuffer sb = new StringBuffer();
        for( final Species species : getSpeciesData().keySet() ) {
            addSpeciesSpecificDomainData( sb, species, html, tax_code_to_id_map, phy );
        }
        return sb;
    }

    private StringBuffer getSpeciesDataInCustomOrder( final boolean html,
                                                      final Map<String, Integer> tax_code_to_id_map,
                                                      final Phylogeny phy ) {
        final StringBuffer sb = new StringBuffer();
        for( final Species order_species : getSpeciesCustomOrder() ) {
            if ( getSpeciesData().keySet().contains( order_species ) ) {
                addSpeciesSpecificDomainData( sb, order_species, html, tax_code_to_id_map, phy );
            }
            else {
                sb.append( PrintableDomainSimilarity.NO_SPECIES );
                sb.append( PrintableDomainSimilarity.SPECIES_SEPARATOR );
            }
        }
        return sb;
    }

    private void init() {
        _detailedness = DomainSimilarityCalculator.Detailedness.PUNCTILIOUS;
    }

    private boolean isTreatAsBinaryComparison() {
        return _treat_as_binary_comparison;
    }

    private StringBuffer toStringBufferDetailedHTML( final Map<String, Integer> tax_code_to_id_map, final Phylogeny phy ) {
        final StringBuffer sb = new StringBuffer();
        sb.append( "<tr>" );
        sb.append( "<td>" );
        sb.append( "<b>" );
        sb.append( "<a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + getDomainId() + "\" target=\"pfam_window\">"
                + getDomainId() + "</a>" );
        sb.append( "</b>" );
        sb.append( "<a name=\"" + getDomainId() + "\">" );
        sb.append( "</td>" );
        sb.append( "<td>" );
        sb.append( "<a href=\"" + SurfacingConstants.GOOGLE_SCHOLAR_SEARCH + getDomainId()
                + "\" target=\"gs_window\">gs</a>" );
        sb.append( "</td>" );
        if ( getMaximalSimilarityScore() > 0 ) {
            sb.append( "<td>" );
            sb.append( ForesterUtil.round( getMeanSimilarityScore(), 3 ) );
            sb.append( "</td>" );
            if ( SurfacingConstants.PRINT_MORE_DOM_SIMILARITY_INFO ) {
                if ( !isTreatAsBinaryComparison() ) {
                    sb.append( "<td>" );
                    sb.append( "(" );
                    sb.append( ForesterUtil.round( getStandardDeviationOfSimilarityScore(), 3 ) );
                    sb.append( ")" );
                    sb.append( "</td>" );
                    sb.append( "<td>" );
                    sb.append( "[" );
                    sb.append( ForesterUtil.round( getMinimalSimilarityScore(), 3 ) );
                    sb.append( "-" );
                    sb.append( ForesterUtil.round( getMaximalSimilarityScore(), 3 ) );
                    sb.append( "]" );
                    sb.append( "</td>" );
                }
            }
        }
        sb.append( "<td>" );
        sb.append( getMaximalDifference() );
        sb.append( "</td>" );
        sb.append( "<td>" );
        if ( isTreatAsBinaryComparison() ) {
            sb.append( getMaximalDifferenceInCounts() );
        }
        else {
            sb.append( Math.abs( getMaximalDifferenceInCounts() ) );
        }
        sb.append( "</td>" );
        if ( !isTreatAsBinaryComparison() ) {
            sb.append( "<td>" );
            sb.append( "<b>" );
            sb.append( getSpeciesData().size() );
            sb.append( "</b>" );
            sb.append( "</td>" );
        }
        if ( ( getSpeciesCustomOrder() == null ) || getSpeciesCustomOrder().isEmpty() ) {
            sb.append( "<td>" );
            sb.append( getSpeciesDataInAlphabeticalOrder( true, tax_code_to_id_map, phy ) );
            sb.append( getDomainDataInAlphabeticalOrder() );
            sb.append( getTaxonomyGroupDistribution( phy ) );
            sb.append( "</td>" );
        }
        else {
            sb.append( "<td>" );
            sb.append( getSpeciesDataInCustomOrder( true, tax_code_to_id_map, phy ) );
            sb.append( getDomainDataInAlphabeticalOrder() );
            sb.append( getTaxonomyGroupDistribution( phy ) );
            sb.append( "</td>" );
        }
        sb.append( "</tr>" );
        return sb;
    }

    private StringBuffer toStringBufferSimpleTabDelimited() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getDomainId() );
        sb.append( "\t" );
        sb.append( getSpeciesDataInAlphabeticalOrder( false, null, null ) );
        sb.append( "\n" );
        return sb;
    }

    public static enum PRINT_OPTION {
        HTML, SIMPLE_TAB_DELIMITED;
    }

    class ValueComparator implements Comparator<String> {

        final private Map<String, Integer> _base;

        public ValueComparator( final Map<String, Integer> base ) {
            _base = base;
        }

        public int compare( final String a, final String b ) {
            if ( _base.get( a ) >= _base.get( b ) ) {
                return -1;
            }
            else {
                return 1;
            } // returning 0 would merge keys
        }
    }
}

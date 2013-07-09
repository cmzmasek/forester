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

import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.species.Species;
import org.forester.surfacing.DomainSimilarityCalculator.Detailedness;
import org.forester.util.ForesterUtil;

public class PrintableDomainSimilarity implements DomainSimilarity {

    final public static String                                           SPECIES_SEPARATOR = "  ";
    final private static char                                            TAB               = '\t';
    final private static int                                             BEFORE            = -1;
    final private static int                                             EQUAL             = 0;
    final private static int                                             AFTER             = 1;
    final private static String                                          NO_SPECIES        = "     ";
    final private double                                                 _min;
    final private double                                                 _max;
    final private double                                                 _mean;
    final private double                                                 _sd;
    final private int                                                    _n;
    private final int                                                    _max_difference_in_counts;
    private final int                                                    _max_difference;
    final private CombinableDomains                                      _combinable_domains;
    final private SortedMap<Species, SpeciesSpecificDomainSimilariyData> _species_data;
    final private DomainSimilaritySortField                              _sort_field;
    private List<Species>                                                _species_order;
    private final boolean                                                _sort_by_species_count_first;
    private DomainSimilarityCalculator.Detailedness                      _detailedness;
    private final boolean                                                _treat_as_binary_comparison;

    /**
     * If go_id_to_term_map not null, detailed GO information is written,
     * only GO ids otherwise.
     * 
     * 
     */
    public PrintableDomainSimilarity( final CombinableDomains combinable_domains,
                                      final double min,
                                      final double max,
                                      final double mean,
                                      final double median,
                                      final double sd,
                                      final int n,
                                      final int max_difference_in_counts,
                                      final int max_difference,
                                      final SortedMap<Species, SpeciesSpecificDomainSimilariyData> species_data,
                                      final DomainSimilaritySortField sort_field,
                                      final boolean sort_by_species_count_first,
                                      final boolean treat_as_binary_comparison ) {
        if ( combinable_domains == null ) {
            throw new IllegalArgumentException( "attempt to use null combinable domains" );
        }
        if ( sort_field == null ) {
            throw new IllegalArgumentException( "attempt to use null sorting" );
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
        _sort_field = sort_field;
        _sort_by_species_count_first = sort_by_species_count_first;
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

    private void addSpeciesSpecificDomainData( final StringBuffer sb,
                                               final Species species,
                                               final boolean html,
                                               final Map<String, Integer> tax_code_to_id_map ) {
        if ( getDetaildness() != DomainSimilarityCalculator.Detailedness.BASIC ) {
            sb.append( "[" );
        }
        if ( html ) {
            sb.append( "<b>" );
            final String tax_code = species.getSpeciesId();
            if ( !ForesterUtil.isEmpty( tax_code )
                    && ( ( tax_code_to_id_map != null ) && tax_code_to_id_map.containsKey( tax_code ) ) ) {
                sb.append( "<a href=\"" + SurfacingConstants.UNIPROT_TAXONOMY_ID_LINK
                        + tax_code_to_id_map.get( tax_code ) + "\" target=\"taxonomy_window\">" + tax_code + "</a>" );
            }
            else {
                sb.append( tax_code );
            }
            sb.append( "</b>" );
        }
        else {
            sb.append( species.getSpeciesId() );
        }
        if ( getDetaildness() != DomainSimilarityCalculator.Detailedness.BASIC ) {
            sb.append( ":" );
            sb.append( getSpeciesData().get( species ).toStringBuffer( getDetaildness(), html ) );
            sb.append( "]" );
        }
        if ( html ) {
            sb.append( "<br>" );
        }
        sb.append( PrintableDomainSimilarity.SPECIES_SEPARATOR );
    }

    private void boldEndIfSortedBy( final DomainSimilaritySortField sort_field, final StringBuffer sb ) {
        if ( getSortField() == sort_field ) {
            sb.append( "</b>" );
        }
    }

    private void boldStartIfSortedBy( final DomainSimilaritySortField sort_field, final StringBuffer sb ) {
        if ( getSortField() == sort_field ) {
            sb.append( "<b>" );
        }
    }

    private int compareByDomainId( final DomainSimilarity other ) {
        return getDomainId().compareTo( other.getDomainId() );
    }

    private int compareBySpeciesCount( final DomainSimilarity domain_similarity ) {
        final int s_this = getSpeciesData().size();
        final int s_other = domain_similarity.getSpeciesData().size();
        if ( s_this < s_other ) {
            return PrintableDomainSimilarity.BEFORE;
        }
        else if ( s_this > s_other ) {
            return PrintableDomainSimilarity.AFTER;
        }
        else {
            return PrintableDomainSimilarity.EQUAL;
        }
    }

    @Override
    public int compareTo( final DomainSimilarity domain_similarity ) {
        if ( this == domain_similarity ) {
            return PrintableDomainSimilarity.EQUAL;
        }
        else if ( domain_similarity == null ) {
            throw new IllegalArgumentException( "attempt to compare " + this.getClass() + " to null" );
        }
        else if ( domain_similarity.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to compare " + this.getClass() + " to "
                    + domain_similarity.getClass() );
        }
        switch ( getSortField() ) {
            case MIN:
                if ( isSortBySpeciesCountFirst() ) {
                    final int i = compareBySpeciesCount( domain_similarity );
                    if ( i != PrintableDomainSimilarity.EQUAL ) {
                        return i;
                    }
                }
                if ( getMinimalSimilarityScore() < domain_similarity.getMinimalSimilarityScore() ) {
                    return PrintableDomainSimilarity.BEFORE;
                }
                else if ( getMinimalSimilarityScore() > domain_similarity.getMinimalSimilarityScore() ) {
                    return PrintableDomainSimilarity.AFTER;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case MAX:
                if ( isSortBySpeciesCountFirst() ) {
                    final int i = compareBySpeciesCount( domain_similarity );
                    if ( i != PrintableDomainSimilarity.EQUAL ) {
                        return i;
                    }
                }
                if ( getMaximalSimilarityScore() < domain_similarity.getMaximalSimilarityScore() ) {
                    return PrintableDomainSimilarity.BEFORE;
                }
                else if ( getMaximalSimilarityScore() > domain_similarity.getMaximalSimilarityScore() ) {
                    return PrintableDomainSimilarity.AFTER;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case MEAN:
                if ( isSortBySpeciesCountFirst() ) {
                    final int i = compareBySpeciesCount( domain_similarity );
                    if ( i != PrintableDomainSimilarity.EQUAL ) {
                        return i;
                    }
                }
                if ( getMeanSimilarityScore() < domain_similarity.getMeanSimilarityScore() ) {
                    return PrintableDomainSimilarity.BEFORE;
                }
                else if ( getMeanSimilarityScore() > domain_similarity.getMeanSimilarityScore() ) {
                    return PrintableDomainSimilarity.AFTER;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case SD:
                if ( isSortBySpeciesCountFirst() ) {
                    final int i = compareBySpeciesCount( domain_similarity );
                    if ( i != PrintableDomainSimilarity.EQUAL ) {
                        return i;
                    }
                }
                if ( getStandardDeviationOfSimilarityScore() < domain_similarity
                        .getStandardDeviationOfSimilarityScore() ) {
                    return PrintableDomainSimilarity.BEFORE;
                }
                else if ( getStandardDeviationOfSimilarityScore() > domain_similarity
                        .getStandardDeviationOfSimilarityScore() ) {
                    return PrintableDomainSimilarity.AFTER;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case MAX_DIFFERENCE:
                if ( isSortBySpeciesCountFirst() ) {
                    final int i = compareBySpeciesCount( domain_similarity );
                    if ( i != PrintableDomainSimilarity.EQUAL ) {
                        return i;
                    }
                }
                if ( getMaximalDifference() > domain_similarity.getMaximalDifference() ) {
                    return PrintableDomainSimilarity.BEFORE;
                }
                else if ( getMaximalDifference() < domain_similarity.getMaximalDifference() ) {
                    return PrintableDomainSimilarity.AFTER;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case ABS_MAX_COUNTS_DIFFERENCE:
                if ( isSortBySpeciesCountFirst() ) {
                    final int i = compareBySpeciesCount( domain_similarity );
                    if ( i != PrintableDomainSimilarity.EQUAL ) {
                        return i;
                    }
                }
                if ( Math.abs( getMaximalDifferenceInCounts() ) > Math.abs( domain_similarity
                        .getMaximalDifferenceInCounts() ) ) {
                    return PrintableDomainSimilarity.BEFORE;
                }
                else if ( Math.abs( getMaximalDifferenceInCounts() ) < Math.abs( domain_similarity
                        .getMaximalDifferenceInCounts() ) ) {
                    return PrintableDomainSimilarity.AFTER;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case MAX_COUNTS_DIFFERENCE:
                if ( getSpeciesData().size() != 2 ) {
                    throw new RuntimeException( "attempt to sort by maximal difference with species not equal to two" );
                }
                if ( isSortBySpeciesCountFirst() ) {
                    final int i = compareBySpeciesCount( domain_similarity );
                    if ( i != PrintableDomainSimilarity.EQUAL ) {
                        return i;
                    }
                }
                if ( getMaximalDifferenceInCounts() > domain_similarity.getMaximalDifferenceInCounts() ) {
                    return PrintableDomainSimilarity.BEFORE;
                }
                else if ( getMaximalDifferenceInCounts() < domain_similarity.getMaximalDifferenceInCounts() ) {
                    return PrintableDomainSimilarity.AFTER;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case SPECIES_COUNT:
                final int i = compareBySpeciesCount( domain_similarity );
                if ( i != PrintableDomainSimilarity.EQUAL ) {
                    return i;
                }
                else {
                    return compareByDomainId( domain_similarity );
                }
            case DOMAIN_ID:
                return compareByDomainId( domain_similarity );
        }
        throw new AssertionError( "Unknown sort method: " + getSortField() );
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

    private CombinableDomains getCombinableDomains() {
        return _combinable_domains;
    }

    private DomainSimilarityCalculator.Detailedness getDetaildness() {
        return _detailedness;
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

    private DomainSimilaritySortField getSortField() {
        return _sort_field;
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
    public SortedMap<Species, SpeciesSpecificDomainSimilariyData> getSpeciesData() {
        return _species_data;
    }

    private StringBuffer getSpeciesDataInAlphabeticalOrder( final boolean html,
                                                            final Map<String, Integer> tax_code_to_id_map ) {
        final StringBuffer sb = new StringBuffer();
        for( final Species species : getSpeciesData().keySet() ) {
            addSpeciesSpecificDomainData( sb, species, html, tax_code_to_id_map );
        }
        return sb;
    }

    private StringBuffer getSpeciesDataInCustomOrder( final boolean html, final Map<String, Integer> tax_code_to_id_map ) {
        final StringBuffer sb = new StringBuffer();
        for( final Species order_species : getSpeciesCustomOrder() ) {
            if ( getSpeciesData().keySet().contains( order_species ) ) {
                addSpeciesSpecificDomainData( sb, order_species, html, tax_code_to_id_map );
            }
            else {
                sb.append( PrintableDomainSimilarity.NO_SPECIES );
                sb.append( PrintableDomainSimilarity.SPECIES_SEPARATOR );
            }
        }
        return sb;
    }

    @Override
    public double getStandardDeviationOfSimilarityScore() {
        return _sd;
    }

    private void init() {
        _detailedness = DomainSimilarityCalculator.Detailedness.PUNCTILIOUS;
    }

    private boolean isSortBySpeciesCountFirst() {
        return _sort_by_species_count_first;
    }

    private boolean isTreatAsBinaryComparison() {
        return _treat_as_binary_comparison;
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
                                        final Map<String, Integer> tax_code_to_id_map ) {
        switch ( print_option ) {
            case SIMPLE_TAB_DELIMITED:
                return toStringBufferSimpleTabDelimited();
            case HTML:
                return toStringBufferDetailedHTML( tax_code_to_id_map );
            default:
                throw new AssertionError( "Unknown print option: " + print_option );
        }
    }

    private StringBuffer toStringBufferDetailedHTML( final Map<String, Integer> tax_code_to_id_map ) {
        final StringBuffer sb = new StringBuffer();
        sb.append( "<tr>" );
        sb.append( "<td>" );
        boldStartIfSortedBy( DomainSimilaritySortField.DOMAIN_ID, sb );
        sb.append( "<a href=\"" + SurfacingConstants.PFAM_FAMILY_ID_LINK + getDomainId() + "\" target=\"pfam_window\">"
                + getDomainId() + "</a>" );
        boldEndIfSortedBy( DomainSimilaritySortField.DOMAIN_ID, sb );
        sb.append( "<a name=\"" + getDomainId() + "\">" );
        sb.append( "</td>" );
        sb.append( "<td>" );
        sb.append( "<a href=\"" + SurfacingConstants.GOOGLE_SCHOLAR_SEARCH + getDomainId()
                + "\" target=\"gs_window\">gs</a>" );
        sb.append( "</td>" );
        sb.append( "<td>" );
        boldStartIfSortedBy( DomainSimilaritySortField.MEAN, sb );
        sb.append( ForesterUtil.round( getMeanSimilarityScore(), 3 ) );
        boldEndIfSortedBy( DomainSimilaritySortField.MEAN, sb );
        sb.append( "</td>" );
        if ( !isTreatAsBinaryComparison() ) {
            sb.append( "<td>" );
            sb.append( "(" );
            boldStartIfSortedBy( DomainSimilaritySortField.SD, sb );
            sb.append( ForesterUtil.round( getStandardDeviationOfSimilarityScore(), 3 ) );
            boldEndIfSortedBy( DomainSimilaritySortField.SD, sb );
            sb.append( ")" );
            sb.append( "</td>" );
            sb.append( "<td>" );
            sb.append( "[" );
            boldStartIfSortedBy( DomainSimilaritySortField.MIN, sb );
            sb.append( ForesterUtil.round( getMinimalSimilarityScore(), 3 ) );
            boldEndIfSortedBy( DomainSimilaritySortField.MIN, sb );
            sb.append( "-" );
            boldStartIfSortedBy( DomainSimilaritySortField.MAX, sb );
            sb.append( ForesterUtil.round( getMaximalSimilarityScore(), 3 ) );
            boldEndIfSortedBy( DomainSimilaritySortField.MAX, sb );
            sb.append( "]" );
            sb.append( "</td>" );
        }
        sb.append( "<td>" );
        boldStartIfSortedBy( DomainSimilaritySortField.MAX_DIFFERENCE, sb );
        sb.append( getMaximalDifference() );
        boldEndIfSortedBy( DomainSimilaritySortField.MAX_DIFFERENCE, sb );
        sb.append( "</td>" );
        sb.append( "<td>" );
        if ( isTreatAsBinaryComparison() ) {
            boldStartIfSortedBy( DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE, sb );
            boldStartIfSortedBy( DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE, sb );
            sb.append( getMaximalDifferenceInCounts() );
            boldEndIfSortedBy( DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE, sb );
            boldStartIfSortedBy( DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE, sb );
        }
        else {
            boldStartIfSortedBy( DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE, sb );
            boldStartIfSortedBy( DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE, sb );
            sb.append( Math.abs( getMaximalDifferenceInCounts() ) );
            boldEndIfSortedBy( DomainSimilaritySortField.ABS_MAX_COUNTS_DIFFERENCE, sb );
            boldStartIfSortedBy( DomainSimilaritySortField.MAX_COUNTS_DIFFERENCE, sb );
        }
        sb.append( "</td>" );
        if ( !isTreatAsBinaryComparison() ) {
            sb.append( "<td>" );
            if ( ( getSortField() == DomainSimilaritySortField.SPECIES_COUNT ) || isSortBySpeciesCountFirst() ) {
                sb.append( "<b>" );
            }
            sb.append( getSpeciesData().size() );
            if ( ( getSortField() == DomainSimilaritySortField.SPECIES_COUNT ) || isSortBySpeciesCountFirst() ) {
                sb.append( "</b>" );
            }
            sb.append( "</td>" );
        }
        if ( ( getSpeciesCustomOrder() == null ) || getSpeciesCustomOrder().isEmpty() ) {
            sb.append( "<td>" );
            sb.append( getSpeciesDataInAlphabeticalOrder( true, tax_code_to_id_map ) );
            sb.append( "</td>" );
        }
        else {
            sb.append( "<td>" );
            sb.append( getSpeciesDataInCustomOrder( true, tax_code_to_id_map ) );
            sb.append( "</td>" );
        }
        sb.append( "</tr>" );
        return sb;
    }

    private StringBuffer toStringBufferSimpleTabDelimited() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getDomainId() );
        switch ( getSortField() ) {
            case MIN:
                sb.append( TAB );
                sb.append( ForesterUtil.round( getMinimalSimilarityScore(), 3 ) );
                break;
            case MAX:
                sb.append( TAB );
                sb.append( ForesterUtil.round( getMaximalSimilarityScore(), 3 ) );
                break;
            case MEAN:
                sb.append( TAB );
                sb.append( ForesterUtil.round( getMeanSimilarityScore(), 3 ) );
                break;
            case SD:
                sb.append( TAB );
                sb.append( ForesterUtil.round( getStandardDeviationOfSimilarityScore(), 3 ) );
                break;
            case MAX_DIFFERENCE:
                sb.append( TAB );
                sb.append( getMaximalDifference() );
            case ABS_MAX_COUNTS_DIFFERENCE:
            case MAX_COUNTS_DIFFERENCE:
                sb.append( TAB );
                if ( isTreatAsBinaryComparison() ) {
                    sb.append( getMaximalDifferenceInCounts() );
                }
                else {
                    sb.append( Math.abs( getMaximalDifferenceInCounts() ) );
                }
                break;
            case SPECIES_COUNT:
                sb.append( TAB );
                sb.append( getSpeciesData().size() );
                break;
            case DOMAIN_ID:
                break;
            default:
                throw new AssertionError( "Unknown sort method: " + getSortField() );
        }
        // ^^     if ( getGoAnnotationOutput() != DomainSimilarityCalculator.GoAnnotationOutput.NONE ) {
        // ^^       sb.append( TAB );
        // ^^       addGoInformation( sb, true, false );
        // ^^   }
        return sb;
    }

    public static enum PRINT_OPTION {
        SIMPLE_TAB_DELIMITED, HTML;
    }
}

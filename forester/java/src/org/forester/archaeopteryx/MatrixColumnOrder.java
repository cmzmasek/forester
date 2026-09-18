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

package org.forester.archaeopteryx;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The order of a heat-map MATRIX's columns -- the per-tab <b>View &rarr; Order Matrix Columns</b> setting.
 * <p>
 * Only the MATRIX columns move, and only among the slots MATRIX columns already occupy: every other column (a colour
 * strip, a symbol, a bar) keeps its exact place. Clustering and frequency need numbers, and the matrix is the block
 * whose column order is what a reader reads patterns from.
 * <p>
 * Every data-driven mode first puts the columns in TABLE order and works from there, so a result never depends on
 * which mode the tab was in before -- switching Alphabetical &rarr; Clustered and Same as Table &rarr; Clustered
 * must give the same clustering.
 * <p>
 * A cell that is not assessed never counts as 0, in any mode: blank is absence of evidence, not evidence of absence.
 * <p>
 * <b>Clustered</b> is complete-linkage hierarchical clustering (Sørensen 1948) of the columns on Euclidean distance --
 * the defaults of R's {@code pheatmap}, {@code heatmap.2} and {@code ComplexHeatmap}, and the clustered heat map of
 * Eisen et al. 1998. It is written to give the SAME column order as {@code hclust(dist(t(m)), method = "complete")$order}
 * in R for a matrix {@code m} whose columns are in table order: the same distance (a pair of strains is used only where
 * both genes were assessed, the squared sum scaled up by {@code tips / used}, R's {@code dist} convention), the same
 * tie-break (the lowest-indexed pair first, clusters labelled by their lowest member index, as R's Fortran
 * {@code hclust} does) and the same leaf order (R's {@code hcass2}: a singleton before a cluster, the lower index first
 * between two singletons, the earlier-formed first between two clusters). One deliberate divergence: two genes no
 * strain was assessed for together have no distance at all; R refuses such input, and here they are treated as
 * maximally distant, so they join last.
 */
final class MatrixColumnOrder {

    /** The modes, in menu order. */
    enum Mode {
        CLUSTERED( "Clustered (co-occurrence)",
                   "Columns whose values agree across the tips sit together: complete-linkage clustering on Euclidean "
                           + "distance, the default of R's heat maps. A tip missing either value is left out of that "
                           + "pair, never read as 0." ),
        TABLE( "Same as Table", "The order the file, or the imported table, lists the columns in." ),
        ALPHABETICAL( "Alphabetical", "By column name, ignoring case." ),
        FREQUENCY( "Frequency",
                   "Highest mean value first, over the tips that have a value -- on 0/1 data, the fraction of tips "
                           + "carrying it." ),
        MANUAL( "Manual",
                "Your own order: moving a field with the arrows in Tools > Annotation Fields selects this, and nothing "
                        + "re-sorts it." );

        private final String _label;
        private final String _tooltip;

        Mode( final String label, final String tooltip ) {
            _label = label;
            _tooltip = tooltip;
        }

        String label() {
            return _label;
        }

        String tooltip() {
            return _tooltip;
        }
    }

    /** A new tab, and Reset to Defaults, cluster the matrix. */
    static final Mode DEFAULT = Mode.CLUSTERED;

    /**
     * {@code specs} with its MATRIX columns re-ordered per {@code mode}, slot-preserving: every non-MATRIX spec keeps
     * its index, and the MATRIX specs are permuted among the indices they occupy. {@link Mode#MANUAL} (and a spec list
     * with fewer than two MATRIX columns) returns the list unchanged.
     */
    static List<AnnotationColumns.ColumnSpec> apply( final List<AnnotationColumns.ColumnSpec> specs, final Mode mode,
                                                     final Phylogeny phy ) {
        if ( ( specs == null ) || ( mode == null ) || ( mode == Mode.MANUAL ) ) {
            return specs;
        }
        final List<Integer> slots = new ArrayList<Integer>();
        final List<AnnotationColumns.ColumnSpec> matrix = new ArrayList<AnnotationColumns.ColumnSpec>();
        final List<String> refs = new ArrayList<String>();
        for( int i = 0; i < specs.size(); ++i ) {
            if ( specs.get( i )._type == AnnotationColumns.Type.MATRIX ) {
                slots.add( i );
                matrix.add( specs.get( i ) );
                refs.add( specs.get( i )._ref );
            }
        }
        if ( matrix.size() < 2 ) {
            return specs;
        }
        final List<String> ordered = order( refs, mode, phy );
        final List<AnnotationColumns.ColumnSpec> out = new ArrayList<AnnotationColumns.ColumnSpec>( specs );
        final boolean[] used = new boolean[ matrix.size() ];
        for( int k = 0; k < ordered.size(); ++k ) {
            for( int m = 0; m < matrix.size(); ++m ) {
                if ( !used[ m ] && matrix.get( m )._ref.equals( ordered.get( k ) ) ) {
                    used[ m ] = true;
                    out.set( slots.get( k ), matrix.get( m ) );
                    break;
                }
            }
        }
        return out;
    }

    /** What an edit of the column list leaves on screen: the specs to show, and the mode the tab is in afterwards. */
    record Resolved( List<AnnotationColumns.ColumnSpec> specs, Mode mode ) {
    }

    /**
     * The Annotation Fields dialog's OK. If the user MOVED a row and the matrix order they left differs from what
     * {@code mode} would produce, their order wins and the tab becomes {@link Mode#MANUAL} -- re-sorting must never
     * silently undo a move. Otherwise {@code mode} is re-applied, which also places any newly added matrix fields. So
     * moving a row and back again, or moving only a colour strip, leaves the tab in its mode.
     */
    static Resolved resolveEdited( final List<AnnotationColumns.ColumnSpec> specs, final boolean moved, final Mode mode,
                                   final Phylogeny phy ) {
        final List<AnnotationColumns.ColumnSpec> by_mode = apply( specs, mode, phy );
        if ( moved && !matrixRefs( specs ).equals( matrixRefs( by_mode ) ) ) {
            return new Resolved( specs, Mode.MANUAL );
        }
        return new Resolved( by_mode, mode );
    }

    /** The refs of the MATRIX columns of {@code specs}, in order (empty for null). */
    static List<String> matrixRefs( final List<AnnotationColumns.ColumnSpec> specs ) {
        final List<String> out = new ArrayList<String>();
        if ( specs != null ) {
            for( final AnnotationColumns.ColumnSpec s : specs ) {
                if ( s._type == AnnotationColumns.Type.MATRIX ) {
                    out.add( s._ref );
                }
            }
        }
        return out;
    }

    /** {@code refs} in {@code mode} order. {@link Mode#MANUAL} returns them as given. */
    static List<String> order( final List<String> refs, final Mode mode, final Phylogeny phy ) {
        if ( ( mode == null ) || ( mode == Mode.MANUAL ) ) {
            return new ArrayList<String>( refs );
        }
        final List<String> table = inSourceOrder( refs, TreePanelUtil.propertyRefsInSourceOrder( phy ) );
        switch ( mode ) {
            case TABLE:
                return table;
            case ALPHABETICAL:
                return alphabetical( table );
            case FREQUENCY:
                return byFrequency( table, phy );
            case CLUSTERED:
                return clustered( table, phy );
            default:
                return table;
        }
    }

    /**
     * {@code refs} re-ordered to follow {@code source_order}. A ref the source order does not name keeps its incoming
     * relative order and goes at the end -- it cannot be placed, so it must not be dropped. That is not hypothetical:
     * an ELEMENT SLOT candidate (taxonomy, sequence, ...) is a colorable field that is not a node property at all, so
     * it never appears in the source order.
     */
    static List<String> inSourceOrder( final List<String> refs, final List<String> source_order ) {
        final List<String> out = new ArrayList<String>();
        for( final String ref : source_order ) {
            if ( refs.contains( ref ) && !out.contains( ref ) ) {
                out.add( ref );
            }
        }
        for( final String ref : refs ) {
            if ( !out.contains( ref ) ) {
                out.add( ref );
            }
        }
        return out;
    }

    /** By the name the column header shows, ignoring case; the raw ref breaks a tie (the chooser's own sort key). */
    static List<String> alphabetical( final List<String> refs ) {
        final List<String> out = new ArrayList<String>( refs );
        Collections.sort( out, ( a, b ) -> {
            final int c = PropertyColorScheme.displayName( a ).compareToIgnoreCase( PropertyColorScheme.displayName( b ) );
            return ( c != 0 ) ? c : a.compareTo( b );
        } );
        return out;
    }

    /**
     * By the mean value over the tips that HAVE a value, highest first -- on 0/1 data exactly the fraction of strains
     * carrying the gene, and just as meaningful for an ordinal certainty or an abundance. A blank is left out, never
     * averaged in as 0. A column with no value at all goes last. Ties keep the incoming (table) order: the sort is
     * stable.
     */
    static List<String> byFrequency( final List<String> refs, final Phylogeny phy ) {
        final Double[][] v = values( refs, phy );
        final double[] mean = new double[ refs.size() ];
        for( int g = 0; g < refs.size(); ++g ) {
            double sum = 0.0;
            int n = 0;
            for( final Double x : v[ g ] ) {
                if ( x != null ) {
                    sum += x;
                    ++n;
                }
            }
            mean[ g ] = ( n > 0 ) ? ( sum / n ) : Double.NEGATIVE_INFINITY;
        }
        final List<Integer> idx = new ArrayList<Integer>();
        for( int g = 0; g < refs.size(); ++g ) {
            idx.add( g );
        }
        Collections.sort( idx, ( a, b ) -> Double.compare( mean[ b ], mean[ a ] ) ); // stable: ties keep table order
        final List<String> out = new ArrayList<String>();
        for( final int g : idx ) {
            out.add( refs.get( g ) );
        }
        return out;
    }

    /** Complete-linkage clustering of the columns on Euclidean distance, in R's {@code hclust} leaf order. */
    static List<String> clustered( final List<String> refs, final Phylogeny phy ) {
        if ( refs.size() < 3 ) {
            return new ArrayList<String>( refs ); // two leaves have one order up to a flip, and R puts index 1 first
        }
        final int[] leaf_order = completeLinkageOrder( distances( values( refs, phy ) ) );
        final List<String> out = new ArrayList<String>();
        for( final int g : leaf_order ) {
            out.add( refs.get( g ) );
        }
        return out;
    }

    /**
     * The numeric value of every column at every tip, {@code [column][tip]}, {@code null} where the tip has no value --
     * the same number the heat-map cell is drawn from. Every tip counts, collapsed or not, like the matrix's shared
     * colour scale: an order that changed when a clade is collapsed would make the columns jump under the reader.
     */
    static Double[][] values( final List<String> refs, final Phylogeny phy ) {
        final List<PhylogenyNode> tips = phy.getExternalNodes();
        final Double[][] v = new Double[ refs.size() ][ tips.size() ];
        for( int g = 0; g < refs.size(); ++g ) {
            for( int t = 0; t < tips.size(); ++t ) {
                v[ g ][ t ] = PropertyColorScheme.parseNumber( PropertyColorScheme.valueFor( tips.get( t ),
                                                                                               refs.get( g ) ) );
            }
        }
        return v;
    }

    /**
     * Euclidean distance between every two columns, R's {@code dist(method = "euclidean")} convention for missing
     * values: only the tips where BOTH columns have a value are used, and the squared sum is scaled up by
     * {@code tips / used}, so a pair assessed on fewer strains is not made to look closer for it. A pair that shares
     * no assessed tip at all has no distance; it is {@code +Infinity} here (R returns NA, and its hclust then refuses).
     */
    static double[][] distances( final Double[][] v ) {
        final int n = v.length;
        final double[][] d = new double[ n ][ n ];
        for( int a = 0; a < n; ++a ) {
            for( int b = a + 1; b < n; ++b ) {
                double sum = 0.0;
                int used = 0;
                final int tips = v[ a ].length;
                for( int t = 0; t < tips; ++t ) {
                    if ( ( v[ a ][ t ] != null ) && ( v[ b ][ t ] != null ) ) {
                        final double dev = v[ a ][ t ] - v[ b ][ t ];
                        sum += dev * dev;
                        ++used;
                    }
                }
                final double dist;
                if ( used == 0 ) {
                    dist = Double.POSITIVE_INFINITY;
                }
                else {
                    dist = Math.sqrt( ( used == tips ) ? sum : ( sum / ( ( double ) used / tips ) ) );
                }
                d[ a ][ b ] = dist;
                d[ b ][ a ] = dist;
            }
        }
        return d;
    }

    /**
     * The leaf order of the complete-linkage dendrogram over the symmetric distance matrix {@code d}, reproducing R's
     * {@code hclust(method = "complete")$order}.
     * <p>
     * Each step merges the closest two clusters; a cluster is labelled by its lowest member index, and a tie goes to
     * the lexicographically lowest {@code (label, label)} pair, which is the pair R's nearest-neighbour scan finds
     * first. The merged cluster's distance to every other cluster is the larger of the two (complete linkage, the
     * Lance-Williams update R uses). The leaf order then follows R's {@code hcass2}: at each merge the left child is a
     * singleton rather than a cluster, the lower index between two singletons, and the earlier-formed between two
     * clusters.
     */
    static int[] completeLinkageOrder( final double[][] d ) {
        final int n = d.length;
        if ( n == 0 ) {
            return new int[ 0 ];
        }
        if ( n == 1 ) {
            return new int[] { 0 };
        }
        final double[][] dist = new double[ n ][];
        for( int i = 0; i < n; ++i ) {
            dist[ i ] = d[ i ].clone();
        }
        final boolean[] active = new boolean[ n ];
        // node id of the cluster currently labelled i: -(i+1) for a singleton, (stage+1) for a merged cluster
        final int[] node = new int[ n ];
        for( int i = 0; i < n; ++i ) {
            active[ i ] = true;
            node[ i ] = -( i + 1 );
        }
        final int[] left = new int[ n - 1 ];
        final int[] right = new int[ n - 1 ];
        for( int stage = 0; stage < ( n - 1 ); ++stage ) {
            int bi = -1;
            int bj = -1;
            double best = Double.NaN;
            for( int i = 0; i < n; ++i ) {
                if ( !active[ i ] ) {
                    continue;
                }
                for( int j = i + 1; j < n; ++j ) {
                    if ( !active[ j ] ) {
                        continue;
                    }
                    // strictly smaller only: the FIRST pair in (i, j) order wins a tie, as in R's scan
                    if ( ( bi < 0 ) || ( dist[ i ][ j ] < best ) ) {
                        best = dist[ i ][ j ];
                        bi = i;
                        bj = j;
                    }
                }
            }
            // R's hcass2 orientation: singleton before cluster; lower index between singletons; earlier stage between
            // clusters. Node ids make all three one rule: singletons are negative, clusters positive by stage.
            final int a = node[ bi ];
            final int b = node[ bj ];
            final int l;
            final int r;
            if ( ( a < 0 ) && ( b < 0 ) ) {
                l = ( -a < -b ) ? a : b;
                r = ( -a < -b ) ? b : a;
            }
            else if ( ( a < 0 ) || ( b < 0 ) ) {
                l = ( a < 0 ) ? a : b;
                r = ( a < 0 ) ? b : a;
            }
            else {
                l = Math.min( a, b );
                r = Math.max( a, b );
            }
            left[ stage ] = l;
            right[ stage ] = r;
            // complete linkage: the merged cluster (kept under the lower label, bi) is as far as its farther half
            for( int k = 0; k < n; ++k ) {
                if ( active[ k ] && ( k != bi ) && ( k != bj ) ) {
                    final double m = Math.max( dist[ bi ][ k ], dist[ bj ][ k ] );
                    dist[ bi ][ k ] = m;
                    dist[ k ][ bi ] = m;
                }
            }
            active[ bj ] = false;
            node[ bi ] = stage + 1;
        }
        // expand from the root, in place, left child before right -- hcass2's IORDER
        final List<Integer> order = new ArrayList<Integer>();
        order.add( left[ n - 2 ] );
        order.add( right[ n - 2 ] );
        for( int stage = n - 2; stage >= 1; --stage ) {
            final int pos = order.indexOf( stage ); // the cluster formed at this stage (1-based node id)
            if ( pos >= 0 ) {
                order.set( pos, left[ stage - 1 ] );
                order.add( pos + 1, right[ stage - 1 ] );
            }
        }
        final int[] out = new int[ n ];
        for( int i = 0; i < n; ++i ) {
            out[ i ] = -order.get( i ) - 1;
        }
        return out;
    }

    private MatrixColumnOrder() {
    }
}

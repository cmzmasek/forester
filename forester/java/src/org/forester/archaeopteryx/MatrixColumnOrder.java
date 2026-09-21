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
import java.util.Arrays;
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
 * <p>
 * <b>Clustered (ignoring shared absence)</b> is the same complete-linkage clustering read from a different distance:
 * the Bray-Curtis dissimilarity (Bray &amp; Curtis 1957), R {@code vegan}'s {@code vegdist} default. Euclidean distance
 * has the double-zero problem -- two genes that are both ABSENT from the same strains are counted as agreeing there --
 * so on a sparse pan-genome the rare genes cluster together merely for being rare. Bray-Curtis drops a tip where both
 * columns are 0 instead of scoring it as agreement; see {@link #brayCurtisDistances}.
 */
final class MatrixColumnOrder {

    /** The modes, in menu order. */
    enum Mode {
        CLUSTERED( "Clustered (co-occurrence)",
                   "Columns whose values agree across the tips sit together: complete-linkage clustering on Euclidean "
                           + "distance, the default of R's heat maps. A tip missing either value is left out of that "
                           + "pair, never read as 0." ),
        CLUSTERED_PRESENCE( "Clustered (ignoring shared absence)",
                            "Columns found in the same tips sit together: a tip where BOTH columns are 0 is left out "
                                    + "of that pair, so two rare genes no longer look alike merely for being rare. "
                                    + "Complete-linkage clustering on the Bray-Curtis dissimilarity, the default of R's "
                                    + "vegdist(). For values that are 0 or more." ),
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
            case CLUSTERED_PRESENCE:
                return clusteredByPresence( table, phy );
            default:
                return table;
        }
    }

    /**
     * {@code refs} re-ordered to follow {@code source_order}. A ref the source order does not name cannot be placed,
     * so it must not be dropped: it goes at the end, and the unplaceable ones are sorted among themselves. That is
     * not hypothetical -- an ELEMENT SLOT candidate (taxonomy, sequence, ...) is a colorable field that is not a
     * node property at all, and {@code propertyRefsInSourceOrder} walks the DISPLAYED tree, so a matrix column
     * whose values all sit outside the current subtree drops out of the source order too.
     * <p>
     * The tail is SORTED rather than left in its incoming order so that this function depends only on its inputs as
     * a SET. Every data-driven mode normalises through here precisely so a result cannot depend on the order the
     * columns happened to be in, and an incoming-order tail broke that promise exactly when it mattered: the column
     * order is computed from the PRE-sort order and the dendrogram that describes it from the POST-sort one, so a
     * single unplaceable ref made the two disagree and the dendrogram silently vanished from a correctly clustered
     * matrix.
     */
    static List<String> inSourceOrder( final List<String> refs, final List<String> source_order ) {
        final List<String> out = new ArrayList<String>();
        for( final String ref : source_order ) {
            if ( refs.contains( ref ) && !out.contains( ref ) ) {
                out.add( ref );
            }
        }
        final List<String> unplaceable = new ArrayList<String>();
        for( final String ref : refs ) {
            if ( !out.contains( ref ) && !unplaceable.contains( ref ) ) {
                unplaceable.add( ref );
            }
        }
        Collections.sort( unplaceable );
        out.addAll( unplaceable );
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
        return leafOrder( refs, completeLinkageOrder( distances( values( refs, phy ) ) ) );
    }

    /**
     * Complete-linkage clustering of the columns on Bray-Curtis dissimilarity, in R's {@code hclust} leaf order --
     * the same dendrogram, read from a distance that ignores the tips where both columns are absent.
     */
    static List<String> clusteredByPresence( final List<String> refs, final Phylogeny phy ) {
        if ( refs.size() < 3 ) {
            return new ArrayList<String>( refs ); // as above: one order up to a flip
        }
        return leafOrder( refs, completeLinkageOrder( brayCurtisDistances( values( refs, phy ) ) ) );
    }

    /** {@code refs} read out in the order the dendrogram's leaves come in. */
    private static List<String> leafOrder( final List<String> refs, final int[] leaf_order ) {
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
     * Bray-Curtis dissimilarity between every two columns: {@code sum|x - y| / sum(x + y)} over the tips where BOTH
     * columns have a value (pairwise deletion, R {@code vegan}'s {@code vegdist(method = "bray", na.rm = TRUE)}).
     * Being a ratio it normalizes itself, so -- unlike Euclidean -- it needs no scale-up for the tips it dropped.
     * <p>
     * A tip where both columns are 0 adds nothing to either sum, so it drops out: the DOUBLE ZERO, which Euclidean
     * reads as agreement, is simply not evidence here. Two genes that are each rare, and never in the same strain, are
     * maximally distant (1) rather than nearly identical. On 0/1 data this is exactly the Sorensen-Dice dissimilarity
     * {@code 1 - 2|A and B| / (|A| + |B|)}.
     * <p>
     * Two divergences from {@code vegdist}, both forced: a pair that shares no assessed tip at all is
     * {@code +Infinity} (vegdist gives NA, and hclust then refuses), and a pair that is 0 at every tip it shares --
     * the undefined 0/0 -- is 0, because those two columns are identical wherever they can be compared (vegdist again
     * gives NaN). The same rule covers the only other way the denominator can fail, negative values cancelling: a pair
     * that cannot be told apart is 0, one that differs over a non-positive total is {@code +Infinity} rather than a
     * negative or NaN distance. Bray-Curtis is meant for values that are 0 or more; vegan warns for negative data too.
     * <p>
     * The result is a DISSIMILARITY, not a metric -- it does not obey the triangle inequality, which complete linkage
     * does not require. On the data it is meant for (values 0 or more) it lies in [0, 1]: 0 for columns that agree
     * wherever both were assessed, 1 for columns that share no tip they are both present at.
     */
    static double[][] brayCurtisDistances( final Double[][] v ) {
        final int n = v.length;
        final double[][] d = new double[ n ][ n ];
        for( int a = 0; a < n; ++a ) {
            for( int b = a + 1; b < n; ++b ) {
                double num = 0.0;
                double den = 0.0;
                int used = 0;
                final int tips = v[ a ].length;
                for( int t = 0; t < tips; ++t ) {
                    if ( ( v[ a ][ t ] != null ) && ( v[ b ][ t ] != null ) ) {
                        num += Math.abs( v[ a ][ t ] - v[ b ][ t ] );
                        den += v[ a ][ t ] + v[ b ][ t ];
                        ++used;
                    }
                }
                final double dist;
                if ( used == 0 ) {
                    dist = Double.POSITIVE_INFINITY;
                }
                else if ( den <= 0.0 ) {
                    dist = ( num == 0.0 ) ? 0.0 : Double.POSITIVE_INFINITY;
                }
                else {
                    dist = num / den;
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
        return completeLinkage( d ).order();
    }

    /**
     * The dendrogram behind a matrix's column order -- what a drawn column dendrogram is made of -- or {@code null}
     * when none may be drawn.
     * <p>
     * {@code refs_in_drawn_order} is the matrix as it is ON SCREEN. The clustering is recomputed from the data and
     * its leaves are then required to BE that order: if they are not, this returns null rather than a dendrogram.
     * That is the whole guard. A dendrogram drawn over a matrix it does not describe would have connectors that
     * cross, a picture asserting a grouping the columns do not have -- and the ways to get there are ordinary: a
     * column dragged by hand (the tab goes MANUAL, but nothing else has to), a figure restored from a file, a cached
     * dendrogram outliving an edit to the data. Checking the answer beats trying to enumerate the causes.
     */
    static Dendrogram dendrogramFor( final List<String> refs_in_drawn_order, final Mode mode, final Phylogeny phy ) {
        if ( ( mode != Mode.CLUSTERED ) && ( mode != Mode.CLUSTERED_PRESENCE ) ) {
            return null; // the other modes did not come from a clustering, so there is no tree behind them
        }
        if ( ( refs_in_drawn_order == null ) || ( refs_in_drawn_order.size() < 3 ) ) {
            return null; // fewer than three columns: clustered() returns them as they are, without clustering
        }
        final List<String> table = inSourceOrder( refs_in_drawn_order, TreePanelUtil.propertyRefsInSourceOrder( phy ) );
        final Double[][] v = values( table, phy );
        final Dendrogram d = completeLinkage( ( mode == Mode.CLUSTERED ) ? distances( v ) : brayCurtisDistances( v ) );
        final List<String> leaves = new ArrayList<String>();
        for( final int g : d.order() ) {
            leaves.add( table.get( g ) );
        }
        return leaves.equals( refs_in_drawn_order ) ? d : null;
    }

    /**
     * The complete-linkage dendrogram of a distance matrix -- the same structure R's {@code hclust} returns, so it can
     * be pinned against it wholesale:
     * <ul>
     * <li>{@code left} / {@code right}: one merge per stage, in R's {@code $merge} node-id convention -- a NEGATIVE
     * value is the singleton of that index (-1 = column 0), a POSITIVE value the cluster formed at that earlier
     * stage. The pair is oriented by R's {@code hcass2} rule (singleton before cluster, lower index between
     * singletons, earlier stage between clusters), which is what makes the leaf order reproducible;</li>
     * <li>{@code height}: the distance each merge happened AT ({@code $height}) -- this is what a drawn dendrogram's
     * bar heights are, and what tells a reader whether a block of columns is one tight cluster or two loose ones.
     * It can be {@code +Infinity} for a pair of columns that share no assessed tip (see {@link #distances}), so
     * anything that draws it has to decide what to do with a merge that has no finite height;</li>
     * <li>{@code order}: the leaves as {@code $order} reads them out, which is the column order itself.</li>
     * </ul>
     */
    record Dendrogram( int[] left, int[] right, double[] height, int[] order ) {

        /** The number of merges: one fewer than the number of columns (0 for a single column). */
        int stages() {
            return height.length;
        }

        // A record whose components are ARRAYS gets equals/hashCode by array IDENTITY and a toString of four
        // "[I@1a2b3c" -- so two structurally identical dendrograms would compare unequal, and a failure message
        // would print nothing a reader could use. Both are spelled out here rather than left as a trap.

        @Override
        public boolean equals( final Object o ) {
            if ( this == o ) {
                return true;
            }
            if ( !( o instanceof Dendrogram ) ) {
                return false;
            }
            final Dendrogram d = ( Dendrogram ) o;
            return Arrays.equals( left, d.left ) && Arrays.equals( right, d.right )
                    && Arrays.equals( height, d.height ) && Arrays.equals( order, d.order );
        }

        @Override
        public int hashCode() {
            return ( ( ( ( 31 * Arrays.hashCode( left ) ) + Arrays.hashCode( right ) ) * 31 )
                    + Arrays.hashCode( height ) ) * 31 + Arrays.hashCode( order );
        }

        @Override
        public String toString() {
            final StringBuilder sb = new StringBuilder( "Dendrogram[order=" );
            sb.append( Arrays.toString( order ) ).append( ", merges=" );
            for( int i = 0; i < stages(); ++i ) {
                sb.append( i > 0 ? ", " : "" ).append( '(' ).append( left[ i ] ).append( ',' ).append( right[ i ] )
                        .append( " @" ).append( height[ i ] ).append( ')' );
            }
            return sb.append( ']' ).toString();
        }
    }

    static Dendrogram completeLinkage( final double[][] d ) {
        final int n = d.length;
        if ( n == 0 ) {
            return new Dendrogram( new int[ 0 ], new int[ 0 ], new double[ 0 ], new int[ 0 ] );
        }
        if ( n == 1 ) {
            return new Dendrogram( new int[ 0 ], new int[ 0 ], new double[ 0 ], new int[] { 0 } );
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
        final double[] height = new double[ n - 1 ];
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
            height[ stage ] = best; // R's $height: the distance the two clusters were apart when they merged
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
        return new Dendrogram( left, right, height, out );
    }

    private MatrixColumnOrder() {
    }
}

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
import java.util.List;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * {@link MatrixColumnOrder}: the five modes of View &gt; Order Matrix Columns.
 * <p>
 * The clustering expectations are R's OWN output, generated deliberately with R 4.5.3 --
 * {@code hclust(dist(t(m)), method = "complete")} on the same small matrices -- not restated here in this test's own
 * words: a second implementation of the rule written by the same hand could agree with the first and guard nothing.
 * The four cases are chosen to reach every branch of the algorithm: no ties, zero-distance ties (duplicate columns),
 * a matrix where EVERY pair ties, and missing values (R's pairwise-deletion scaling). A wider cross-check against R --
 * the 40-gene demo, the pinned contract fixture, random, binary and block-missing matrices -- was run by hand; this
 * test keeps the small cases that exercise each rule, so the suite needs no R.
 */
public final class MatrixColumnOrderTest {

    private static final double EPS = 1.0e-6;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MatrixColumnOrder: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return testDistancesMatchR() && testMissingIsNotZero() && testNoOverlap() && testClusteredOrderMatchesR()
                && testLinkageIsComplete()
                && testClusteredIgnoresIncomingOrder() && testFrequency() && testAlphabetical() && testTableOrder()
                && testManualAndEdges() && testApplyIsSlotPreserving() && testDialogResolution() && testDeterministic();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [MatrixColumnOrderTest] " + msg );
        return false;
    }

    // ---- fixtures --------------------------------------------------------------------------------------------------

    /**
     * A star tree whose tips carry the columns of {@code m} ({@code m[strain][gene]}, {@code null} = not assessed) as
     * {@code meta:<gene>} properties, stated in {@code genes} order -- the table order.
     */
    private static Phylogeny tree( final String[] genes, final Integer[][] m ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int s = 0; s < m.length; ++s ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "s" + s );
            final PropertiesList pl = new PropertiesList();
            for( int g = 0; g < genes.length; ++g ) {
                if ( m[ s ][ g ] != null ) {
                    pl.addProperty( new Property( "meta:" + genes[ g ], String.valueOf( m[ s ][ g ] ), "", "xsd:integer",
                                                  AppliesTo.NODE ) );
                }
            }
            tip.getNodeData().setProperties( pl );
            root.addAsChild( tip );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static List<String> refs( final String... genes ) {
        final List<String> out = new ArrayList<String>();
        for( final String g : genes ) {
            out.add( "meta:" + g );
        }
        return out;
    }

    /** {@code strains x genes} from per-gene columns, for readability. */
    private static Integer[][] columns( final Integer[]... cols ) {
        final Integer[][] m = new Integer[ cols[ 0 ].length ][ cols.length ];
        for( int g = 0; g < cols.length; ++g ) {
            for( int s = 0; s < cols[ g ].length; ++s ) {
                m[ s ][ g ] = cols[ g ][ s ];
            }
        }
        return m;
    }

    /** R: cbind(g1=c(0,0,1,4,4,3), g2=c(0,1,1,4,3,3), g3=c(4,4,4,0,0,1), g4=c(4,3,4,1,0,0), g5=c(2,2,2,2,2,2)). */
    private static final String[]  NT_GENES = { "g1", "g2", "g3", "g4", "g5" };
    private static final Integer[][] NT     = columns( new Integer[] { 0, 0, 1, 4, 4, 3 },
                                                       new Integer[] { 0, 1, 1, 4, 3, 3 },
                                                       new Integer[] { 4, 4, 4, 0, 0, 1 },
                                                       new Integer[] { 4, 3, 4, 1, 0, 0 },
                                                       new Integer[] { 2, 2, 2, 2, 2, 2 } );

    // ---- distance -------------------------------------------------------------------------------------------------

    /** R's dist() of the no-ties matrix, in R's own (column-major lower triangle) order. */
    private static boolean testDistancesMatchR() {
        final double[] r = { 1.414214, 8.774964, 8.246211, 4.242641, 7.937254, 7.483315, 3.464102, 1.732051, 4.582576,
                4.242641 };
        final double[][] d = MatrixColumnOrder.distances( MatrixColumnOrder.values( refs( NT_GENES ),
                                                                                    tree( NT_GENES, NT ) ) );
        int k = 0;
        for( int a = 0; a < 5; ++a ) {
            for( int b = a + 1; b < 5; ++b ) {
                if ( Math.abs( d[ a ][ b ] - r[ k ] ) > EPS ) {
                    return fail( "distance " + NT_GENES[ a ] + "-" + NT_GENES[ b ] + " should be R's " + r[ k ] + ", got "
                            + d[ a ][ b ] );
                }
                if ( d[ a ][ b ] != d[ b ][ a ] ) {
                    return fail( "the distance matrix must be symmetric" );
                }
                ++k;
            }
        }
        return true;
    }

    /**
     * Missing is not zero. X has BLANKS where Z has zeros. Pairwise deletion compares X and Z only where both were
     * assessed (they agree: 0), and X with the all-zero Y only on those two strains, scaled up by 4/2 -- R gives 8.
     * Filling the blanks with 0 would give 5.657 instead, so this value can only come from the right rule.
     */
    private static boolean testMissingIsNotZero() {
        final String[] genes = { "X", "Z", "Y" };
        final Phylogeny phy = tree( genes, columns( new Integer[] { 4, 4, null, null }, new Integer[] { 4, 4, 0, 0 },
                                                    new Integer[] { 0, 0, 0, 0 } ) );
        final double[][] d = MatrixColumnOrder.distances( MatrixColumnOrder.values( refs( genes ), phy ) );
        if ( Math.abs( d[ 0 ][ 1 ] ) > EPS ) {
            return fail( "X and Z agree wherever both were assessed, so their distance is 0 (R), got " + d[ 0 ][ 1 ] );
        }
        if ( Math.abs( d[ 0 ][ 2 ] - 8.0 ) > EPS ) {
            return fail( "X-Y must be R's 8.000000 (pairwise deletion, scaled 4/2); 5.657 would mean a blank was read as"
                    + " 0; got " + d[ 0 ][ 2 ] );
        }
        if ( Math.abs( d[ 1 ][ 2 ] - 5.656854 ) > EPS ) {
            return fail( "Z-Y (no blanks) must be R's 5.656854, got " + d[ 1 ][ 2 ] );
        }
        if ( !MatrixColumnOrder.clustered( refs( genes ), phy ).equals( refs( "Y", "X", "Z" ) ) ) {
            return fail( "the missing-value case must cluster in R's order 3 1 2 (Y X Z), got "
                    + MatrixColumnOrder.clustered( refs( genes ), phy ) );
        }
        return true;
    }

    /** Two genes no strain was assessed for together have NO distance: +Infinity (R: NA), so they join last. */
    private static boolean testNoOverlap() {
        final String[] genes = { "P", "Q", "Q2" };
        final Phylogeny phy = tree( genes, columns( new Integer[] { 1, 1, null, null }, new Integer[] { null, null, 2, 2 },
                                                    new Integer[] { null, null, 2, 3 } ) );
        final double[][] d = MatrixColumnOrder.distances( MatrixColumnOrder.values( refs( genes ), phy ) );
        if ( !Double.isInfinite( d[ 0 ][ 1 ] ) || !Double.isInfinite( d[ 0 ][ 2 ] ) ) {
            return fail( "a pair with no jointly assessed strain must be +Infinity, got " + d[ 0 ][ 1 ] + ", "
                    + d[ 0 ][ 2 ] );
        }
        final List<String> o = MatrixColumnOrder.clustered( refs( genes ), phy );
        if ( ( o.size() != 3 ) || !o.containsAll( refs( genes ) ) ) {
            return fail( "an unmeasurable pair must still leave every column in the order, got " + o );
        }
        if ( Math.abs( o.indexOf( "meta:Q" ) - o.indexOf( "meta:Q2" ) ) != 1 ) {
            return fail( "the two comparable genes must join first and sit together, got " + o );
        }
        return true;
    }

    // ---- clustered order ------------------------------------------------------------------------------------------

    private static boolean testClusteredOrderMatchesR() {
        // no ties -- R: 3 4 5 1 2
        final List<String> nt = MatrixColumnOrder.clustered( refs( NT_GENES ), tree( NT_GENES, NT ) );
        if ( !nt.equals( refs( "g3", "g4", "g5", "g1", "g2" ) ) ) {
            return fail( "no-ties case must match R's order 3 4 5 1 2, got " + nt );
        }
        // duplicate columns (zero-distance ties) -- R: cbind(a=c(1,2,3), b=c(1,2,3), c=c(3,2,1), d=c(1,2,3)) -> 3 4 1 2
        final String[] dup = { "a", "b", "c", "d" };
        final List<String> du = MatrixColumnOrder.clustered( refs( dup ),
                tree( dup, columns( new Integer[] { 1, 2, 3 }, new Integer[] { 1, 2, 3 }, new Integer[] { 3, 2, 1 },
                                    new Integer[] { 1, 2, 3 } ) ) );
        if ( !du.equals( refs( "c", "d", "a", "b" ) ) ) {
            return fail( "zero-distance ties must break as R does (order 3 4 1 2), got " + du );
        }
        // every pair equidistant -- R: diag(4) -> 4 3 1 2. The tie-break AND the leaf-order rules decide everything.
        final String[] eq = { "e1", "e2", "e3", "e4" };
        final List<String> ed = MatrixColumnOrder.clustered( refs( eq ),
                tree( eq, columns( new Integer[] { 1, 0, 0, 0 }, new Integer[] { 0, 1, 0, 0 }, new Integer[] { 0, 0, 1, 0 },
                                   new Integer[] { 0, 0, 0, 1 } ) ) );
        if ( !ed.equals( refs( "e4", "e3", "e1", "e2" ) ) ) {
            return fail( "an all-ties matrix must follow R's tie-break and hcass2 leaf order (4 3 1 2), got " + ed );
        }
        return true;
    }

    /**
     * The LINKAGE itself. On the small cases above a wrong linkage passes -- on 3-5 columns the dendrograms of
     * complete, average and single linkage coincide, so sabotaging the linkage survived them all. This 6x6 matrix was
     * found by searching random matrices in R (seed 20260918) for one with NO tied distances on which complete,
     * average, single and McQuitty linkage give four DIFFERENT orders -- so it isolates the linkage and nothing else.
     * R: complete 6 3 5 4 1 2; average 6 1 4 2 3 5; single 1 6 4 2 3 5; mcquitty 6 4 1 2 3 5.
     */
    private static boolean testLinkageIsComplete() {
        final String[] genes = { "h1", "h2", "h3", "h4", "h5", "h6" };
        final Phylogeny phy = tree( genes, columns( new Integer[] { 4, 4, 3, 3, 3, 2 }, new Integer[] { 1, 1, 4, 3, 3, 0 },
                                                    new Integer[] { 1, 3, 4, 1, 0, 0 }, new Integer[] { 0, 4, 4, 1, 1, 4 },
                                                    new Integer[] { 2, 0, 2, 0, 0, 1 },
                                                    new Integer[] { 4, 1, 0, 0, 2, 4 } ) );
        final List<String> o = MatrixColumnOrder.clustered( refs( genes ), phy );
        if ( !o.equals( refs( "h6", "h3", "h5", "h4", "h1", "h2" ) ) ) {
            return fail( "complete linkage must give R's order 6 3 5 4 1 2 (average/single/McQuitty each differ), got "
                    + o );
        }
        return true;
    }

    /**
     * Every data-driven mode starts from TABLE order, so the clustering does not depend on the order the columns
     * arrived in (the previous mode's). Fed alphabetically, the all-ties matrix above must still give R's order.
     */
    private static boolean testClusteredIgnoresIncomingOrder() {
        final String[] eq = { "e1", "e2", "e3", "e4" };
        final Phylogeny phy = tree( eq, columns( new Integer[] { 1, 0, 0, 0 }, new Integer[] { 0, 1, 0, 0 },
                                                 new Integer[] { 0, 0, 1, 0 }, new Integer[] { 0, 0, 0, 1 } ) );
        final List<String> scrambled = refs( "e3", "e1", "e4", "e2" );
        final List<String> o = MatrixColumnOrder.order( scrambled, MatrixColumnOrder.Mode.CLUSTERED, phy );
        if ( !o.equals( refs( "e4", "e3", "e1", "e2" ) ) ) {
            return fail( "CLUSTERED must not depend on the incoming order (the previous mode's), got " + o );
        }
        return true;
    }

    // ---- frequency, alphabetical, table ---------------------------------------------------------------------------

    /**
     * Mean over the ASSESSED strains, highest first. "hi" is 4,blank,blank: mean 4 beats "mid" 4,4,3 (3.67) -- if a
     * blank were averaged in as 0, "hi" would fall to 1.33 and come LAST. "none" has no value at all: last. "tie_a"
     * and "tie_b" share a mean and must keep their table order.
     */
    private static boolean testFrequency() {
        final String[] genes = { "none", "tie_a", "mid", "tie_b", "hi" };
        final Phylogeny phy = tree( genes, columns( new Integer[] { null, null, null }, new Integer[] { 2, 2, 2 },
                                                    new Integer[] { 4, 4, 3 }, new Integer[] { 1, 3, 2 },
                                                    new Integer[] { 4, null, null } ) );
        final List<String> o = MatrixColumnOrder.order( refs( genes ), MatrixColumnOrder.Mode.FREQUENCY, phy );
        if ( !o.equals( refs( "hi", "mid", "tie_a", "tie_b", "none" ) ) ) {
            return fail( "FREQUENCY must be mean-over-assessed, descending, ties in table order, no-value last; got "
                    + o );
        }
        return true;
    }

    /**
     * By the displayed name, ignoring case. displayName capitalises a property's first letter, so "aac6" is shown as
     * "Aac6" and cannot tell a case-sensitive sort from the right one (the first version of this test used it, and a
     * case-sensitive sort survived). "ibc" shows as "Ibc": ignoring case it sorts BEFORE "IS26" (b &lt; s), while an
     * ASCII sort puts "IS26" first (S &lt; b) -- so this pair decides it.
     */
    private static boolean testAlphabetical() {
        final String[] genes = { "rpoB", "IS26", "ibc", "Tet" };
        final Phylogeny phy = tree( genes, columns( new Integer[] { 1 }, new Integer[] { 1 }, new Integer[] { 1 },
                                                    new Integer[] { 1 } ) );
        final List<String> o = MatrixColumnOrder.order( refs( genes ), MatrixColumnOrder.Mode.ALPHABETICAL, phy );
        if ( !o.equals( refs( "ibc", "IS26", "rpoB", "Tet" ) ) ) {
            return fail( "ALPHABETICAL must ignore case (Ibc before IS26), got " + o );
        }
        return true;
    }

    /** The source's own order, whatever order the columns are currently in; an unplaceable ref is kept, at the end. */
    private static boolean testTableOrder() {
        final String[] genes = { "zeta", "alpha", "mu" };
        final Phylogeny phy = tree( genes, columns( new Integer[] { 1 }, new Integer[] { 2 }, new Integer[] { 3 } ) );
        final List<String> incoming = new ArrayList<String>( refs( "mu", "alpha", "zeta" ) );
        incoming.add( 1, "tax:scientific_name" ); // an element-slot ref: a column, but never a node property
        final List<String> o = MatrixColumnOrder.order( incoming, MatrixColumnOrder.Mode.TABLE, phy );
        final List<String> want = new ArrayList<String>( refs( "zeta", "alpha", "mu" ) );
        want.add( "tax:scientific_name" );
        if ( !o.equals( want ) ) {
            return fail( "TABLE must restore the source order and keep an unplaceable ref at the end, got " + o );
        }
        return true;
    }

    // ---- manual, edges, apply ---------------------------------------------------------------------------------------

    private static boolean testManualAndEdges() {
        final Phylogeny phy = tree( NT_GENES, NT );
        final List<String> mine = refs( "g5", "g1", "g4", "g2", "g3" );
        if ( !MatrixColumnOrder.order( mine, MatrixColumnOrder.Mode.MANUAL, phy ).equals( mine ) ) {
            return fail( "MANUAL must keep the user's order exactly" );
        }
        for( final MatrixColumnOrder.Mode m : MatrixColumnOrder.Mode.values() ) {
            if ( !MatrixColumnOrder.order( new ArrayList<String>(), m, phy ).isEmpty() ) {
                return fail( m + " of no columns must be empty" );
            }
            if ( !MatrixColumnOrder.order( refs( "g2" ), m, phy ).equals( refs( "g2" ) ) ) {
                return fail( m + " of one column must be that column" );
            }
        }
        if ( MatrixColumnOrder.DEFAULT != MatrixColumnOrder.Mode.CLUSTERED ) {
            return fail( "the default mode must be CLUSTERED" );
        }
        if ( MatrixColumnOrder.completeLinkageOrder( new double[ 0 ][ 0 ] ).length != 0 ) {
            return fail( "an empty distance matrix has an empty order" );
        }
        return true;
    }

    /**
     * Only the MATRIX columns move, and only among the slots they occupy: a colour strip at index 1 and a symbol at
     * index 3 must stay exactly there. A spec list with fewer than two MATRIX columns, MANUAL, and null are unchanged.
     */
    private static boolean testApplyIsSlotPreserving() {
        final Phylogeny phy = tree( NT_GENES, NT );
        final AnnotationColumns.ColumnSpec strip = new AnnotationColumns.ColumnSpec( "data:host",
                                                                                     AnnotationColumns.Type.COLOR_STRIP );
        final AnnotationColumns.ColumnSpec sym = new AnnotationColumns.ColumnSpec( "data:flag",
                                                                                   AnnotationColumns.Type.SYMBOL );
        final List<AnnotationColumns.ColumnSpec> specs = new ArrayList<AnnotationColumns.ColumnSpec>();
        specs.add( mx( "g1" ) );
        specs.add( strip );
        specs.add( mx( "g2" ) );
        specs.add( sym );
        specs.add( mx( "g3" ) );
        specs.add( mx( "g4" ) );
        specs.add( mx( "g5" ) );
        final List<AnnotationColumns.ColumnSpec> out = MatrixColumnOrder.apply( specs, MatrixColumnOrder.Mode.CLUSTERED,
                                                                                phy );
        if ( ( out.size() != 7 ) || ( out.get( 1 ) != strip ) || ( out.get( 3 ) != sym ) ) {
            return fail( "non-MATRIX columns must keep their exact slots (1 and 3)" );
        }
        final List<String> matrix = new ArrayList<String>();
        for( final int slot : new int[] { 0, 2, 4, 5, 6 } ) {
            if ( out.get( slot )._type != AnnotationColumns.Type.MATRIX ) {
                return fail( "slot " + slot + " must still hold a MATRIX column" );
            }
            matrix.add( out.get( slot )._ref );
        }
        if ( !matrix.equals( refs( "g3", "g4", "g5", "g1", "g2" ) ) ) {
            return fail( "the MATRIX slots must hold R's clustered order g3 g4 g5 g1 g2, got " + matrix );
        }
        if ( specs.get( 0 )._ref.equals( "meta:g3" ) ) {
            return fail( "apply must not mutate the list it was given" );
        }
        if ( MatrixColumnOrder.apply( specs, MatrixColumnOrder.Mode.MANUAL, phy ) != specs ) {
            return fail( "MANUAL must return the list unchanged" );
        }
        final List<AnnotationColumns.ColumnSpec> one = Arrays.asList( strip, mx( "g1" ), sym );
        if ( MatrixColumnOrder.apply( one, MatrixColumnOrder.Mode.ALPHABETICAL, phy ) != one ) {
            return fail( "fewer than two MATRIX columns: nothing to order, list unchanged" );
        }
        if ( MatrixColumnOrder.apply( null, MatrixColumnOrder.Mode.CLUSTERED, phy ) != null ) {
            return fail( "null specs stay null" );
        }
        return true;
    }

    /**
     * The Annotation Fields dialog's OK ({@link MatrixColumnOrder#resolveEdited}). Each case is paired with the one
     * that differs only by the thing under test, and must come out the OTHER way:
     * a move that changes the matrix order -> Manual, the user's order kept; the same list with NO move -> the mode
     * re-sorts it; a move that leaves the matrix order as the mode has it (moved and back, or only a colour strip
     * moved) -> the mode stays.
     */
    private static boolean testDialogResolution() {
        final Phylogeny phy = tree( NT_GENES, NT );
        final AnnotationColumns.ColumnSpec strip = new AnnotationColumns.ColumnSpec( "data:host",
                                                                                     AnnotationColumns.Type.COLOR_STRIP );
        // the clustered order is g3 g4 g5 g1 g2; the user drags g1 to the front
        final List<AnnotationColumns.ColumnSpec> dragged = Arrays.asList( mx( "g1" ), mx( "g3" ), mx( "g4" ), mx( "g5" ),
                                                                          mx( "g2" ), strip );
        final MatrixColumnOrder.Resolved moved = MatrixColumnOrder.resolveEdited( dragged, true,
                MatrixColumnOrder.Mode.CLUSTERED, phy );
        if ( ( moved.mode() != MatrixColumnOrder.Mode.MANUAL ) || !moved.specs().equals( dragged ) ) {
            return fail( "a move that changes the matrix order must win and make the tab Manual, got " + moved.mode() );
        }
        final MatrixColumnOrder.Resolved not_moved = MatrixColumnOrder.resolveEdited( dragged, false,
                MatrixColumnOrder.Mode.CLUSTERED, phy );
        if ( ( not_moved.mode() != MatrixColumnOrder.Mode.CLUSTERED )
                || !MatrixColumnOrder.matrixRefs( not_moved.specs() ).equals( refs( "g3", "g4", "g5", "g1", "g2" ) ) ) {
            return fail( "the SAME list without a move must be re-sorted by the mode (e.g. newly added fields), got "
                    + not_moved.mode() + " " + MatrixColumnOrder.matrixRefs( not_moved.specs() ) );
        }
        // moved, but the matrix order is exactly the mode's (moved and back, or only the strip moved)
        final List<AnnotationColumns.ColumnSpec> strip_moved = Arrays.asList( strip, mx( "g3" ), mx( "g4" ), mx( "g5" ),
                                                                              mx( "g1" ), mx( "g2" ) );
        final MatrixColumnOrder.Resolved harmless = MatrixColumnOrder.resolveEdited( strip_moved, true,
                MatrixColumnOrder.Mode.CLUSTERED, phy );
        if ( harmless.mode() != MatrixColumnOrder.Mode.CLUSTERED ) {
            return fail( "a move that leaves the matrix order as the mode has it must NOT flip the tab to Manual" );
        }
        final MatrixColumnOrder.Resolved manual = MatrixColumnOrder.resolveEdited( dragged, false,
                MatrixColumnOrder.Mode.MANUAL, phy );
        if ( ( manual.mode() != MatrixColumnOrder.Mode.MANUAL ) || !manual.specs().equals( dragged ) ) {
            return fail( "a Manual tab keeps the dialog's order even without a move" );
        }
        return true;
    }

    private static boolean testDeterministic() {
        final Phylogeny phy = tree( NT_GENES, NT );
        for( final MatrixColumnOrder.Mode m : MatrixColumnOrder.Mode.values() ) {
            final List<String> a = MatrixColumnOrder.order( refs( NT_GENES ), m, phy );
            final List<String> b = MatrixColumnOrder.order( refs( NT_GENES ), m, phy );
            if ( !a.equals( b ) ) {
                return fail( m + " must be deterministic: " + a + " vs " + b );
            }
        }
        return true;
    }

    private static AnnotationColumns.ColumnSpec mx( final String gene ) {
        return new AnnotationColumns.ColumnSpec( "meta:" + gene, AnnotationColumns.Type.MATRIX );
    }
}

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
 * {@link MatrixColumnOrder}: the six modes of View &gt; Order Matrix Columns.
 * <p>
 * The clustering expectations are R's OWN output, generated deliberately with R 4.5.3 --
 * {@code hclust(dist(t(m)), method = "complete")} for Euclidean and {@code vegan::vegdist(t(m), method = "bray")}
 * (vegan 2.7-2) for Bray-Curtis, on the same small matrices -- not restated here in this test's own
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
                && testManualAndEdges() && testApplyIsSlotPreserving() && testDialogResolution() && testDeterministic()
                && testBrayCurtisMatchesVegan() && testDoubleZeroIsNotAgreement() && testBrayCurtisPairwiseDeletion()
                && testBrayCurtisIsSorensenOnBinary() && testBrayCurtisEdges();
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
        // the same for the other clustered mode: every pair of these columns shares no presence, so Bray-Curtis ties
        // them all at 1.0 and only the tie-break and leaf order are left -- R's vegdist + hclust give 4 3 1 2 too
        final List<String> b = MatrixColumnOrder.order( scrambled, MatrixColumnOrder.Mode.CLUSTERED_PRESENCE, phy );
        if ( !b.equals( refs( "e4", "e3", "e1", "e2" ) ) ) {
            return fail( "CLUSTERED_PRESENCE must not depend on the incoming order either, got " + b );
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

    // ---- Bray-Curtis: clustering that ignores shared absence ---------------------------------------------------------

    /**
     * The double-zero fixture, 6 genes over 8 strains, built so that the two Clustered modes MUST disagree: two core
     * genes (h1, h2, present everywhere), two clade genes on disjoint halves (h3, h4), and two RARE genes each present
     * in one strain only -- and in DIFFERENT strains (h5 in the first, h6 in the last). h5 and h6 share no strain they
     * are both present in, and agree only by being absent.
     * <p>
     * R: {@code cbind(h1=..., h2=..., h3=..., h4=..., h5=..., h6=...)} of the columns below.
     */
    private static final String[]    DZ_GENES = { "h1", "h2", "h3", "h4", "h5", "h6" };
    private static final Integer[][] DZ       = columns( new Integer[] { 4, 4, 4, 4, 4, 4, 4, 4 },
                                                         new Integer[] { 4, 4, 3, 4, 4, 4, 4, 4 },
                                                         new Integer[] { 4, 4, 4, 4, 0, 0, 0, 0 },
                                                         new Integer[] { 0, 0, 0, 0, 4, 4, 4, 4 },
                                                         new Integer[] { 4, 0, 0, 0, 0, 0, 0, 0 },
                                                         new Integer[] { 0, 0, 0, 0, 0, 0, 0, 4 } );

    /** Every pair of the double-zero fixture, as R's {@code vegdist(t(m), method = "bray")} gives it. */
    private static boolean testBrayCurtisMatchesVegan() {
        final double[][] r = {
                { 0.0, 0.015873015873, 0.333333333333, 0.333333333333, 0.777777777778, 0.777777777778 },
                { 0.015873015873, 0.0, 0.361702127660, 0.319148936170, 0.771428571429, 0.771428571429 },
                { 0.333333333333, 0.361702127660, 0.0, 1.0, 0.6, 1.0 },
                { 0.333333333333, 0.319148936170, 1.0, 0.0, 1.0, 0.6 },
                { 0.777777777778, 0.771428571429, 0.6, 1.0, 0.0, 1.0 },
                { 0.777777777778, 0.771428571429, 1.0, 0.6, 1.0, 0.0 } };
        final double[][] d = MatrixColumnOrder
                .brayCurtisDistances( MatrixColumnOrder.values( refs( DZ_GENES ), tree( DZ_GENES, DZ ) ) );
        for( int a = 0; a < 6; ++a ) {
            for( int b = 0; b < 6; ++b ) {
                if ( Math.abs( d[ a ][ b ] - r[ a ][ b ] ) > EPS ) {
                    return fail( "Bray-Curtis " + DZ_GENES[ a ] + "-" + DZ_GENES[ b ] + " should be vegdist's "
                            + r[ a ][ b ] + ", got " + d[ a ][ b ] );
                }
            }
            if ( d[ a ][ a ] != 0.0 ) {
                return fail( "a column is at distance 0 from itself, got " + d[ a ][ a ] );
            }
        }
        return true;
    }

    /**
     * What the mode is FOR. On the same fixture the two distances disagree about the two rare genes, in opposite
     * directions: Euclidean counts the six strains where both h5 and h6 are absent as agreement, so they become each
     * other's NEAREST column (5.657, against 6.93 and more to everything else); Bray-Curtis drops those strains, so
     * they are each other's FARTHEST (1.0 -- they are never present together) and each one's nearest is the clade gene
     * it actually occurs with (0.6). The whole orders differ too: R gives {@code 4 1 2 3 5 6} for Euclidean and
     * {@code 3 5 6 4 1 2} for Bray-Curtis.
     */
    private static boolean testDoubleZeroIsNotAgreement() {
        final Phylogeny phy = tree( DZ_GENES, DZ );
        final Double[][] v = MatrixColumnOrder.values( refs( DZ_GENES ), phy );
        if ( nearest( MatrixColumnOrder.distances( v ), 4 ) != 5 ) {
            return fail( "Euclidean must make the two rare genes each other's nearest (the double-zero problem this "
                    + "mode exists for); if it no longer does, the fixture stopped isolating the thing under test" );
        }
        final double[][] bc = MatrixColumnOrder.brayCurtisDistances( v );
        if ( nearest( bc, 4 ) != 2 ) {
            return fail( "Bray-Curtis must put h5 nearest the clade gene it occurs with (h3), got index "
                    + nearest( bc, 4 ) );
        }
        if ( nearest( bc, 5 ) != 3 ) {
            return fail( "Bray-Curtis must put h6 nearest h4, got index " + nearest( bc, 5 ) );
        }
        if ( ( bc[ 4 ][ 5 ] != 1.0 ) || ( bc[ 3 ][ 2 ] != 1.0 ) ) {
            return fail( "columns never present in the same tip are maximally distant (1.0), got " + bc[ 4 ][ 5 ] );
        }
        final List<String> euclid = MatrixColumnOrder.order( refs( DZ_GENES ), MatrixColumnOrder.Mode.CLUSTERED, phy );
        if ( !euclid.equals( refs( "h4", "h1", "h2", "h3", "h5", "h6" ) ) ) {
            return fail( "CLUSTERED must match R's Euclidean order 4 1 2 3 5 6, got " + euclid );
        }
        final List<String> bray = MatrixColumnOrder.order( refs( DZ_GENES ), MatrixColumnOrder.Mode.CLUSTERED_PRESENCE,
                                                           phy );
        if ( !bray.equals( refs( "h3", "h5", "h6", "h4", "h1", "h2" ) ) ) {
            return fail( "CLUSTERED_PRESENCE must match R's Bray-Curtis order 3 5 6 4 1 2, got " + bray );
        }
        return true;
    }

    /** The index of the column closest to {@code g}, itself excluded. */
    private static int nearest( final double[][] d, final int g ) {
        int best = -1;
        for( int k = 0; k < d.length; ++k ) {
            if ( ( k != g ) && ( ( best < 0 ) || ( d[ g ][ k ] < d[ g ][ best ] ) ) ) {
                best = k;
            }
        }
        return best;
    }

    /**
     * Missing cells are deleted pairwise (vegan's {@code na.rm = TRUE}), not read as 0 and not scaled up: being a
     * ratio, Bray-Curtis needs no scaling. The neighbouring case is the SAME matrix with its four blanks filled with
     * 0, which R makes a different number on five of the six pairs -- so these values can only come from the right
     * rule. (The two orders happen to coincide on this matrix, which is why the VALUES are what is pinned.)
     */
    private static boolean testBrayCurtisPairwiseDeletion() {
        final String[] genes = { "g1", "g2", "g3", "g4" };
        final Integer[][] m = columns( new Integer[] { 4, 4, 0, 2, null, 3 }, new Integer[] { 4, null, 3, 0, 1, 2 },
                                       new Integer[] { 0, 0, null, 1, 4, 2 }, new Integer[] { 2, 0, 4, null, 3, 1 } );
        final double[][] d = MatrixColumnOrder
                .brayCurtisDistances( MatrixColumnOrder.values( refs( genes ), tree( genes, m ) ) );
        // vegdist(t(m), method = "bray", na.rm = TRUE); the same matrix with 0 for NA gives the second column
        final double[][] pairs = { { 0, 1, 0.333333333333, 0.478260869565 }, { 0, 2, 0.625, 0.7 },
                { 0, 3, 0.666666666667, 0.739130434783 }, { 1, 2, 0.571428571429, 0.647058823529 },
                { 2, 3, 0.333333333333, 0.529411764706 } };
        for( final double[] p : pairs ) {
            final double got = d[ ( int ) p[ 0 ] ][ ( int ) p[ 1 ] ];
            if ( Math.abs( got - p[ 2 ] ) > EPS ) {
                return fail( "pairwise deletion: " + genes[ ( int ) p[ 0 ] ] + "-" + genes[ ( int ) p[ 1 ] ]
                        + " must be vegdist's " + p[ 2 ] + " (filling the blanks with 0 would give " + p[ 3 ]
                        + "), got " + got );
            }
        }
        if ( Math.abs( d[ 1 ][ 3 ] - 0.3 ) > EPS ) {
            return fail( "g2-g4, the one pair zero-filling does NOT change, must still be 0.3, got " + d[ 1 ][ 3 ] );
        }
        return true;
    }

    /**
     * On 0/1 data Bray-Curtis IS the Sorensen-Dice dissimilarity, {@code 1 - 2|A and B| / (|A| + |B|)} -- the reason
     * it is the right distance for presence/absence. Checked twice over: against vegdist's numbers for this matrix
     * (R 4.5.3, {@code set.seed(7)}), and against the Sorensen formula computed here from SET COUNTS, which is a
     * different formula rather than this test restating the one under test.
     */
    private static boolean testBrayCurtisIsSorensenOnBinary() {
        final String[] genes = { "b1", "b2", "b3", "b4", "b5" };
        final Integer[][] bin = columns( new Integer[] { 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 },
                                         new Integer[] { 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1 },
                                         new Integer[] { 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1 },
                                         new Integer[] { 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1 },
                                         new Integer[] { 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0 } );
        final double[][] d = MatrixColumnOrder
                .brayCurtisDistances( MatrixColumnOrder.values( refs( genes ), tree( genes, bin ) ) );
        final double[] vegan = { 0.75, 0.75, 1.0, 0.555555555556, 0.4, 0.454545454545, 0.636363636364, 0.272727272727,
                0.454545454545, 0.5 };
        int k = 0;
        for( int a = 0; a < 5; ++a ) {
            for( int b = a + 1; b < 5; ++b ) {
                if ( Math.abs( d[ a ][ b ] - vegan[ k ] ) > EPS ) {
                    return fail( "binary " + genes[ a ] + "-" + genes[ b ] + " must be vegdist's " + vegan[ k ]
                            + ", got " + d[ a ][ b ] );
                }
                int shared = 0;
                int in_a = 0;
                int in_b = 0;
                for( int t = 0; t < bin.length; ++t ) {
                    in_a += bin[ t ][ a ];
                    in_b += bin[ t ][ b ];
                    shared += ( ( bin[ t ][ a ] == 1 ) && ( bin[ t ][ b ] == 1 ) ) ? 1 : 0;
                }
                final double sorensen = 1.0 - ( ( 2.0 * shared ) / ( in_a + in_b ) );
                if ( Math.abs( d[ a ][ b ] - sorensen ) > EPS ) {
                    return fail( "on 0/1 data Bray-Curtis must equal the Sorensen-Dice dissimilarity " + sorensen
                            + " for " + genes[ a ] + "-" + genes[ b ] + ", got " + d[ a ][ b ] );
                }
                ++k;
            }
        }
        return true;
    }

    /**
     * The three cases vegdist cannot answer, each decided here so clustering never sees a NaN: no jointly assessed
     * tip at all (vegdist: NA) is +Infinity, so the pair joins last; two columns that are 0 at every tip they share
     * (vegdist: NaN, "empty rows") are identical wherever they can be compared, so 0; and a pair whose values are
     * signed and cancel to a non-positive total -- Bray-Curtis is meant for values that are 0 or more -- is
     * +Infinity if the columns differ at all, never a negative distance or a NaN.
     */
    private static boolean testBrayCurtisEdges() {
        final String[] genes = { "P", "Q", "Q2" };
        final Phylogeny phy = tree( genes, columns( new Integer[] { 1, 1, null, null },
                                                    new Integer[] { null, null, 2, 2 },
                                                    new Integer[] { null, null, 2, 3 } ) );
        final double[][] d = MatrixColumnOrder.brayCurtisDistances( MatrixColumnOrder.values( refs( genes ), phy ) );
        if ( !Double.isInfinite( d[ 0 ][ 1 ] ) || !Double.isInfinite( d[ 0 ][ 2 ] ) ) {
            return fail( "a pair with no jointly assessed tip must be +Infinity, got " + d[ 0 ][ 1 ] );
        }
        final List<String> o = MatrixColumnOrder.clusteredByPresence( refs( genes ), phy );
        if ( ( o.size() != 3 ) || !o.containsAll( refs( genes ) ) ) {
            return fail( "an unmeasurable pair must still leave every column in the order, got " + o );
        }
        final String[] zero = { "z1", "z2", "z3" };
        final double[][] dz = MatrixColumnOrder.brayCurtisDistances( MatrixColumnOrder
                .values( refs( zero ), tree( zero, columns( new Integer[] { 0, 0, 0 }, new Integer[] { 0, 0, 0 },
                                                            new Integer[] { 0, 1, 0 } ) ) ) );
        if ( dz[ 0 ][ 1 ] != 0.0 ) {
            return fail( "two columns that are 0 at every shared tip are identical there: distance 0, got "
                    + dz[ 0 ][ 1 ] );
        }
        if ( dz[ 0 ][ 2 ] != 1.0 ) {
            return fail( "an all-zero column against one with a single 1 shares no presence: 1.0, got " + dz[ 0 ][ 2 ] );
        }
        final String[] neg = { "n1", "n2" };
        final double[][] dn = MatrixColumnOrder.brayCurtisDistances( MatrixColumnOrder
                .values( refs( neg ), tree( neg, columns( new Integer[] { 2, -2 }, new Integer[] { -2, 2 } ) ) ) );
        if ( !Double.isInfinite( dn[ 0 ][ 1 ] ) ) {
            return fail( "a non-positive total with columns that differ must be +Infinity, never negative or NaN, got "
                    + dn[ 0 ][ 1 ] );
        }
        return true;
    }

    private static AnnotationColumns.ColumnSpec mx( final String gene ) {
        return new AnnotationColumns.ColumnSpec( "meta:" + gene, AnnotationColumns.Type.MATRIX );
    }
}

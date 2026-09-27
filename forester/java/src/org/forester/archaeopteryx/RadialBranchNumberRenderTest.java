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

import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Where the support and branch-length NUMBERS land in the two radial layouts (circular and unrooted).
 * <p>
 * Three properties, all measured off the rendered pixels rather than off the placement arithmetic:
 * <ol>
 * <li><b>Nothing is drawn across its own branch.</b> The headline defect: the numbers used to ride the branch at a
 * fixed +/-2 px BASELINE offset, flipped in sign on the far half of the fan to keep them on the same visual side --
 * but a baseline is the BOTTOM of the text, so the flipped sign put the whole ascent back across the branch. 6 of 6
 * circular branches were struck through by their own spoke.</li>
 * <li><b>The two numbers sit on OPPOSITE sides</b> of the branch, the same way round on every branch -- the
 * rectangular layouts' convention (length above, support below), carried into the fan.</li>
 * <li><b>Both are actually drawn.</b> Deliberately paired with (1), because (1) alone is satisfied by drawing
 * nothing at all -- and that is not a hypothetical: claiming the two sides as two separate boxes made a branch's
 * own support refuse its own length (the occupancy map is axis-aligned in DEVICE space, so on a diagonal branch the
 * two bounds overlap although the boxes do not), and every support value on the canvas disappeared while "nothing
 * overprints" stayed green.</li>
 * </ol>
 * And a fourth, on a DENSE fan: <b>a number never lies across a branch that is not its own</b>. In the unrooted
 * layout sibling branches leave one node a few degrees apart, so a number that outreaches its short branch lay over
 * its sibling's or back over the parent's (3 of 16 drawn numbers on the bat phylogeny, 4 of 25 on the animal tree of
 * life, more in circular). A branch is now an obstacle like another number is ({@link BranchObstacles}), under the
 * same Auto-hide switch -- so the fan is checked both ways: with auto-hide on, no number ink lands on tree ink and
 * some numbers are refused for exactly that reason while others are still drawn; with auto-hide off, the refusals
 * stop and the ink lands on the tree again, which proves both the gate and that the fixture really is crowded.
 * <p>
 * Each number is isolated by rendering the same layout with it alone: the pixel difference against a render with
 * neither is exactly that number's ink, so nothing here depends on guessing a colour. Headful; a green no-op when
 * headless.
 */
public final class RadialBranchNumberRenderTest {

    /** Ink further than this from a branch midpoint belongs to some other branch. */
    private static final int    RADIUS       = 30;
    /** A tree pixel counts as INK when any channel is below this: a line's core, not its antialiasing fringe. A
     *  number placed flush against its own branch can share one outermost fringe pixel with it (measured on the
     *  96-spoke fan: tree #f1f1f1 under number #dddddd, one pixel), which is not "drawn across a branch". */
    private static final int    INK_MAX      = 0xE8;
    /** Minimum white the ink must leave between itself and the branch line it labels. */
    private static final double MIN_CLEAR_PX = 1.0;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RadialBranchNumberRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            return radialNumbersOk();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "threw: " + e );
        }
    }

    private static boolean radialNumbersOk() throws Exception {
        final Phylogeny phy = fanTree();
        final Configuration conf = new Configuration();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "radialnum" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final MainFrame frame = mf[ 0 ];
            try {
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                // the support SYMBOL rides the branch line itself; it is a separate feature with its own occupancy
                // map, and its ink would be counted as the numbers' if it were drawn here
                frame.getOptions().setSupportVisualization( Options.SUPPORT_VISUALIZATION.NONE );
                tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                frame.showWhole();
                tp.getControlPanel().setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                for( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] {
                        PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                    tp.setPhylogenyGraphicsType( gt );
                    checkLayout( tp, phy, gt, ok );
                }
                // (4) the dense fan: numbers that would cross a sibling's or the parent's branch
                final Phylogeny dense = denseStar();
                tp.setTree( dense );
                tp.recalculateMaxDistanceToRoot();
                for( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] {
                        PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                    tp.setPhylogenyGraphicsType( gt );
                    checkDenseFan( tp, gt, ok, 900 );
                }
                // (5) short inner legs in CIRCULAR: a number on a leg a few px long reaches inward past the parent's
                // radius, where the parent's ARC runs -- an arc is a branch line too. Circular only: the unrooted
                // layout has no arcs, and this tree is not crowded there.
                tp.setTree( shortLegs() );
                tp.recalculateMaxDistanceToRoot();
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                checkDenseFan( tp, PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, ok, 900 );
                // (6) the demo pair, measured where the paint is: the file that says it fits must draw every number
                // and refuse none; the one that says it crowds must refuse some, keep some, and cross nothing
                checkDemoPair( tp, ok );
                checkInkBand( tp, ok );
                // (the search-hit policy -- a hit's number is NOT exempt -- is pinned in CrowdedBranchDataTest, in
                // the rectangular layout: here every refused number is refused by a BRANCH first, so a fixture on
                // this fan could not tell an exempting funnel from a faithful one)
            }
            catch ( final Exception e ) {
                e.printStackTrace();
                fail( ok, "threw: " + e );
            }
        } );
        return ok[ 0 ];
    }

    private static void checkLayout( final TreePanel tp, final Phylogeny phy, final PHYLOGENY_GRAPHICS_TYPE gt,
                                     final boolean[] ok ) {
        final int w = 1300; // big enough that the branch midpoints separate -- asserted below, not assumed
        final BufferedImage none = render( tp, w, false, false );
        final BufferedImage len_only = render( tp, w, false, true );
        final BufferedImage sup_only = render( tp, w, true, false );
        final BufferedImage both = render( tp, w, true, true );
        if ( ( none.getWidth() != both.getWidth() ) || ( none.getHeight() != both.getHeight() )
                || ( len_only.getWidth() != both.getWidth() ) || ( sup_only.getWidth() != both.getWidth() ) ) {
            fail( ok, gt + ": the four renders must share a size (turning a number on must not move the tree)" );
            return;
        }

        // (1) no number may be painted over ink that was already on the canvas -- which, with blank tip labels and
        // no support symbols, is the tree itself. This is the defect, and it is whole-image: a number drawn across
        // ANY branch counts, not only across its own.
        int drawn = 0;
        for( int y = 0; y < both.getHeight(); ++y ) {
            for( int x = 0; x < both.getWidth(); ++x ) {
                if ( both.getRGB( x, y ) != none.getRGB( x, y ) ) {
                    ++drawn;
                }
            }
        }
        final int over_ink = inkOnInk( both, none );
        if ( over_ink > 0 ) {
            fail( ok, gt + ": a branch number is drawn over the tree (" + over_ink + " of " + drawn
                    + " number pixels land on ink)" );
        }
        if ( drawn < 200 ) {
            fail( ok, gt + ": the branch numbers must be drawn at all (" + drawn + " number pixels)" );
        }

        // Fixture precondition: ink is attributed to the nearest branch midpoint, so the midpoints must be further
        // apart than the radius that attribution uses -- otherwise one branch's numbers could answer for another's
        // and every per-branch check below would pin nothing.
        final List<PhylogenyNode> branches = internalBranches( phy );
        if ( branches.size() < 6 ) {
            fail( ok, gt + ": the fixture must offer several branches in different directions (" + branches.size()
                    + ")" );
            return;
        }
        for( int i = 0; i < branches.size(); ++i ) {
            for( int j = i + 1; j < branches.size(); ++j ) {
                final double[] a = midAndAngle( branches.get( i ), phy, gt );
                final double[] b = midAndAngle( branches.get( j ), phy, gt );
                if ( Math.hypot( a[ 0 ] - b[ 0 ], a[ 1 ] - b[ 1 ] ) <= ( 2 * RADIUS ) ) {
                    fail( ok, gt + ": fixture too crowded to attribute ink to a branch (two midpoints "
                            + Math.round( Math.hypot( a[ 0 ] - b[ 0 ], a[ 1 ] - b[ 1 ] ) ) + " px apart)" );
                    return;
                }
            }
        }

        // (2) + (3): per branch, on which side of its own line does each number's ink sit, and how close does it
        // come? Measured as a SIGNED perpendicular offset from the branch line, in the frame of the branch's own
        // OUTWARD direction -- so "the same side" means the same side going round the fan, which is the property
        // the far-half flip exists to keep and the thing a fixed baseline offset cannot express.
        int len_side = 0, sup_side = 0;
        boolean seen = false;
        for( final PhylogenyNode node : branches ) {
            final double[] ma = midAndAngle( node, phy, gt );
            final double[] len = perpendicularSpan( len_only, none, ma );
            final double[] sup = perpendicularSpan( sup_only, none, ma );
            if ( ( len == null ) || ( sup == null ) ) {
                fail( ok, gt + ": branch at " + Math.round( Math.toDegrees( ma[ 2 ] ) )
                        + " deg must carry BOTH its length and its support (length ink=" + ( len != null )
                        + ", support ink=" + ( sup != null ) + ")" );
                continue;
            }
            // each number wholly on one side: its nearest and furthest edge agree in sign
            if ( ( Math.signum( len[ 0 ] ) != Math.signum( len[ 1 ] ) )
                    || ( Math.signum( sup[ 0 ] ) != Math.signum( sup[ 1 ] ) ) ) {
                fail( ok, gt + ": branch at " + Math.round( Math.toDegrees( ma[ 2 ] ) )
                        + " deg has a number straddling its own line (length spans " + len[ 0 ] + ".." + len[ 1 ]
                        + ", support spans " + sup[ 0 ] + ".." + sup[ 1 ] + ")" );
                continue;
            }
            final int ls = (int) Math.signum( len[ 0 ] ), ss = (int) Math.signum( sup[ 0 ] );
            if ( ls == ss ) {
                fail( ok, gt + ": the length and the support must take OPPOSITE sides of the branch (both on side "
                        + ls + " at " + Math.round( Math.toDegrees( ma[ 2 ] ) ) + " deg)" );
            }
            if ( !seen ) {
                len_side = ls;
                sup_side = ss;
                seen = true;
            }
            else if ( ( ls != len_side ) || ( ss != sup_side ) ) {
                fail( ok, gt + ": every branch must put its length on the same side going round the fan (branch at "
                        + Math.round( Math.toDegrees( ma[ 2 ] ) ) + " deg has length on side " + ls + ", the first on "
                        + len_side + ")" );
            }
            final double clear = Math.min( Math.abs( len[ 0 ] ), Math.abs( sup[ 0 ] ) );
            if ( clear < MIN_CLEAR_PX ) {
                fail( ok, gt + ": a number touches its branch at " + Math.round( Math.toDegrees( ma[ 2 ] ) )
                        + " deg (nearest ink " + clear + " px off the line)" );
            }
        }

        // (4) the space RESERVED is the space the ink uses. Checked with one number shown at a time, which is the
        // case that can tell them apart: with both drawn the pair's box straddles the branch and is very nearly
        // centred on it anyway, so a reservation that ignored the offset would look right. With one, the box sits
        // wholly on that number's side, and a claim centred on the branch would reserve the EMPTY side and leave
        // the ink unprotected -- crowding would then let a neighbour draw over it.
        renderRecording( tp, w, false, true, ok, gt, "branch length", branches, phy, none );
        renderRecording( tp, w, true, false, ok, gt, "support", branches, phy, none );
    }

    /**
     * On a fan too dense for its numbers: with auto-hide ON no number ink may land on tree ink, at least one number
     * must have been refused for a branch, and some must still be drawn (the rule hides what crosses, not
     * everything); with auto-hide OFF nothing is refused and the ink lands on the tree again.
     */
    private static void checkDenseFan( final TreePanel tp, final PHYLOGENY_GRAPHICS_TYPE gt, final boolean[] ok,
                                       final int w ) {
        final ControlPanel cp = tp.getControlPanel();
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        final BufferedImage none_on = render( tp, w, false, false );
        final BufferedImage both_on = render( tp, w, true, true );
        final int refused_on = tp.numbersRefusedForBranchesForTest();
        final int drawn_on = tp.numbersDrawnForTest();
        final int over_on = inkOnInk( both_on, none_on );
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, false );
        final BufferedImage none_off = render( tp, w, false, false );
        final BufferedImage both_off = render( tp, w, true, true );
        final int refused_off = tp.numbersRefusedForBranchesForTest();
        final int over_off = inkOnInk( both_off, none_off );
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        if ( over_off == 0 ) {
            fail( ok, gt + ": precondition -- with auto-hide off the dense fan's numbers must land on its branches"
                    + " (0 number pixels on ink), so the fixture is not dense enough to show anything" );
            return;
        }
        if ( refused_off != 0 ) {
            fail( ok, gt + ": with auto-hide off nothing may be refused for a branch (" + refused_off + " were)" );
        }
        if ( over_on != 0 ) {
            fail( ok, gt + ": with auto-hide on a number is still drawn across a branch (" + over_on
                    + " number pixels on tree ink; " + over_off + " with auto-hide off) at "
                    + inkOnInkWhere( both_on, none_on ) );
        }
        if ( refused_on == 0 ) {
            fail( ok, gt + ": the dense fan must make the rule refuse at least one number for a branch" );
        }
        if ( drawn_on == 0 ) {
            fail( ok, gt + ": the rule must hide what crosses, not everything (0 numbers drawn, " + refused_on
                    + " refused)" );
        }
    }

    /** The first few such pixels, with the tree ink's colour there -- a solid line reads very differently from a
     *  faint antialiasing fringe. */
    private static String inkOnInkWhere( final BufferedImage with, final BufferedImage without ) {
        final StringBuilder sb = new StringBuilder();
        int n = 0;
        for( int y = 0; ( y < with.getHeight() ) && ( n < 5 ); ++y ) {
            for( int x = 0; ( x < with.getWidth() ) && ( n < 5 ); ++x ) {
                if ( ( with.getRGB( x, y ) != without.getRGB( x, y ) ) && isInk( without.getRGB( x, y ) ) ) {
                    sb.append( '(' ).append( x ).append( ',' ).append( y ).append( " tree=#" )
                            .append( Integer.toHexString( without.getRGB( x, y ) & 0xFFFFFF ) ).append( " number=#" )
                            .append( Integer.toHexString( with.getRGB( x, y ) & 0xFFFFFF ) ).append( ") " );
                    ++n;
                }
            }
        }
        return sb.toString();
    }

    /** The shipped demo pair (forester/demo/radial-numbers-fit.xml / -crowded.xml), in the UNROOTED layout. */
    /**
     * The INK band a number offers the branch-obstacle test: a digits-only number's stops at the digits' cap height
     * and the baseline (the line box carries air above and below), a number that is NOT digits only -- a support
     * label with its standard deviation, "88(0.05)" -- gets the whole line box, because a bracket or a parenthesis
     * reaches into that air. The second branch no fixture had entered until 2026-09-27 (a review find), so a
     * digitsOnly that said yes to a parenthesis would have passed everything. Read off the recorded claims:
     * [8] top, [9] bottom, [11] ink top, [12] ink bottom, all in the branch's frame.
     */
    private static void checkInkBand( final TreePanel tp, final boolean[] ok ) {
        final Phylogeny phy = fanTree();
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getBranchData().isHasConfidences() ) {
                for( final Confidence c : n.getBranchData().getConfidences() ) {
                    c.setStandardDeviation( 0.05 );
                }
            }
        }
        tp.setTree( phy );
        tp.recalculateMaxDistanceToRoot();
        tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        tp.getControlPanel().setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        final List<PhylogenyNode> branches = internalBranches( phy );
        try {
            tp.getOptions().setShowConfidenceStddev( false );
            render( tp, 900, true, true );
            int trimmed = 0, recorded = 0;
            for( final PhylogenyNode node : branches ) {
                final float[] c = tp.radialNumberClaimForTest( node );
                if ( c == null ) {
                    continue;
                }
                ++recorded;
                if ( ( c[ 11 ] > ( c[ 8 ] + 0.5 ) ) && ( c[ 12 ] < ( c[ 9 ] - 0.5 ) ) ) {
                    ++trimmed;
                }
            }
            if ( ( recorded == 0 ) || ( trimmed != recorded ) ) {
                fail( ok, "digits-only numbers must offer an ink band trimmed to cap height and baseline on BOTH sides; "
                        + trimmed + " of " + recorded + " recorded claims did" );
            }
            tp.getOptions().setShowConfidenceStddev( true );
            render( tp, 900, true, true );
            int full = 0;
            recorded = 0;
            for( final PhylogenyNode node : branches ) {
                final float[] c = tp.radialNumberClaimForTest( node );
                if ( c == null ) {
                    continue;
                }
                ++recorded;
                if ( ( Math.abs( c[ 11 ] - c[ 8 ] ) < 0.01 ) || ( Math.abs( c[ 12 ] - c[ 9 ] ) < 0.01 ) ) {
                    ++full; // the support's side reaches the line box's edge
                }
            }
            if ( ( recorded == 0 ) || ( full != recorded ) ) {
                fail( ok, "a support label with its standard deviation, 88(0.05), is not digits only and must offer the "
                        + "WHOLE line box on its side; " + full + " of " + recorded + " recorded claims did" );
            }
        }
        finally {
            tp.getOptions().setShowConfidenceStddev( false );
        }
    }

    private static void checkDemoPair( final TreePanel tp, final boolean[] ok ) {
        final Phylogeny fit = demo( "radial-numbers-fit.xml", ok );
        final Phylogeny crowded = demo( "radial-numbers-crowded.xml", ok );
        if ( ( fit == null ) || ( crowded == null ) ) {
            return;
        }
        final ControlPanel cp = tp.getControlPanel();
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        tp.setTree( fit );
        tp.recalculateMaxDistanceToRoot();
        final int w = 900;
        BufferedImage none = render( tp, w, false, false );
        BufferedImage both = render( tp, w, true, true );
        if ( tp.numbersRefusedForBranchesForTest() != 0 ) {
            fail( ok, "radial-numbers-fit.xml claims every number has room, but " + tp.numbersRefusedForBranchesForTest()
                    + " were refused for a branch" );
        }
        if ( tp.numbersDrawnForTest() != fit.getNumberOfExternalNodes() ) {
            fail( ok, "radial-numbers-fit.xml claims none of its " + fit.getNumberOfExternalNodes()
                    + " numbers is dropped, but " + tp.numbersDrawnForTest() + " were drawn" );
        }
        if ( inkOnInk( both, none ) != 0 ) {
            fail( ok, "radial-numbers-fit.xml draws a number across a branch (" + inkOnInk( both, none ) + " px)" );
        }
        tp.setTree( crowded );
        tp.recalculateMaxDistanceToRoot();
        // at the radial diameter of the 900x600 window the demo's description quotes (min(728, 513) of the panel
        // inside it): on a larger canvas the same 32 spokes have room, and a fan that crowds nowhere shows nothing
        checkDenseFan( tp, PHYLOGENY_GRAPHICS_TYPE.UNROOTED, ok, 513 ); // refuses some, keeps some, crosses nothing;
                                                                        // with auto-hide off the crossings appear
        // the exact split is QUOTED -- the demo catalogue row and the file's own description both say "16 of the 32
        // are drawn, 12 dropped for a spoke and 4 more for another number" -- so it is asserted, not printed: a
        // quoted number that only a note reports rots without anyone noticing (a review find, 2026-09-27)
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        render( tp, 513, true, true );
        final int drawn = tp.numbersDrawnForTest();
        final int for_branch = tp.numbersRefusedForBranchesForTest();
        final int for_number = tp.numbersSuppressedForTest() - for_branch - tp.numbersRefusedForLabelsForTest();
        if ( ( drawn != 16 ) || ( for_branch != 12 ) || ( for_number != 4 ) ) {
            fail( ok, "radial-numbers-crowded.xml at the 900x600 window: the README row and the generator's description "
                    + "quote 16 drawn, 12 dropped for a spoke and 4 for another number; got " + drawn + " / " + for_branch
                    + " / " + for_number + " -- if the rule or the font changed on purpose, update BOTH sentences and "
                    + "regenerate the demo (demosAreFreshOk insists)" );
        }
    }

    private static Phylogeny demo( final String name, final boolean[] ok ) {
        final java.io.File file = new java.io.File( System.getProperty( "user.dir" ), "forester/demo/" + name );
        if ( !file.exists() ) {
            fail( ok, "demo tree missing: " + file.getAbsolutePath() );
            return null;
        }
        try {
            return org.forester.phylogeny.factories.ParserBasedPhylogenyFactory.getInstance().create( file,
                    org.forester.io.parsers.phyloxml.PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            fail( ok, "could not read " + name + ": " + e );
            return null;
        }
    }

    /** How many pixels {@code with} changes where {@code without} already had ink (see {@link #INK_MAX}). */
    private static int inkOnInk( final BufferedImage with, final BufferedImage without ) {
        int n = 0;
        for( int y = 0; y < with.getHeight(); ++y ) {
            for( int x = 0; x < with.getWidth(); ++x ) {
                if ( ( with.getRGB( x, y ) != without.getRGB( x, y ) ) && isInk( without.getRGB( x, y ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean isInk( final int rgb ) {
        return ( ( ( rgb >> 16 ) & 0xFF ) < INK_MAX ) || ( ( ( rgb >> 8 ) & 0xFF ) < INK_MAX ) || ( ( rgb & 0xFF ) < INK_MAX );
    }

    /**
     * A star of 96 short tips, blank names, a length on every branch: in both radial layouts the spokes at a
     * number's midpoint are closer than the number is tall, so some numbers cannot be placed without lying across
     * a neighbour -- and some can.
     */
    private static Phylogeny denseStar() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 96; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "" );
            tip.setDistanceToParent( 0.25 + ( ( i % 3 ) * 0.05 ) );
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.recalculateNumberOfExternalDescendants( false );
        return phy;
    }

    /**
     * Inner nodes hanging on legs a few px long behind long tip branches, each with a support value so its pair is
     * two numbers tall: in the circular phylogram that pair reaches inward past the parent's radius and across the
     * parent's arc, while the tips' own numbers, far out on the ring, have room.
     */
    private static Phylogeny shortLegs() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int k = 0; k < 4; ++k ) {
            final PhylogenyNode a = new PhylogenyNode();
            a.setDistanceToParent( 0.3 );
            a.getBranchData().addConfidence( new Confidence( 90 + k, "bootstrap" ) );
            for( int j = 0; j < 2; ++j ) {
                final PhylogenyNode b = new PhylogenyNode();
                b.setDistanceToParent( 0.01 );
                b.getBranchData().addConfidence( new Confidence( 80 + j, "bootstrap" ) );
                for( int i = 0; i < 4; ++i ) {
                    final PhylogenyNode tip = new PhylogenyNode();
                    tip.setName( "" );
                    tip.setDistanceToParent( 0.4 );
                    b.addAsChild( tip );
                }
                a.addAsChild( b );
            }
            root.addAsChild( a );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.recalculateNumberOfExternalDescendants( false );
        return phy;
    }

    /** Renders with exactly one of the two numbers on and asserts, per branch, that the box it reserved contains
     *  the ink it drew. {@code none} is the same layout with neither number, so the difference is that number. */
    private static void renderRecording( final TreePanel tp, final int w, final boolean support, final boolean length,
                                         final boolean[] ok, final PHYLOGENY_GRAPHICS_TYPE gt, final String what,
                                         final List<PhylogenyNode> branches, final Phylogeny phy,
                                         final BufferedImage none ) {
        final BufferedImage img = render( tp, w, support, length );
        for( final PhylogenyNode node : branches ) {
            final float[] claim = tp.radialNumberClaimForTest( node );
            if ( claim == null ) {
                fail( ok, gt + ": no " + what + " reservation was recorded for the branch to node " + node.getId() );
                continue;
            }
            final double[] ma = midAndAngle( node, phy, gt );
            final int[] ink = inkBounds( img, none, ma );
            if ( ink == null ) {
                fail( ok, gt + ": the " + what + " drew nothing on the branch at "
                        + Math.round( Math.toDegrees( ma[ 2 ] ) ) + " deg" );
                continue;
            }
            // a pixel of slack: antialiasing writes just outside the glyph box the reservation is built from
            if ( ( ink[ 0 ] < ( claim[ 0 ] - 1.5 ) ) || ( ink[ 1 ] < ( claim[ 1 ] - 1.5 ) )
                    || ( ink[ 2 ] > ( claim[ 0 ] + claim[ 2 ] + 1.5 ) )
                    || ( ink[ 3 ] > ( claim[ 1 ] + claim[ 3 ] + 1.5 ) ) ) {
                fail( ok, gt + ": the " + what + " at " + Math.round( Math.toDegrees( ma[ 2 ] ) )
                        + " deg is drawn outside the space it reserved (ink " + ink[ 0 ] + "," + ink[ 1 ] + ".."
                        + ink[ 2 ] + "," + ink[ 3 ] + " vs reserved " + claim[ 0 ] + "," + claim[ 1 ] + " "
                        + claim[ 2 ] + "x" + claim[ 3 ] + ")" );
            }
        }
    }

    /** Device bounding box {@code {x0,y0,x1,y1}} of the ink {@code with} adds near this branch's midpoint. */
    private static int[] inkBounds( final BufferedImage with, final BufferedImage without,
                                    final double[] mid_and_angle ) {
        final double mx = mid_and_angle[ 0 ], my = mid_and_angle[ 1 ];
        int x0 = Integer.MAX_VALUE, y0 = Integer.MAX_VALUE, x1 = Integer.MIN_VALUE, y1 = Integer.MIN_VALUE;
        final int lx = Math.max( 0, (int) ( mx - RADIUS ) ), hx = Math.min( with.getWidth() - 1,
                                                                            (int) ( mx + RADIUS ) );
        final int ly = Math.max( 0, (int) ( my - RADIUS ) ), hy = Math.min( with.getHeight() - 1,
                                                                            (int) ( my + RADIUS ) );
        for( int y = ly; y <= hy; ++y ) {
            for( int x = lx; x <= hx; ++x ) {
                if ( ( with.getRGB( x, y ) != without.getRGB( x, y ) )
                        && ( Math.hypot( ( x + 0.5 ) - mx, ( y + 0.5 ) - my ) <= RADIUS ) ) {
                    x0 = Math.min( x0, x );
                    y0 = Math.min( y0, y );
                    x1 = Math.max( x1, x );
                    y1 = Math.max( y1, y );
                }
            }
        }
        return ( x1 < 0 ) ? null : new int[] { x0, y0, x1, y1 };
    }

    /**
     * The signed perpendicular distances of the nearest and furthest ink of ONE number from its branch line, or
     * null when that number drew nothing near this branch. Positive is the side 90 deg clockwise of the branch's
     * outward direction (device axes, y down); the sign itself is arbitrary, its CONSISTENCY across the fan is the
     * point. {@code with} is a render carrying just this one number, {@code without} the same layout carrying
     * neither, so the difference between them is exactly this number's ink.
     */
    private static double[] perpendicularSpan( final BufferedImage with, final BufferedImage without,
                                               final double[] mid_and_angle ) {
        final double mx = mid_and_angle[ 0 ], my = mid_and_angle[ 1 ], angle = mid_and_angle[ 2 ];
        final double ux = Math.cos( angle ), uy = Math.sin( angle ); // along the branch, outward
        double nearest = Double.MAX_VALUE, furthest = 0;
        boolean any = false;
        final int x0 = Math.max( 0, (int) ( mx - RADIUS ) ), x1 = Math.min( with.getWidth() - 1, (int) ( mx + RADIUS ) );
        final int y0 = Math.max( 0, (int) ( my - RADIUS ) ), y1 = Math.min( with.getHeight() - 1,
                                                                            (int) ( my + RADIUS ) );
        for( int y = y0; y <= y1; ++y ) {
            for( int x = x0; x <= x1; ++x ) {
                if ( with.getRGB( x, y ) == without.getRGB( x, y ) ) {
                    continue;
                }
                final double dx = ( x + 0.5 ) - mx, dy = ( y + 0.5 ) - my;
                if ( Math.hypot( dx, dy ) > RADIUS ) {
                    continue;
                }
                final double perp = ( dx * uy ) - ( dy * ux ); // the branch-normal component
                if ( !any || ( Math.abs( perp ) < Math.abs( nearest ) ) ) {
                    nearest = perp;
                }
                if ( !any || ( Math.abs( perp ) > Math.abs( furthest ) ) ) {
                    furthest = perp;
                }
                any = true;
            }
        }
        return any ? new double[] { nearest, furthest } : null;
    }

    /** {@code {mid x, mid y, branch direction}} for a node's incoming branch, as each radial painter draws it: a
     *  radial leg along the node's own spoke in CIRCULAR, a straight parent-&gt;node line in UNROOTED. */
    private static double[] midAndAngle( final PhylogenyNode node, final Phylogeny phy,
                                         final PHYLOGENY_GRAPHICS_TYPE gt ) {
        final PhylogenyNode parent = node.getParent();
        if ( gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            final PhylogenyNode root = phy.getRoot();
            final double angle = Math.atan2( node.getYcoord() - root.getYcoord(),
                                             node.getXcoord() - root.getXcoord() );
            final double pr = Math.hypot( parent.getXcoord() - root.getXcoord(),
                                          parent.getYcoord() - root.getYcoord() );
            final double ix = root.getXcoord() + ( Math.cos( angle ) * pr );
            final double iy = root.getYcoord() + ( Math.sin( angle ) * pr );
            return new double[] { ( node.getXcoord() + ix ) / 2.0, ( node.getYcoord() + iy ) / 2.0, angle };
        }
        return new double[] { ( node.getXcoord() + parent.getXcoord() ) / 2.0,
                              ( node.getYcoord() + parent.getYcoord() ) / 2.0,
                              Math.atan2( node.getYcoord() - parent.getYcoord(),
                                          node.getXcoord() - parent.getXcoord() ) };
    }

    private static List<PhylogenyNode> internalBranches( final Phylogeny phy ) {
        final List<PhylogenyNode> out = new ArrayList<PhylogenyNode>();
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isRoot() && !n.isExternal() ) {
                out.add( n );
            }
        }
        return out;
    }

    private static BufferedImage render( final TreePanel tp, final int w, final boolean support,
                                         final boolean length ) {
        tp.setShows( DisplayOption.WRITE_CONFIDENCE_VALUES, support );
        tp.setShows( DisplayOption.WRITE_BRANCH_LENGTH_VALUES, length );
        tp.setRecordRadialNumberClaimsForTest( true ); // a fresh recording of THIS paint's reservations
        tp.setSize( w, w );
        tp.fitRadialTo( w, w ); // the radial canvas is its own square, not the panel's -- setSize alone leaves it
        tp.calcParametersForPainting( w, w );
        tp.resetPreferredSize(); // the tree lays out at its PREFERRED size; render there so node coords map 1:1
        final int pw = (int) Math.ceil( tp.getPreferredSize().getWidth() );
        final int ph = (int) Math.ceil( tp.getPreferredSize().getHeight() );
        final BufferedImage img = new BufferedImage( pw, ph, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
        final ExportTheme theme = ExportTheme.applyIf( tp, true ); // a white background, so "ink" is "not white"
        try {
            tp.paintPhylogeny( g, false, true, pw, ph, 0, 0 );
        }
        finally {
            theme.restore();
            g.dispose();
        }
        return img;
    }

    /**
     * A balanced 8-tip tree whose branches fan out in every direction, with a length and a confidence on each
     * internal branch and BLANK tip names -- the only ink near a branch midpoint is then the tree and the two
     * numbers. Internal branches are made the longest so their midpoints stay well clear of one another.
     */
    private static Phylogeny fanTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 4; ++i ) {
            final PhylogenyNode inner = new PhylogenyNode();
            inner.setDistanceToParent( 0.4 );
            inner.getBranchData().addConfidence( new Confidence( 88 + i, "bootstrap" ) );
            for( int j = 0; j < 2; ++j ) {
                final PhylogenyNode inner2 = new PhylogenyNode();
                inner2.setDistanceToParent( 0.35 );
                inner2.getBranchData().addConfidence( new Confidence( 70 + j, "bootstrap" ) );
                for( int k = 0; k < 2; ++k ) {
                    final PhylogenyNode tip = new PhylogenyNode();
                    tip.setName( "" );
                    tip.setDistanceToParent( 0.3 );
                    inner2.addAsChild( tip );
                }
                inner.addAsChild( inner2 );
            }
            root.addAsChild( inner );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.recalculateNumberOfExternalDescendants( false );
        return phy;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [RadialBranchNumberRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [RadialBranchNumberRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private RadialBranchNumberRenderTest() {
    }
}

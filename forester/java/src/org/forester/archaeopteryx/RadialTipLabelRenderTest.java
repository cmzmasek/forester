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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Tip labels in the two radial layouts: a label is drawn only where no label already drawn this pass is in its way
 * (first come, under "Auto-hide Labels"), a found node's label always, and the branch numbers keep clear of the
 * labels that were drawn.
 * <p>
 * Why this exists: the unrooted layout never thinned its tip labels at all, and circular hid every k-th by index --
 * a proxy that hid labels which did not overlap and kept ones that did (measured on the bat phylogeny at 1100 x 850:
 * 66 overlapping pairs among 34 flat labels unrooted, 21 circular; along the spoke, 40 and 0). Everything here is
 * read off the boxes the paint recorded for the labels it DREW ({@code labelBoxesForTest}), tested with the same
 * oriented-box overlap the rule uses -- and, for the numbers, off the pixels. Headful; a green no-op when headless.
 */
public final class RadialTipLabelRenderTest {

    private static final int TIPS = 96;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RadialTipLabelRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            return labelsOk();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "threw: " + e );
        }
    }

    private static boolean labelsOk() throws Exception {
        final Phylogeny phy = starOfTwoLengths();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy },
                new Configuration(), "radiallabels" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final MainFrame frame = mf[ 0 ];
            try {
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                final ControlPanel cp = tp.getControlPanel();
                frame.getOptions().setSupportVisualization( Options.SUPPORT_VISUALIZATION.NONE );
                cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                frame.showWhole();
                tp.setRecordLabelBoxesForTest( true );

                // (1) UNROOTED, labels along the spoke: the short tips crowd, the long ones do not
                frame.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                crowdedFan( tp, cp, "UNROOTED radial", ok, true );
                inkInsideBoxes( tp, "UNROOTED radial", ok );

                // (2) the found node: a label that the rule would hide is drawn all the same, and claims its place
                foundNode( tp, frame, phy, "UNROOTED radial", 900, ok );

                // (3) CIRCULAR along the spoke: every tip label on the ring at even spacing -- nothing overlaps,
                // nothing is hidden -- as a cladogram AND as a phylogram, because every circular phylogram carries
                // its labels to the outer ring (Christian, 2026-09-27; before that the phylogram left the short
                // tips' labels at a third of the radius, where they crowded: 72 drawn, 24 hidden)
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                for( final Options.PHYLOGENY_DISPLAY_TYPE shape : new Options.PHYLOGENY_DISPLAY_TYPE[] {
                        Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM, Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM } ) {
                    cp.setTreeDisplayType( shape );
                    render( tp, 900, false, true );
                    if ( ( tp.labelsHiddenForTest() != 0 ) || ( tp.labelsDrawnForTest() != TIPS ) ) {
                        fail( ok, "CIRCULAR radial " + shape + ": on the ring at even spacing every label has room (drawn "
                                + tp.labelsDrawnForTest() + " hidden " + tp.labelsHiddenForTest() + " of " + TIPS + ")" );
                    }
                    if ( overlappingPairs( tp.labelBoxesForTest(), -1 ) != 0 ) {
                        fail( ok, "CIRCULAR radial " + shape + ": the drawn labels must not overlap" );
                    }
                }
                // (3a) the same two guards CIRCULAR along the spoke: every label pixel inside its recorded box (the ring
                //      anchor is not the tip anchor), and the found-node exemption at 450 px, where the ring hides half
                //      the names. Both had been asked of unrooted only -- the layout that had the problem first
                cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                inkInsideBoxes( tp, "CIRCULAR radial", ok );
                foundNode( tp, frame, phy, "CIRCULAR radial", 450, ok );

                // (4) CIRCULAR lying flat, as a phylogram: flat labels stack on each other, and the rule hides the
                // ones that would -- exactly, so the rest are clear
                frame.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.HORIZONTAL );
                crowdedFan( tp, cp, "CIRCULAR horizontal", ok, false );
                inkInsideBoxes( tp, "CIRCULAR horizontal", ok );
                frame.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );

                // (4a) rotating the fan: WHICH names survive a turn is promised by neither program, but the rule is
                //      a pure function of the geometry on screen, and no rotation overprints (joint with JS, 2026-09-27)
                rotationOk( tp, cp, ok );

                // (4b) a REFUSED name is never recorded, so it blocks nothing after it. On this uniform ring at 450 px
                //      each name collides only with its two neighbours, and first come keeps every OTHER one -- exactly
                //      half. Recording refused names instead CHAINS: each refused name blocks the next, which is refused
                //      and recorded in turn, and only the FIRST name survives (measured under that mutant: 1 of 96). The
                //      mutant survived every other assertion here (2026-09-27), the property having lived in
                //      OrientedOccupancy.claim() until the caller took over asking and placing
                frame.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
                cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                render( tp, 450, false, true );
                if ( tp.labelsDrawnForTest() != ( TIPS / 2 ) ) {
                    fail( ok, "CIRCULAR radial @450 px: on a uniform ring where only neighbours collide, first come keeps "
                            + "exactly every other name, " + ( TIPS / 2 ) + " of " + TIPS + "; got " + tp.labelsDrawnForTest()
                            + " drawn, " + tp.labelsHiddenForTest() + " hidden (one survivor = refused names were recorded and chained)" );
                }

                // (4c) a tip with nothing to draw -- a taxonomy that is "shown" but renders empty, no name, no image --
                //      is neither drawn nor hidden (a review find, 2026-09-27: it was counted drawn, box of zero width)
                emptyTaxonomyDrawsNothing( tp, cp, ok );

                // (4d) a COLLAPSED clade's radial label claims its space like a tip's (a review find, 2026-09-27: it was
                //      painted with a bare drawString, invisible to the rule, and names were granted the same space)
                collapsedCladeLabel( tp, cp, ok );
                tp.setTree( phy );
                tp.recalculateMaxDistanceToRoot();

                // (5) the numbers keep clear of the labels: with auto-hide on, no number ink lands on label ink
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, false );
                final BufferedImage lab_off = render( tp, 900, false, true );
                final BufferedImage both_off = render( tp, 900, true, true );
                final int over_off = inkOnInk( both_off, lab_off );
                cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                final BufferedImage lab_on = render( tp, 900, false, true );
                final BufferedImage both_on = render( tp, 900, true, true );
                final int over_on = inkOnInk( both_on, lab_on );
                if ( over_off == 0 ) {
                    fail( ok, "precondition -- with auto-hide off some number must land on a label or a branch here, "
                            + "or the fixture cannot show the numbers giving way (0 px)" );
                }
                if ( over_on != 0 ) {
                    fail( ok, "with auto-hide on a number is drawn across a label or a branch (" + over_on
                            + " number pixels on ink; " + over_off + " with auto-hide off)" );
                }
                if ( tp.numbersRefusedForLabelsForTest() == 0 ) {
                    fail( ok, "the fixture must make at least one number give way to a LABEL (none refused for one)" );
                }

                // (9) the shipped demo pair, at the radial diameter of the 1100x850 window the README quotes
                // (min(928, 763) of the panel inside it): the file that says every name fits must hide none; the
                // one that says it crowds must hide some, keep some, and overlap nothing
                frame.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
                cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                final Phylogeny fit_demo = demo( "radial-labels-fit.xml", ok );
                if ( fit_demo != null ) {
                    tp.setTree( fit_demo );
                    tp.recalculateMaxDistanceToRoot();
                    render( tp, 763, false, true );
                    if ( ( tp.labelsHiddenForTest() != 0 ) || ( tp.labelsDrawnForTest() != fit_demo.getNumberOfExternalNodes() ) ) {
                        fail( ok, "radial-labels-fit.xml claims every name has room, got drawn " + tp.labelsDrawnForTest()
                                + " hidden " + tp.labelsHiddenForTest() + " of " + fit_demo.getNumberOfExternalNodes() );
                    }
                }
                final Phylogeny crowded_demo = demo( "radial-labels-crowded.xml", ok );
                if ( crowded_demo != null ) {
                    tp.setTree( crowded_demo );
                    tp.recalculateMaxDistanceToRoot();
                    render( tp, 763, false, true );
                    refusalsCaused( tp, "radial-labels-crowded.xml at the 1100x850 window", ok );
                    if ( ( tp.labelsHiddenForTest() == 0 ) || ( tp.labelsDrawnForTest() == 0 )
                            || ( overlappingPairs( tp.labelBoxesForTest(), -1 ) != 0 ) ) {
                        fail( ok, "radial-labels-crowded.xml must hide some names, keep some and overlap none; got drawn "
                                + tp.labelsDrawnForTest() + " hidden " + tp.labelsHiddenForTest() + " overlapping "
                                + overlappingPairs( tp.labelBoxesForTest(), -1 ) );
                    }
                    // the exact counts are QUOTED -- the demo catalogue row and the file's own description (written by
                    // DemoTreeGenerator) both say "72 of the 96 names are drawn and 24 hidden" -- so they are asserted,
                    // not printed: a note here let a mutant that chained refused names read 49 / 47 and pass (2026-09-27)
                    if ( ( tp.labelsDrawnForTest() != 72 ) || ( tp.labelsHiddenForTest() != 24 ) ) {
                        fail( ok, "radial-labels-crowded.xml at the 1100x850 window: the README row and the generator's "
                                + "description quote 72 drawn / 24 hidden of 96; got " + tp.labelsDrawnForTest() + " / "
                                + tp.labelsHiddenForTest() + " -- if the rule or the font changed on purpose, update BOTH "
                                + "sentences and regenerate the demo (demosAreFreshOk insists)" );
                    }
                }

                // (10) what a tip draws off its label goes with the label: on a fan too dense for its names, a tip
                // whose name the rule hides draws no domain architecture either -- the rectangular layout hides a
                // thinned row's label and track together, and an architecture without its name is identifiable only
                // by its spoke (archaeopteryx.js's find on their side, 2026-09-27). Labels switched off hide nothing
                // by the rule, so every architecture still draws then.
                final Phylogeny domained = domainedStar();
                tp.setTree( domained );
                tp.recalculateMaxDistanceToRoot();
                cp.setShowDomainArchitecturesForTest( true );
                try {
                    frame.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
                    tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                    cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                    render( tp, 900, false, true );
                    if ( tp.labelsHiddenForTest() == 0 ) {
                        fail( ok, "UNROOTED domains: precondition -- the domained fan must hide some names" );
                    }
                    else if ( tp.radialDomainsDrawnForTest() != tp.labelsDrawnForTest() ) {
                        fail( ok, "UNROOTED domains: an architecture goes with its name -- " + tp.radialDomainsDrawnForTest()
                                + " architectures drawn for " + tp.labelsDrawnForTest() + " names drawn ("
                                + tp.labelsHiddenForTest() + " hidden)" );
                    }
                    cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, false );
                    render( tp, 900, false, true );
                    if ( tp.radialDomainsDrawnForTest() != TIPS ) {
                        fail( ok, "UNROOTED domains: with auto-hide off every architecture draws (" + tp.radialDomainsDrawnForTest()
                                + " of " + TIPS + ")" );
                    }
                    cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                    render( tp, 900, false, false ); // names OFF: nothing hidden by the rule, every architecture drawn
                    if ( tp.radialDomainsDrawnForTest() != TIPS ) {
                        fail( ok, "UNROOTED domains: with names off nothing is hidden by the rule, so every architecture "
                                + "draws (" + tp.radialDomainsDrawnForTest() + " of " + TIPS + ")" );
                    }
                    tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    render( tp, 900, false, true );
                    if ( tp.radialDomainsDrawnForTest() != tp.labelsDrawnForTest() ) {
                        fail( ok, "CIRCULAR domains: an architecture goes with its name -- " + tp.radialDomainsDrawnForTest()
                                + " architectures for " + tp.labelsDrawnForTest() + " names" );
                    }
                }
                finally {
                    cp.setShowDomainArchitecturesForTest( false );
                }

                // (8) INTERNAL-node labels (clade names) follow the same rule, after every tip: a clade name never
                // displaces a tip name, the larger clade's name wins where two meet, and a found clade is always
                // drawn (Christian, 2026-09-26)
                cladeLabels( tp, cp, ok );

                // (6) the shipped default: the radial layouts open with labels along the spoke (Christian,
                // 2026-09-25; archaeopteryx.js has the same default, so the two agree by contract)
                if ( Options.createInstance().getNodeLabelDirection() != NODE_LABEL_DIRECTION.RADIAL ) {
                    fail( ok, "a fresh Options must default the node-label direction to RADIAL, got "
                            + Options.createInstance().getNodeLabelDirection() );
                }

                // (7) a real tree, the one the rule was measured on: the bat phylogeny in UNROOTED with radial
                // labels -- 40 overlapping pairs among 34 labels without the rule; with it, some hidden, some kept,
                // none overlapping
                final Phylogeny bat = demo( "bat-phylogeny.xml", ok );
                if ( bat != null ) {
                    tp.setTree( bat );
                    tp.recalculateMaxDistanceToRoot();
                    tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                    cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
                    render( tp, 900, false, true );
                    refusalsCaused( tp, "bat-phylogeny.xml UNROOTED radial", ok );
                    final int pairs = overlappingPairs( tp.labelBoxesForTest(), -1 );
                    if ( ( tp.labelsHiddenForTest() == 0 ) || ( tp.labelsDrawnForTest() == 0 ) || ( pairs != 0 ) ) {
                        fail( ok, "bat-phylogeny.xml UNROOTED radial: expected some labels hidden, some drawn and none "
                                + "overlapping; got drawn " + tp.labelsDrawnForTest() + " hidden "
                                + tp.labelsHiddenForTest() + " overlapping pairs " + pairs );
                    }
                }
            }
            catch ( final Exception e ) {
                e.printStackTrace();
                fail( ok, "threw: " + e );
            }
        } );
        return ok[ 0 ];
    }

    /**
     * On a fan too dense for its labels: auto-hide ON hides some, keeps some, and the kept ones never overlap;
     * auto-hide OFF draws every one, and they do overlap (so the fixture is really crowded). {@code all_fit_expected}
     * is not asserted -- it only names the layout in the messages.
     */
    private static void crowdedFan( final TreePanel tp, final ControlPanel cp, final String what, final boolean[] ok,
                                    final boolean unrooted ) {
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        render( tp, 900, false, true );
        final List<double[]> on = tp.labelBoxesForTest();
        final int drawn_on = tp.labelsDrawnForTest(), hidden_on = tp.labelsHiddenForTest();
        final int pairs_on = overlappingPairs( on, -1 );
        refusalsCaused( tp, what + " (auto-hide on)", ok );
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, false );
        render( tp, 900, false, true );
        final int drawn_off = tp.labelsDrawnForTest(), hidden_off = tp.labelsHiddenForTest();
        final int pairs_off = overlappingPairs( tp.labelBoxesForTest(), -1 );
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        if ( ( drawn_off != TIPS ) || ( hidden_off != 0 ) ) {
            fail( ok, what + ": with auto-hide off every label is drawn (drawn " + drawn_off + " hidden " + hidden_off
                    + " of " + TIPS + ")" );
        }
        if ( pairs_off == 0 ) {
            fail( ok, what + ": precondition -- with auto-hide off the labels must overlap (0 pairs), so this fan "
                    + "is not crowded and shows nothing" );
            return;
        }
        if ( pairs_on != 0 ) {
            fail( ok, what + ": with auto-hide on, drawn labels still overlap (" + pairs_on + " pairs; " + pairs_off
                    + " with it off)" );
        }
        if ( hidden_on == 0 ) {
            fail( ok, what + ": the crowded fan must make the rule hide at least one label" );
        }
        if ( drawn_on == 0 ) {
            fail( ok, what + ": the rule must hide what overlaps, not everything (0 drawn, " + hidden_on + " hidden)" );
        }
        if ( ( drawn_on + hidden_on ) != TIPS ) {
            fail( ok, what + ": drawn + hidden must account for every tip (" + drawn_on + " + " + hidden_on + " != "
                    + TIPS + ")" );
        }
    }

    /**
     * Every pixel of label ink lies inside the box some drawn label recorded (grown by a pixel for antialiasing):
     * the box the rule reasons with is where the text actually is. Independent of the rule's own arithmetic -- a
     * box recorded at the wrong angle, or centred off its text, leaves ink outside every box and fails here even
     * though the recorded boxes would still agree with each other.
     */
    private static void inkInsideBoxes( final TreePanel tp, final String what, final boolean[] ok ) {
        // The layout of a radial view depends on the labels (their reserve sets the radius), so a render without
        // them is a DIFFERENT tree, and the difference is the whole picture. Lay out WITH the names, then paint
        // without them: same geometry, and the difference between the two paints is exactly the label ink.
        final BufferedImage without = renderLaidOutWithLabelsPaintedWithout( tp, 900 );
        final BufferedImage with = render( tp, 900, false, true );
        final List<double[]> boxes = tp.labelBoxesForTest();
        if ( ( without.getWidth() != with.getWidth() ) || ( without.getHeight() != with.getHeight() ) ) {
            fail( ok, what + ": the two renders must share a size for the label ink to be isolated" );
            return;
        }
        // A label's dotted LEADER (circular: from a short tip out to its anchor on the ring) appears exactly when
        // the label does, so it is in this difference too, and it lies outside the label's box by design: skip
        // pixels on any tip->anchor segment, exactly, rather than by colour.
        final List<double[]> leaders = new java.util.ArrayList<double[]>();
        for( final PhylogenyNode tip : tp.getPhylogeny().getExternalNodes() ) {
            final java.awt.geom.Point2D.Double a = tp.circularLabelAnchorForTest( tip );
            if ( Math.hypot( a.x - tip.getXcoord(), a.y - tip.getYcoord() ) >= 1 ) {
                leaders.add( new double[] { tip.getXcoord(), tip.getYcoord(), a.x, a.y } );
            }
        }
        int ink = 0, outside = 0;
        for( int y = 0; y < with.getHeight(); ++y ) {
            for( int x = 0; x < with.getWidth(); ++x ) {
                if ( with.getRGB( x, y ) == without.getRGB( x, y ) ) {
                    continue;
                }
                if ( onALeader( x + 0.5, y + 0.5, leaders ) ) {
                    continue;
                }
                ++ink;
                boolean in = false;
                for( final double[] b : boxes ) {
                    final double dx = ( x + 0.5 ) - b[ 0 ], dy = ( y + 0.5 ) - b[ 1 ];
                    final double u = ( dx * Math.cos( b[ 4 ] ) ) + ( dy * Math.sin( b[ 4 ] ) );
                    final double v = ( -dx * Math.sin( b[ 4 ] ) ) + ( dy * Math.cos( b[ 4 ] ) );
                    if ( ( Math.abs( u ) <= ( b[ 2 ] + 1.5 ) ) && ( Math.abs( v ) <= ( b[ 3 ] + 1.5 ) ) ) {
                        in = true;
                        break;
                    }
                }
                if ( !in ) {
                    ++outside;
                }
            }
        }
        if ( ink < 500 ) {
            fail( ok, what + ": precondition -- the labels must put ink on the canvas (" + ink + " px)" );
        }
        if ( outside > 0 ) {
            fail( ok, what + ": " + outside + " of " + ink + " label pixels lie outside every recorded label box -- "
                    + "the boxes are not where the text is" );
        }
    }

    /** Whether the point lies within 1.5 px of any tip->anchor leader segment. */
    private static boolean onALeader( final double px, final double py, final List<double[]> leaders ) {
        for( final double[] s : leaders ) {
            final double dx = s[ 2 ] - s[ 0 ], dy = s[ 3 ] - s[ 1 ];
            final double len2 = ( dx * dx ) + ( dy * dy );
            double t = ( ( ( px - s[ 0 ] ) * dx ) + ( ( py - s[ 1 ] ) * dy ) ) / len2;
            t = Math.max( 0, Math.min( 1, t ) );
            final double ex = s[ 0 ] + ( t * dx ), ey = s[ 1 ] + ( t * dy );
            if ( Math.hypot( px - ex, py - ey ) <= 2.5 ) { // the dotted stroke is wider than a pixel
                return true;
            }
        }
        return false;
    }

    /** Lays the tree out with the names shown, then paints it with them hidden -- see {@link #inkInsideBoxes}. */
    private static BufferedImage renderLaidOutWithLabelsPaintedWithout( final TreePanel tp, final int w ) {
        tp.setShows( DisplayOption.WRITE_BRANCH_LENGTH_VALUES, false );
        tp.setShows( DisplayOption.SHOW_NODE_NAMES, true );
        tp.setSize( w, w );
        tp.fitRadialTo( w, w );
        tp.calcParametersForPainting( w, w );
        tp.resetPreferredSize();
        final int pw = (int) Math.ceil( tp.getPreferredSize().getWidth() );
        final int ph = (int) Math.ceil( tp.getPreferredSize().getHeight() );
        tp.setShows( DisplayOption.SHOW_NODE_NAMES, false );
        final BufferedImage img = new BufferedImage( pw, ph, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
        final ExportTheme theme = ExportTheme.applyIf( tp, true );
        try {
            tp.paintPhylogeny( g, false, true, pw, ph, 0, 0 );
        }
        finally {
            theme.restore();
            g.dispose();
            tp.setShows( DisplayOption.SHOW_NODE_NAMES, true );
        }
        return img;
    }

    private static boolean overlaps( final double[] a, final double[] b ) {
        return OrientedOccupancy.overlap( a[ 0 ], a[ 1 ], a[ 2 ], a[ 3 ], Math.cos( a[ 4 ] ), Math.sin( a[ 4 ] ),
                b[ 0 ], b[ 1 ], b[ 2 ], b[ 3 ], Math.cos( b[ 4 ] ), Math.sin( b[ 4 ] ) );
    }

    /** The clade-label half of the rule, on a bush of named clades, in both radial layouts. */
    private static void cladeLabels( final TreePanel tp, final ControlPanel cp, final boolean[] ok ) {
        final Phylogeny bush = bushOfNamedClades();
        tp.setTree( bush );
        tp.recalculateMaxDistanceToRoot();
        tp.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        for( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.UNROOTED,
                PHYLOGENY_GRAPHICS_TYPE.CIRCULAR } ) {
            tp.setPhylogenyGraphicsType( gt );
            final String what = gt + " clade labels";
            // tips alone
            tp.setShowInternalDataForThisTab( false );
            render( tp, 900, false, true );
            final Set<Long> tips_alone = drawnIds( tp.labelBoxesForTest() );
            // tips and clade names
            tp.setShowInternalDataForThisTab( true );
            render( tp, 900, false, true );
            final List<double[]> boxes = tp.labelBoxesForTest();
            refusalsCaused( tp, what, ok );
            final Set<Long> tips_with = new HashSet<Long>();
            final List<double[]> internal = new java.util.ArrayList<double[]>();
            for( final double[] b : boxes ) {
                final PhylogenyNode n = bush.getNode( (long) b[ 5 ] );
                if ( n.isExternal() ) {
                    tips_with.add( Long.valueOf( n.getId() ) );
                }
                else {
                    internal.add( b );
                }
            }
            if ( !tips_with.equals( tips_alone ) ) {
                fail( ok, what + ": a clade name must never displace a tip name -- the drawn tip labels changed when "
                        + "the clade labels were switched on (" + tips_alone.size() + " -> " + tips_with.size() + ")" );
            }
            if ( ( tp.internalLabelsHiddenForTest() == 0 ) || ( tp.internalLabelsDrawnForTest() == 0 ) ) {
                fail( ok, what + ": the bush must make the rule hide some clade names and keep some (drawn "
                        + tp.internalLabelsDrawnForTest() + " hidden " + tp.internalLabelsHiddenForTest() + ")" );
            }
            if ( overlappingPairs( boxes, -1 ) != 0 ) {
                fail( ok, what + ": drawn labels (tips and clades together) must not overlap ("
                        + overlappingPairs( boxes, -1 ) + " pairs)" );
            }
            ancestorsFirst( internal, bush, what, ok, true );
            // the order the rule OFFERS clade names in: larger clade first (tip count), shallower among equals. The
            // bush has a size inversion for this -- superclade R (4 tips, depth 1) against P and Q (5 tips, depth
            // 3): R must be offered AFTER them, which a depth-first order would get backwards.
            final List<Long> offered = tp.internalLabelOrderForTest();
            boolean inversion_seen = false;
            for( int i = 1; i < offered.size(); ++i ) {
                final PhylogenyNode a = bush.getNode( offered.get( i - 1 ) ), b = bush.getNode( offered.get( i ) );
                final int ta = a.getNumberOfExternalNodes(), tb = b.getNumberOfExternalNodes();
                if ( ( ta < tb ) || ( ( ta == tb ) && ( depth( a ) > depth( b ) ) ) ) {
                    fail( ok, what + ": clade names must be offered larger clade first, shallower among equals; got "
                            + a.getName() + " (" + ta + " tips, depth " + depth( a ) + ") before " + b.getName() + " ("
                            + tb + " tips, depth " + depth( b ) + ")" );
                    break;
                }
                if ( depth( a ) > depth( b ) ) {
                    inversion_seen = true; // a deeper clade offered before a shallower one: size outranked depth
                }
            }
            if ( !inversion_seen ) {
                fail( ok, what + ": precondition -- the bush must offer some deeper, larger clade before a shallower "
                        + "smaller one, or size-before-depth was never exercised" );
            }
            // a found clade is drawn although the rule would hide it
            PhylogenyNode hidden_clade = null;
            for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = bush.iteratorPreorder(); it
                    .hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( !n.isRoot() && !n.isExternal() && !drawn( boxes, n.getId() ) ) {
                    hidden_clade = n;
                    break;
                }
            }
            if ( hidden_clade != null ) {
                final Set<Long> found = new HashSet<Long>();
                found.add( Long.valueOf( hidden_clade.getId() ) );
                tp.setFoundNodes0( found );
                render( tp, 900, false, true );
                tp.setFoundNodes0( null );
                if ( !drawn( tp.labelBoxesForTest(), hidden_clade.getId() ) ) {
                    fail( ok, what + ": a FOUND clade's label must be drawn although the rule would hide it" );
                }
            }
            // the numbers keep clear of CLADE names too. In circular the tip names sit on the ring, outside every
            // number, so a clade name inside it is the only label a number can meet there -- and on this bush none
            // does (measured: 0 of 54 drawn numbers gave way to a clade name at 900 px, 7 to a branch), so the
            // label half of the precondition is asked of UNROOTED only, where the bush makes numbers yield to clade
            // names; circular asserts the outcome alone. No demo carries both drawn numbers and clade names in a
            // radial layout (measured over five), so this bush is the fixture.
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, false );
            final BufferedImage clade_off = render( tp, 900, false, true );
            final BufferedImage clade_both_off = render( tp, 900, true, true );
            final int clade_over_off = inkOnInk( clade_both_off, clade_off );
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
            final BufferedImage clade_on = render( tp, 900, false, true );
            final BufferedImage clade_both_on = render( tp, 900, true, true );
            final int clade_over_on = inkOnInk( clade_both_on, clade_on );
            if ( ( clade_over_off == 0 ) || ( ( gt == PHYLOGENY_GRAPHICS_TYPE.UNROOTED )
                    && ( tp.numbersRefusedForLabelsForTest() == 0 ) ) ) {
                fail( ok, what + ": precondition -- with auto-hide off some number must land on a clade name or a branch ("
                        + clade_over_off + " px) and, unrooted, with it on at least one number must give way to a LABEL ("
                        + tp.numbersRefusedForLabelsForTest() + " did; " + tp.numbersRefusedForBranchesForTest()
                        + " to a branch, " + tp.numbersDrawnForTest() + " drawn)" );
            }
            if ( clade_over_on != 0 ) {
                fail( ok, what + ": with auto-hide on a number is drawn across a clade name or a branch (" + clade_over_on
                        + " number pixels on ink; " + clade_over_off + " with auto-hide off)" );
            }
            render( tp, 900, false, true ); // numbers off again for what follows
            // the recorded ORDER is per pass, like the boxes: a paint that queues NO clade label (Show Internal Data
            // off again) must not go on reporting the previous paint's -- painted under the SAME recording, so the
            // list spans both passes (a review find, 2026-09-27: the empty-queue early return skipped the clear)
            tp.setShowInternalDataForThisTab( false );
            paint( tp, 900 );
            if ( !tp.internalLabelOrderForTest().isEmpty() ) {
                fail( ok, what + ": a paint that offers no clade label must record an EMPTY order, got the previous "
                        + "paint's " + tp.internalLabelOrderForTest().size() + " entries" );
            }
            tp.setShowInternalDataForThisTab( true );
            // auto-hide off: every clade label drawn
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, false );
            render( tp, 900, false, true );
            if ( tp.internalLabelsHiddenForTest() != 0 ) {
                fail( ok, what + ": with auto-hide off no clade label is hidden (" + tp.internalLabelsHiddenForTest()
                        + " were)" );
            }
            // the numbers give way to clade names as they do to tip names: painted after every label
            final BufferedImage lab_off = render( tp, 900, false, true );
            final BufferedImage both_off = render( tp, 900, true, true );
            final int over_off = inkOnInk( both_off, lab_off );
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
            final BufferedImage lab_on = render( tp, 900, false, true );
            final BufferedImage both_on = render( tp, 900, true, true );
            final int over_on = inkOnInk( both_on, lab_on );
            if ( over_off == 0 ) {
                fail( ok, what + ": precondition -- with auto-hide off some number must land on a label here (0 px)" );
            }
            else if ( over_on != 0 ) {
                fail( ok, what + ": with auto-hide on a number is drawn across a label (" + over_on + " px; "
                        + over_off + " with it off)" );
            }
        }
        // the bat phylogeny, whose family and suborder names sat across its tips: some are now hidden, none overlap
        final Phylogeny bat = demo( "bat-phylogeny.xml", ok );
        if ( bat != null ) {
            tp.setTree( bat );
            tp.recalculateMaxDistanceToRoot();
            tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            tp.setShowInternalDataForThisTab( true );
            render( tp, 900, false, true );
            if ( ( tp.internalLabelsHiddenForTest() == 0 ) || ( overlappingPairs( tp.labelBoxesForTest(), -1 ) != 0 ) ) {
                fail( ok, "bat-phylogeny.xml UNROOTED: expected some clade labels hidden and no drawn label overlapping; got "
                        + "internal drawn " + tp.internalLabelsDrawnForTest() + " hidden " + tp.internalLabelsHiddenForTest()
                        + " overlapping pairs " + overlappingPairs( tp.labelBoxesForTest(), -1 ) );
            }
            final List<double[]> bat_internal = new java.util.ArrayList<double[]>();
            for( final double[] b : tp.labelBoxesForTest() ) {
                if ( !bat.getNode( (long) b[ 5 ] ).isExternal() ) {
                    bat_internal.add( b );
                }
            }
            ancestorsFirst( bat_internal, bat, "bat-phylogeny.xml UNROOTED", ok, false );
        }
    }

    /** Larger clade first: no drawn clade label is painted before a drawn ancestor's. {@code require_pair} makes
     *  the check refuse to pass vacuously: on a fixture built for it, at least one drawn label must have a drawn
     *  ancestor, or the order was never examined. */
    private static void ancestorsFirst( final List<double[]> internal, final Phylogeny phy, final String what,
                                        final boolean[] ok, final boolean require_pair ) {
        int pairs = 0;
        for( int i = 0; i < internal.size(); ++i ) {
            final PhylogenyNode n = phy.getNode( (long) internal.get( i )[ 5 ] );
            for( int j = 0; j < internal.size(); ++j ) {
                if ( ( i != j ) && isAncestor( phy.getNode( (long) internal.get( j )[ 5 ] ), n ) ) {
                    ++pairs;
                    if ( j > i ) {
                        fail( ok, what + ": clade labels must claim largest clade first, but a node's label was "
                                + "painted before its ancestor's" );
                        return;
                    }
                }
            }
        }
        if ( require_pair && ( pairs == 0 ) ) {
            fail( ok, what + ": precondition -- no drawn clade label has a drawn ancestor, so the order was never "
                    + "examined (" + internal.size() + " clade labels drawn)" );
        }
    }

    private static int depth( final PhylogenyNode n ) {
        int d = 0;
        for( PhylogenyNode p = n; !p.isRoot(); p = p.getParent() ) {
            ++d;
        }
        return d;
    }

    private static boolean isAncestor( final PhylogenyNode ancestor, final PhylogenyNode node ) {
        for( PhylogenyNode p = node.getParent(); p != null; p = p.getParent() ) {
            if ( p == ancestor ) {
                return true;
            }
        }
        return false;
    }

    /**
     * Rotating the fan (JS's line, 2026-09-27, measured here first). WHICH names survive a turn is NOT promised by
     * either program: on this uniform 96-star at 260 px a 1.28 rad turn keeps 24 either way, a DIFFERENT 24 (the
     * auto-fit excluded -- same diameter; the irregular trees swap nothing). Two things ARE promised, and are what a
     * user can hold us to: the rule is a pure function of the geometry on screen -- rotating BACK restores exactly
     * the names kept before, nothing is carried from one pass to the next -- and the rotated view does not overprint.
     * Both at the crowded 260 px and at 619 px, one pixel under the size where every name fits the ring (95 kept:
     * the last name in preorder meets the first, which was placed first).
     */
    private static void rotationOk( final TreePanel tp, final ControlPanel cp, final boolean[] ok ) {
        final double a0 = tp.getStartingAngle();
        final NODE_LABEL_DIRECTION dir0 = tp.getOptions().getNodeLabelDirection();
        try {
            tp.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
            cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
            for( final int size : new int[] { 260, 619 } ) {
                render( tp, size, false, true );
                final Set<Long> kept0 = drawnIds( tp.labelBoxesForTest() );
                if ( kept0.isEmpty() || ( kept0.size() == TIPS ) ) {
                    fail( ok, "CIRCULAR rotation @" + size + " px: precondition -- the ring must hide some names and keep "
                            + "some, kept " + kept0.size() + " of " + TIPS );
                }
                for( final double rot : new double[] { 1.28, Math.PI } ) {
                    tp.setStartingAngle( a0 + rot );
                    paint( tp, size ); // no refit: the same diameter, a rigid turn
                    refusalsCaused( tp, "CIRCULAR rotation @" + size + " px by " + rot + " rad", ok );
                    final int pairs = overlappingPairs( tp.labelBoxesForTest(), -1 );
                    if ( ( pairs != 0 ) || ( tp.labelsDrawnForTest() == 0 ) ) {
                        fail( ok, "CIRCULAR rotation @" + size + " px by " + rot + " rad: the drawn names must not overlap ("
                                + pairs + " pairs among " + tp.labelsDrawnForTest() + " drawn)" );
                    }
                    tp.setStartingAngle( a0 );
                    paint( tp, size );
                    final Set<Long> back = drawnIds( tp.labelBoxesForTest() );
                    if ( !back.equals( kept0 ) ) {
                        fail( ok, "CIRCULAR rotation @" + size + " px: turning back must restore exactly the names kept "
                                + "before -- the rule is a function of the geometry on screen, not of how the fan got "
                                + "there (" + kept0.size() + " before, " + back.size() + " after, "
                                + differing( kept0, back ) + " differ)" );
                    }
                }
            }
        }
        finally {
            tp.setStartingAngle( a0 );
            tp.getOptions().setNodeLabelDirection( dir0 );
        }
    }

    /**
     * The tooltip's own sentence, "nothing is hidden unless something already drawn is really in its way", as an
     * invariant: every box the rule REFUSED this pass overlaps a box it DREW (a name, a clade name or an image). A
     * refusal with clear space round it is a hiding without a cause -- an index rule, or a refused box recorded and
     * blocking the one after it (which chains until one name of 96 survives). It asks of each refusal "was this
     * caused by something on the screen", which is true of a refusal wherever the asking lives, so it pins the
     * property in the CALLER without anyone deciding to put it there (archaeopteryx.js's form, 2026-09-27).
     */
    private static void refusalsCaused( final TreePanel tp, final String what, final boolean[] ok ) {
        final List<double[]> drawn = tp.labelBoxesForTest();
        final List<double[]> refused = tp.refusedLabelBoxesForTest();
        int uncaused = 0;
        for( final double[] r : refused ) {
            boolean touches = false;
            for( final double[] d : drawn ) {
                if ( overlaps( r, d ) ) {
                    touches = true;
                    break;
                }
            }
            if ( !touches ) {
                ++uncaused;
            }
        }
        if ( uncaused != 0 ) {
            fail( ok, what + ": " + uncaused + " of " + refused.size() + " refused labels touch NOTHING that was drawn -- "
                    + "a label is hidden only where something already drawn is really in its way" );
        }
    }

    private static int differing( final Set<Long> a, final Set<Long> b ) {
        int n = 0;
        for( final Long id : a ) {
            if ( !b.contains( id ) ) {
                ++n;
            }
        }
        for( final Long id : b ) {
            if ( !a.contains( id ) ) {
                ++n;
            }
        }
        return n;
    }

    /** A tip whose taxonomy is shown but renders EMPTY, with no name and no image, draws nothing -- and so counts as
     *  neither drawn nor hidden, and records no box. */
    private static void emptyTaxonomyDrawsNothing( final TreePanel tp, final ControlPanel cp, final boolean[] ok ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < TIPS; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "" );
            tip.setDistanceToParent( 0.5 );
            tip.getNodeData().setTaxonomy( new org.forester.phylogeny.data.Taxonomy() ); // present, but empty
            root.addAsChild( tip );
        }
        final Phylogeny star = new Phylogeny();
        star.setRoot( root );
        star.setRooted( true );
        star.recalculateNumberOfExternalDescendants( false );
        tp.setTree( star );
        tp.recalculateMaxDistanceToRoot();
        final boolean sci0 = cp.isCheckboxSelected( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES );
        final boolean com0 = cp.isCheckboxSelected( DisplayOption.SHOW_TAXONOMY_COMMON_NAMES );
        try {
            cp.setCheckbox( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, true );
            cp.setCheckbox( DisplayOption.SHOW_TAXONOMY_COMMON_NAMES, true );
            for( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.UNROOTED,
                    PHYLOGENY_GRAPHICS_TYPE.CIRCULAR } ) {
                tp.setPhylogenyGraphicsType( gt );
                render( tp, 450, false, true );
                if ( ( tp.labelsDrawnForTest() != 0 ) || ( tp.labelsHiddenForTest() != 0 )
                        || !tp.labelBoxesForTest().isEmpty() ) {
                    fail( ok, gt + ": a tip with an EMPTY taxonomy, no name and no image draws nothing, so none is drawn or "
                            + "hidden (drawn " + tp.labelsDrawnForTest() + " hidden " + tp.labelsHiddenForTest() + ", "
                            + tp.labelBoxesForTest().size() + " boxes recorded)" );
                }
            }
        }
        finally {
            cp.setCheckbox( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, sci0 );
            cp.setCheckbox( DisplayOption.SHOW_TAXONOMY_COMMON_NAMES, com0 );
        }
    }

    /**
     * A collapsed clade's radial label goes through the rule as a tip's does -- both halves of it. With the clade
     * painted FIRST (the root's first child) it PLACES: nothing drawn after it overlaps it, and where it crowds a
     * neighbouring name that name is refused for it. With the clade painted LAST it ASKS: where the neighbours' names
     * already stand it is refused, and a search hit inside it draws it regardless. Painted with a bare drawString it
     * did neither, and names were granted the same space (a review find, 2026-09-27). The first fixture, the clade
     * bush, did not exercise the claim: its ten-tip wedges held every neighbour clear (measured).
     */
    private static void collapsedCladeLabel( final TreePanel tp, final ControlPanel cp, final boolean[] ok ) {
        tp.getOptions().setNodeLabelDirection( NODE_LABEL_DIRECTION.RADIAL );
        cp.setCheckbox( DisplayOption.DYNAMICALLY_HIDE_DATA, true );
        tp.setShowInternalDataForThisTab( false );
        for( final boolean clade_first : new boolean[] { true, false } ) {
            final Phylogeny star = starWithCollapsibleClade( clade_first );
            final PhylogenyNode x = star.getRoot().getChildNode( clade_first ? 0 : ( star.getRoot().getNumberOfDescendants() - 1 ) );
            x.setCollapse( true );
            tp.setTree( star );
            tp.recalculateMaxDistanceToRoot();
            tp.updateSetOfCollapsedExternalNodes();
            boolean caused_somewhere = false, hidden_somewhere = false;
            int drawn_somewhere = 0;
            for( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.UNROOTED,
                    PHYLOGENY_GRAPHICS_TYPE.CIRCULAR } ) {
                tp.setPhylogenyGraphicsType( gt );
                for( final int size : new int[] { 450, 300 } ) {
                    final String what = gt + " collapsed clade painted " + ( clade_first ? "first" : "last" ) + " @" + size + " px";
                    render( tp, size, false, true );
                    final double[] collapsed = boxOf( tp.labelBoxesForTest(), x.getId() );
                    final int pairs = overlappingPairs( tp.labelBoxesForTest(), -1 );
                    if ( pairs != 0 ) {
                        fail( ok, what + ": drawn names, the collapsed clade's included, must not overlap (" + pairs + " pairs)" );
                    }
                    refusalsCaused( tp, what, ok );
                    if ( collapsed == null ) {
                        if ( boxOf( tp.refusedLabelBoxesForTest(), x.getId() ) == null ) {
                            fail( ok, what + ": the collapsed clade's label was neither drawn nor refused -- it is not going "
                                    + "through the rule" );
                            continue;
                        }
                        hidden_somewhere = true;
                        // the found exemption: a hit inside the collapsed clade draws its label regardless
                        final Set<Long> found = new HashSet<Long>();
                        found.add( Long.valueOf( x.getChildNode( 0 ).getId() ) );
                        tp.setFoundNodes0( found );
                        try {
                            render( tp, size, false, true );
                            if ( boxOf( tp.labelBoxesForTest(), x.getId() ) == null ) {
                                fail( ok, what + ": a collapsed clade holding a search hit must draw its label although the "
                                        + "rule would hide it" );
                            }
                        }
                        finally {
                            tp.setFoundNodes0( null );
                        }
                        continue;
                    }
                    ++drawn_somewhere;
                    for( final double[] r : tp.refusedLabelBoxesForTest() ) {
                        caused_somewhere |= overlaps( r, collapsed );
                    }
                }
            }
            if ( clade_first && ( ( drawn_somewhere == 0 ) || !caused_somewhere ) ) {
                fail( ok, "precondition, clade painted FIRST: its label must be drawn and some later name refused FOR it "
                        + "(drawn in " + drawn_somewhere + " cases, caused a refusal: " + caused_somewhere + ")" );
            }
            if ( !clade_first && !hidden_somewhere ) {
                fail( ok, "precondition, clade painted LAST: in some layout and size the neighbours' names must already "
                        + "stand where its label would go, so that it is refused -- or the ASK is never exercised" );
            }
            x.setCollapse( false );
            tp.updateSetOfCollapsedExternalNodes();
        }
    }

    private static double[] boxOf( final List<double[]> boxes, final long id ) {
        for( final double[] b : boxes ) {
            if ( ( (long) b[ 5 ] ) == id ) {
                return b;
            }
        }
        return null;
    }

    /** Sixty named tips of alternating lengths plus a two-tip clade on short branches with a long name -- the clade to
     *  collapse in {@link #collapsedCladeLabel} -- as the root's FIRST child (painted before every tip) or its LAST. */
    private static Phylogeny starWithCollapsibleClade( final boolean clade_first ) {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode clade = new PhylogenyNode();
        clade.setName( "Collapsed_clade_with_a_long_name" );
        clade.setDistanceToParent( 0.1 );
        for( int i = 0; i < 2; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "inner_" + i );
            tip.setDistanceToParent( 0.1 );
            clade.addAsChild( tip );
        }
        if ( clade_first ) {
            root.addAsChild( clade );
        }
        for( int i = 0; i < 60; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "lineage_" + i + "_named" );
            tip.setDistanceToParent( ( ( i % 2 ) == 0 ) ? 0.2 : 0.6 );
            root.addAsChild( tip );
        }
        if ( !clade_first ) {
            root.addAsChild( clade );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.recalculateNumberOfExternalDescendants( false );
        return phy;
    }

    private static Set<Long> drawnIds( final List<double[]> boxes ) {
        final Set<Long> ids = new HashSet<Long>();
        for( final double[] b : boxes ) {
            ids.add( Long.valueOf( (long) b[ 5 ] ) );
        }
        return ids;
    }

    /**
     * Four named superclades, each with a SINGLE child clade X (a node with one child is legal phyloXML), X with two
     * clades P and Q of five tips each, every root-to-tip path the same length so all tips sit on the ring in
     * circular. A node with one child takes exactly its child's direction in BOTH layouts (unrooted: the child gets
     * the whole wedge; circular: a node's angle is the mean of its children's -- of its CHILDREN, not its tips,
     * which is why an earlier fixture that leaned on a big clade being "nearly collinear" with its parent collided
     * in unrooted and not in circular). So X's long name, starting just past the superclade's node, runs along the
     * superclade's long name and meets it: the larger clade's name, claimed first, wins and X's is hidden. P and Q
     * fan out either side and fit, and they are the drawn descendants the order check needs. Tip labels are never
     * in the way here: a clade name never displaces one, and on the ring they never crowd a clade name (a first
     * version relied on that and hid nothing in circular). A fifth superclade R has FOUR SHORT tips right past its
     * node, so that in unrooted its long name runs into its own tips' labels: under the rule (tips first) R's name
     * is the one hidden; painted inline before its subtree, as unrooted's painter would without the queue, it would
     * hide two tip names instead -- which is what the "drawn tip set is identical" check is for, and without R it
     * had nothing to see (a mutant painting clade labels inline survived the fixture without R). The superclade
     * names sit nearest the centre and are
     * kept SHORT so they end before the short tips' labels begin (a long one was hidden in circular by a tip label
     * 13 px off its line, which left the order check with no ancestor/descendant pair to look at, and a mutant
     * that dropped the depth sort survived it); they are the ancestors the order check needs, and it now insists
     * on having at least one such pair. (A first version gave every clade short tips, and the rule hid all
     * eight -- correctly.)
     */
    private static Phylogeny bushOfNamedClades() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int sc = 0; sc < 4; ++sc ) {
            final PhylogenyNode superclade = new PhylogenyNode();
            superclade.setName( "Superclade_" + sc ); // SHORT: in unrooted a name runs along its branch, and a long
                                                       // one here would cover X, P and Q downstream and hide all three
            superclade.setDistanceToParent( 0.15 );
            final PhylogenyNode x = new PhylogenyNode(); // the single child: same direction as its parent
            x.setName( "Clade_X_" + sc + "_with_a_long_family_name" );
            x.setDistanceToParent( 0.05 );
            for( int c = 0; c < 2; ++c ) {
                final PhylogenyNode clade = new PhylogenyNode();
                clade.setName( ( ( c == 0 ) ? "P_" : "Q_" ) + sc );
                clade.setDistanceToParent( 0.3 ); // far enough out that P's and Q's names clear the superclade's:
                                                  // at 0.1 they began where an 80 px superclade name ends (a graze
                                                  // the rule rightly refused, and the order check went vacuous)
                for( int i = 0; i < 5; ++i ) {
                    final PhylogenyNode tip = new PhylogenyNode();
                    tip.setName( "tip_" + sc + "_" + c + "_" + i + "_named" );
                    tip.setDistanceToParent( 0.35 ); // every path 0.15 + 0.05 + 0.3 + 0.35
                    clade.addAsChild( tip );
                }
                x.addAsChild( clade );
            }
            superclade.addAsChild( x );
            root.addAsChild( superclade );
        }
        final PhylogenyNode r = new PhylogenyNode();
        r.setName( "Superclade_R_with_short_tips_and_a_long_name" );
        r.setDistanceToParent( 0.15 );
        for( int i = 0; i < 4; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "tip_R_" + i + "_named" );
            tip.setDistanceToParent( 0.05 );
            r.addAsChild( tip );
        }
        root.addAsChild( r );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.recalculateNumberOfExternalDescendants( false );
        return phy;
    }

    /** The two-length star again, every tip carrying a two-domain architecture. */
    private static Phylogeny domainedStar() {
        final Phylogeny phy = starOfTwoLengths();
        for( final PhylogenyNode tip : phy.getExternalNodes() ) {
            final org.forester.phylogeny.data.DomainArchitecture da = new org.forester.phylogeny.data.DomainArchitecture();
            da.addDomain( new org.forester.phylogeny.data.ProteinDomain( "PF_A", 10, 80, 1e-6 ) );
            da.addDomain( new org.forester.phylogeny.data.ProteinDomain( "PF_B", 110, 190, 1e-6 ) );
            da.setTotalLength( 200 );
            final org.forester.phylogeny.data.Sequence seq = new org.forester.phylogeny.data.Sequence();
            seq.setDomainArchitecture( da );
            tip.getNodeData().setSequence( seq );
        }
        return phy;
    }

    /** A tip whose label the last paint did not draw, or -1. */
    /**
     * The found node: a label that the rule would hide is drawn all the same, and claims its place; with "Bold Found
     * Labels" on its box is measured bold. Run in BOTH radial layouts -- the ring is a different anchor from the tip,
     * and a guard written for the layout that had the problem does not grow on its own when a second layout
     * acquires it (archaeopteryx.js's find on their side, 2026-09-27). {@code size} must crowd the fan: 900 px for
     * the unrooted star, 450 px for its ring (half the names hidden).
     */
    private static void foundNode( final TreePanel tp, final MainFrame frame, final Phylogeny phy, final String what,
                                   final int size, final boolean[] ok ) {
        try {
            render( tp, size, false, true ); // auto-hide on again (crowdedFan ends on an auto-hide-off paint)
            final long hidden_id = aHiddenTip( tp, phy );
            if ( hidden_id < 0 ) {
                fail( ok, what + ": precondition -- some tip's label must have been hidden" );
            }
            else {
                final Set<Long> found = new HashSet<Long>();
                found.add( Long.valueOf( hidden_id ) );
                tp.setFoundNodes0( found );
                render( tp, size, false, true );
                if ( !drawn( tp.labelBoxesForTest(), hidden_id ) ) {
                    fail( ok, what + ": a FOUND node's label must be drawn although the rule would hide it" );
                }
                final int pairs = overlappingPairs( tp.labelBoxesForTest(), hidden_id );
                if ( pairs != 0 ) {
                    fail( ok, what + ": with a found label placed, every OTHER pair must still be clear ("
                            + pairs + " overlap)" );
                }
                // the found label claims its place: whatever is drawn AFTER it keeps clear of it (what was drawn
                // before it may not, and that is intended -- the hit is drawn on top)
                final List<double[]> boxes = tp.labelBoxesForTest();
                int found_at = -1;
                for( int i = 0; i < boxes.size(); ++i ) {
                    if ( ( (long) boxes.get( i )[ 5 ] ) == hidden_id ) {
                        found_at = i;
                    }
                }
                for( int i = found_at + 1; ( found_at >= 0 ) && ( i < boxes.size() ); ++i ) {
                    if ( overlaps( boxes.get( found_at ), boxes.get( i ) ) ) {
                        fail( ok, what + ": a label drawn after the found one lies across it -- the found "
                                + "label must claim its place" );
                        break;
                    }
                }
                // the BOLD branch of the reservation: with "Bold Found Labels" on, the hit's name is painted
                // bold and its box must be measured bold too. Until this ran, that branch was never entered by
                // any test (the option is off by default), so a box built from the regular font would have
                // passed everything -- archaeopteryx.js found the same hole on their side, 2026-09-26.
                final double plain_hw = ( found_at >= 0 ) ? boxes.get( found_at )[ 2 ] : -1;
                frame.getOptions().setBoldFoundLabels( true );
                try {
                    render( tp, size, false, true );
                    double bold_hw = -1;
                    for( final double[] b : tp.labelBoxesForTest() ) {
                        if ( ( (long) b[ 5 ] ) == hidden_id ) {
                            bold_hw = b[ 2 ];
                        }
                    }
                    if ( bold_hw <= plain_hw ) {
                        fail( ok, "precondition -- a bold hit's label must reserve a wider box than the plain one ("
                                + bold_hw + " vs " + plain_hw + "), or the bold branch was never entered" );
                    }
                    inkInsideBoxes( tp, what + ", bold search hit", ok );
                }
                finally {
                    frame.getOptions().setBoldFoundLabels( false );
                }
                tp.setFoundNodes0( null );
            }
        }
        finally {
            tp.setFoundNodes0( null );
        }
    }

    private static long aHiddenTip( final TreePanel tp, final Phylogeny phy ) {
        for( final PhylogenyNode tip : phy.getExternalNodes() ) {
            if ( !drawn( tp.labelBoxesForTest(), tip.getId() ) ) {
                return tip.getId();
            }
        }
        return -1;
    }

    private static boolean drawn( final List<double[]> boxes, final long id ) {
        for( final double[] b : boxes ) {
            if ( ( (long) b[ 5 ] ) == id ) {
                return true;
            }
        }
        return false;
    }

    /** Pairs of drawn labels whose oriented boxes overlap, leaving out every pair that includes {@code except}. */
    private static int overlappingPairs( final List<double[]> boxes, final long except ) {
        int n = 0;
        for( int i = 0; i < boxes.size(); ++i ) {
            final double[] a = boxes.get( i );
            if ( ( (long) a[ 5 ] ) == except ) {
                continue;
            }
            for( int j = i + 1; j < boxes.size(); ++j ) {
                final double[] b = boxes.get( j );
                if ( ( (long) b[ 5 ] ) == except ) {
                    continue;
                }
                if ( OrientedOccupancy.overlap( a[ 0 ], a[ 1 ], a[ 2 ], a[ 3 ], Math.cos( a[ 4 ] ), Math.sin( a[ 4 ] ),
                        b[ 0 ], b[ 1 ], b[ 2 ], b[ 3 ], Math.cos( b[ 4 ] ), Math.sin( b[ 4 ] ) ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** How many pixels {@code with} changes where {@code without} already had ink (a channel under 0xE8: a line's
     *  core or a glyph, not an antialiasing fringe). */
    private static int inkOnInk( final BufferedImage with, final BufferedImage without ) {
        int n = 0;
        for( int y = 0; y < with.getHeight(); ++y ) {
            for( int x = 0; x < with.getWidth(); ++x ) {
                final int rgb = without.getRGB( x, y );
                final boolean ink = ( ( ( rgb >> 16 ) & 0xFF ) < 0xE8 ) || ( ( ( rgb >> 8 ) & 0xFF ) < 0xE8 )
                        || ( ( rgb & 0xFF ) < 0xE8 );
                if ( ink && ( with.getRGB( x, y ) != rgb ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static BufferedImage render( final TreePanel tp, final int w, final boolean numbers,
                                         final boolean labels ) {
        tp.setShows( DisplayOption.WRITE_BRANCH_LENGTH_VALUES, numbers );
        tp.setShows( DisplayOption.SHOW_NODE_NAMES, labels );
        tp.setRecordLabelBoxesForTest( true ); // a fresh recording of THIS paint
        return paint( tp, w );
    }

    /** One more paint at {@code w} under the CURRENT recording -- so a test can watch what a second pass does to
     *  the lists the first one filled. */
    private static BufferedImage paint( final TreePanel tp, final int w ) {
        tp.setSize( w, w );
        tp.fitRadialTo( w, w );
        tp.calcParametersForPainting( w, w );
        tp.resetPreferredSize();
        final int pw = (int) Math.ceil( tp.getPreferredSize().getWidth() );
        final int ph = (int) Math.ceil( tp.getPreferredSize().getHeight() );
        final BufferedImage img = new BufferedImage( pw, ph, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
        final ExportTheme theme = ExportTheme.applyIf( tp, true );
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
     * A star of 96 named tips whose lengths alternate short and long: in the unrooted layout the short tips sit at
     * a third of the fan's radius, where the spokes are closer than a label is tall, while the long tips at the rim
     * have room -- so the rule must hide some and keep some. The long spokes' branch-length numbers sit at half
     * their length, level with the short tips' labels on the neighbouring spokes, which is where a number meets a
     * label.
     */
    private static Phylogeny starOfTwoLengths() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < TIPS; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "tip_" + i + "_named" );
            tip.setDistanceToParent( ( ( i % 2 ) == 0 ) ? 0.2 : 0.6 );
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.recalculateNumberOfExternalDescendants( false );
        return phy;
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

    private static boolean fail( final String msg ) {
        System.out.println( "  [RadialTipLabelRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [RadialTipLabelRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private RadialTipLabelRenderTest() {
    }
}

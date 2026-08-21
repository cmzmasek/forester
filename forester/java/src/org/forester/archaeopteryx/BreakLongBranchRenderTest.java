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
import java.io.File;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders forester/demo/long-branch-break.xml (an ingroup with ~0.1-0.3 branches + one ~4.0 outgroup) as a phylogram
 * with "Break Long Branches" OFF then ON, and asserts: (1) the depth scale grows (the informative part reclaims the
 * width the outgroup consumed) -- an ingroup tip is drawn much farther from the root; and (2) a break glyph appears
 * across the outgroup branch (ink above/below its horizontal segment, where there is none when off). Also covers the
 * ALIGNED phylogram, preferred-size (no ballooning), a long INTERNAL branch, and RADIAL parity -- the circular AND
 * unrooted phylograms decompress the ingroup + draw the rotated "//" glyph on the capped spoke. Headful; a green
 * no-op when headless. Dogfoods the demo.
 */
public final class BreakLongBranchRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "BreakLongBranchRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return breakRenderOk();
    }

    private static boolean breakRenderOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/long-branch-break.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "break" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    for( final PhylogenyNode tip : phy.getExternalNodes() ) {
                        tip.setName( "" ); // blank labels so the only ON-vs-OFF pixel change is the layout / glyph
                    }
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    final int w = 900, h = 620;
                    frame.showWhole();

                    // (1) decompression: capping the outgroup shrinks the effective tree height, so the depth scale
                    // (xCorrectionFactor) grows -- the informative part reclaims the width the outgroup consumed.
                    tp.setSize( w, h );
                    o.setBreakLongBranches( false );
                    tp.calcParametersForPainting( w, h );
                    final float corr_off = tp.getXcorrectionFactor();
                    o.setBreakLongBranches( true );
                    tp.calcParametersForPainting( w, h );
                    final float corr_on = tp.getXcorrectionFactor();
                    if ( corr_on <= ( corr_off * 1.8f ) ) {
                        fail( ok, "Break Long Branches must grow the depth scale (corr on=" + corr_on + " off="
                                + corr_off + ")" );
                    }

                    // (2) the break glyph, isolated without relying on node-coord<->pixel mapping: render the SAME
                    // capped layout twice -- once normally (the outlier branch drawn broken -> glyph) and once with the
                    // glyph suppressed (the PAINT_BREAK_GLYPH test seam) so the layout is byte-identical. The pixel diff
                    // is exactly the glyph, so it also mutation-verifies: drop the dispatch and with_glyph loses it too.
                    final BufferedImage with_glyph = renderAt( tp, w, h );
                    final BufferedImage no_glyph;
                    TreePanel.PAINT_BREAK_GLYPH = false;
                    try {
                        no_glyph = renderAt( tp, w, h );
                    }
                    finally {
                        TreePanel.PAINT_BREAK_GLYPH = true; // always restore the production seam, even on an exception
                    }
                    final int diff = diffPixels( with_glyph, no_glyph );
                    if ( diff < 20 ) {
                        fail( ok, "a break glyph must be drawn across a capped branch (glyph pixel diff=" + diff + ")" );
                    }

                    // (3) while capping, the full-width scale AXIS is suppressed (it would be positioned via the UNCAPPED
                    // max distance and can't represent a truncated branch) but the small scale BAR is kept.
                    o.setBreakLongBranches( false );
                    final boolean axis_applies_off = tp.scaleAxisAppliesToLayout();
                    o.setBreakLongBranches( true );
                    final boolean axis_applies_on = tp.scaleAxisAppliesToLayout();
                    if ( !axis_applies_off || axis_applies_on ) {
                        fail( ok, "the scale axis must apply on this phylogram but be suppressed while capping (off="
                                + axis_applies_off + " on=" + axis_applies_on + ")" );
                    }
                    // the scale BAR still draws while capping: with break on, toggling "Show Scale" changes the image
                    // (isShowScale is chrome, layout-identical), so the diff isolates the bar. If it were suppressed the
                    // two renders would be identical.
                    o.setShowScale( false );
                    final BufferedImage bar_off = renderAt( tp, w, h );
                    o.setShowScale( true );
                    final BufferedImage bar_on = renderAt( tp, w, h );
                    o.setShowScale( false );
                    if ( ( bar_off.getWidth() == bar_on.getWidth() ) && ( bar_off.getHeight() == bar_on.getHeight() )
                            && ( diffPixels( bar_off, bar_on ) < 10 ) ) {
                        fail( ok, "the scale bar must still draw while capping (Show-Scale on-vs-off diff="
                                + diffPixels( bar_off, bar_on ) + ")" );
                    }
                    // the VERTICAL scale-axis reserve is suppressed while capping too (mirrors the horizontal
                    // scaleAxisBottomReserve): with the axis not drawn, no phantom side band should be reserved.
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
                    o.setShowScaleAxis( true );
                    o.setBreakLongBranches( false );
                    final int vreserve_off = tp.verticalScaleAxisReserve();
                    o.setBreakLongBranches( true );
                    final int vreserve_on = tp.verticalScaleAxisReserve();
                    if ( ( vreserve_off <= 0 ) || ( vreserve_on != 0 ) ) {
                        fail( ok, "the vertical scale-axis reserve must apply off but be 0 while capping (off="
                                + vreserve_off + " on=" + vreserve_on + ")" );
                    }
                    o.setShowScaleAxis( false );
                    o.setBreakLongBranches( false );
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

                    // (4) an ALIGNED "A" phylogram is capped too: the depth scale decompresses AND the aligned label
                    // column lands ON-canvas (anchored to the CAPPED extent, not off-canvas at the uncapped outlier).
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM );
                    o.setBreakLongBranches( false );
                    tp.calcParametersForPainting( w, h );
                    tp.resetPreferredSize();
                    final float aligned_off = tp.getXcorrectionFactor();
                    o.setBreakLongBranches( true );
                    tp.calcParametersForPainting( w, h );
                    final float aligned_on = tp.getXcorrectionFactor();
                    if ( aligned_on <= ( aligned_off * 1.8f ) ) {
                        fail( ok, "an aligned phylogram must be capped too (corr off=" + aligned_off + " on="
                                + aligned_on + ")" );
                    }
                    // the aligned label column sits ~ (drawn extent)*corr + a small leader gap. With the capped-extent
                    // anchor col_x/corr ~= cappedHeight(1.6) + small; anchored to the uncapped outlier it is ~= maxDist
                    // (4.0) -- the tree squashed into the left, labels far to the right. Discriminate at 2.8x corr.
                    final double col_x = tp.alignedLabelColumnX();
                    if ( col_x > ( aligned_on * 2.8 ) ) {
                        fail( ok, "the aligned label column must track the capped tree extent, not the outlier (col_x="
                                + col_x + ", corr=" + aligned_on + ", ratio=" + ( col_x / aligned_on ) + ")" );
                    }
                    o.setBreakLongBranches( false );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );

                    // (5) the re-fit-on-toggle guard is option-INDEPENDENT (so turning the option OFF also re-fits and
                    // the uncapped layout is restored): true for this unaligned phylogram, false for a cladogram.
                    if ( !tp.breakLongBranchesRelevantToLayout() ) {
                        fail( ok, "an unaligned branch-length phylogram must be relevant to the re-fit-on-toggle" );
                    }
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
                    if ( tp.breakLongBranchesRelevantToLayout() ) {
                        fail( ok, "a cladogram must not be relevant to the re-fit-on-toggle" );
                    }
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );

                    // (6) the PREFERRED SIZE (scroll extent + overview) uses the DRAWN (capped) height -- corr * capped
                    // height -- not corr * the UNCAPPED height. The latter balloons the extent (a capped branch makes
                    // corr large while the uncapped height stays huge -> a hugely oversized scroll extent -> clipped
                    // labels / wrong overview). So break-on must not blow the preferred width far past break-off.
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    o.setBreakLongBranches( false );
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    tp.resetPreferredSize();
                    final double pref_w_off = tp.getPreferredSize().getWidth();
                    o.setBreakLongBranches( true );
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    tp.resetPreferredSize();
                    final double pref_w_on = tp.getPreferredSize().getWidth();
                    if ( pref_w_on > ( pref_w_off * 1.5 ) ) {
                        fail( ok, "the preferred size must use the reduced (capped) lengths, not balloon (off="
                                + pref_w_off + " on=" + pref_w_on + ")" );
                    }
                    o.setBreakLongBranches( false );

                    // (7) a long INTERNAL branch (a whole clade behind it), not just a tip -- the fails_8_20 case that
                    // surfaced the preferred-size ballooning. Capping it must DECOMPRESS (corr grows = the internal
                    // branch is capped, pulling its clade in) AND not balloon the preferred size.
                    tp.setTree( internalLongBranchTree() );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    o.setBreakLongBranches( false );
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    tp.resetPreferredSize();
                    final float icorr_off = tp.getXcorrectionFactor();
                    final double ipref_off = tp.getPreferredSize().getWidth();
                    o.setBreakLongBranches( true );
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    tp.resetPreferredSize();
                    final float icorr_on = tp.getXcorrectionFactor();
                    final double ipref_on = tp.getPreferredSize().getWidth();
                    if ( icorr_on <= ( icorr_off * 1.8f ) ) {
                        fail( ok, "a long INTERNAL branch must be capped -- the depth scale must decompress (corr off="
                                + icorr_off + " on=" + icorr_on + ")" );
                    }
                    if ( ipref_on > ( ipref_off * 1.5 ) ) {
                        fail( ok, "a long internal branch must not balloon the preferred size (off=" + ipref_off
                                + " on=" + ipref_on + ")" );
                    }
                    o.setBreakLongBranches( false );

                    // (8) RADIAL parity: the circular AND unrooted phylograms honour Break Long Branches. A tree with a
                    // long OUTGROUP tip (20.0) behind an ingroup of ~0.3 branches (median 0.3 -> cap 2.4): capping fans
                    // the ingroup OUT (its radius from the centre grows) while the outgroup spoke is pulled IN, and a
                    // "//" glyph is drawn on the capped spoke. Radii are taken from the ROOT node, which sits at the
                    // layout centre in both radial layouts. Both assertions together mutation-verify the two mechanisms:
                    // the absolute-radius growth catches the normalizer/urt-factor cap; the ingroup/outgroup ratio
                    // growth catches the per-spoke length cap.
                    for( final Options.PHYLOGENY_GRAPHICS_TYPE gt : new Options.PHYLOGENY_GRAPHICS_TYPE[] {
                            Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                        final Phylogeny rphy = radialBreakTree();
                        tp.setTree( rphy );
                        tp.recalculateMaxDistanceToRoot(); // production pairs setTree with this; radial depth reads it
                        tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                        tp.setPhylogenyGraphicsType( gt );
                        final PhylogenyNode a = rphy.getNode( "A" ); // a deep ingroup tip
                        final PhylogenyNode out = rphy.getNode( "OUT" ); // the long outgroup tip
                        final boolean circular = ( gt == Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                        o.setBreakLongBranches( false );
                        if ( circular ? tp.breakLongBranchesActiveCircular() : tp.breakLongBranchesActiveUnrooted() ) {
                            fail( ok, gt + ": radial capping must be inactive while the option is off" );
                        }
                        renderAt( tp, w, h ); // side effect: lays out the node coordinates at this size
                        final double a_off = radiusFromRoot( a, rphy );
                        final double out_off = radiusFromRoot( out, rphy );
                        o.setBreakLongBranches( true );
                        if ( !( circular ? tp.breakLongBranchesActiveCircular()
                                         : tp.breakLongBranchesActiveUnrooted() ) ) {
                            fail( ok, gt + ": radial capping must be active while the option is on (phylogram + lengths)" );
                        }
                        renderAt( tp, w, h );
                        final double a_on = radiusFromRoot( a, rphy );
                        final double out_on = radiusFromRoot( out, rphy );
                        if ( a_on <= ( a_off * 2.0 ) ) {
                            fail( ok, gt + ": the ingroup must decompress radially (A radius off=" + a_off + " on="
                                    + a_on + ")" );
                        }
                        final double ratio_off = ( out_off > 0 ) ? ( a_off / out_off ) : 0;
                        final double ratio_on = ( out_on > 0 ) ? ( a_on / out_on ) : 0;
                        if ( ratio_on <= ( ratio_off * 2.0 ) ) {
                            fail( ok, gt + ": the outgroup spoke must be capped (ingroup/outgroup radius ratio off="
                                    + ratio_off + " on=" + ratio_on + ")" );
                        }
                        // the "//" break glyph renders on the capped radial spoke -- isolated via the PAINT_BREAK_GLYPH
                        // seam (layout identical), so the diff is exactly the glyph (and mutation-verifies its dispatch).
                        final BufferedImage rg_with = renderAt( tp, w, h );
                        final BufferedImage rg_without;
                        TreePanel.PAINT_BREAK_GLYPH = false;
                        try {
                            rg_without = renderAt( tp, w, h );
                        }
                        finally {
                            TreePanel.PAINT_BREAK_GLYPH = true;
                        }
                        final int rdiff = diffPixels( rg_with, rg_without );
                        if ( rdiff < 20 ) {
                            fail( ok, gt + ": a break glyph must be drawn on the capped radial spoke (diff=" + rdiff
                                    + ")" );
                        }
                        if ( !tp.breakLongBranchesRelevantToRadialLayout() ) {
                            fail( ok, gt + ": a radial branch-length phylogram must be relevant to the radial re-fit" );
                        }
                        o.setBreakLongBranches( false );
                    }

                    // (9) the UNROOTED OVERVIEW thumbnail caps too (mirrors paintUnrooted). The export render path skips
                    // the overview, so drive a SCREEN paint (paintComponent) with the tree larger than the viewport so
                    // the overview turns on, and read the thumbnail coordinates (XSecondary, taken from the root, which
                    // sits at the overview centre). With break on the outlier's overview spoke is pulled IN (the
                    // thumbnail ingroup/outgroup radius ratio grows) -- else the thumbnail shows the outlier dominating
                    // while the main view caps it (the fossils_8_20-class overview inconsistency, radial edition).
                    {
                        final Phylogeny rphy = radialBreakTree();
                        tp.setTree( rphy );
                        tp.recalculateMaxDistanceToRoot();
                        tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                        final PhylogenyNode a = rphy.getNode( "A" );
                        final PhylogenyNode out = rphy.getNode( "OUT" );
                        o.setShowOverview( true );
                        final int ow = 1400, oh = 1400;
                        frame.showWhole();
                        tp.setSize( ow, oh ); // larger than the frame viewport (updateOvSizes turns the overview on)
                        o.setBreakLongBranches( false );
                        tp.calcParametersForPainting( ow, oh ); // runs updateOvSizes -> isOvOn, then draws the thumbnail
                        paintScreen( tp, ow, oh );
                        if ( !tp.isOvOn() ) {
                            fail( ok, "precondition: the overview must be ON so the unrooted thumbnail is exercised" );
                        }
                        final double ov_ratio_off = ovRatio( a, out, rphy );
                        o.setBreakLongBranches( true );
                        tp.calcParametersForPainting( ow, oh );
                        paintScreen( tp, ow, oh );
                        final double ov_ratio_on = ovRatio( a, out, rphy );
                        if ( ov_ratio_on <= ( ov_ratio_off * 2.0 ) ) {
                            fail( ok, "the unrooted overview must cap the outlier too (thumbnail ingroup/outgroup ratio "
                                    + "off=" + ov_ratio_off + " on=" + ov_ratio_on + ")" );
                        }
                        o.setShowOverview( false );
                        o.setBreakLongBranches( false );
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** ((Homo:0.18,Pan:0.15):0.2, (Mus:0.25,Rat:0.2):10.0, out:0.5) -- a long INTERNAL branch (10.0) behind a clade
     *  (the fails_8_20 shape). median {..} = 0.2 -> cap 1.6, so the 10.0 internal branch caps. */
    private static Phylogeny internalLongBranchTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode c1 = branchChild( root, null, 0.2 );
        branchChild( c1, "Homo", 0.18 );
        branchChild( c1, "Pan", 0.15 );
        final PhylogenyNode c2 = branchChild( root, null, 10.0 ); // the long internal branch
        branchChild( c2, "Mus", 0.25 );
        branchChild( c2, "Rat", 0.2 );
        branchChild( root, "out", 0.5 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /** ((A:0.3,B:0.3):0.3, OUT:20.0) -- a long OUTGROUP tip behind an ingroup of ~0.3 branches. median {0.3,0.3,0.3,20}
     *  = 0.3 -> cap 2.4, so the 20.0 outgroup spoke caps hard (the radial decompression is dramatic, ~8x). */
    private static Phylogeny radialBreakTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode c = branchChild( root, null, 0.3 );
        branchChild( c, "A", 0.3 );
        branchChild( c, "B", 0.3 );
        branchChild( root, "OUT", 20.0 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    /** Pixel distance from a node to the tree's root node -- the root sits at the layout centre in both radial
     *  layouts (circular: fraction 0 -> centre; unrooted: getWidth()/2, getHeight()/2), so this is the node's radius. */
    private static double radiusFromRoot( final PhylogenyNode n, final Phylogeny phy ) {
        final PhylogenyNode root = phy.getRoot();
        final double dx = n.getXcoord() - root.getXcoord();
        final double dy = n.getYcoord() - root.getYcoord();
        return Math.sqrt( ( dx * dx ) + ( dy * dy ) );
    }

    /** Ingroup/outgroup radius ratio in the OVERVIEW thumbnail (from the XSecondary coords set by paintUnrootedLite;
     *  the root's XSecondary is the overview centre). Grows when the outlier spoke is capped in the thumbnail. */
    private static double ovRatio( final PhylogenyNode ingroup, final PhylogenyNode outgroup, final Phylogeny phy ) {
        final PhylogenyNode root = phy.getRoot();
        final double ai = Math.hypot( ingroup.getXSecondary() - root.getXSecondary(),
                ingroup.getYSecondary() - root.getYSecondary() );
        final double ao = Math.hypot( outgroup.getXSecondary() - root.getXSecondary(),
                outgroup.getYSecondary() - root.getYSecondary() );
        return ( ao > 0 ) ? ( ai / ao ) : 0;
    }

    /** Drive the on-screen paint (paintComponent) so the overview thumbnail is drawn (the export path skips it). */
    private static void paintScreen( final TreePanel tp, final int w, final int h ) {
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            tp.paintComponent( g );
        }
        finally {
            g.dispose();
        }
    }

    private static PhylogenyNode branchChild( final PhylogenyNode parent, final String name, final double dtp ) {
        final PhylogenyNode n = new PhylogenyNode();
        if ( name != null ) {
            n.setName( name );
        }
        n.setDistanceToParent( dtp );
        parent.addAsChild( n );
        return n;
    }

    /** Lay the tree out at w x h and paint it directly into an image in the white export theme, leaving the node
     *  coordinates at THAT size (unlike renderPhylogenyToImage, whose trailing showWhole re-fits them to the panel). */
    private static BufferedImage renderAt( final TreePanel tp, final int w, final int h ) {
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        tp.resetPreferredSize(); // the tree lays out at its PREFERRED size; render there so node coords map 1:1
        final int pw = (int) Math.ceil( tp.getPreferredSize().getWidth() );
        final int ph = (int) Math.ceil( tp.getPreferredSize().getHeight() );
        final BufferedImage img = new BufferedImage( pw, ph, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON );
        final ExportTheme theme = ExportTheme.applyIf( tp, true ); // predictable white background for the pixel checks
        try {
            tp.paintPhylogeny( g, false, true, pw, ph, 0, 0 );
        }
        finally {
            theme.restore();
            g.dispose();
        }
        return img;
    }

    /** Number of pixels that differ between two same-size images. */
    private static int diffPixels( final BufferedImage a, final BufferedImage b ) {
        int n = 0;
        for( int y = 0; y < a.getHeight(); ++y ) {
            for( int x = 0; x < a.getWidth(); ++x ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [BreakLongBranchRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [BreakLongBranchRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private BreakLongBranchRenderTest() {
    }
}

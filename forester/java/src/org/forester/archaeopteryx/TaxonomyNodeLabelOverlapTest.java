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

import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Guards that an external node's italic scientific name and the node-data label drawn after it do NOT overlap.
 * Both label segments must start from the SAME horizontal offset ({@link TreePanel#labelSegmentStartX}); when they
 * differed (taxonomy +3, node data +2) the node data began a pixel inside the taxonomy box, so an italic scientific
 * name's right overhang bled into the following node name.
 *
 * Part 1 (always runs) pins the shared label-start offset the two segments use. Part 2 (needs a display) renders
 * [italic scientific name][node name] with the two segments forced to distinct colors and asserts the rightmost
 * taxonomy pixel is strictly left of the leftmost node-name pixel. Antialiasing is OFF so each glyph pixel is
 * exactly one segment's color: this catches glyph-BODY overlap (the offset regression this guards), though a
 * hair-thin anti-aliased overhang alone would round away -- the offset fix is what keeps the bodies apart.
 * Guarded to a no-op when headless.
 */
public final class TaxonomyNodeLabelOverlapTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TaxonomyNodeLabelOverlap: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = testAbutInvariant();
        if ( !GraphicsEnvironment.isHeadless() ) {
            ok &= testNoRenderedOverlap( "Homo sapiens" );
            ok &= testNoRenderedOverlap( "Drosophila melanogasterf" ); // ends in an italic 'f': largest right overhang
            ok &= testInternalGapMatchesExternal();
        }
        return ok;
    }

    /**
     * Pins the shared label-start offset that both segments rely on. labelSegmentStartX is affine in prior_width
     * (slope 1), so the node-data start (prior = taxonomy end) equals the taxonomy start plus the taxonomy width --
     * the abutment the taxonomy and node-data paths get by calling this ONE helper. This part pins the helper's
     * concrete offset and arithmetic (a change to LABEL_GAP_AFTER_NODE_SHAPE or the formula fails here); it does
     * NOT by itself prove the paint methods call the helper -- testNoRenderedOverlap (headful) guards the painted
     * result, i.e. that paintTaxonomy and paintNodeData actually abut.
     */
    private static boolean testAbutInvariant() {
        boolean ok = true;
        final float base = 137.4f;
        final int half = 3;
        final float prior = 11f; // x already consumed before the taxonomy segment
        final float taxo_w = 58f;
        final float taxo_start = TreePanel.labelSegmentStartX( base, half, prior );
        final float data_start = TreePanel.labelSegmentStartX( base, half, prior + taxo_w );
        if ( data_start != taxo_start + taxo_w ) {
            ok = fail( "node data must begin exactly where the taxonomy label ends: taxoStart=" + taxo_start
                    + " + taxoW=" + taxo_w + " != dataStart=" + data_start );
        }
        // Pin the concrete offset: node-shape half + a 2px gap, with no prior segment. This is the value the OLD
        // bug got wrong (taxonomy used +3, node data +2); both must now be this same +2.
        if ( TreePanel.labelSegmentStartX( base, half, 0f ) != ( base + half + 2f ) ) {
            ok = fail( "the shared label-start offset must be (node-shape half + 2px), got "
                    + TreePanel.labelSegmentStartX( base, half, 0f ) + " != " + ( base + half + 2f ) );
        }
        if ( !( taxo_start > base + half ) ) {
            ok = fail( "a label segment must start beyond the node shape (base+half=" + ( base + half )
                    + "), got " + taxo_start );
        }
        return ok;
    }

    private static boolean testNoRenderedOverlap( final String sci ) {
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( sci ) }, conf, "ov" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    // show the scientific name + the sequence name (the two labels whose overlap is under test)
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, true );
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_SEQ_NAMES, true );
                    final Options o = frame.getOptions();
                    o.setAntialiasPrint( false ); // crisp single-color glyphs, so a pixel is exactly one segment's color
                    o.setUseItalicScientificNames( true );
                    o.setGraphicsExportWhiteBackground( false ); // keep our forced colors (no light-theme override)
                    // Markers deliberately magenta / cyan, NOT pure red/green/blue/yellow: a scheme's own
                    // found-node and duplication-box colors are 0xFF0000 / 0x00FF00 / 0x0000FF etc., so a
                    // red/blue marker could be confused with such an element; magenta and cyan are in no scheme.
                    final int tax_rgb = 0xFF00FF, seq_rgb = 0x00FFFF;
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.TAXONOMY, new Color( tax_rgb ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.SEQUENCE, new Color( seq_rgb ) );
                    tp.getTreeColorSet().setColorSchema( 0 ); // load the just-set scheme-0 colors into the active fields
                    frame.showWhole();
                    final int w = 900, h = 400;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    int tax_max_x = -1, seq_min_x = Integer.MAX_VALUE;
                    for( int y = 0; y < img.getHeight(); ++y ) {
                        for( int x = 0; x < img.getWidth(); ++x ) {
                            final int p = img.getRGB( x, y ) & 0xFFFFFF;
                            if ( ( p == tax_rgb ) && ( x > tax_max_x ) ) {
                                tax_max_x = x;
                            }
                            if ( ( p == seq_rgb ) && ( x < seq_min_x ) ) {
                                seq_min_x = x;
                            }
                        }
                    }
                    // both segments must actually have been painted, else the test would pass vacuously
                    if ( ( tax_max_x < 0 ) || ( seq_min_x == Integer.MAX_VALUE ) ) {
                        fail( ok, "'" + sci + "': expected both a taxonomy and a node-name segment to be drawn "
                                + "(taxMaxX=" + tax_max_x + ", seqMinX=" + seq_min_x + ")" );
                    }
                    else if ( seq_min_x <= tax_max_x ) {
                        fail( ok, "'" + sci + "': the node name overlaps the italic scientific name -- leftmost "
                                + "node-name pixel " + seq_min_x + " is not right of rightmost taxonomy pixel "
                                + tax_max_x );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "'" + sci + "': unexpected " + t );
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

    /**
     * The taxonomy-&gt;node-data gap of an INTERNAL node's above-the-branch label must be no wider than an EXTERNAL
     * node's, and must NOT overlap. The internal path once added INTERNAL_LABEL_SEGMENT_GAP on top of the trailing
     * space the taxonomy already carries (a double gap), so internal labels read too spread out. Comparing internal
     * TO external (same font, same name) cancels font/rasterizer differences, so the gap DIFFERENCE is the robust,
     * platform-independent signal. Also run with the worst-case italic-'f' overhang, to confirm gap=0 did not tip
     * the internal label into overlap (the 5px used to buffer it).
     */
    private static boolean testInternalGapMatchesExternal() {
        FontResources.registerBundledFonts(); // pin the bundled font so glyph metrics are reproducible
        boolean ok = internalGapOk( "Homo sapiens" );
        ok &= internalGapOk( "Drosophila melanogasterf" ); // italic 'f': the largest right overhang
        return ok;
    }

    private static boolean internalGapOk( final String sci ) {
        try {
            final Configuration conf = new Configuration();
            // root -> [ A(non-root internal) -> (C,D), B ]; every node carries a scientific name + a sequence name
            final Phylogeny phy = new Phylogeny();
            final PhylogenyNode root = new PhylogenyNode();
            final PhylogenyNode a = withData( "Aint", sci );
            final PhylogenyNode c = withData( "Cext", sci );
            a.addAsChild( c );
            a.addAsChild( withData( "Dext", sci ) );
            root.addAsChild( a );
            root.addAsChild( withData( "Bext", sci ) );
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "gap" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, true );
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_SEQ_NAMES, true );
                    final Options o = frame.getOptions();
                    o.setAntialiasPrint( false );
                    o.setUseItalicScientificNames( true );
                    o.setGraphicsExportWhiteBackground( false );
                    final int tax = 0xFF00FF, seq = 0x00FFFF;
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.TAXONOMY, new Color( tax ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.SEQUENCE, new Color( seq ) );
                    tp.getTreeColorSet().setColorSchema( 0 );
                    frame.showWhole();
                    final int w = 1100, h = 500;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    // band derived from the ACTUAL label-font height: the internal label sits ABOVE its node, so a
                    // fixed symmetric band would clip its top on a large/HiDPI font and skew the measurement; a whole
                    // font-height +margin above/below covers the label yet stays well inside the ~80px row spacing
                    final int band = tp.getMainPanel().getTreeFontSet().getFontMetricsLarge().getHeight() + 4;
                    final int internal_gap = columnGap( img, a.getYcoord(), band, tax, seq ); // A is the internal node
                    final int external_gap = columnGap( img, c.getYcoord(), band, tax, seq ); // C is an external leaf
                    if ( ( internal_gap == ABSENT ) || ( external_gap == ABSENT ) ) {
                        fail( ok, "'" + sci + "': internal/external label not both drawn (internal=" + internal_gap
                                + " external=" + external_gap + ")" );
                    }
                    else if ( internal_gap <= 0 ) {
                        fail( ok, "'" + sci + "': internal taxonomy/name OVERLAP (gap=" + internal_gap + ")" );
                    }
                    // the internal gap must not exceed the external one by more than a couple px (the old bug added a
                    // whole extra INTERNAL_LABEL_SEGMENT_GAP -> visibly too large)
                    else if ( internal_gap > external_gap + 3 ) {
                        fail( ok, "'" + sci + "': internal taxonomy->name gap (" + internal_gap
                                + ") is too large vs external (" + external_gap + ")" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "'" + sci + "': internal-gap unexpected " + t );
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

    /** Sentinel for {@link #columnGap}: one of the two colors was absent from the window. */
    private static final int ABSENT = Integer.MIN_VALUE;

    /**
     * The horizontal gap (px) in a +/-{@code band}px window around {@code yc} between the rightmost taxonomy-colored
     * ink (tax_max) and the GLOBAL leftmost node-data-colored ink (seq_min): {@code seq_min - tax_max}. For a
     * [taxonomy][data] label the data is entirely right of the taxonomy, so a positive result is clear space and a
     * result &lt;= 0 means the two segments OVERLAP (this is why seq_min is the global leftmost, not restricted to
     * the right of tax_max -- otherwise an overlap could never be detected). Returns {@link #ABSENT} if either color
     * is missing from the window.
     */
    private static int columnGap( final BufferedImage img, final float yc, final int band, final int tax,
                                  final int seq ) {
        final int y0 = Math.max( 0, (int) yc - band );
        final int y1 = Math.min( img.getHeight(), (int) yc + band );
        int tax_max = -1;
        int seq_min = Integer.MAX_VALUE;
        for( int y = y0; y < y1; ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int p = img.getRGB( x, y ) & 0xFFFFFF;
                if ( p == tax ) {
                    tax_max = Math.max( tax_max, x );
                }
                else if ( p == seq ) {
                    seq_min = Math.min( seq_min, x );
                }
            }
        }
        return ( ( tax_max < 0 ) || ( seq_min == Integer.MAX_VALUE ) ) ? ABSENT : ( seq_min - tax_max );
    }

    private static PhylogenyNode withData( final String name, final String sci ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        n.getNodeData().setTaxonomy( t );
        final Sequence s = new Sequence();
        s.setName( "SEQ" );
        n.getNodeData().addSequence( s );
        return n;
    }

    /** Cladogram (tips aligned at one x) of leaves that each carry an italic scientific name and a node name. */
    private static Phylogeny tree( final String sci ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 3; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "n" + i );
            final Taxonomy t = new Taxonomy();
            t.setScientificName( sci );
            leaf.getNodeData().setTaxonomy( t );
            final Sequence s = new Sequence();
            s.setName( "ABCDEF" );
            leaf.getNodeData().addSequence( s );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [TaxonomyNodeLabelOverlapTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [TaxonomyNodeLabelOverlapTest] " + msg );
        ok[ 0 ] = false;
    }

    private TaxonomyNodeLabelOverlapTest() {
    }
}

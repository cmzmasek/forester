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
            conf.setDisplayTaxonomyScientificNames( true );
            conf.setDisplaySequenceNames( true );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( sci ) }, conf, "ov" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
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

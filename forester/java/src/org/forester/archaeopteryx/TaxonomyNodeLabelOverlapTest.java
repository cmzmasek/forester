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
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.JPanel;
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
 *
 * At HIGH RESOLUTION (a 2x Retina screen, the 4x raster export, 300 dpi) a label part is drawn at its
 * near-fractional width while the component's FontMetrics measure it at 1x: an italic name ran up to 3.5 px past
 * its measure and swallowed the space before the node data (the Peribunyaviridae 4x PNG, 2026-10-29), and a
 * right-aligned internal label ran into its branch. {@code testDrawnAdvanceWidth} (always runs) pins the measure
 * against INK; {@code testHiResLabelsMatchScaleOne} (headful) renders through the export at each scale against 1x.
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
        // after the checks above: they look for exact pure-colour pixels, which a 1x italic in the bundled font
        // never has (text is antialiased here even with antialiasing off), and the next line registers that font
        ok &= testDrawnAdvanceWidth();
        if ( !GraphicsEnvironment.isHeadless() ) {
            for( final double s : HI_RES_SCALES ) {
                ok &= testHiResLabelsMatchScaleOne( s );
            }
            ok &= testRadialOneSpaceAfterName( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            ok &= testRadialOneSpaceAfterName( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        }
        return ok;
    }

    /** A Retina screen, the default 4x raster export, a 300 dpi fixed-size export. */
    private static final double[] HI_RES_SCALES   = { 2.0, 4.0, 300.0 / 72.0 };
    /** Tip names from the Peribunyaviridae tree whose 4x PNG lost the space after the name (2026-10-29). Termeil
     *  is the control: its 1x advance is no shorter than its drawn one at 4x, so it kept its gap even then. */
    private static final String[]  TIP_NAMES       = { "Guaroa virus", "Yacaaba virus", "Shamonda virus",
            "Akabane virus", "Termeil virus" };
    private static final String    INTERNAL_NAME   = "Orthobunyavirus akabaneense";
    /** Roman letters whose 1x advance is short of the drawn one at 2x / 4x / 300 dpi in Source Sans 3 12 (about
     *  +5 px in all); the roman words of real metadata ("Peru", "Australia") drift the other way, so they could not
     *  show a short measure of the node-data segment. */
    private static final String    INTERNAL_DATA   = "SDGTWZSDGT";
    private static final String    TIP_DATA        = "Peru";

    private static Font bundledFont( final int style ) {
        FontResources.registerBundledFonts();
        return new Font( FontResources.DEFAULT_FIGURE_FONT, style, 12 );
    }

    /**
     * {@link TreePanel#drawnAdvanceWidth} is what a label part ADVANCES when drawn: measured by INK (the left edge
     * of a bar drawn after the text in the same call, minus the bar drawn alone), never by asking the font again.
     * At 1x it equals the component's integer metrics exactly, so the screen at 1x and a scale-1 raster are
     * unchanged. At 2x / 4x / 300 dpi it covers the drawn advance and rounds up by less than a pixel, where the
     * component metrics fall short (the precondition: by 2 px or more for some name, the space they swallowed).
     */
    private static boolean testDrawnAdvanceWidth() {
        boolean ok = true;
        final Font it = bundledFont( Font.ITALIC );
        if ( !FontResources.DEFAULT_FIGURE_FONT.equals( it.getFamily() ) ) {
            return fail( "the bundled " + FontResources.DEFAULT_FIGURE_FONT + " is not registered (got "
                    + it.getFamily() + "): the measured fixtures below would be meaningless" );
        }
        final FontMetrics component = new JPanel().getFontMetrics( it );
        for( final boolean aa : new boolean[] { false, true } ) {
            for( final String s : TIP_NAMES ) {
                final Graphics2D g = scaledGraphics( 1.0, aa );
                final int w = TreePanel.drawnAdvanceWidth( g, it, s );
                g.dispose();
                if ( w != component.stringWidth( s ) ) {
                    ok = fail( "at 1x (aa=" + aa + ") '" + s + "' must measure as the component metrics do: "
                            + w + " != " + component.stringWidth( s ) );
                }
                final double ink = inkAdvance( it, s, 1.0, aa );
                if ( Math.abs( ink - w ) > 1e-9 ) {
                    ok = fail( "at 1x (aa=" + aa + ") '" + s + "' advances " + ink + " px, measured " + w );
                }
            }
        }
        for( final double scale : HI_RES_SCALES ) {
            double worst_short = 0;
            for( final String s : TIP_NAMES ) {
                final double ink = inkAdvance( it, s, scale, false );
                final Graphics2D g = scaledGraphics( scale, false );
                final int w = TreePanel.drawnAdvanceWidth( g, it, s );
                g.dispose();
                if ( ( w < ink - 1e-6 ) || ( w >= ink + 1.0 + ( 1.0 / scale ) ) ) {
                    ok = fail( "at " + scale + "x '" + s + "' advances " + ink + " px by ink, measured " + w
                            + " (must cover it, by less than a pixel)" );
                }
                worst_short = Math.max( worst_short, ink - component.stringWidth( s ) );
            }
            if ( worst_short < 2.0 ) {
                ok = fail( "precondition at " + scale + "x: no fixture name is drawn 2 px wider than the 1x metrics"
                        + " say (worst " + worst_short + "), so nothing here would show a short measure" );
            }
        }
        return ok;
    }

    private static Graphics2D scaledGraphics( final double scale, final boolean aa ) {
        final Graphics2D g = new BufferedImage( 4, 4, BufferedImage.TYPE_INT_RGB ).createGraphics();
        g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, aa ? RenderingHints.VALUE_TEXT_ANTIALIAS_ON
                : RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
        g.scale( scale, scale );
        return g;
    }

    /** The advance of {@code s} as drawString lays it out, from ink: a bar drawn after s in the SAME call stands
     *  exactly one advance right of the bar drawn alone. User-space px, to 1/scale. */
    private static double inkAdvance( final Font f, final String s, final double scale, final boolean aa ) {
        return barLeft( f, s + "|", scale, aa ) - barLeft( f, "|", scale, aa );
    }

    private static double barLeft( final Font f, final String text, final double scale, final boolean aa ) {
        final int w = (int) Math.ceil( 300 * scale ), h = (int) Math.ceil( 30 * scale );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        g.setRenderingHint( RenderingHints.KEY_TEXT_ANTIALIASING, aa ? RenderingHints.VALUE_TEXT_ANTIALIAS_ON
                : RenderingHints.VALUE_TEXT_ANTIALIAS_OFF );
        g.scale( scale, scale );
        g.setColor( Color.BLACK );
        g.setFont( f );
        g.drawString( text, 10f, 22f );
        g.dispose();
        // the bar is the rightmost ink: walk left from its right edge to the first blank column
        int x = w - 1;
        while ( ( x >= 0 ) && !inkColumn( img, x ) ) {
            --x;
        }
        while ( ( x > 0 ) && inkColumn( img, x - 1 ) ) {
            --x;
        }
        return x / scale;
    }

    private static boolean inkColumn( final BufferedImage img, final int x ) {
        for( int y = 0; y < img.getHeight(); ++y ) {
            if ( ( img.getRGB( x, y ) & 0xFF ) < 128 ) {
                return true;
            }
        }
        return false;
    }

    /**
     * A raster at {@code scale} places every label segment where the 1x raster does, scaled -- through the app's
     * own export ({@link AptxUtil#renderPhylogenyToImage}). Two placements the 1x component metrics broke at high
     * resolution (the 4x PNG of 2026-10-29): (1) a tip's node data starts after the italic name's space, not inside
     * it -- the clear gap between the two inks, in user px, may not shrink by more than a pixel; (2) an internal
     * label right-aligned to its node (taxonomy, then node data) ends where it ends at 1x -- not past it into the
     * branch. The precondition: the fixture's names are drawn at least 2 px wider than the 1x metrics say at this
     * scale, so the old measure would fail both.
     */
    private static boolean testHiResLabelsMatchScaleOne( final double scale ) {
        try {
            final Configuration conf = new Configuration();
            final Phylogeny phy = new Phylogeny();
            final PhylogenyNode root = new PhylogenyNode();
            final PhylogenyNode a = withData( "A", INTERNAL_NAME, INTERNAL_DATA );
            final PhylogenyNode[] tips = new PhylogenyNode[ TIP_NAMES.length ];
            for( int i = 0; i < TIP_NAMES.length; ++i ) {
                tips[ i ] = withData( "t" + i, TIP_NAMES[ i ], TIP_DATA );
                ( i < 2 ? a : root ).addAsChild( tips[ i ] );
            }
            root.addAsChild( a );
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "hires" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, true );
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_SEQ_NAMES, true );
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_NODE_NAMES, false );
                    final Options o = frame.getOptions();
                    o.setAntialiasExport( false ); // a pixel is exactly one segment's colour
                    o.setUseItalicScientificNames( true );
                    o.setInternalLabelsAboveBranch( true );
                    o.setGraphicsExportWhiteBackground( false );
                    final int tax = 0xFF00FF, seq = 0x00FFFF;
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.TAXONOMY, new Color( tax ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.SEQUENCE, new Color( seq ) );
                    tp.getTreeColorSet().setColorSchema( 0 );
                    // the measured fixture is Source Sans 3 at 12: at 11 every one of these names is drawn NARROWER
                    // than its 1x metrics, and the test would see nothing
                    frame.getMainPanel().getTreeFontSet().setBaseFont( bundledFont( Font.PLAIN ) );
                    frame.showWhole();
                    final int w = 1100, h = 700;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final Font large = frame.getMainPanel().getTreeFontSet().getLargeFont();
                    if ( !FontResources.DEFAULT_FIGURE_FONT.equals( large.getFamily() ) || ( large.getSize() != 12 ) ) {
                        fail( ok, scale + "x: the tree font did not stay " + FontResources.DEFAULT_FIGURE_FONT
                                + " 12 (got " + large + ")" );
                        return;
                    }
                    // precondition: at this scale the old measure is short by 2 px or more, for a tip name, the
                    // internal name and the internal node data alike
                    final Font italic = large.deriveFont( large.getStyle() | Font.ITALIC );
                    final Graphics2D probe = scaledGraphics( scale, false );
                    double tip_short = 0;
                    for( final String s : TIP_NAMES ) {
                        tip_short = Math.max( tip_short, TreePanel.drawnAdvanceWidth( probe, italic, s )
                                - tp.getFontMetrics( italic ).stringWidth( s ) );
                    }
                    final int internal_short = TreePanel.drawnAdvanceWidth( probe, italic, INTERNAL_NAME )
                            - tp.getFontMetrics( italic ).stringWidth( INTERNAL_NAME );
                    final int data_short = TreePanel.drawnAdvanceWidth( probe, large, INTERNAL_DATA )
                            - tp.getFontMetrics( large ).stringWidth( INTERNAL_DATA );
                    probe.dispose();
                    if ( ( tip_short < 2 ) || ( internal_short < 2 ) || ( data_short < 2 ) ) {
                        fail( ok, "precondition at " + scale + "x: the 1x metrics must fall 2 px short (tip names "
                                + tip_short + ", internal name " + internal_short + ", internal data "
                                + data_short + ")" );
                        return;
                    }
                    final BufferedImage one = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final BufferedImage hi = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, scale, false );
                    final double actual_scale = (double) hi.getWidth() / one.getWidth();
                    if ( Math.abs( actual_scale - scale ) > 0.01 ) {
                        fail( ok, "the export was not rendered at " + scale + "x (got " + actual_scale + ")" );
                        return;
                    }
                    final int band = large.getSize() + 8;
                    for( int i = 0; i < tips.length; ++i ) {
                        final double g1 = gapUser( one, 1.0, tips[ i ].getYcoord(), band, tax, seq );
                        final double gs = gapUser( hi, scale, tips[ i ].getYcoord(), band, tax, seq );
                        if ( Double.isNaN( g1 ) || Double.isNaN( gs ) ) {
                            fail( ok, scale + "x '" + TIP_NAMES[ i ] + "': name and node data not both drawn" );
                        }
                        else if ( gs < g1 - 1.0 ) {
                            fail( ok, scale + "x '" + TIP_NAMES[ i ] + "': the gap after the italic name shrank from "
                                    + g1 + " px at 1x to " + gs + " px" );
                        }
                    }
                    // the internal label sits ABOVE its branch: a band from a font height above the node to just
                    // below it, clear of the tip rows (about 60 px away)
                    for( final int colour : new int[] { tax, seq } ) {
                        final double e1 = rightmostInkUser( one, 1.0, a.getYcoord(), band, colour );
                        final double es = rightmostInkUser( hi, scale, a.getYcoord(), band, colour );
                        final String what = ( colour == tax ) ? "name" : "node data";
                        if ( Double.isNaN( e1 ) || Double.isNaN( es ) ) {
                            fail( ok, scale + "x: the internal label's " + what + " was not drawn" );
                        }
                        else if ( es > e1 + 1.0 ) {
                            fail( ok, scale + "x: the internal label's " + what + " ends at " + es
                                    + " px, past where it ends at 1x (" + e1 + "), toward its node" );
                        }
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, scale + "x: unexpected " + t );
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
     * Circular and unrooted put ONE space between the scientific name and the node data, as the rectangular label
     * does: the taxonomy label ends with its own separator, so the node-data string drops its leading one. Read off
     * the box each drawn label claims (its width is the width of exactly what is drawn: the taxonomy parts, italic
     * name in italics, then the node-data string). A tip WITHOUT a taxonomy keeps its leading space, so its label
     * stands where it stood -- the deliberate non-change. Until 2026-10-29 the radial label had two spaces there;
     * at high resolution the short measure of the italic name ate about one of them, which hid it.
     */
    private static boolean testRadialOneSpaceAfterName( final Options.PHYLOGENY_GRAPHICS_TYPE type ) {
        try {
            final Configuration conf = new Configuration();
            final Phylogeny phy = new Phylogeny();
            final PhylogenyNode root = new PhylogenyNode();
            for( final String sci : TIP_NAMES ) {
                root.addAsChild( withData( "x", sci, TIP_DATA ) );
            }
            for( int i = 0; i < 2; ++i ) { // the controls: node data alone
                final PhylogenyNode n = new PhylogenyNode();
                final Sequence seq = new Sequence();
                seq.setName( TIP_DATA + i );
                n.getNodeData().addSequence( seq );
                root.addAsChild( n );
            }
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "radial" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, true );
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_SEQ_NAMES, true );
                    tp.getControlPanel().setCheckbox( DisplayOption.SHOW_NODE_NAMES, false );
                    frame.getOptions().setUseItalicScientificNames( true );
                    frame.getMainPanel().getTreeFontSet().setBaseFont( bundledFont( Font.PLAIN ) );
                    tp.setPhylogenyGraphicsType( type );
                    final int w = 1100, h = 1100;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    tp.setRecordLabelBoxesForTest( true );
                    AptxUtil.renderPhylogenyToImage( w, h, tp, frame.getOptions(), false, 1, false );
                    final java.util.List<double[]> boxes = tp.labelBoxesForTest();
                    tp.setRecordLabelBoxesForTest( false );
                    final Font large = frame.getMainPanel().getTreeFontSet().getLargeFont();
                    final FontMetrics roman = tp.getFontMetrics( large );
                    final FontMetrics italic = tp.getFontMetrics( large.deriveFont( large.getStyle() | Font.ITALIC ) );
                    if ( ( boxes == null ) || ( boxes.size() != phy.getNumberOfExternalNodes() ) ) {
                        fail( ok, type + ": every tip label must be drawn (drew "
                                + ( boxes == null ? "none" : boxes.size() ) + " of " + phy.getNumberOfExternalNodes()
                                + ")" );
                        return;
                    }
                    int with_taxonomy = 0, without = 0;
                    for( final double[] b : boxes ) {
                        final PhylogenyNode n = phy.getNode( (long) b[ 5 ] );
                        final String data = n.getNodeData().getSequence().getName();
                        final int expected;
                        if ( n.getNodeData().isHasTaxonomy() ) {
                            ++with_taxonomy;
                            final String sci = n.getNodeData().getTaxonomy().getScientificName();
                            expected = italic.stringWidth( sci ) + roman.stringWidth( " " ) + roman.stringWidth( data );
                        }
                        else {
                            ++without;
                            expected = roman.stringWidth( " " + data );
                        }
                        final double drawn = 2.0 * b[ 2 ];
                        if ( Math.abs( drawn - expected ) > 1e-6 ) {
                            fail( ok, type + " '" + ( n.getNodeData().isHasTaxonomy()
                                    ? n.getNodeData().getTaxonomy().getScientificName() + " " : "" ) + data
                                    + "': the label is " + drawn + " px wide, " + expected + " px with "
                                    + ( n.getNodeData().isHasTaxonomy() ? "one space after the name"
                                            : "its leading space" ) );
                        }
                    }
                    if ( ( with_taxonomy != TIP_NAMES.length ) || ( without != 2 ) ) {
                        fail( ok, type + ": expected " + TIP_NAMES.length + " named tips and 2 controls, read "
                                + with_taxonomy + " and " + without );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, type + ": unexpected " + t );
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
     * Whether {@code rgb} is CORE ink of {@code colour} (magenta or cyan) on white: at least half covered. Text is
     * antialiased on macOS even with antialiasing off, and a thin 1x italic stroke covers no pixel fully, so an
     * exact-colour test finds nothing there; the half-coverage edge is where the glyph's edge is at every scale.
     */
    private static boolean isInk( final int rgb, final int colour ) {
        final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
        if ( b < 200 ) {
            return false;
        }
        return ( colour == 0xFF00FF ) ? ( ( r >= 200 ) && ( g <= 127 ) ) : ( ( g >= 200 ) && ( r <= 127 ) );
    }

    /** The clear gap, in USER px, between the rightmost {@code tax} ink and the leftmost {@code seq} ink in a band
     *  of +/-{@code band} user px around {@code yc} on an image rendered at {@code scale}; NaN when either is
     *  absent. The leftmost seq ink is global, so an overlap comes out negative. */
    private static double gapUser( final BufferedImage img, final double scale, final float yc, final int band,
                                   final int tax, final int seq ) {
        final int y0 = Math.max( 0, (int) Math.floor( ( yc - band ) * scale ) );
        final int y1 = Math.min( img.getHeight(), (int) Math.ceil( ( yc + band ) * scale ) );
        int tax_max = -1, seq_min = Integer.MAX_VALUE;
        for( int y = y0; y < y1; ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int p = img.getRGB( x, y ) & 0xFFFFFF;
                if ( isInk( p, tax ) ) {
                    tax_max = Math.max( tax_max, x );
                }
                else if ( isInk( p, seq ) ) {
                    seq_min = Math.min( seq_min, x );
                }
            }
        }
        return ( ( tax_max < 0 ) || ( seq_min == Integer.MAX_VALUE ) ) ? Double.NaN
                : ( ( seq_min - ( tax_max + 1 ) ) / scale );
    }

    /** The right edge of {@code colour}'s ink in the band from {@code band} px above {@code yc} to 2 px below it,
     *  in USER px (the column after the last inked one); NaN when absent. */
    private static double rightmostInkUser( final BufferedImage img, final double scale, final float yc,
                                            final int band, final int colour ) {
        final int y0 = Math.max( 0, (int) Math.floor( ( yc - band ) * scale ) );
        final int y1 = Math.min( img.getHeight(), (int) Math.ceil( ( yc + 2 ) * scale ) );
        int max = -1;
        for( int y = y0; y < y1; ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( isInk( img.getRGB( x, y ) & 0xFFFFFF, colour ) && ( x > max ) ) {
                    max = x;
                }
            }
        }
        return ( max < 0 ) ? Double.NaN : ( ( max + 1 ) / scale );
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
                    o.setAntialiasExport( false ); // crisp single-color glyphs, so a pixel is exactly one segment's color
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
                    o.setAntialiasExport( false );
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
        return withData( name, sci, "SEQ" );
    }

    private static PhylogenyNode withData( final String name, final String sci, final String seq_name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        n.getNodeData().setTaxonomy( t );
        final Sequence s = new Sequence();
        s.setName( seq_name );
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

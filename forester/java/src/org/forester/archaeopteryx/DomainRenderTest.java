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
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.HashSet;
import java.util.Set;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the domain-architecture demo in the new FLAT style and exercises the modernised palette, the glow toggle, and
 * the three label modes -- in particular the E-value-cutoff-aware, draggable domain legend. Headful; a green no-op when
 * headless. Dogfoods forester/demo/domain-architectures.xml.
 */
public final class DomainRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DomainRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return domainRenderOk();
    }

    private static boolean domainRenderOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/domain-architectures.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "dom" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( true );
                    tp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    if ( !tp.getControlPanel().isShowDomainArchitectures() ) {
                        fail( ok, "the demo tree should auto-enable Domain Architectures" );
                    }
                    final int w = 820, h = 460;

                    // 1. the modern palette gives each domain a DISTINCT colour (was a name-hash). The demo has 6
                    //    distinct domains: SH3, SH2, Pkinase, PH, SAM, Ig.
                    final String[] doms = { "SH3", "SH2", "Pkinase", "PH", "SAM", "Ig" };
                    final Set<Integer> colours = new HashSet<>();
                    for ( final String d : doms ) {
                        colours.add( RenderableDomainArchitecture.colorFor( d ).getRGB() );
                    }
                    if ( colours.size() < doms.length ) {
                        fail( ok, "the domain palette must give each of the 6 domains a distinct colour, got "
                                + colours.size() );
                    }

                    // 2. the flat domain boxes render (many colourful pixels)
                    o.setDomainLabelMode( Options.DOMAIN_LABEL_MODE.ON_DOMAINS );
                    o.setShowDomainGlow( false );
                    tp.setDomainEvalueThresholdExpForTest( -3 );
                    final BufferedImage img_off = render( tp, o, w, h );
                    if ( countColorful( img_off ) < 400 ) {
                        fail( ok, "the flat domain boxes must paint many colourful pixels, got " + countColorful( img_off ) );
                    }

                    // 3. glow ON visibly changes the render (adds the hue halos) -- a pixel diff over glow OFF
                    o.setShowDomainGlow( true );
                    final BufferedImage img_glow = render( tp, o, w, h );
                    final int diff = diffPixels( img_off, img_glow );
                    if ( diff < 100 ) {
                        fail( ok, "domain glow must visibly change the render (halos), pixel diff " + diff );
                    }
                    o.setShowDomainGlow( false );

                    // 4. the LEGEND label mode draws a draggable, hit-testable domain legend; ON_DOMAINS draws none
                    final Rectangle bounds = new Rectangle( 0, 0, w, h );
                    tp.setSize( w, h );
                    o.setDomainLabelMode( Options.DOMAIN_LABEL_MODE.ON_DOMAINS );
                    drawLegend( tp, w, h, bounds );
                    if ( tp.getDomainLegendBoundsForTest() != null ) {
                        fail( ok, "no domain legend must draw in ON_DOMAINS mode" );
                    }
                    o.setDomainLabelMode( Options.DOMAIN_LABEL_MODE.LEGEND );
                    tp.setDomainEvalueThresholdExpForTest( -3 );
                    drawLegend( tp, w, h, bounds );
                    final Rectangle lb = tp.getDomainLegendBoundsForTest();
                    if ( lb == null ) {
                        fail( ok, "the domain legend must draw in LEGEND mode" );
                    }
                    else {
                        if ( !tp.isOnDomainLegend( mouseAt( tp, lb.x + ( lb.width / 2 ), lb.y + ( lb.height / 2 ) ) ) ) {
                            fail( ok, "isOnDomainLegend must be true inside the legend box" );
                        }
                        if ( tp.isOnDomainLegend( mouseAt( tp, lb.x - 30, lb.y + ( lb.height / 2 ) ) ) ) {
                            fail( ok, "isOnDomainLegend must be false outside the legend box" );
                        }
                        if ( !tp.isOnAnyLegend( mouseAt( tp, lb.x + ( lb.width / 2 ), lb.y + ( lb.height / 2 ) ) ) ) {
                            fail( ok, "isOnAnyLegend must be true over the domain legend" );
                        }
                    }

                    // 4b. turning OFF the domain display removes the legend (no orphan legend with no boxes)
                    tp.getControlPanel().setShowDomainArchitecturesForTest( false );
                    drawLegend( tp, w, h, bounds );
                    if ( tp.getDomainLegendBoundsForTest() != null ) {
                        fail( ok, "no domain legend must draw when Show Domain Architectures is off" );
                    }
                    tp.getControlPanel().setShowDomainArchitecturesForTest( true );

                    // 4c. in a RADIAL layout the domain legend draws ONLY when labels are RADIAL -- then domain boxes
                    //     ride the spoke, so the legend is informative. Under horizontal labels the boxes (and the
                    //     legend) are suppressed, since a spoke-riding bar would clash with the upright labels.
                    final Options.PHYLOGENY_GRAPHICS_TYPE prev_type = tp.getPhylogenyGraphicsType();
                    final Options.NODE_LABEL_DIRECTION prev_dir = o.getNodeLabelDirection();
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    o.setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL );
                    drawLegend( tp, w, h, bounds );
                    if ( tp.getDomainLegendBoundsForTest() == null ) {
                        fail( ok, "the domain legend must draw in a circular layout with RADIAL labels" );
                    }
                    o.setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.HORIZONTAL );
                    drawLegend( tp, w, h, bounds );
                    if ( tp.getDomainLegendBoundsForTest() != null ) {
                        fail( ok, "no domain legend in a circular layout with HORIZONTAL labels (domains suppressed)" );
                    }
                    o.setNodeLabelDirection( prev_dir );
                    tp.setPhylogenyGraphicsType( prev_type );

                    // 4d. enabling domains in a radial layout auto-enables "Radial Labels" (option 2): a spoke-riding
                    //     domain bar needs radial labels, so a rectangular tree does NOT flip, but a circular one does.
                    o.setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.HORIZONTAL );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    frame.enableRadialLabelsIfDomainsInRadialLayout();
                    if ( o.getNodeLabelDirection() != Options.NODE_LABEL_DIRECTION.HORIZONTAL ) {
                        fail( ok, "a rectangular layout must NOT auto-enable radial labels" );
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    frame.enableRadialLabelsIfDomainsInRadialLayout();
                    if ( o.getNodeLabelDirection() != Options.NODE_LABEL_DIRECTION.RADIAL ) {
                        fail( ok, "a circular layout with domains on must auto-enable radial labels" );
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    o.setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.HORIZONTAL );

                    // 5. the legend is E-value-cutoff AWARE: below every domain's e-value (1e-6) NOTHING passes and the
                    //    legend is empty (bounds null); a permissive threshold draws it.
                    tp.setDomainEvalueThresholdExpForTest( -12 ); // 1e-12 < 1e-6 -> nothing passes
                    drawLegend( tp, w, h, bounds );
                    if ( tp.getDomainLegendBoundsForTest() != null ) {
                        fail( ok, "the domain legend must be empty when no domain passes the E-value cutoff" );
                    }
                    tp.setDomainEvalueThresholdExpForTest( 0 ); // 1 -> everything passes
                    drawLegend( tp, w, h, bounds );
                    if ( tp.getDomainLegendBoundsForTest() == null ) {
                        fail( ok, "the domain legend must draw when the E-value cutoff admits domains" );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "threw: " + t );
                }
            } );
            SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            return fail( "threw: " + t );
        }
    }

    private static BufferedImage render( final TreePanel tp, final Options o, final int w, final int h ) {
        tp.setSize( w, h );
        tp.calcParametersForPainting( w, h );
        return AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
    }

    private static void drawLegend( final TreePanel tp, final int w, final int h, final Rectangle bounds ) {
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.drawDomainLegendForTest( g, bounds, true );
        g.dispose();
    }

    private static MouseEvent mouseAt( final TreePanel tp, final int x, final int y ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), 0, x, y, 1, false );
    }

    private static int diffPixels( final BufferedImage a, final BufferedImage b ) {
        int n = 0;
        final int w = Math.min( a.getWidth(), b.getWidth() ), h = Math.min( a.getHeight(), b.getHeight() );
        for ( int y = 0; y < h; ++y ) {
            for ( int x = 0; x < w; ++x ) {
                final int p = a.getRGB( x, y ), q = b.getRGB( x, y );
                final int dr = Math.abs( ( ( p >> 16 ) & 0xFF ) - ( ( q >> 16 ) & 0xFF ) );
                final int dg = Math.abs( ( ( p >> 8 ) & 0xFF ) - ( ( q >> 8 ) & 0xFF ) );
                final int db = Math.abs( ( p & 0xFF ) - ( q & 0xFF ) );
                if ( ( dr + dg + db ) > 24 ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static int countColorful( final BufferedImage img ) {
        int n = 0;
        for ( int y = 0; y < img.getHeight(); ++y ) {
            for ( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) ), min = Math.min( r, Math.min( g, b ) );
                if ( ( ( max - min ) > 40 ) && ( max > 90 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [DomainRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [DomainRenderTest] " + msg );
        ok[ 0 ] = false;
    }
}

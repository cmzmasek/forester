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
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Tests "Size by property" ({@link PropertySizeScale}). Headless: the value-&gt;diameter mapping is
 * area-proportional and monotonic, a value-less node keeps the base size, and {@link PropertyColorScheme#numericRefs}
 * offers only numeric columns. Headful: a rendered tree really draws a bigger tip dot for a bigger value.
 */
public final class PropertySizeScaleTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PropertySizeScale: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = mathOk() && numericRefsOk() && legendMathOk() && demoTreeOk();
        if ( !GraphicsEnvironment.isHeadless() ) {
            ok &= sizeRenderOk();
            ok &= sizeLegendRenderOk();
            ok &= legendInteractionOk();
        }
        return ok;
    }

    /** The size legend's INTERACTION wiring: hit-test (inside / outside / stale-after-off), isOnAnyLegend, and -- via
     *  the REAL MouseListener -- that a click where the top-most size legend OVERLAPS the color legend resets the SIZE
     *  legend (not the color one). Covers isOnSizeLegend/isOnAnyLegend/handleSizeLegendClick + the startLegendDrag/
     *  dragLegend size branch + the mouseClicked size-first dispatch. */
    private static boolean legendInteractionOk() {
        try {
            final Phylogeny phy = treeWith( "data:abundance", "1", "50", "100" );
            final List<PhylogenyNode> leaves = phy.getExternalNodes();
            addProp( leaves.get( 0 ), "data:host", "Human" ); // a categorical prop so Color-by has a legend too
            addProp( leaves.get( 1 ), "data:host", "Avian" );
            addProp( leaves.get( 2 ), "data:host", "Swine" );
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "sizeint" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final int w = 700, h = 400;
                    tp.setSize( w, h ); // so getVisibleRect() (used by dragLegend) is (0,0,w,h)
                    final Rectangle bounds = new Rectangle( 0, 0, w, h );
                    final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );

                    // hit-test on the size legend alone
                    tp.setSizeByPropertyRef( "data:abundance" );
                    drawBoth( tp, img, bounds ); // Color-by off here -> only the size legend draws
                    final Rectangle sb = tp.getSizeLegendBounds();
                    if ( sb == null ) {
                        fail( ok, "size legend bounds must be recorded when drawn draggable" );
                    }
                    else {
                        if ( !tp.isOnSizeLegend( mouseAt( tp, sb.x + ( sb.width / 2 ), sb.y + ( sb.height / 2 ) ) ) ) {
                            fail( ok, "isOnSizeLegend must be true inside the box" );
                        }
                        if ( tp.isOnSizeLegend( mouseAt( tp, sb.x - 40, sb.y + ( sb.height / 2 ) ) ) ) {
                            fail( ok, "isOnSizeLegend must be false outside the box" );
                        }
                        if ( !tp.isOnAnyLegend( mouseAt( tp, sb.x + ( sb.width / 2 ), sb.y + ( sb.height / 2 ) ) ) ) {
                            fail( ok, "isOnAnyLegend must be true over the size legend" );
                        }
                        // stale-bounds guard: turning Size-by off makes the hit-test false at the SAME point
                        final MouseEvent center = mouseAt( tp, sb.x + ( sb.width / 2 ), sb.y + ( sb.height / 2 ) );
                        tp.setSizeByPropertyRef( null );
                        if ( tp.isOnSizeLegend( center ) ) {
                            fail( ok, "after Size-by off, isOnSizeLegend must be false (guarded by isSizeByProperty)" );
                        }
                        if ( tp.getSizeLegendBounds() != null ) {
                            fail( ok, "turning Size-by off must clear the stale size-legend bounds" );
                        }
                    }

                    // overlap dispatch: a click on the on-top size legend resets IT, not the color legend underneath
                    tp.setColorByPropertyRef( "data:host" );
                    tp.setSizeByPropertyRef( "data:abundance" );
                    drawBoth( tp, img, bounds );
                    final Rectangle cb = tp.getPropertyLegendBounds();
                    final Rectangle sb0 = tp.getSizeLegendBounds();
                    if ( ( cb == null ) || ( sb0 == null ) ) {
                        fail( ok, "both legends must draw for the overlap test" );
                    }
                    else {
                        // drag the size legend onto the color legend so the boxes overlap
                        tp.startLegendDrag( mouseAt( tp, sb0.x + ( sb0.width / 2 ), sb0.y + ( sb0.height / 2 ) ) );
                        tp.dragLegend( mouseAt( tp, cb.x + ( cb.width / 2 ), cb.y + ( cb.height / 2 ) ) );
                        drawBoth( tp, img, bounds );
                        final Rectangle sbOver = tp.getSizeLegendBounds();
                        final Rectangle cbNow = tp.getPropertyLegendBounds();
                        if ( !sbOver.intersects( cbNow ) ) {
                            fail( ok, "precondition: the dragged size legend should overlap the color legend" );
                        }
                        else {
                            final int ox = Math.max( sbOver.x, cbNow.x ) + 2;
                            final int oy = Math.max( sbOver.y, cbNow.y ) + 2;
                            new MouseListener( tp ).mouseClicked( mouseAtClicks( tp, ox, oy, 2 ) ); // real dispatch
                            drawBoth( tp, img, bounds );
                            if ( tp.getSizeLegendBounds().intersects( tp.getPropertyLegendBounds() ) ) {
                                fail( ok, "double-clicking the on-top size legend must RESET it (clear the overlap)" );
                            }
                            if ( !tp.getPropertyLegendBounds().equals( cb ) ) {
                                fail( ok, "the color legend must NOT move when the size legend is clicked in overlap" );
                            }
                        }
                    }
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

    private static void drawBoth( final TreePanel tp, final BufferedImage img, final Rectangle bounds ) {
        final Graphics2D g = img.createGraphics();
        tp.drawLegendForTest( g, bounds, true );      // color/rank legend -> _property_legend_bounds (if active)
        tp.drawSizeLegendForTest( g, bounds, true );  // size legend -> _size_legend_bounds
        g.dispose();
    }

    private static MouseEvent mouseAt( final TreePanel tp, final int x, final int y ) {
        return mouseAtClicks( tp, x, y, 1 );
    }

    private static MouseEvent mouseAtClicks( final TreePanel tp, final int x, final int y, final int clicks ) {
        return new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), 0, x, y, clicks, false );
    }

    /** The committed demo tree (forester/demo/size-by-property.xml) must really demonstrate the feature: it parses,
     *  offers the numeric ref, and the largest-value tip scales up from the smallest -- so the demo cannot silently
     *  rot (see DemoTreeGenerator / the demo gallery). */
    private static boolean demoTreeOk() {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/size-by-property.xml" );
        if ( !file.exists() ) {
            return fail( "demo tree missing: " + file.getAbsolutePath() );
        }
        try {
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance().create( file, parser );
            if ( ( parser.getErrorCount() > 0 ) || ( phys.length != 1 ) || phys[ 0 ].isEmpty() ) {
                return fail( "demo tree did not parse into one non-empty tree" );
            }
            final Phylogeny phy = phys[ 0 ];
            if ( !PropertyColorScheme.numericRefs( phy ).contains( "data:read_count" ) ) {
                return fail( "demo tree must offer numeric 'data:read_count' to size by" );
            }
            final PropertySizeScale scale = new PropertySizeScale( phy, "data:read_count" );
            if ( scale.isEmpty() ) {
                return fail( "demo tree size scale must not be empty" );
            }
            final float base = 12f;
            float d_min = -1, d_max = -1;
            for( final PhylogenyNode leaf : phy.getExternalNodes() ) {
                if ( "A/HongKong/1997".equals( leaf.getName() ) ) { // read_count 120 = the smallest
                    d_min = scale.diameterFor( leaf, base );
                }
                else if ( "A/Anhui/2013".equals( leaf.getName() ) ) { // read_count 88000 = the largest
                    d_max = scale.diameterFor( leaf, base );
                }
            }
            if ( ( d_min < 0 ) || ( d_max <= d_min ) ) {
                return fail( "demo tree's max-value tip must scale larger than the min-value tip (min=" + d_min
                        + " max=" + d_max + ")" );
            }
        }
        catch ( final Exception e ) {
            return fail( "demo tree could not be read: " + e.getMessage() );
        }
        return true;
    }

    private static boolean mathOk() {
        boolean ok = true;
        final float base = 10f;
        // pure diameter(): area-proportional, clamped
        ok &= approx( "diameter t=0", PropertySizeScale.diameter( 0.0, base ), base ); // MIN_FACTOR=1 -> base
        ok &= approx( "diameter t=1", PropertySizeScale.diameter( 1.0, base ), 3 * base ); // MAX_FACTOR=3
        ok &= approx( "diameter clamp<0", PropertySizeScale.diameter( -5, base ), base );
        ok &= approx( "diameter clamp>1", PropertySizeScale.diameter( 9, base ), 3 * base );
        // AREA is linear in t: at t=0.5 the AREA is the midpoint of the endpoints' areas (NOT the diameter -- a
        // linear-diameter map would give 20; area-proportional gives sqrt(500) ~= 22.36)
        final double mid = PropertySizeScale.diameter( 0.5, base );
        ok &= approx( "area-proportional midpoint", (float) ( mid * mid ),
                      (float) ( ( ( base * base ) + ( 9 * base * base ) ) / 2.0 ) );
        if ( Math.abs( mid - 20.0 ) < 0.5 ) {
            ok = fail( "t=0.5 diameter " + mid + " looks LINEAR (should be area-proportional ~22.36)" );
        }

        // a scale over real leaf values: min->base, max->3x, monotonic; a value-less leaf -> base
        final Phylogeny phy = treeWith( "data:v", "0", "5", "10", null );
        final List<PhylogenyNode> leaves = phy.getExternalNodes();
        final PropertySizeScale scale = new PropertySizeScale( phy, "data:v" );
        if ( scale.isEmpty() ) {
            ok = fail( "scale over numeric values must not be empty" );
        }
        final float d0 = scale.diameterFor( leaves.get( 0 ), base ); // value 0 (min)
        final float d1 = scale.diameterFor( leaves.get( 1 ), base ); // value 5 (mid)
        final float d2 = scale.diameterFor( leaves.get( 2 ), base ); // value 10 (max)
        final float dn = scale.diameterFor( leaves.get( 3 ), base ); // no value
        ok &= approx( "min value -> base", d0, base );
        ok &= approx( "max value -> 3x base", d2, 3 * base );
        if ( !( ( d0 < d1 ) && ( d1 < d2 ) ) ) {
            ok = fail( "diameter must increase with value: " + d0 + " < " + d1 + " < " + d2 );
        }
        if ( dn != base ) {
            ok = fail( "a node with no value must keep the base size, got " + dn );
        }
        // hasValueFor tells the painter which tips actually carry a value (only those get a size dot)
        if ( !scale.hasValueFor( leaves.get( 0 ) ) ) {
            ok = fail( "a tip with a numeric value must report hasValueFor == true" );
        }
        if ( scale.hasValueFor( leaves.get( 3 ) ) ) {
            ok = fail( "a value-less tip must report hasValueFor == false" );
        }
        // min/max value + value-driven diameter (used by the size legend to draw sample dots)
        if ( ( scale.getMinValue() != 0 ) || ( scale.getMaxValue() != 10 ) ) {
            ok = fail( "min/max value must reflect the data (got " + scale.getMinValue() + ".." + scale.getMaxValue() + ")" );
        }
        ok &= approx( "diameterForValue min", scale.diameterForValue( 0, base ), base );
        ok &= approx( "diameterForValue max", scale.diameterForValue( 10, base ), 3 * base );
        final float dmv = scale.diameterForValue( 5, base );
        if ( !( ( scale.diameterForValue( 0, base ) < dmv ) && ( dmv < scale.diameterForValue( 10, base ) ) ) ) {
            ok = fail( "diameterForValue must increase with value" );
        }

        // all-equal values: no spread -> every node at the base size (and numericRefs won't offer it, see below)
        final Phylogeny flat = treeWith( "data:v", "5", "5", "5" );
        final PropertySizeScale flat_scale = new PropertySizeScale( flat, "data:v" );
        if ( flat_scale.diameterFor( flat.getExternalNodes().get( 0 ), base ) != base ) {
            ok = fail( "an all-equal property must render every node at the base size" );
        }
        // no numeric values at all -> empty
        if ( !new PropertySizeScale( treeWith( "data:v", "cat", "dog" ), "data:v" ).isEmpty() ) {
            ok = fail( "a non-numeric property must yield an empty size scale" );
        }
        return ok;
    }

    /** The size-legend's pure helpers: sample values (min/mid/max, or one when flat) and value formatting. */
    private static boolean legendMathOk() {
        boolean ok = true;
        final double[] s = PropertySizeScale.sampleValues( 10, 90 );
        if ( ( s.length != 3 ) || ( s[ 0 ] != 10 ) || ( s[ 1 ] != 50 ) || ( s[ 2 ] != 90 ) ) {
            ok = fail( "samples must be min/mid/max, got " + java.util.Arrays.toString( s ) );
        }
        final double[] flat = PropertySizeScale.sampleValues( 7, 7 );
        if ( ( flat.length != 1 ) || ( flat[ 0 ] != 7 ) ) {
            ok = fail( "a no-spread range must yield a single sample value" );
        }
        if ( PropertySizeScale.sampleValues( 5, 3 ).length != 1 ) { // defensive: min > max (empty scale)
            ok = fail( "min > max must yield a single sample value" );
        }
        if ( !"120".equals( PropertySizeScale.formatValue( 120.0 ) ) ) {
            ok = fail( "a whole number must format as an integer, got " + PropertySizeScale.formatValue( 120.0 ) );
        }
        if ( !"5.3".equals( PropertySizeScale.formatValue( 5.3 ) ) ) {
            ok = fail( "a fraction must format with decimals (US-locale), got " + PropertySizeScale.formatValue( 5.3 ) );
        }
        // a large fractional (e.g. a year midpoint) keeps its integer part -- not collapsed by significant figures
        if ( !"2016.5".equals( PropertySizeScale.formatValue( 2016.5 ) ) ) {
            ok = fail( "a large fraction must keep its value, got " + PropertySizeScale.formatValue( 2016.5 ) );
        }
        // SMALL magnitudes (a 0..1 property: p-values / pident / posterior) must NOT collapse to "0"/duplicate labels
        if ( !"0.005".equals( PropertySizeScale.formatValue( 0.005 ) ) ) {
            ok = fail( "a small value must keep ~3 sig digits, got " + PropertySizeScale.formatValue( 0.005 ) );
        }
        // the min/mid/max of a small range must be three DISTINCT, non-zero labels (the old 0.## rounded them to
        // "0"/"0.01"/"0.01" -- the key lied about which dot was which)
        final double[] tiny = PropertySizeScale.sampleValues( 0.001, 0.009 );
        final String lo = PropertySizeScale.formatValue( tiny[ 0 ] );
        final String mid = PropertySizeScale.formatValue( tiny[ 1 ] );
        final String hi = PropertySizeScale.formatValue( tiny[ 2 ] );
        if ( lo.equals( mid ) || mid.equals( hi ) || "0".equals( lo ) ) {
            ok = fail( "a small range must yield three distinct non-zero labels, got " + lo + "/" + mid + "/" + hi );
        }
        return ok;
    }

    private static boolean numericRefsOk() {
        boolean ok = true;
        if ( !PropertyColorScheme.numericRefs( treeWith( "data:year", "2015", "2016", "2020", "2024" ) )
                .contains( "data:year" ) ) {
            ok = fail( "a numeric property must appear in numericRefs" );
        }
        if ( PropertyColorScheme.numericRefs( treeWith( "repseq:host", "cat", "cat", "dog", "dog" ) )
                .contains( "repseq:host" ) ) {
            ok = fail( "a categorical property must NOT appear in numericRefs" );
        }
        return ok;
    }

    private static boolean sizeRenderOk() {
        try {
            final Phylogeny phy = treeWith( "data:abundance", "1", "50", "100", null );
            final List<PhylogenyNode> leaves = phy.getExternalNodes();
            // node names are painted in the SEQUENCE color too, which we force to cyan below to measure the dot;
            // blank them so the ONLY cyan ink is the property dot itself (else dotWidth would measure label ink,
            // and a no-op sizing regression could pass on the glyph-width difference of the labels)
            for( final PhylogenyNode leaf : leaves ) {
                leaf.setName( "" );
            }
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "size" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setAntialiasPrint( false );
                    o.setGraphicsExportWhiteBackground( false );
                    final int dot = 0x00FFFF; // the size-only dot is drawn in the SEQUENCE color -- force it to cyan
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.SEQUENCE, new Color( dot ) );
                    tp.getTreeColorSet().setColorSchema( 0 );
                    tp.setSizeByPropertyRef( "data:abundance" );
                    if ( !tp.isSizeByProperty() ) {
                        fail( ok, "Size by should be active after setSizeByPropertyRef" );
                    }
                    frame.showWhole();
                    final int w = 700, h = 400;
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );
                    final BufferedImage img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int small = dotWidth( img, leaves.get( 0 ).getXcoord(), leaves.get( 0 ).getYcoord(), dot ); // abundance 1 (min)
                    final int big = dotWidth( img, leaves.get( 2 ).getXcoord(), leaves.get( 2 ).getYcoord(), dot );   // abundance 100 (max)
                    final int none = dotWidth( img, leaves.get( 3 ).getXcoord(), leaves.get( 3 ).getYcoord(), dot );  // no value -> no size dot
                    if ( ( small < 0 ) || ( big < 0 ) ) {
                        fail( ok, "expected a tip dot at both the min and max node (small=" + small + " big=" + big + ")" );
                    }
                    else if ( big <= small ) {
                        fail( ok, "the max-value tip dot (" + big + "px) must be larger than the min-value one ("
                                + small + "px)" );
                    }
                    if ( none >= 0 ) {
                        fail( ok, "a value-less tip must draw NO size dot in size-only mode, got width " + none );
                    }
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

    /** Diameter (px) of the {@code color} dot at the tip ({@code xc},{@code yc}): the longest vertical run of that
     *  color in the node's OWN x-column. The tip label is drawn to the RIGHT of the node and the size legend sits at
     *  a figure corner, so neither crosses this column -- the measurement is the dot alone, wherever the legend is.
     *  -1 if none. */
    private static int dotWidth( final BufferedImage img, final float xc, final float yc, final int color ) {
        final int x = Math.round( xc );
        if ( ( x < 0 ) || ( x >= img.getWidth() ) ) {
            return -1;
        }
        final int y0 = Math.max( 0, (int) yc - 20 );
        final int y1 = Math.min( img.getHeight(), (int) yc + 20 );
        int run = 0, best = -1;
        for( int y = y0; y < y1; ++y ) {
            if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == color ) {
                if ( ++run > best ) {
                    best = run;
                }
            }
            else {
                run = 0;
            }
        }
        return best;
    }

    /** The size legend renders decodably AND coexists with the color legend: with BOTH Color-by and Size-by on, the
     *  two legend boxes are both drawn and do NOT overlap, and the size legend's max-value sample dot is clearly
     *  larger than the base dot (so a reader can decode size -> value). */
    private static boolean sizeLegendRenderOk() {
        try {
            final Phylogeny phy = treeWith( "data:abundance", "1", "50", "100", "80" );
            final List<PhylogenyNode> leaves = phy.getExternalNodes();
            final String[] hosts = { "Human", "Avian", "Swine", "Avian" }; // a categorical prop to Color-by
            for( int i = 0; i < leaves.size(); i++ ) {
                addProp( leaves.get( i ), "data:host", hosts[ i ] );
            }
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "sizelegend" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final int cyan = 0x00FFFF; // the size dots + legend ink are drawn in the SEQUENCE color
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.SEQUENCE, new Color( cyan ) );
                    tp.getTreeColorSet().setColorSchema( 0 );
                    tp.setColorByPropertyRef( "data:host" );     // color legend (categorical) -> top-right
                    tp.setSizeByPropertyRef( "data:abundance" ); // size legend -> bottom-right (a color legend is present)
                    if ( !tp.isColorByProperty() || !tp.isSizeByProperty() ) {
                        fail( ok, "both Color by and Size by should be active" );
                    }
                    final float base = tp.baseDotSize(); // single source shared with the tip renderer + the legend
                    final int w = 700, h = 400;
                    final Rectangle bounds = new Rectangle( 0, 0, w, h );
                    final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                    final Graphics2D g = img.createGraphics();
                    tp.drawLegendForTest( g, bounds, true );     // color legend (sets _property_legend_bounds)
                    tp.drawSizeLegendForTest( g, bounds, true );  // size legend (sets _size_legend_bounds)
                    g.dispose();
                    final Rectangle color_bounds = tp.getPropertyLegendBounds();
                    final Rectangle size_bounds = tp.getSizeLegendBounds();
                    if ( ( color_bounds == null ) || ( size_bounds == null ) ) {
                        fail( ok, "both legends must draw (color=" + color_bounds + " size=" + size_bounds + ")" );
                    }
                    else if ( color_bounds.intersects( size_bounds ) ) {
                        fail( ok, "the size legend must not overlap the color legend by default" );
                    }
                    if ( size_bounds != null ) {
                        // measure INSIDE the 1px border (else the box outline's full-width top/bottom edges dominate);
                        // the widest neutral run must be the max-value sample dot, clearly bigger than the base dot
                        final int widest = widestRun( img, cyan, size_bounds.y + 2, ( size_bounds.y + size_bounds.height ) - 2,
                                size_bounds.x + 2, ( size_bounds.x + size_bounds.width ) - 2 );
                        if ( widest < Math.round( 2.5f * base ) ) {
                            fail( ok, "the size legend's max sample dot must scale up (~" + Math.round( 3 * base )
                                    + "px), widest neutral run was only " + widest + "px" );
                        }
                    }
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

    /** Widest contiguous horizontal run (px) of {@code color} within [x0,x1) x [y0,y1). */
    private static int widestRun( final BufferedImage img, final int color, int y0, int y1, int x0, int x1 ) {
        y0 = Math.max( 0, y0 );
        y1 = Math.min( img.getHeight(), y1 );
        x0 = Math.max( 0, x0 );
        x1 = Math.min( img.getWidth(), x1 );
        int best = 0;
        for( int y = y0; y < y1; ++y ) {
            int run = 0;
            for( int x = x0; x < x1; ++x ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == color ) {
                    if ( ++run > best ) {
                        best = run;
                    }
                }
                else {
                    run = 0;
                }
            }
        }
        return best;
    }

    private static void addProp( final PhylogenyNode n, final String ref, final String value ) {
        PropertiesList pl = n.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            n.getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( ref, value, "", "xsd:string", AppliesTo.NODE ) );
    }

    private static Phylogeny treeWith( final String ref, final String... values ) {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        int i = 0;
        for( final String v : values ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "leaf" + ( i++ ) );
            if ( v != null ) {
                final PropertiesList pl = new PropertiesList();
                pl.addProperty( new Property( ref, v, "", "xsd:string", AppliesTo.NODE ) );
                n.getNodeData().setProperties( pl );
            }
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean approx( final String name, final float actual, final float expected ) {
        if ( Math.abs( actual - expected ) > 0.5f ) {
            System.out.println( "  [PropertySizeScaleTest] " + name + ": expected ~" + expected + ", got " + actual );
            return false;
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [PropertySizeScaleTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [PropertySizeScaleTest] " + msg );
        ok[ 0 ] = false;
    }

    private PropertySizeScaleTest() {
    }
}

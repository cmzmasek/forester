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
        boolean ok = mathOk() && numericRefsOk() && demoTreeOk();
        if ( !GraphicsEnvironment.isHeadless() ) {
            ok &= sizeRenderOk();
        }
        return ok;
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
                    final int small = dotWidth( img, leaves.get( 0 ).getYcoord(), dot ); // abundance 1 (min)
                    final int big = dotWidth( img, leaves.get( 2 ).getYcoord(), dot );   // abundance 100 (max)
                    final int none = dotWidth( img, leaves.get( 3 ).getYcoord(), dot );  // no value -> no size dot
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

    /** Width (px) of the {@code color} ink in a +/-20px band around {@code yc}, or -1 if none. */
    private static int dotWidth( final BufferedImage img, final float yc, final int color ) {
        final int y0 = Math.max( 0, (int) yc - 20 );
        final int y1 = Math.min( img.getHeight(), (int) yc + 20 );
        int min_x = Integer.MAX_VALUE, max_x = -1;
        for( int y = y0; y < y1; ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) == color ) {
                    min_x = Math.min( min_x, x );
                    max_x = Math.max( max_x, x );
                }
            }
        }
        return ( max_x < 0 ) ? -1 : ( ( max_x - min_x ) + 1 );
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

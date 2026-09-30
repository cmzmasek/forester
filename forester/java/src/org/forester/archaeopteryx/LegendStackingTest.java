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

import java.awt.Point;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.Options.TREE_ORIENTATION;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Legends placed at their DEFAULT spots never open on top of each other ({@link TreePanel}'s
 * claimDefaultLegendSpot and sharedLegendHoldsTopRight). The bugs it guards against: in circular and unrooted the
 * ancestral-pie legend went to the top right ON the Color-by legend (a rule written when Color-by was suppressed
 * radially); Size-by + pies without Color-by shared the top right in every layout; the Size-by legend, the
 * internal-taxonomy key and the domain legend all took the same bottom-right spot. Also pinned: the first legend at
 * a corner keeps it exactly; a dragged legend neither claims nor yields; a pass forgets the last one's spots.
 */
public final class LegendStackingTest {

    private static final int INSET = 10; // TreePanel.LEGEND_EDGE_INSET: the gap kept between stacked legends

    private static final Object[][] LAYOUTS = {
            { "rectangular root-left", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_LEFT },
            { "rectangular root-top", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_TOP },
            { "rectangular root-bottom", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_BOTTOM },
            { "circular", PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, TREE_ORIENTATION.ROOT_LEFT },
            { "unrooted", PHYLOGENY_GRAPHICS_TYPE.UNROOTED, TREE_ORIENTATION.ROOT_LEFT } };

    public static void main( final String[] args ) {
        System.out.println( test() ? "LegendStackingTest: OK." : "LegendStackingTest: FAILED." );
        System.exit( 0 );
    }

    public static boolean test() {
        try {
            return pieTreeOk() & domainTreeOk();
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected " + e );
        }
    }

    /** Color-by / Size-by / pies on the pie demo, in every layout. */
    private static boolean pieTreeOk() throws Exception {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/ancestral-pie-charts.xml" );
        if ( !file.exists() ) {
            return fail( "demo tree missing: " + file.getAbsolutePath() );
        }
        final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        addTipProperties( phy );
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "stacking" ) );
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                for ( final Object[] layout : LAYOUTS ) {
                    final String name = (String) layout[ 0 ];
                    // (1) Color-by + pies: the pies go to the bottom LEFT in EVERY layout, clear of the Color-by key
                    set( tp, "beast:location", null, "location" );
                    paint( mf[ 0 ], tp, layout );
                    final Rectangle color = tp.getPropertyLegendBounds();
                    final Rectangle pies = tp.getAncestralPieLegendBounds();
                    if ( ( color == null ) || ( pies == null ) ) {
                        ok[ 0 ] = fail( name + ": Color-by + pies must draw both legends, got " + color + " / " + pies );
                        continue;
                    }
                    final Rectangle view = tp.getVisibleRect();
                    if ( color.intersects( pies ) ) {
                        ok[ 0 ] = fail( name + ": the pie legend " + pies + " opens on the Color-by legend " + color );
                    }
                    if ( ( pies.x != ( view.x + INSET ) ) || ( ( pies.y + pies.height ) != ( view.y + view.height - INSET ) ) ) {
                        ok[ 0 ] = fail( name + ": beside a Color-by legend the pie legend belongs at the bottom left, got "
                                + pies + " in " + view );
                    }
                    // (2) Size-by + pies, no Color-by: both want the top right; the first keeps it, the second stacks
                    set( tp, null, "data:sz", "location" );
                    paint( mf[ 0 ], tp, layout );
                    final Rectangle size = tp.getSizeLegendBounds();
                    final Rectangle pies2 = tp.getAncestralPieLegendBounds();
                    if ( ( size == null ) || ( pies2 == null ) || ( tp.getPropertyLegendBounds() != null ) ) {
                        ok[ 0 ] = fail( name + ": Size-by + pies must draw exactly those two, got " + size + " / " + pies2 );
                        continue;
                    }
                    if ( ( size.y != ( view.y + INSET ) ) || ( ( size.x + size.width ) != ( view.x + view.width - INSET ) ) ) {
                        ok[ 0 ] = fail( name + ": the first legend at the top right must keep the corner, got " + size );
                    }
                    if ( pies2.y < ( size.y + size.height + INSET ) ) {
                        ok[ 0 ] = fail( name + ": the pie legend " + pies2 + " must stack under the Size-by legend " + size );
                    }
                    if ( ( pies2.x + pies2.width ) != ( size.x + size.width ) ) {
                        ok[ 0 ] = fail( name + ": a stacked legend stays flush with its corner's edge, got " + pies2 );
                    }
                    // (3) a pass forgets the last pass's spots: painting again moves nothing
                    paint( mf[ 0 ], tp, layout );
                    if ( !size.equals( tp.getSizeLegendBounds() ) || !pies2.equals( tp.getAncestralPieLegendBounds() ) ) {
                        ok[ 0 ] = fail( name + ": a second paint moved the legends: " + tp.getSizeLegendBounds() + " / "
                                + tp.getAncestralPieLegendBounds() );
                    }
                    // (4) a DRAGGED legend neither claims nor yields: the Size-by legend dragged ONTO the top-right
                    // corner (flush with the view's top) stays exactly there, and the pie legend still takes its
                    // plain corner -- the user placed the one, and the other does not dodge it
                    tp.setSizeLegendOffsetForTest( new Point( view.width - size.width, 0 ) );
                    paint( mf[ 0 ], tp, layout );
                    final Rectangle dragged = tp.getSizeLegendBounds();
                    final Rectangle pies3 = tp.getAncestralPieLegendBounds();
                    if ( ( dragged == null ) || ( dragged.y != view.y ) ) {
                        ok[ 0 ] = fail( name + ": a dragged legend stays where it was dragged, got " + dragged );
                    }
                    else if ( !dragged.intersects( pies3 ) ) {
                        ok[ 0 ] = fail( name + ": precondition: the dragged legend must sit on the pies' corner, got "
                                + dragged + " / " + pies3 );
                    }
                    if ( ( pies3 == null ) || ( pies3.y != ( view.y + INSET ) ) ) {
                        ok[ 0 ] = fail( name + ": beside a DRAGGED Size-by legend the pies take the plain corner, got "
                                + pies3 );
                    }
                    tp.setSizeLegendOffsetForTest( null );
                }
                // (6) a SHORT window: the Color-by key (top right) and the Size-by key (bottom right, beside it) meet.
                // Stacking up would push the Size-by key off the top, so it goes sideways -- left of the Color-by key,
                // inside the view; and in a view too small for that too, it keeps its plain corner (overlap beats
                // out of sight)
                set( tp, "beast:location", "data:sz", null );
                final java.awt.Dimension frame_size = ( (JFrame) mf[ 0 ] ).getSize();
                for ( final Object[] layout : LAYOUTS ) {
                    final String name = (String) layout[ 0 ];
                    paint( mf[ 0 ], tp, layout );
                    final int color_h = tp.getPropertyLegendBounds().height;
                    final int size_h = tp.getSizeLegendBounds().height;
                    // a view in which the two keys, each at its corner, are closer than the gap
                    shrinkViewTo( mf[ 0 ], tp, ( color_h + size_h + INSET + INSET ) - 4, frame_size );
                    paint( mf[ 0 ], tp, layout );
                    final Rectangle view = tp.getVisibleRect();
                    final Rectangle color = tp.getPropertyLegendBounds();
                    final Rectangle size = tp.getSizeLegendBounds();
                    if ( ( view.height >= ( color_h + size_h + INSET + INSET + INSET ) )
                            || ( view.height < ( Math.max( color_h, size_h ) + INSET + INSET ) ) ) {
                        ok[ 0 ] = fail( name + ": precondition: a view where the two corners meet but one key fits, got "
                                + view );
                        ( (JFrame) mf[ 0 ] ).setSize( frame_size );
                        ( (JFrame) mf[ 0 ] ).validate();
                        continue;
                    }
                    if ( !view.contains( size ) || !view.contains( color ) ) {
                        ok[ 0 ] = fail( name + ": in a short window both keys stay in view, got " + color + " / " + size
                                + " in " + view );
                    }
                    final Rectangle apart = new Rectangle( color );
                    apart.grow( INSET - 1, INSET - 1 );
                    if ( apart.intersects( size ) ) {
                        ok[ 0 ] = fail( name + ": in a short window the keys are still kept apart, got " + color + " / "
                                + size );
                    }
                    if ( ( size.x + size.width ) > color.x ) {
                        ok[ 0 ] = fail( name + ": the Size-by key goes LEFT of the Color-by key, got " + size + " / " + color );
                    }
                    // too small for either way: the plain bottom-right corner
                    shrinkViewTo( mf[ 0 ], tp, size_h + INSET + INSET + 2, frame_size );
                    ( (JFrame) mf[ 0 ] ).setSize( 260, ( (JFrame) mf[ 0 ] ).getHeight() );
                    ( (JFrame) mf[ 0 ] ).validate();
                    paint( mf[ 0 ], tp, layout );
                    final Rectangle tiny = tp.getVisibleRect();
                    final Rectangle size2 = tp.getSizeLegendBounds();
                    if ( tiny.width >= ( tp.getPropertyLegendBounds().width + size2.width + ( 3 * INSET ) ) ) {
                        ok[ 0 ] = fail( name + ": precondition: a view too narrow for the keys side by side, got " + tiny );
                    }
                    // the plain corner, clamped into the view as the corner itself is
                    else if ( ( size2.x != Math.max( tiny.x, ( tiny.x + tiny.width ) - size2.width - INSET ) )
                            || ( size2.y != Math.max( tiny.y, ( tiny.y + tiny.height ) - size2.height - INSET ) ) ) {
                        ok[ 0 ] = fail( name + ": with no room either way a key keeps its plain corner, got " + size2
                                + " in " + tiny );
                    }
                    ( (JFrame) mf[ 0 ] ).setSize( frame_size );
                    ( (JFrame) mf[ 0 ] ).validate();
                }
                set( tp, null, null, null );
                // (5) outside a paint pass nothing is claimed: a legend drawn on its own lands on its plain corner
                final BufferedImage scratch = new BufferedImage( 400, 300, BufferedImage.TYPE_INT_ARGB );
                tp.setAncestralPieTrait( "location" );
                for ( int i = 0; i < 2; ++i ) {
                    tp.drawAncestralPieLegendForTest( scratch.createGraphics(), new Rectangle( 0, 0, 400, 300 ), true,
                                                      false );
                    final Rectangle b = tp.getAncestralPieLegendBounds();
                    if ( ( b == null ) || ( b.y != INSET ) ) {
                        ok[ 0 ] = fail( "a legend drawn outside a paint pass must take its plain corner, draw " + i
                                + " got " + b );
                    }
                }
                tp.setAncestralPieTrait( null );
            } );
        }
        finally {
            dispose( mf );
        }
        return ok[ 0 ];
    }

    /** Three legends that all want the bottom right: Size-by (beside a Color-by key), the internal-taxonomy key and
     *  the domain legend. Domain boxes draw radially only with radial labels, so rectangular root-left and circular
     *  (radial labels) are the layouts where all three are drawn. */
    private static boolean domainTreeOk() throws Exception {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/domain-architectures.xml" );
        if ( !file.exists() ) {
            return fail( "demo tree missing: " + file.getAbsolutePath() );
        }
        final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        addTipProperties( phy );
        int k = 0;
        for ( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() ) {
                final Taxonomy t = new Taxonomy();
                t.setScientificName( "Order" + ( k++ ) );
                t.setRank( "order" );
                n.getNodeData().setTaxonomy( t );
            }
        }
        final MainFrame[] mf = new MainFrame[ 1 ];
        final boolean[] ok = { true };
        final Options.NODE_LABEL_DIRECTION[] saved_dir = { null };
        try {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, new Configuration(), "stacking-br" ) );
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final Options o = mf[ 0 ].getOptions();
                saved_dir[ 0 ] = o.getNodeLabelDirection();
                o.setShowInternalTaxonomyKey( true );
                o.setDomainLabelMode( Options.DOMAIN_LABEL_MODE.LEGEND );
                mf[ 0 ].setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL );
                set( tp, "data:grp", "data:sz", null );
                for ( final Object[] layout : new Object[][] { LAYOUTS[ 0 ], LAYOUTS[ 3 ] } ) {
                    final String name = (String) layout[ 0 ];
                    paint( mf[ 0 ], tp, layout );
                    final Rectangle view = tp.getVisibleRect();
                    final List<Rectangle> b = Arrays.asList( tp.getSizeLegendBounds(),
                                                             tp.getInternalTaxaKeyBoundsForTest(),
                                                             tp.getDomainLegendBoundsForTest() );
                    final String[] names = { "Size-by", "internal taxa", "domain" };
                    if ( b.contains( null ) ) {
                        ok[ 0 ] = fail( name + ": all three bottom-right legends must be drawn, got " + b );
                        continue;
                    }
                    final Rectangle size = b.get( 0 );
                    if ( ( ( size.y + size.height ) != ( view.y + view.height - INSET ) )
                            || ( ( size.x + size.width ) != ( view.x + view.width - INSET ) ) ) {
                        ok[ 0 ] = fail( name + ": the first legend at the bottom right must keep the corner, got " + size );
                    }
                    for ( int i = 0; i < 3; ++i ) {
                        for ( int j = i + 1; j < 3; ++j ) {
                            final Rectangle apart = new Rectangle( b.get( i ) );
                            apart.grow( INSET - 1, INSET - 1 );
                            if ( apart.intersects( b.get( j ) ) ) {
                                ok[ 0 ] = fail( name + ": the " + names[ i ] + " legend " + b.get( i ) + " and the "
                                        + names[ j ] + " legend " + b.get( j ) + " are not kept apart" );
                            }
                        }
                    }
                    // stacked UP from the bottom corner, in paint order
                    if ( !( ( b.get( 1 ).y < size.y ) && ( b.get( 2 ).y < b.get( 1 ).y ) ) ) {
                        ok[ 0 ] = fail( name + ": bottom-right legends stack upward in paint order, got " + b );
                    }
                }
                o.setShowInternalTaxonomyKey( false );
                o.setDomainLabelMode( Options.DOMAIN_LABEL_MODE.NONE );
                mf[ 0 ].setNodeLabelDirection( saved_dir[ 0 ] );
                set( tp, null, null, null );
            } );
        }
        finally {
            dispose( mf );
        }
        return ok[ 0 ];
    }

    private static void addTipProperties( final Phylogeny phy ) {
        int i = 1;
        for ( final PhylogenyNode leaf : phy.getExternalNodes() ) {
            PropertiesList pl = leaf.getNodeData().getProperties();
            if ( pl == null ) {
                pl = new PropertiesList();
                leaf.getNodeData().setProperties( pl );
            }
            pl.addProperty( new Property( "data:sz", Integer.toString( i ), "", "xsd:decimal",
                                          Property.AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:grp", ( ( i++ % 2 ) == 0 ) ? "even" : "odd", "", "xsd:string",
                                          Property.AppliesTo.NODE ) );
        }
    }

    private static void set( final TreePanel tp, final String color, final String size, final String pies ) {
        tp.setColorByPropertyRef( color );
        tp.setSizeByPropertyRef( size );
        tp.setAncestralPieTrait( pies );
    }

    private static void paint( final MainFrame frame, final TreePanel tp, final Object[] layout ) {
        tp.setPhylogenyGraphicsType( (PHYLOGENY_GRAPHICS_TYPE) layout[ 1 ] );
        tp.setTreeOrientation( (TREE_ORIENTATION) layout[ 2 ] );
        frame.showWhole();
        final BufferedImage img = new BufferedImage( Math.max( 1, tp.getWidth() ), Math.max( 1, tp.getHeight() ),
                                                     BufferedImage.TYPE_INT_RGB );
        tp.printAll( img.getGraphics() );
    }

    /** Shrinks the frame until the tree panel's view is {@code h} tall (the panel sits in a scroll pane, so only
     *  the frame's size moves the viewport). */
    private static void shrinkViewTo( final MainFrame frame, final TreePanel tp, final int h,
                                      final java.awt.Dimension frame_size ) {
        final JFrame f = (JFrame) frame;
        f.setSize( frame_size );
        f.validate();
        final int chrome = f.getHeight() - tp.getVisibleRect().height;
        f.setSize( frame_size.width, chrome + h );
        f.validate();
    }

    private static void dispose( final MainFrame[] mf ) throws Exception {
        SwingUtilities.invokeAndWait( () -> {
            if ( mf[ 0 ] != null ) {
                ( (JFrame) mf[ 0 ] ).dispose();
            }
        } );
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [LegendStackingTest] " + message );
        return false;
    }

    private LegendStackingTest() {
    }
}

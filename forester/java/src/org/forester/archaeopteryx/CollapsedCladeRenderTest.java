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
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * A collapsed clade on screen, as Archaeopteryx.js draws it (archaeopteryx.js 407f204): it takes 1 + log2(tips)/4 rows
 * in the layout (the gap to a neighbour is the mean of the weights), its wedge runs from the node to the clade's
 * nearest and farthest tips (one step in a cladogram), it is filled in the dominant Color-by colour, outlined and
 * counted while a search hits inside and filled in the hit colour when all its tips are hits, stays bright while the
 * rest dims, and is named by its node name, else a 95% Color-by value, else the tips' common name prefix. Rectangular
 * (root left and top), circular and unrooted. Headful.
 * <p>
 * The tree: ((a1:1,a2:1)A:1,(((c1:0.5,c2:2)C1:0.5,(c3:1,c4:1.5)C2:1)C:1,(d1:1,d2:1)D:1)CD:1,e:1). C holds 4 tips,
 * 1.0 / 2.5 / 2.0 / 2.5 from C, so it takes 1.5 rows and weighs 2; with it collapsed the rows are a1 a2 C d1 d2 e = 7.
 */
public final class CollapsedCladeRenderTest {

    private static final int W = 900;
    private static final int H = 700;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CollapsedCladeRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString(
                    "((a1:1,a2:1)A:1,(((c1:0.5,c2:2)C1:0.5,(c3:1,c4:1.5)C2:1)C:1,(d1:1,d2:1)D:1)CD:1,e:1)root" );
            for( final String t : new String[] { "a1", "a2", "d2", "e" } ) {
                lineage( phy.getNode( t ), "Delta" );
            }
            lineage( phy.getNode( "d1" ), "Omicron" ); // Omicron is also on screen, so it has a colour there
            for( final String t : new String[] { "c1", "c2", "c3", "c4" } ) {
                lineage( phy.getNode( t ), "Omicron" );
            }
            SwingUtilities.invokeAndWait( () -> {
                mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(), "collapsed" );
                mf[ 0 ].setSize( W, H );
            } );
            final String[] failure = { null };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    exercise( mf[ 0 ], phy );
                }
                catch ( final AssertionError e ) {
                    failure[ 0 ] = e.getMessage();
                }
            } );
            if ( failure[ 0 ] != null ) {
                System.out.println( "  [CollapsedCladeRenderTest] " + failure[ 0 ] );
                return false;
            }
            return true;
        }
        catch ( final Throwable e ) {
            e.printStackTrace( System.out );
            return false;
        }
        finally {
            if ( mf[ 0 ] != null ) {
                try {
                    SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
                }
                catch ( final Exception e ) {
                    // nothing left to release
                }
            }
        }
    }

    private static void exercise( final MainFrame mf, final Phylogeny phy ) {
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = tp.getControlPanel();
        tp.setSize( W, H );
        tp.getOptions().setShowOverview( false );
        cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
        final PhylogenyNode a1 = phy.getNode( "a1" ), a2 = phy.getNode( "a2" ), d1 = phy.getNode( "d1" );
        final PhylogenyNode c = phy.getNode( "C" ), A = phy.getNode( "A" );
        tp.collapse( c );
        ck( c.isCollapse(), "precondition: C collapses" );
        layout( tp, cp );

        // ---- rows: C weighs 2 (1.5 rows); the gap to each neighbour is the mean of the weights --------------------
        ck( near( tp.rowWeight( c ), 2.0 ), "a collapsed 4-tip clade weighs 2 rows, got " + tp.rowWeight( c ) );
        ck( near( tp.rowWeight( phy.getRoot() ), 7.0 ), "the tree takes 7 rows (6 displayed leaves, C counting 2), got "
                + tp.rowWeight( phy.getRoot() ) );
        final double tip_gap = a2.getYcoord() - a1.getYcoord();
        ck( tip_gap > 1, "precondition: the tips are laid out apart, gap " + tip_gap );
        ck( near( c.getYcoord() - a2.getYcoord(), 1.5 * tip_gap, 0.02 * tip_gap ),
            "a tip and a collapsed 4-tip clade sit (1 + 2) / 2 = 1.5 tip gaps apart: " + ( c.getYcoord() - a2.getYcoord() )
                    + " vs " + ( 1.5 * tip_gap ) );
        ck( near( d1.getYcoord() - c.getYcoord(), 1.5 * tip_gap, 0.02 * tip_gap ), "...on both sides" );

        // ---- the wedge: apex at C, upper edge to the nearest tip (1.0), lower to the farthest (2.5) --------------
        double[] w = tp.collapsedWedgeRectangular( c );
        final double unit = a1.getXcoord() - A.getXcoord(); // a branch of length 1
        ck( near( w[ 0 ], c.getXcoord() ) && near( w[ 1 ], c.getYcoord() ), "the apex is the clade's node" );
        ck( near( w[ 2 ] - w[ 0 ], 1.0 * unit, 0.5 ), "the upper edge reaches the nearest tip (1.0 away): " + ( w[ 2 ] - w[ 0 ] )
                + " vs " + unit );
        ck( near( w[ 3 ] - w[ 0 ], 2.5 * unit, 0.5 ), "the lower edge reaches the farthest tip (2.5 away): " + ( w[ 3 ] - w[ 0 ] )
                + " vs " + ( 2.5 * unit ) );
        ck( near( w[ 4 ], Math.max( 6, 1.5 * tip_gap * 0.82 ) / 2, 0.01 ), "the wedge is 1.5 rows x 0.82 tall: half "
                + w[ 4 ] + " vs " + ( Math.max( 6, 1.5 * tip_gap * 0.82 ) / 2 ) );

        // ---- the label: node name, else a 95% Color-by value, else the tips' common name prefix --------------------
        ck( "C · 4 tips".equals( tp.collapsedLookForTest( c ).label ), "named by its node: "
                + tp.collapsedLookForTest( c ).label );
        c.setName( "" );
        tp.setColorByPropertyRef( null ); // the tree may have opened coloured by lineage
        ck( "4 tips".equals( tp.collapsedLookForTest( c ).label ), "no name, no Color-by, no long prefix: just the count, got "
                + tp.collapsedLookForTest( c ).label );
        tp.setColorByPropertyRef( "data:lineage" );
        ck( tp.isColorByProperty(), "precondition: Color-by host is on" );
        ck( "Omicron · 4 tips".equals( tp.collapsedLookForTest( c ).label ), "all 4 tips Omicron names the clade, got "
                + tp.collapsedLookForTest( c ).label );
        lineage( phy.getNode( "c4" ), "Delta" );
        ck( "4 tips".equals( tp.collapsedLookForTest( c ).label ), "3 of 4 (75%) is under 95%: no name, got "
                + tp.collapsedLookForTest( c ).label );
        final String[] tips = { "c1", "c2", "c3", "c4" };
        for( int i = 0; i < tips.length; ++i ) {
            phy.getNode( tips[ i ] ).setName( "SARS_CoV_2/human/x" + i );
        }
        ck( "SARS_CoV_2/human · 4 tips".equals( tp.collapsedLookForTest( c ).label ),
            "then the tips' common name prefix, its trailing separator dropped, got " + tp.collapsedLookForTest( c ).label );

        // ---- colours: the dominant Color-by colour, 22% fill, 90% outline ------------------------------------------
        Color omicron = tp.getPropertyBasedColor( phy.getNode( "d1" ) );
        TreePanel.CollapsedLook look = tp.collapsedLookForTest( c );
        ck( sameRgb( look.fill, omicron ), "3 Omicron + 1 Delta fills in the Omicron colour, got " + look.fill + " vs " + omicron );
        ck( look.fill.getAlpha() == 56 && look.stroke.getAlpha() == 230, "fill 22% / outline 90% opaque, got "
                + look.fill.getAlpha() + " / " + look.stroke.getAlpha() );
        ck( look.stroke_width == 1f && !look.full, "no hit: a 1 px outline" );
        ck( sameRgb( look.ink, tp.getTreeColorSet().getSequenceColor() ), "the label in the tip label colour" );
        // the paint follows: the wedge's middle takes the fill over the background
        shot( tp );
        w = tp.collapsedWedgeRectangular( c );
        final BufferedImage img = shot( tp );
        final int cx = (int) Math.round( ( w[ 0 ] + w[ 2 ] + w[ 3 ] ) / 3 ), cy = (int) Math.round( w[ 1 ] );
        final Color bg = tp.getTreeColorSet().getBackgroundColor();
        final Color px = new Color( img.getRGB( cx, cy ) );
        final Color expected = blend( omicron, bg, 0.22 );
        ck( distance( px, expected ) < distance( px, bg ) && distance( px, expected ) < 40,
            "the wedge's middle is painted in the fill (" + px + ", expected ~" + expected + ", background " + bg + ")" );
        tp.setColorByPropertyRef( null );
        ck( sameRgb( tp.collapsedLookForTest( c ).fill, tp.getTreeColorSet().getBranchColor() ),
            "without Color-by the wedge takes the branch colour" );
        // a value that only the hidden tips carry has no legend row on screen, but it keeps its REMEMBERED colour, so
        // the wedge keeps the colour its tips wore (Christian, 2026-09-13, option (b))
        lineage( phy.getNode( "d1" ), "Delta" );
        lineage( phy.getNode( "SARS_CoV_2/human/x3" ), "Omicron" );
        tp.setColorByPropertyRef( "data:lineage" );
        ck( !tp.getPropertyColorScheme().getValueColors().containsKey( "Omicron" ),
            "precondition: Omicron is only hidden, so the legend has no row for it" );
        ck( sameRgb( tp.collapsedLookForTest( c ).fill, omicron ),
            "a value carried only by hidden tips keeps its remembered colour, got " + tp.collapsedLookForTest( c ).fill
                    + " vs " + omicron );
        // ...also a value that was NEVER on screen: after Reset forgets the memory, re-choosing the field remembers
        // every value the whole tree carries, so the hidden-only value still gets a colour of its own
        tp.resetColorStateToDefaults();
        tp.setColorByPropertyRef( "data:lineage" );
        final Color delta = tp.getPropertyColorScheme().getValueColors().get( "Delta" );
        final Color never = tp.collapsedLookForTest( c ).fill;
        ck( ( delta != null ) && !sameRgb( never, delta ) && !sameRgb( never, tp.getTreeColorSet().getBranchColor() ),
            "a never-shown hidden value gets its own remembered colour (not Delta's, not the branch colour), got "
                    + never + " (Delta " + delta + ")" );
        ck( !tp.getPropertyColorScheme().getValueColors().containsKey( "Omicron" ), "...still without a legend row" );
        lineage( phy.getNode( "d1" ), "Omicron" );
        tp.resetColorStateToDefaults();
        tp.setColorByPropertyRef( "data:lineage" );
        omicron = tp.getPropertyBasedColor( phy.getNode( "d1" ) ); // the fresh assignment follows the tips on screen

        // ---- a search hit inside: outline + count; all tips hits: filled, bold, in the hit colour -----------------
        final Color found = tp.getTreeColorSet().getFoundColor0();
        tp.setFoundNodes0( ids( phy, "SARS_CoV_2/human/x0" ) );
        look = tp.collapsedLookForTest( c );
        ck( look.label.endsWith( " [1/4]" ), "one hit inside counts [1/4], got " + look.label );
        ck( sameRgb( look.stroke, found ) && ( look.stroke_width == 1.5f ), "outlined 1.5 px in the hit colour" );
        ck( sameRgb( look.fill, omicron ) && !look.full, "a partial hit keeps the Color-by fill" );
        tp.setFoundNodes0( ids( phy, "SARS_CoV_2/human/x0", "SARS_CoV_2/human/x1", "SARS_CoV_2/human/x2",
                                "SARS_CoV_2/human/x3" ) );
        look = tp.collapsedLookForTest( c );
        ck( look.full && sameRgb( look.fill, found ) && ( look.fill.getAlpha() == 115 ),
            "every tip a hit: filled 45% in the hit colour, got " + look.fill + " alpha " + look.fill.getAlpha() );
        ck( sameRgb( look.ink, found ) && look.label.endsWith( " [4/4]" ), "...and labelled in it" );

        // ---- Dim Non-Matches: the clade dims without a hit, stays bright holding one, and a hidden hit engages it --
        tp.getOptions().setDimNonMatches( true );
        tp.setFoundNodes0( ids( phy, "e" ) );
        shot( tp );
        ck( tp.hasVisibleFoundNodeForTest(), "precondition: a drawn hit engages dimming" );
        ck( !sameRgb( tp.collapsedLookForTest( c ).ink, tp.getTreeColorSet().getSequenceColor() ),
            "a collapsed clade holding no hit dims with the rest" );
        tp.setFoundNodes0( ids( phy, "SARS_CoV_2/human/x2" ) );
        shot( tp );
        ck( tp.hasVisibleFoundNodeForTest(), "a hit hidden inside a collapsed clade counts as on screen (the rest dims)" );
        ck( sameRgb( tp.collapsedLookForTest( c ).ink, tp.getTreeColorSet().getSequenceColor() ),
            "the clade holding the hit stays bright" );
        tp.getOptions().setDimNonMatches( false );
        tp.setFoundNodes0( null );

        // ---- root on top: the same wedge, drawn through the rotation ------------------------------------------------
        tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_TOP );
        layout( tp, cp );
        w = tp.collapsedWedgeRectangular( c );
        final BufferedImage top = shot( tp );
        final Point2D.Double dev = tp.screenPoint( ( w[ 0 ] + w[ 2 ] + w[ 3 ] ) / 3, w[ 1 ] );
        ck( distance( new Color( top.getRGB( (int) Math.round( dev.x ), (int) Math.round( dev.y ) ) ), bg ) > 8,
            "root on top: the wedge is painted where the rotation puts it" );
        tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );

        // ---- cladogram: the wedge is one step deep and ends on the tip column ---------------------------------------
        cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
        layout( tp, cp );
        w = tp.collapsedWedgeRectangular( c );
        ck( near( w[ 2 ], w[ 3 ] ), "a cladogram wedge's two edges end together" );
        if ( !tp.getOptions().getCladogramType().equals( Options.CLADOGRAM_TYPE.NON_LINED_UP ) ) {
            ck( near( w[ 3 ], a1.getXcoord(), 0.5 ), "the lined-up cladogram wedge ends on the tip column: " + w[ 3 ] + " vs "
                    + a1.getXcoord() );
            ck( c.getXcoord() < a1.getXcoord() - 1, "...its apex one step inside it" );
        }
        cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );

        // ---- circular: rows as angles, the wedge along the spoke ----------------------------------------------------
        radial( tp, cp, Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        final double row = ( 2 * Math.PI ) / 7;
        final double ga = tp.circularAngleForTest( a2 ) - tp.circularAngleForTest( a1 );
        final double gc = tp.circularAngleForTest( c ) - tp.circularAngleForTest( a2 );
        final double gd = tp.circularAngleForTest( d1 ) - tp.circularAngleForTest( c );
        ck( near( ga, row, 1e-6 ), "circular: two tips sit one row apart (2pi / 7), got " + ga );
        ck( near( gc, 1.5 * row, 1e-6 ) && near( gd, 1.5 * row, 1e-6 ),
            "circular: a tip and the collapsed clade sit 1.5 rows apart, got " + gc + " / " + gd );
        double[] r = tp.collapsedWedgeRadialForTest( c );
        ck( ( r[ 0 ] > 1 ) && near( r[ 1 ] / r[ 0 ], 2.5, 0.02 ), "circular phylogram: the wedge reaches the nearest (1.0) and "
                + "farthest (2.5) tips along the spoke, got " + r[ 0 ] + " / " + r[ 1 ] );
        cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
        radial( tp, cp, Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        r = tp.collapsedWedgeRadialForTest( c );
        ck( near( r[ 0 ], r[ 1 ] ) && ( r[ 0 ] > 1 ), "circular cladogram: one step, both edges to the ring" );
        cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );

        // ---- unrooted: the same wedge along the spoke ---------------------------------------------------------------
        radial( tp, cp, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
        r = tp.collapsedWedgeRadialForTest( c );
        ck( ( r[ 0 ] > 1 ) && near( r[ 1 ] / r[ 0 ], 2.5, 0.02 ), "unrooted: nearest / farthest along the spoke, got "
                + r[ 0 ] + " / " + r[ 1 ] );

        // ---- opening it restores plain rows --------------------------------------------------------------------------
        radial( tp, cp, Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        tp.collapse( c );
        ck( !c.isCollapse(), "precondition: C opens again" );
        layout( tp, cp );
        ck( near( tp.rowWeight( phy.getRoot() ), 9.0 ), "opened: 9 tips, 9 rows, got " + tp.rowWeight( phy.getRoot() ) );
    }

    private static void layout( final TreePanel tp, final ControlPanel cp ) {
        cp.displayedPhylogenyMightHaveChanged( true );
        tp.calcParametersForPainting( W, H );
        shot( tp );
    }

    private static void radial( final TreePanel tp, final ControlPanel cp, final Options.PHYLOGENY_GRAPHICS_TYPE type ) {
        tp.getOptions().setPhylogenyGraphicsType( type );
        tp.setPhylogenyGraphicsType( type );
        layout( tp, cp );
    }

    private static BufferedImage shot( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.printAll( g );
        g.dispose();
        return img;
    }

    private static void lineage( final PhylogenyNode n, final String value ) {
        final PropertiesList props = new PropertiesList();
        props.addProperty( new Property( "data:lineage", value, "", "xsd:string", AppliesTo.NODE ) );
        n.getNodeData().setProperties( props );
    }

    private static Set<Long> ids( final Phylogeny phy, final String... names ) {
        final Set<Long> out = new HashSet<>();
        for( final String n : names ) {
            out.add( phy.getNode( n ).getId() );
        }
        return out;
    }

    private static boolean sameRgb( final Color a, final Color b ) {
        return ( a != null ) && ( b != null ) && ( a.getRed() == b.getRed() ) && ( a.getGreen() == b.getGreen() )
                && ( a.getBlue() == b.getBlue() );
    }

    private static Color blend( final Color fg, final Color bg, final double alpha ) {
        return new Color( (int) Math.round( ( fg.getRed() * alpha ) + ( bg.getRed() * ( 1 - alpha ) ) ),
                          (int) Math.round( ( fg.getGreen() * alpha ) + ( bg.getGreen() * ( 1 - alpha ) ) ),
                          (int) Math.round( ( fg.getBlue() * alpha ) + ( bg.getBlue() * ( 1 - alpha ) ) ) );
    }

    private static double distance( final Color a, final Color b ) {
        return Math.abs( a.getRed() - b.getRed() ) + Math.abs( a.getGreen() - b.getGreen() )
                + Math.abs( a.getBlue() - b.getBlue() );
    }

    private static boolean near( final double a, final double b ) {
        return near( a, b, 1e-3 );
    }

    private static boolean near( final double a, final double b, final double tol ) {
        return Math.abs( a - b ) <= tol;
    }

    private static void ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            throw new AssertionError( msg );
        }
    }

    private CollapsedCladeRenderTest() {
    }
}

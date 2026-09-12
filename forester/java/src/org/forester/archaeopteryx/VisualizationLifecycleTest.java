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
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * The visualization-selection LIFECYCLE on a real tab -- JS-authoritative (forester.js 79dd9de, Christian 2026-09-12).
 * A VIEW (entering or leaving a subtree) never changes what "Color by" offers or what is chosen; it only re-summarizes
 * the chosen field over the tips on screen. An EDIT (a deletion, an undo) re-derives the candidates but KEEPS a chosen
 * field while it still carries a value anywhere. And -- desktop only, since Archaeopteryx.js cannot edit node data -- a
 * node-data write is an edit too, so the colouring follows the new values instead of going stale.
 */
public final class VisualizationLifecycleTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "VisualizationLifecycle: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // needs a display
        }
        try {
            return viewNeverReclassifies() && editKeepsChosenField() && nodeDataWriteRecolours();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** A one-value clade stays coloured, the menu is the tree's, and a field absent from the view keeps its legend. */
    private static boolean viewNeverReclassifies() throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree() }, new Configuration(), "vislife1" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            final ControlPanel cp = tp.getControlPanel();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setColorByPropertyRef( "x:Genus" );
            cp.populateColorByPropertyBox();
            final List<String> menu = new ArrayList<String>( cp.colorByPropertyRefs() );
            check( ok, "Genus and Sparse offered on the whole tree: " + menu,
                   menu.contains( "x:Genus" ) && menu.contains( "x:Sparse" ) );
            // into clade A, where every tip is one genus: CLASSIFYING the view would refuse Genus outright
            tp.subTree( tp.getPhylogeny().getRoot().getChildNode( 0 ) );
            check( ok, "now in a subtree", tp.isCurrentTreeIsSubtree() );
            check( ok, "the chosen field stands", "x:Genus".equals( tp.getColorByPropertyRef() ) );
            check( ok, "the one-value clade is COLOURED, with one legend row", tp.isColorByProperty()
                    && ( tp.getPropertyColorScheme().getValueColors().size() == 1 ) );
            cp.populateColorByPropertyBox(); // what a tab switch does while a subtree is displayed
            check( ok, "the menu is the TREE's, not the view's: " + cp.colorByPropertyRefs(),
                   menu.equals( cp.colorByPropertyRefs() ) );
            // a chosen field with NO value in the view stays chosen; nothing is coloured, the legend is "no value"
            tp.setColorByPropertyRef( "x:Sparse" );
            check( ok, "a field absent from the view stays chosen", "x:Sparse".equals( tp.getColorByPropertyRef() ) );
            check( ok, "nothing is coloured", !tp.isColorByProperty() );
            check( ok, "but its legend is there, the no-value row for all six tips", tp.hasColorByPropertyLegend()
                    && ( tp.getPropertyColorScheme().missingCount() == 6 ) );
            paintLegend( tp ); // the empty-scheme legend path must draw without throwing
            tp.superTree();
            check( ok, "back on the tree, the field colours again",
                   tp.isColorByProperty() && "x:Sparse".equals( tp.getColorByPropertyRef() ) );
            check( ok, "and the menu never moved: " + cp.colorByPropertyRefs(), menu.equals( cp.colorByPropertyRefs() ) );
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    /** A deletion that leaves one host keeps the chosen Host, appended and flagged; an undo makes it ordinary again. */
    private static boolean editKeepsChosenField() throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree() }, new Configuration(), "vislife2" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            final ControlPanel cp = tp.getControlPanel();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setColorByPropertyRef( "x:Host" );
            cp.populateColorByPropertyBox();
            // delete clade A: only "duck" is left, which the rules would now refuse
            tp.deleteNodeOrSubtreeConfirmed( tp.getPhylogeny().getRoot().getChildNode( 0 ), false );
            final List<String> menu = cp.colorByPropertyRefs();
            check( ok, "Host is still chosen", "x:Host".equals( tp.getColorByPropertyRef() ) );
            check( ok, "Host is KEPT in the menu, last: " + menu,
                   !menu.isEmpty() && "x:Host".equals( menu.get( menu.size() - 1 ) ) );
            check( ok, "and flagged kept", ( tp.visualizationCandidate( "x:Host" ) != null )
                    && tp.visualizationCandidate( "x:Host" )._kept );
            check( ok, "an unchosen field refused by the edit is gone: " + menu, !menu.contains( "x:Genus" ) );
            check( ok, "the remaining tips are still coloured", tp.isColorByProperty() && tp.getPropertyColorScheme()
                    .getValueColors().keySet().equals( new java.util.HashSet<String>( Arrays.asList( "Duck" ) ) ) );
            // an undo is an edit too: the restored tree offers Host as an ordinary candidate
            tp.undo();
            check( ok, "after the undo Host is still chosen", "x:Host".equals( tp.getColorByPropertyRef() ) );
            check( ok, "and an ordinary candidate again", ( tp.visualizationCandidate( "x:Host" ) != null )
                    && !tp.visualizationCandidate( "x:Host" )._kept );
            check( ok, "Genus is back in the menu: " + cp.colorByPropertyRefs(),
                   cp.colorByPropertyRefs().contains( "x:Genus" ) );
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    /**
     * DESKTOP ONLY: a node-data write is an EDIT. Every write path ends in {@code setEdited(true)} (the node-data
     * window's refresh does), which re-derives the candidates and rebuilds the scheme -- before that, a tip given a NEW
     * value kept no colour and the legend went stale until something else rebuilt it.
     */
    private static boolean nodeDataWriteRecolours() throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { tree() }, new Configuration(), "vislife3" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            final ControlPanel cp = tp.getControlPanel();
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setColorByPropertyRef( "x:Zone" );
            cp.populateColorByPropertyBox();
            check( ok, "two zones before the write", tp.getPropertyColorScheme().getValueColors().size() == 2 );
            final PhylogenyNode a0 = tip( tp, "a0" );
            setValue( a0, "x:Zone", "E" );
            tp.setEdited( true );
            check( ok, "the new zone is COLOURED, not left grey", tp.getPropertyColorScheme().colorFor( a0 ) != null );
            check( ok, "and the legend has three rows: " + tp.getPropertyColorScheme().getValueColors().keySet(),
                   tp.getPropertyColorScheme().getValueColors().size() == 3 );
            // removing the chosen field from every tip drops the choice: there is nothing left to colour by
            for( final PhylogenyNode n : tp.getPhylogeny().getExternalNodes() ) {
                setValue( n, "x:Zone", null );
            }
            tp.setEdited( true );
            check( ok, "a chosen field with no value left falls back to none", tp.getColorByPropertyRef() == null );
            check( ok, "and leaves the menu: " + cp.colorByPropertyRefs(), !cp.colorByPropertyRefs().contains( "x:Zone" ) );
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    private static void paintLegend( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( 800, 600, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        try {
            tp.drawLegendForTest( g, new Rectangle( 0, 0, 800, 600 ), true );
        }
        finally {
            g.dispose();
        }
    }

    /** Two clades of six tips. Genus is one value per clade; Host is cat/dog in A and duck in B; Zone alternates;
     *  Sparse is carried only in clade B. */
    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int c = 0; c < 2; ++c ) {
            final PhylogenyNode clade = new PhylogenyNode();
            for( int i = 0; i < 6; ++i ) {
                final PhylogenyNode n = new PhylogenyNode();
                n.setName( ( ( c == 0 ) ? "a" : "b" ) + i );
                final PropertiesList pl = new PropertiesList();
                pl.addProperty( prop( "x:Genus", ( c == 0 ) ? "Mastadenovirus" : "Aviadenovirus" ) );
                pl.addProperty( prop( "x:Host", ( c == 0 ) ? ( ( ( i % 2 ) == 0 ) ? "cat" : "dog" ) : "duck" ) );
                pl.addProperty( prop( "x:Zone", ( ( i % 2 ) == 0 ) ? "N" : "S" ) );
                if ( c == 1 ) {
                    pl.addProperty( prop( "x:Sparse", ( ( i % 2 ) == 0 ) ? "p" : "q" ) );
                }
                n.getNodeData().setProperties( pl );
                clade.addAsChild( n );
            }
            root.addAsChild( clade );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static Property prop( final String ref, final String value ) {
        return new Property( ref, value, "", "xsd:string", AppliesTo.NODE );
    }

    /** Replaces (or, with a null value, removes) the node's property {@code ref}. */
    private static void setValue( final PhylogenyNode n, final String ref, final String value ) {
        final PropertiesList pl = new PropertiesList();
        for( final Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( !ref.equals( p.getRef() ) ) {
                pl.addProperty( p );
            }
        }
        if ( value != null ) {
            pl.addProperty( prop( ref, value ) );
        }
        n.getNodeData().setProperties( pl );
    }

    private static PhylogenyNode tip( final TreePanel tp, final String name ) {
        for( final PhylogenyNode n : tp.getPhylogeny().getExternalNodes() ) {
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        throw new IllegalStateException( "no tip " + name );
    }

    private static void check( final boolean[] ok, final String what, final boolean cond ) {
        if ( !cond ) {
            System.out.println( "  [VisualizationLifecycleTest] FAILED: " + what );
            ok[ 0 ] = false;
        }
    }

    private VisualizationLifecycleTest() {
        // not instantiable
    }
}

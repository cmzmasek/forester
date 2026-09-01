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

import java.awt.GraphicsEnvironment;
import java.util.Arrays;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * A composed figure survives being saved and reopened.
 * <p>
 * Before this, only three things travelled with a tree, and none of them was the figure: the annotation columns,
 * the clade marks, colour-by and the properties shown in the labels were all lost on save/reload. The round trip
 * below is the whole point -- capture what is drawn, put it on the tree, read it back, and get the same figure.
 */
public final class FigureSpecTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FigureSpec: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return codec() && parsing() && roundTrip() && perTabRestore();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [FigureSpecTest] " + msg );
        return false;
    }

    /** A phyloXML property value collapses whitespace on reload, so anything stored in one must be escaped. */
    private static boolean codec() {
        for( final String s : new String[] { "data:host", "a b", "  lead and trail  ", "semi;colon", "eq=uals",
                                             "pipe|bar", "tilde~x", "back\\slash", "tab\there", "" } ) {
            final String round = PropertyTextCodec.unesc( PropertyTextCodec.esc( s ) );
            if ( !s.equals( round ) ) {
                return fail( "codec round trip failed for [" + s + "] -> [" + round + "]" );
            }
        }
        // the separators must not survive escaping, or a value containing one would split the record
        final String esc = PropertyTextCodec.esc( "a;b=c|d~e f" );
        for( final char c : new char[] { ';', '=', '|', '~', ' ' } ) {
            if ( esc.indexOf( c ) >= 0 ) {
                return fail( "'" + c + "' must not survive escaping: " + esc );
            }
        }
        return true;
    }

    private static boolean parsing() {
        if ( FigureSpec.parse( null ) != null ) {
            return fail( "no property -> no figure" );
        }
        if ( FigureSpec.parse( "" ) != null ) {
            return fail( "an empty property -> no figure" );
        }
        // a version this build does not know must yield NO figure rather than a misread one
        if ( FigureSpec.parse( "v99;layout=CIRCULAR" ) != null ) {
            return fail( "an unknown version must not be parsed as if it were v1" );
        }
        // an unknown KEY is ignored, so a figure written by a newer build still opens
        final FigureSpec f = FigureSpec.parse( "v1;layout=CIRCULAR;someFutureThing=42" );
        if ( ( f == null ) || !"CIRCULAR".equals( f.get( "layout" ) ) ) {
            return fail( "a newer build's extra key must be ignored, not fatal" );
        }
        return true;
    }

    /** The real thing: compose a figure, write it to the tree, read it back into a fresh panel. */
    private static boolean roundTrip() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final boolean[] ok = { true };
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { tree() }, new Configuration(), "figspec" ) );
            final Phylogeny[] saved = new Phylogeny[ 1 ];
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                // compose a figure: a layout, an overlay, a label choice, and a display toggle
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                tp.setAnnotationColumns( Arrays.asList(
                        new AnnotationColumns.ColumnSpec( "data:host", AnnotationColumns.Type.COLOR_STRIP ) ) );
                tp.setLabelPropertyRefs( Arrays.asList( "data:host" ) );
                tp.setColorByPropertyRef( "data:host" );
                tp.setShows( DisplayOption.SHOW_TAX_RANK, true );
                tp.setShows( DisplayOption.SHOW_NODE_NAMES, false );
                tp.syncFigureToTree(); // the production seam: what every phyloXML save calls
                saved[ 0 ] = tp.getPhylogeny().copy(); // as a save/reload would hand it back
            } );
            // the figure must actually be ON the tree, in the aptx: namespace so it stays out of the user's way
            final FigureSpec read = FigureSpec.readFrom( saved[ 0 ] );
            if ( read == null ) {
                ok[ 0 ] = fail( "the figure was not stored on the tree" );
                return false;
            }
            if ( !TreePanelUtil.isInternalPropertyRef( FigureSpec.FIGURE_REF ) ) {
                ok[ 0 ] = fail( "the figure property must be internal (aptx:), never shown as user data" );
            }
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );

            // ...and opening that tree afresh must reproduce the figure
            final MainFrame[] mf2 = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf2[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { saved[ 0 ] }, new Configuration(), "figspec2" ) );
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf2[ 0 ].getMainPanel().getCurrentTreePanel();
                if ( tp.getPhylogenyGraphicsType() != Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                    ok[ 0 ] = fail( "the layout was not restored, got " + tp.getPhylogenyGraphicsType() );
                }
                if ( ( tp.getAnnotationColumnSpecs() == null ) || ( tp.getAnnotationColumnSpecs().size() != 1 )
                        || !"data:host".equals( tp.getAnnotationColumnSpecs().get( 0 )._ref ) ) {
                    ok[ 0 ] = fail( "the annotation column was not restored: " + tp.getAnnotationColumnSpecs() );
                }
                if ( ( tp.getLabelPropertyRefs() == null )
                        || !tp.getLabelPropertyRefs().equals( Arrays.asList( "data:host" ) ) ) {
                    ok[ 0 ] = fail( "the label fields were not restored: " + tp.getLabelPropertyRefs() );
                }
                if ( !"data:host".equals( tp.getColorByPropertyRef() ) ) {
                    ok[ 0 ] = fail( "colour-by was not restored: " + tp.getColorByPropertyRef() );
                }
                if ( !tp.shows( DisplayOption.SHOW_TAX_RANK ) || tp.shows( DisplayOption.SHOW_NODE_NAMES ) ) {
                    ok[ 0 ] = fail( "the display toggles were not restored -- which labels are drawn IS the figure" );
                }
                // and Clear All Overlays takes the overlays off again without touching the layout
                FigureSpec.overlaysOff().applyTo( tp );
                if ( ( tp.getAnnotationColumnSpecs() != null ) && !tp.getAnnotationColumnSpecs().isEmpty() ) {
                    ok[ 0 ] = fail( "clearing overlays must remove the annotation columns" );
                }
                if ( tp.getColorByPropertyRef() != null ) {
                    ok[ 0 ] = fail( "clearing overlays must remove colour-by" );
                }
                if ( tp.getPhylogenyGraphicsType() != Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                    ok[ 0 ] = fail( "clearing OVERLAYS must not change the layout" );
                }
                if ( !tp.shows( DisplayOption.SHOW_TAX_RANK ) ) {
                    ok[ 0 ] = fail( "clearing OVERLAYS must not strip the labels" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf2[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * A figure belongs to ITS tab. Both halves of that are easy to get wrong, and both were: the
     * phylogram/cladogram choice is stored per TAB INDEX, so restoring a figure before its tab was selected wrote
     * it onto whichever tab happened to be in front, and capturing a figure read the front tab's choice for every
     * tab -- which is exactly what "Save All" does with several tabs open.
     */
    private static boolean perTabRestore() {
        if (GraphicsEnvironment.isHeadless()) {
            return true;
        }
        try {
            final boolean[] ok = { true };
            final Phylogeny first = tree();
            final Phylogeny second = tree();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { first }, new Configuration(), "figspec3" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainPanel mp = mf[ 0 ].getMainPanel();
                // give the FIRST tab an opinion that differs from both the defaults and the figure below
                mp.getControlPanel().setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM );
                mp.getCurrentTreePanel().setShows( DisplayOption.SHOW_TAX_RANK, false );
                FigureSpec.writeToTree( second,
                                        FigureSpec.parse( "v1;displaytype=ALIGNED_PHYLOGRAM;show.SHOW_TAX_RANK=true" ) );
                mp.addPhylogenyInNewTab( second, new Configuration(), "second", null );

                if ( mp.getControlPanel().treeDisplayTypeAt( 0 ) != Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM ) {
                    ok[ 0 ] = fail( "opening a figure in a NEW tab must not change the first tab's display type, got "
                            + mp.getControlPanel().treeDisplayTypeAt( 0 ) );
                }
                if ( mp.getTreePanels().get( 0 ).shows( DisplayOption.SHOW_TAX_RANK ) ) {
                    ok[ 0 ] = fail( "...nor the first tab's labels" );
                }
                if ( mp.getControlPanel().treeDisplayTypeAt( 1 ) != Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM ) {
                    ok[ 0 ] = fail( "the new tab must take the figure's display type, got "
                            + mp.getControlPanel().treeDisplayTypeAt( 1 ) );
                }
                if ( !mp.getTreePanels().get( 1 ).shows( DisplayOption.SHOW_TAX_RANK ) ) {
                    ok[ 0 ] = fail( "the new tab must take the figure's labels" );
                }
                if ( !mp.getControlPanel().isCheckboxSelected( DisplayOption.SHOW_TAX_RANK ) ) {
                    ok[ 0 ] = fail( "the checkboxes must show the restored figure -- the new tab is the current one" );
                }
                // ...and "Save All": every tab stamps its OWN figure while only one of them is in front
                mp.getTreePanels().get( 0 ).syncFigureToTree();
                mp.getTreePanels().get( 1 ).syncFigureToTree();
                if ( !"UNALIGNED_PHYLOGRAM".equals( FigureSpec.readFrom( first ).get( "displaytype" ) )
                        || !"ALIGNED_PHYLOGRAM".equals( FigureSpec.readFrom( second ).get( "displaytype" ) ) ) {
                    ok[ 0 ] = fail( "each tab must save its own display type, got "
                            + FigureSpec.readFrom( first ).get( "displaytype" ) + " and "
                            + FigureSpec.readFrom( second ).get( "displaytype" ) );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String host : new String[] { "cat", "dog", "cat", "fish" } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "t_" + host + n.getId() );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:host", host, "", "xsd:string", AppliesTo.NODE ) );
            n.getNodeData().setProperties( pl );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private FigureSpecTest() {
    }
}

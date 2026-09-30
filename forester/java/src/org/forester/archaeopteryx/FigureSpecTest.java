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
        return codec() && parsing() && storage() && roundTrip() && perTabRestore();
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

    /** A file as 0.11.117 to 0.11.172 wrote it: the figure on the ROOT CLADE, not under {@code <phylogeny>}. */
    private static final String LEGACY_ROOT_CLADE_FIGURE = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            + "<phyloxml xmlns=\"http://www.phyloxml.org\">\n"
            + "<phylogeny rooted=\"true\">\n"
            + "<clade>\n"
            + "<property ref=\"aptx:figure\" datatype=\"xsd:string\" applies_to=\"phylogeny\">"
            + "v1;layout=CIRCULAR;columns=data:host\\sCOLOR_STRIP\\sCIRCLE\\sfalse</property>\n"
            + "<clade><name>a</name>"
            + "<property ref=\"data:host\" datatype=\"xsd:string\" applies_to=\"node\">cat</property></clade>\n"
            + "<clade><name>b</name>"
            + "<property ref=\"data:host\" datatype=\"xsd:string\" applies_to=\"node\">dog</property></clade>\n"
            + "</clade>\n"
            + "</phylogeny>\n"
            + "</phyloxml>\n";

    /**
     * WHERE the figure lives: a direct child of {@code <phylogeny>}, like every other piece of per-tree app state
     * (Christian, 2026-09-30) -- never on the root clade, where 0.11.117 to 0.11.172 put it. A HARD BREAK: a figure
     * on a clade is ignored, and saving deletes it, so a re-saved old file carries no figure that is never read.
     */
    private static boolean storage() {
        try {
            // written: at the phylogeny level, and nowhere on the nodes
            final Phylogeny phy = tree();
            FigureSpec.writeToTree( phy, FigureSpec.parse( "v1;layout=CIRCULAR" ) );
            if ( ( figuresAtPhylogenyLevel( phy ) != 1 ) || ( figuresOnNodes( phy ) != 0 ) ) {
                return fail( "a written figure must be ONE phylogeny-level property and none on a node, got "
                        + figuresAtPhylogenyLevel( phy ) + " / " + figuresOnNodes( phy ) );
            }
            // ...and in the saved TEXT, a direct child of <phylogeny> (asked of the XML itself, not of our reader)
            final StringBuffer xml = new org.forester.io.writers.PhylogenyWriter().toPhyloXML( phy, 0 );
            final String parents = figureParentsInXml( xml.toString() );
            if ( !"phylogeny".equals( parents ) ) {
                return fail( "in the file the figure must sit directly under <phylogeny>, found under: " + parents );
            }
            final Phylogeny reread = parseXsdValidating( xml );
            if ( ( FigureSpec.readFrom( reread ) == null )
                    || !"CIRCULAR".equals( FigureSpec.readFrom( reread ).get( "layout" ) ) ) {
                return fail( "a figure saved under <phylogeny> must read back" );
            }
            // an old file: its figure on the root clade is IGNORED -- a well-formed v1 value, so only WHERE it sits
            // can be the reason (parse() of the same value must give a figure)
            final Phylogeny legacy = parseXsdValidating( new StringBuffer( LEGACY_ROOT_CLADE_FIGURE ) );
            if ( ( figuresAtPhylogenyLevel( legacy ) != 0 ) || ( figuresOnNodes( legacy ) != 1 ) ) {
                return fail( "fixture: the legacy file must carry its figure on the root clade only" );
            }
            if ( FigureSpec.parse( "v1;layout=CIRCULAR;columns=data:host\\sCOLOR_STRIP\\sCIRCLE\\sfalse" ) == null ) {
                return fail( "fixture: the legacy figure's value must itself be a readable figure" );
            }
            if ( FigureSpec.readFrom( legacy ) != null ) {
                return fail( "a figure on the root clade (0.11.117 to 0.11.172) must be ignored, not read" );
            }
            // ...and saving a figure onto it leaves exactly one, under <phylogeny>: the dead clade copy is deleted
            FigureSpec.writeToTree( legacy, FigureSpec.parse( "v1;layout=RECTANGULAR" ) );
            if ( ( figuresAtPhylogenyLevel( legacy ) != 1 ) || ( figuresOnNodes( legacy ) != 0 ) ) {
                return fail( "saving must delete an old clade copy, got " + figuresAtPhylogenyLevel( legacy )
                        + " at the phylogeny level / " + figuresOnNodes( legacy ) + " on nodes" );
            }
            final String saved = figureParentsInXml( new org.forester.io.writers.PhylogenyWriter()
                    .toPhyloXML( legacy, 0 ).toString() );
            if ( !"phylogeny".equals( saved ) ) {
                return fail( "a re-saved old file must hold the figure under <phylogeny> only, found: " + saved );
            }
            // both present (a hand-edited file): the phylogeny level wins
            final Phylogeny both = parseXsdValidating( new StringBuffer( LEGACY_ROOT_CLADE_FIGURE ) );
            both.setProperties( new PropertiesList() );
            both.getProperties().addProperty( new Property( FigureSpec.FIGURE_REF, "v1;layout=RECTANGULAR", "",
                                                            "xsd:string", AppliesTo.PHYLOGENY ) );
            if ( !"RECTANGULAR".equals( FigureSpec.readFrom( both ).get( "layout" ) ) ) {
                return fail( "with a figure in both places the one under <phylogeny> must win" );
            }
            // no figure, or an empty one: removed from both places
            FigureSpec.writeToTree( both, null );
            if ( ( FigureSpec.readFrom( both ) != null ) || ( figuresAtPhylogenyLevel( both ) != 0 )
                    || ( figuresOnNodes( both ) != 0 ) ) {
                return fail( "writing no figure must remove it from BOTH places" );
            }
            // an EMPTY figure (not null: capture() of no panel) must not be written as a bare "v1"
            final FigureSpec empty = FigureSpec.capture( null );
            if ( ( empty == null ) || !empty.isEmpty() ) {
                return fail( "fixture: capture(null) must give an empty, non-null figure" );
            }
            FigureSpec.writeToTree( phy, empty );
            if ( ( figuresAtPhylogenyLevel( phy ) != 0 ) || ( figuresOnNodes( phy ) != 0 ) ) {
                return fail( "an empty figure must not be written" );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "storage: " + e );
        }
    }

    private static Phylogeny parseXsdValidating( final StringBuffer xml ) throws java.io.IOException {
        final org.forester.io.parsers.phyloxml.PhyloXmlParser p = org.forester.io.parsers.phyloxml.PhyloXmlParser
                .createPhyloXmlParserXsdValidating();
        p.setSource( xml );
        return p.parse()[ 0 ];
    }

    private static int figuresAtPhylogenyLevel( final Phylogeny phy ) {
        int n = 0;
        if ( phy.getProperties() != null ) {
            for( final Property p : phy.getProperties().getProperties() ) {
                if ( FigureSpec.FIGURE_REF.equals( p.getRef() ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static int figuresOnNodes( final Phylogeny phy ) {
        int n = 0;
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PropertiesList pl = it.next().getNodeData().getProperties();
            if ( pl != null ) {
                for( final Property p : pl.getProperties() ) {
                    if ( FigureSpec.FIGURE_REF.equals( p.getRef() ) ) {
                        ++n;
                    }
                }
            }
        }
        return n;
    }

    /** The parent element of every {@code aptx:figure} property in the XML text, comma-separated in document order. */
    private static String figureParentsInXml( final String xml ) throws Exception {
        final javax.xml.parsers.DocumentBuilderFactory f = javax.xml.parsers.DocumentBuilderFactory.newInstance();
        f.setNamespaceAware( true );
        final org.w3c.dom.NodeList props = f.newDocumentBuilder()
                .parse( new org.xml.sax.InputSource( new java.io.StringReader( xml ) ) )
                .getElementsByTagNameNS( "*", "property" );
        final StringBuilder sb = new StringBuilder();
        for( int i = 0; i < props.getLength(); ++i ) {
            final org.w3c.dom.Element e = (org.w3c.dom.Element) props.item( i );
            if ( FigureSpec.FIGURE_REF.equals( e.getAttribute( "ref" ) ) ) {
                if ( sb.length() > 0 ) {
                    sb.append( ',' );
                }
                sb.append( e.getParentNode().getLocalName() );
            }
        }
        return sb.toString();
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
                // the radial layouts' label direction is part of a figure (a review find, 2026-09-27: the default
                // became RADIAL, and a spec that could not pin it would render differently across versions)
                final TreePanel tp1 = mp.getTreePanels().get( 1 );
                tp1.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                FigureSpec.parse( "v1;labeldirection=HORIZONTAL" ).applyTo( tp1 );
                if ( tp1.getOptions().getNodeLabelDirection() != Options.NODE_LABEL_DIRECTION.HORIZONTAL ) {
                    ok[ 0 ] = fail( "a figure must be able to pin flat labels, got "
                            + tp1.getOptions().getNodeLabelDirection() );
                }
                // ...through the MENU, not the option alone: updateOptions rewrites the option from the checkbox on
                // every menu action, so a direction set behind the checkbox's back lasted exactly until the next
                // click anywhere in the menus (a second review find, 2026-09-27)
                if ( mf[ 0 ]._label_direction_cbmi.isSelected() ) {
                    ok[ 0 ] = fail( "a figure's flat labels must uncheck the Radial Labels menu item" );
                }
                mf[ 0 ].updateOptions( mf[ 0 ].getOptions() );
                if ( tp1.getOptions().getNodeLabelDirection() != Options.NODE_LABEL_DIRECTION.HORIZONTAL ) {
                    ok[ 0 ] = fail( "a figure's label direction must survive the next menu action (updateOptions), got "
                            + tp1.getOptions().getNodeLabelDirection() );
                }
                // a RECTANGULAR figure carries no direction -- it draws none, and applying one would reset the
                // user's frame-wide preference on opening; a radial figure carries it
                if ( FigureSpec.capture( tp1 ).get( "labeldirection" ) != null ) {
                    ok[ 0 ] = fail( "a rectangular figure must not carry a label direction, got "
                            + FigureSpec.capture( tp1 ).get( "labeldirection" ) );
                }
                tp1.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                if ( !"HORIZONTAL".equals( FigureSpec.capture( tp1 ).get( "labeldirection" ) ) ) {
                    ok[ 0 ] = fail( "a captured radial figure must carry the label direction, got "
                            + FigureSpec.capture( tp1 ).get( "labeldirection" ) );
                }
                FigureSpec.parse( "v1;labeldirection=RADIAL" ).applyTo( tp1 );
                if ( ( tp1.getOptions().getNodeLabelDirection() != Options.NODE_LABEL_DIRECTION.RADIAL )
                        || !mf[ 0 ]._label_direction_cbmi.isSelected() ) {
                    ok[ 0 ] = fail( "...and radial ones, menu item included" );
                }
                tp1.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
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

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

import java.awt.Component;
import java.awt.Container;
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.AbstractButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods.NDF;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * "Metadata", one word for the fields a node carries (Christian, 2026-10-01, in both programs, "everywhere, but only
 * what users see"): phyloXML calls them properties, the earlier menus called them annotations, and users know
 * neither. Machine-facing names (the figure keys labelprops= / show.SHOW_PROPERTIES, enum and class names, the
 * phyloXML &lt;property&gt; element, demo FILE names) deliberately keep the old word.
 * <p>
 * The pinned strings are the ones no other test reads; the live sweep walks every menu item, control and tooltip of a
 * real window and fails on the old words in the field sense, so a new control cannot bring "property" back unnoticed.
 * The window parts are headful; a green no-op when headless.
 */
public final class MetadataWordingTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MetadataWording: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( !pinnedOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return liveWindowOk();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [MetadataWordingTest] " + msg );
        return false;
    }

    /** The old words, in the sense of the fields a node carries. */
    private static final Pattern OLD = Pattern.compile( "(?i)(annotation field|annotation column|import annotation"
            + "|re-import annotation|annotation table|annotation source|annotation import|\\bpropert(y|ies)\\b)" );

    /** What may still say it, and why: the Tree Properties WINDOW (the tree's own characteristics, kept in both
     *  programs), and the one deliberate "(its phyloXML properties: ...)" that tells a phyloXML user where the
     *  metadata comes from. Removed before the text is matched. */
    private static final String[] ALLOWED = { "Tree Properties", "Tree properties", "tree properties",
            "its phyloXML properties" };

    static String offending( final String text ) {
        if ( text == null ) {
            return null;
        }
        String t = text;
        for( final String a : ALLOWED ) {
            t = t.replace( a, "" );
        }
        return OLD.matcher( t ).find() ? text : null;
    }

    private static boolean pinnedOk() {
        if ( !"Metadata".equals( DisplayOption.SHOW_PROPERTIES.title() ) ) {
            return fail( "the Display Data checkbox must read \"Metadata\", got \"" + DisplayOption.SHOW_PROPERTIES.title()
                    + "\"" );
        }
        // the node window's section, read through the code that names it (writeTo reports the sections it changed):
        // SEC_PROPERTIES itself is a compile-time constant, so comparing it here would compare this class's own copy
        final PhylogenyNode n = tree().getRoot().getChildNode( 0 );
        final NodeDataDraft base = NodeDataDraft.from( n );
        final NodeDataDraft edited = base.copy();
        edited.properties.get( 0 ).value = "lion";
        final String sections = edited.writeTo( n, base ).toString();
        if ( !"[Metadata]".equals( sections ) ) {
            return fail( "the node window's section must read \"Metadata\", got " + sections );
        }
        if ( !"Any Metadata Field".equals( SearchField.ofNdf( NDF.Properties ).label() ) ) {
            return fail( "the search scope must read \"Any Metadata Field\", got \""
                    + SearchField.ofNdf( NDF.Properties ).label() + "\"" );
        }
        for( final SearchField f : SearchField.stringMenuFields() ) {
            if ( offending( f.label() ) != null ) {
                return fail( "a search scope still says it: \"" + f.label() + "\"" );
            }
        }
        // the filter itself: it must flag the old words and pass what is allowed, or the sweep below proves nothing
        if ( ( offending( "Import Annotations (CSV/TSV)..." ) == null ) || ( offending( "color by a node property" ) == null )
                || ( offending( "Annotation Fields…" ) == null ) || ( offending( "Properties" ) == null ) ) {
            return fail( "the filter must flag the old words" );
        }
        if ( ( offending( "Tree Properties…" ) != null ) || ( offending( "Metadata Fields…" ) != null )
                || ( offending( "Read [&...] Annotations (BEAST, MrBayes)" ) != null )
                || ( offending( "Show the metadata a node carries (its phyloXML properties: the ref)" ) != null ) ) {
            return fail( "the filter must pass the new words and the deliberate exceptions" );
        }
        return true;
    }

    private static boolean liveWindowOk() {
        final boolean[] ok = { true };
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { tree() }, new Configuration(), "metadata-wording" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final List<String> seen = new ArrayList<>();
                    final JMenuBar bar = mf[ 0 ].getJMenuBar();
                    for( int i = 0; i < bar.getMenuCount(); ++i ) {
                        collect( bar.getMenu( i ), seen );
                    }
                    collect( mf[ 0 ].getContentPane(), seen );
                    // the sweep must have REACHED what was renamed, or a clean result means nothing
                    for( final String must : new String[] { "Metadata Fields…", "Import Metadata (CSV/TSV)...",
                            "Import Metadata from URL...", "Re-import Metadata", "Metadata" } ) {
                        if ( !seen.contains( must ) ) {
                            ok[ 0 ] = fail( "fixture: the sweep never reached \"" + must + "\"" );
                        }
                    }
                    for( final String s : seen ) {
                        if ( offending( s ) != null ) {
                            ok[ 0 ] = fail( "a user-visible text still says it: \"" + s + "\"" );
                        }
                    }
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return fail( "unexpected: " + e );
        }
        return ok[ 0 ];
    }

    /** Every text and tooltip a user can see in {@code c}: menus and their submenus, buttons, check boxes, labels. */
    private static void collect( final Component c, final List<String> seen ) {
        if ( c instanceof AbstractButton ) {
            add( seen, ( (AbstractButton) c ).getText() );
        }
        if ( c instanceof JLabel ) {
            add( seen, ( (JLabel) c ).getText() );
        }
        if ( c instanceof JComponent ) {
            add( seen, ( (JComponent) c ).getToolTipText() );
        }
        if ( c instanceof JMenu ) {
            for( final Component m : ( (JMenu) c ).getMenuComponents() ) {
                collect( m, seen );
            }
        }
        if ( c instanceof Container ) {
            for( final Component k : ( (Container) c ).getComponents() ) {
                collect( k, seen );
            }
        }
    }

    private static void add( final List<String> seen, final String s ) {
        if ( ( s != null ) && !s.isEmpty() ) {
            seen.add( s );
        }
    }

    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String host : new String[] { "cat", "dog", "fish" } ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "t_" + host );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:host", host, "", "xsd:string", AppliesTo.NODE ) );
            n.getNodeData().setProperties( pl );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private MetadataWordingTest() {
    }
}

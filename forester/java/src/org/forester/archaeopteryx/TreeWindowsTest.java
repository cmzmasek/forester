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
import java.awt.Rectangle;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.TreeText.Format;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Tests for the two per-tab tree windows inside a real {@link MainFrameApplication} (needs a display):
 * <ul>
 * <li>{@link TreePropertiesFrame}: opens once per tab (a second open brings the same window forward), shows the
 * file / tip subtitle, gets the dirty bullet, is marked stale by an edit and re-reads the tree (keeping unwritten
 * edits across an undo that swaps the tree), writes through the chrome, and is closed with its tab;</li>
 * <li>{@link TreeTextFrame}: shows the format the menu item asked for, switches formats in place, tints markup
 * and leaves names plain, wraps by default only for the one-line formats, finds and steps through hits (wrapping
 * around, "no matches" in red), regenerates after the tree changed, and is closed with its tab;</li>
 * <li>the View menu carries "Tree Properties…" (Cmd-I) and the three "as ..." items.</li>
 * </ul>
 */
public final class TreeWindowsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeWindows: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = TreeFactsTest.fixture();
            phy.setName( "Fixture tree" );
            final Phylogeny other = TreePropertiesEditTest.threeLevel( "Other" );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy, other }, new Configuration(), "fixture.xml" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    menu( ok, mf[ 0 ] );
                    propertiesFrame( ok, mf[ 0 ] );
                    textFrame( ok, mf[ 0 ] );
                    closedWithTab( ok, mf[ 0 ] );
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    ok[ 0 ] = false;
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void menu( final boolean[] ok, final MainFrame mf ) {
        JMenu view = null;
        for( int i = 0; i < mf.getJMenuBar().getMenuCount(); ++i ) {
            final JMenu m = mf.getJMenuBar().getMenu( i ); // null for a non-menu bar component (e.g. glue)
            if ( ( m != null ) && "View".equals( m.getText() ) ) {
                view = m;
            }
        }
        check( ok, "View menu", view != null );
        final StringBuilder labels = new StringBuilder();
        JMenuItem props = null;
        for( int i = 0; ( view != null ) && ( i < view.getItemCount() ); ++i ) {
            final JMenuItem it = view.getItem( i );
            if ( it != null ) {
                labels.append( it.getText() ).append( '|' );
                if ( MainFrame.TREE_PROPERTIES_LABEL.equals( it.getText() ) ) {
                    props = it;
                }
            }
        }
        check( ok, "menu starts with Tree Properties then the three text formats: " + labels,
               labels.toString().startsWith( "Tree Properties…|as phyloXML|as Newick|as Nexus|" ) );
        check( ok, "no old items", !labels.toString().contains( "Basic Tree Information" )
                && !labels.toString().contains( "Edit Tree Name" ) );
        final int mask = java.awt.Toolkit.getDefaultToolkit().getMenuShortcutKeyMaskEx();
        check( ok, "Cmd-I", ( props != null ) && KeyStroke.getKeyStroke( java.awt.event.KeyEvent.VK_I, mask )
                .equals( props.getAccelerator() ) );
        check( ok, "tooltips", ( props != null ) && ( props.getToolTipText() != null ) );
    }

    private static void propertiesFrame( final boolean[] ok, final MainFrame mf ) throws Exception {
        mf.getMainPanel().getTabbedPane().setSelectedIndex( 0 ); // the fixture tree (the last tab opens selected)
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
        check( ok, "nothing open", tp.treePropertiesFrameForTest() == null );
        mf.showTreeProperties();
        final TreePropertiesFrame f = tp.treePropertiesFrameForTest();
        check( ok, "opened", ( f != null ) && f.isDisplayable() );
        if ( f == null ) {
            return;
        }
        check( ok, "title: " + f.getTitle(), "Tree Properties: Fixture tree".equals( f.getTitle() ) );
        check( ok, "subtitle names the file and tips: " + f.getForm().headerSubtitleForTest(),
               f.getForm().headerSubtitleForTest().startsWith( tp.getTreeFile().getName() + " · 5 tips" ) );
        final Rectangle usable = GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
        check( ok, "fits the usable screen: " + f.getBounds(), usable.contains( f.getBounds() ) );
        check( ok, "time axis section (the panel provides it; the fixture has dates)",
               f.getForm().hasSectionForTest( TreeFacts.TIME_AXIS ) );
        // opening again: the same window
        mf.showTreeProperties();
        check( ok, "same window on re-open", tp.treePropertiesFrameForTest() == f );
        // an edit elsewhere marks it stale; the coalesced refresh is pending, and re-reads on demand
        tp.setEdited( true );
        check( ok, "stale after an edit", f.isRefreshPendingForTest() );
        f.rebindNow();
        check( ok, "refresh done", !f.isRefreshPendingForTest() );
        check( ok, "subtitle notes unsaved changes: " + f.getForm().headerSubtitleForTest(),
               f.getForm().headerSubtitleForTest().endsWith( "unsaved changes" ) );
        // typing makes it dirty (bullet + mark); writing through the chrome renames tree, tab and title
        f.getForm().setTextForTest( TreePropertiesDraft.NAME, "Renamed tree" );
        check( ok, "dirty bullet", f.getTitle().startsWith( "• " ) && f.isDocumentModifiedMark() );
        check( ok, "status", "Unsaved changes".equals( f.statusTextForTest() ) );
        check( ok, "written", f.writeNow() );
        check( ok, "tree renamed", "Renamed tree".equals( tp.getPhylogeny().getName() ) );
        check( ok, "tab renamed", "Renamed tree".equals( mf.getMainPanel().getTabbedPane().getTitleAt( 0 ) ) );
        check( ok, "title follows: " + f.getTitle(), "Tree Properties: Renamed tree".equals( f.getTitle() ) );
        check( ok, "undoable", tp.canUndo() );
        // an unwritten edit survives an undo that swaps the tree underneath; the baseline is the restored tree
        f.getForm().setTextForTest( TreePropertiesDraft.DESCRIPTION, "unwritten" );
        mf.undo();
        check( ok, "undo restored the name", "Fixture tree".equals( tp.getPhylogeny().getName() ) );
        check( ok, "window still open after undo", f.isDisplayable() && ( tp.treePropertiesFrameForTest() == f ) );
        f.rebindNow();
        check( ok, "unwritten edit kept", "unwritten".equals( f.getForm().collect().description ) );
        check( ok, "baseline is the restored tree", "Fixture tree".equals( f.getForm().baseline().name ) );
        check( ok, "still dirty", f.isDirty() );
        check( ok, "title back: " + f.getTitle(), f.getTitle().endsWith( "Tree Properties: Fixture tree" ) );
        check( ok, "write goes to the restored tree", f.writeNow() && "unwritten".equals( tp.getPhylogeny().getDescription() ) );
        // close releases the slot; a re-open makes a new window
        f.close();
        check( ok, "slot released", tp.treePropertiesFrameForTest() == null );
        check( ok, "disposed", !f.isDisplayable() );
        mf.showTreeProperties();
        check( ok, "re-opened as a new window", ( tp.treePropertiesFrameForTest() != null )
                && ( tp.treePropertiesFrameForTest() != f ) );
        tp.treePropertiesFrameForTest().close();
    }

    private static void textFrame( final boolean[] ok, final MainFrame mf ) {
        mf.getMainPanel().getTabbedPane().setSelectedIndex( 0 );
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
        check( ok, "nothing open", tp.treeTextFrameForTest() == null );
        mf.viewAsText( Format.NEWICK );
        final TreeTextFrame f = tp.treeTextFrameForTest();
        check( ok, "opened", ( f != null ) && f.isDisplayable() );
        if ( f == null ) {
            return;
        }
        check( ok, "on Newick", f.format() == Format.NEWICK );
        check( ok, "title: " + f.getTitle(), f.getTitle().startsWith( "Tree as Newick: " ) );
        final String nwk = TreeText.render( tp.getPhylogeny(), Format.NEWICK,
                                            tp.getOptions().getNhConversionSupportValueStyle() );
        check( ok, "text is the rendering", nwk.equals( f.textForTest() ) );
        check( ok, "subtitle: " + f.subtitleForTest(), f.subtitleForTest().startsWith( "Newick · " )
                && f.subtitleForTest().contains( "1 line" ) );
        check( ok, "Newick wraps by default", f.isWrapForTest() );
        // tint: a bracket is muted, a name is plain
        final Color muted = FormWidgets.mutedColor();
        check( ok, "bracket muted", muted.equals( f.colorAtForTest( nwk.indexOf( '(' ) ) ) );
        final int name_at = nwk.indexOf( "B:" );
        check( ok, "name plain", ( name_at >= 0 ) && !muted.equals( f.colorAtForTest( name_at ) ) );
        // the menu's other items switch the SAME window
        mf.viewAsText( Format.PHYLOXML );
        check( ok, "same window", tp.treeTextFrameForTest() == f );
        check( ok, "switched", f.format() == Format.PHYLOXML );
        check( ok, "phyloXML does not wrap", !f.isWrapForTest() );
        check( ok, "phyloXML text", f.textForTest().contains( "<phyloxml" ) );
        f.setWrapForTest( true );
        check( ok, "wrap toggled", f.isWrapForTest() );
        mf.viewAsText( Format.NEWICK );
        check( ok, "wrap remembered per format", f.isWrapForTest() );
        mf.viewAsText( Format.PHYLOXML );
        check( ok, "... and per format", f.isWrapForTest() );
        // find: hits, stepping with wrap-around, no matches
        f.setFindTextForTest( "clade" );
        final int hits = f.hitCountForTest();
        check( ok, "hits found: " + hits, hits >= 2 );
        check( ok, "first hit selected", f.currentHitForTest() == 0 );
        f.stepForTest( 1 );
        check( ok, "next", f.currentHitForTest() == 1 );
        f.stepForTest( -1 );
        f.stepForTest( -1 );
        check( ok, "previous wraps around", f.currentHitForTest() == hits - 1 );
        f.stepForTest( 1 );
        check( ok, "next wraps around", f.currentHitForTest() == 0 );
        check( ok, "count label", f.findFieldForTest().getText().equals( "clade" ) );
        f.setFindTextForTest( "zzzz-not-there" );
        check( ok, "no matches", f.hitCountForTest() == 0 );
        f.setFindTextForTest( "CLADE" );
        check( ok, "case-insensitive", f.hitCountForTest() == hits );
        f.setFindTextForTest( "" );
        check( ok, "cleared", f.hitCountForTest() == 0 );
        // the tree changed: stale, then regenerated on demand
        final PhylogenyNode tip = tp.getPhylogeny().getFirstExternalNode();
        tip.setName( "RENAMED_TIP" );
        tp.setEdited( true );
        check( ok, "stale", f.isStaleForTest() );
        check( ok, "old text until refreshed", !f.textForTest().contains( "RENAMED_TIP" ) );
        f.refreshNow();
        check( ok, "regenerated", f.textForTest().contains( "RENAMED_TIP" ) && !f.isStaleForTest() );
        mf.viewAsText( Format.NEXUS );
        check( ok, "other format regenerated too", f.textForTest().contains( "RENAMED_TIP" ) );
        f.close();
        check( ok, "slot released", ( tp.treeTextFrameForTest() == null ) && !f.isDisplayable() );
    }

    private static void closedWithTab( final boolean[] ok, final MainFrame mf ) {
        final MainFrameApplication app = (MainFrameApplication) mf;
        mf.getMainPanel().getTabbedPane().setSelectedIndex( 1 ); // the "Other" tree
        final TreePanel tp = mf.getMainPanel().getCurrentTreePanel();
        mf.showTreeProperties();
        mf.viewAsText( Format.NEXUS );
        final TreePropertiesFrame pf = tp.treePropertiesFrameForTest();
        final TreeTextFrame tf = tp.treeTextFrameForTest();
        check( ok, "both open", ( pf != null ) && ( tf != null ) );
        tp.showNodeEditFrame( tp.getPhylogeny().getFirstExternalNode() );
        check( ok, "a node window too", tp.openNodeFrameCountForTest() == 1 );
        app.closeTabAt( 1 );
        check( ok, "tab closed", mf.getMainPanel().getTabbedPane().getTabCount() == 1 );
        check( ok, "node window closed with the tab", tp.openNodeFrameCountForTest() == 0 );
        check( ok, "properties window closed with the tab", ( pf != null ) && !pf.isDisplayable() );
        check( ok, "text window closed with the tab", ( tf != null ) && !tf.isDisplayable() );
    }

    private static void check( final boolean[] ok, final String what, final boolean condition ) {
        if ( !condition ) {
            System.out.println( "  [TreeWindowsTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private TreeWindowsTest() {
    }
}

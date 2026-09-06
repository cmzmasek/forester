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
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Headful (skip-when-headless) tests for {@link NodeFrame}, the node window: title and status follow the form's
 * dirty/valid state, the Write button is enabled only when there is something valid to write, a write clears the
 * dirty state, the VIEW window has no Write button, and the window fits the usable screen.
 */
public final class NodeFrameTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeFrame: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // needs a display
        }
        try {
            final boolean[] ok = { true };
            final PhylogenyNode n = NodeDataDraftTest.richNode();
            final NodeFrame[] f = new NodeFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> f[ 0 ] = new NodeFrame( n, null, 0, NodeDataForm.Mode.EDIT ) );
            SwingUtilities.invokeAndWait( () -> {
                final NodeFrame fr = f[ 0 ];
                check( ok, "edit title", "Edit Node: BRCA1 clade".equals( fr.getTitle() ) );
                check( ok, "clean on open", !fr.isDirty() );
                check( ok, "write disabled while clean", !fr.writeButtonForTest().isEnabled() );
                check( ok, "no document-modified mark", !fr.isDocumentModifiedMark() );
                check( ok, "write is the default button", fr.getRootPane().getDefaultButton() == fr.writeButtonForTest() );
                final Rectangle usable = GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
                check( ok, "fits the usable screen: " + fr.getBounds() + " vs " + usable,
                       usable.contains( fr.getBounds() ) );
                fr.getForm().setTextForTest( NodeDataDraft.NAME, "renamed" );
                check( ok, "dirty title gets the dot", fr.getTitle().startsWith( "• " ) );
                check( ok, "document-modified mark", fr.isDocumentModifiedMark() );
                check( ok, "status says unsaved", "Unsaved changes".equals( fr.statusTextForTest() ) );
                check( ok, "write enabled", fr.writeButtonForTest().isEnabled() );
                fr.getForm().setTextForTest( NodeDataDraft.BRANCH_LENGTH, "zzz" );
                check( ok, "status shows the problem", fr.statusTextForTest().contains( "Branch length" ) );
                check( ok, "write disabled while invalid", !fr.writeButtonForTest().isEnabled() );
                check( ok, "writeNow refuses", !fr.writeNow() );
                fr.getForm().setTextForTest( NodeDataDraft.BRANCH_LENGTH, "0.5" );
                check( ok, "writeNow writes", fr.writeNow() );
                check( ok, "node renamed", "renamed".equals( n.getName() ) );
                check( ok, "clean title after write", "Edit Node: BRCA1 clade".equals( fr.getTitle() ) );
                check( ok, "status confirms", "Written to the tree.".equals( fr.statusTextForTest() ) );
                check( ok, "write disabled again", !fr.writeButtonForTest().isEnabled() );
                fr.requestClose(); // clean -> closes without asking
                check( ok, "closed", !fr.isDisplayable() );
            } );
            SwingUtilities.invokeAndWait( () -> f[ 0 ] = new NodeFrame( n, null, 0, NodeDataForm.Mode.VIEW ) );
            SwingUtilities.invokeAndWait( () -> {
                final NodeFrame fr = f[ 0 ];
                check( ok, "view title", "Node: renamed".equals( fr.getTitle() ) );
                check( ok, "view has no write button", fr.writeButtonForTest() == null );
                check( ok, "view never dirty", !fr.isDirty() );
                fr.close();
                check( ok, "view closed", !fr.isDisplayable() );
            } );
            return ok[ 0 ] && slotsReleasedByIdentity();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * The tree panel tracks its open node windows in a compacted array. A window must release ITS slot when it
     * closes, whatever order the windows close in -- the old code released "the slot I was opened into", which
     * goes stale as soon as an earlier window closes and the array shifts.
     */
    private static boolean slotsReleasedByIdentity() throws Exception {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( final String name : new String[] { "a", "b", "c" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.1 );
            root.addAsChild( tip );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "slots" ) );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
            final List<NodeFrame> frames = new ArrayList<>();
            for( final PhylogenyNode tip : tp.getPhylogeny().getRoot().getDescendants() ) {
                tp.showNodeEditFrame( tip );
                for( final java.awt.Window w : java.awt.Window.getWindows() ) {
                    if ( ( w instanceof NodeFrame ) && w.isDisplayable() && !frames.contains( w ) ) {
                        frames.add( (NodeFrame) w );
                    }
                }
            }
            check( ok, "three windows tracked", ( tp.openNodeFrameCountForTest() == 3 ) && ( frames.size() == 3 ) );
            frames.get( 0 ).close(); // the FIRST one: the other two shift down a slot
            check( ok, "two left", tp.openNodeFrameCountForTest() == 2 );
            frames.get( 1 ).close(); // now in slot 0, but opened into slot 1: the old code nulled slot 1 -- i.e.
                                     // dropped the still-open THIRD window and kept tracking this closed one
            check( ok, "one left", tp.openNodeFrameCountForTest() == 1 );
            check( ok, "the open window is the one still tracked", tp.isNodeFrameTrackedForTest( frames.get( 2 ) ) );
            check( ok, "the closed window is not", !tp.isNodeFrameTrackedForTest( frames.get( 1 ) ) );
            frames.get( 2 ).close();
            check( ok, "none left", tp.openNodeFrameCountForTest() == 0 );
            frames.get( 2 ).close(); // closing twice is harmless
            check( ok, "still none", tp.openNodeFrameCountForTest() == 0 );
            ( (JFrame) mf[ 0 ] ).dispose();
        } );
        return ok[ 0 ];
    }

    private static void check( final boolean[] ok, final String what, final boolean condition ) {
        if ( !condition ) {
            System.out.println( "  [NodeFrameTest] " + what );
            ok[ 0 ] = false;
        }
    }

    private NodeFrameTest() {
        // not instantiable
    }
}

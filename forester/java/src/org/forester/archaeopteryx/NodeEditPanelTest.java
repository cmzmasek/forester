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

import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Tests that the "Edit Node Data" dialog ({@link NodeEditPanel}) lets the user edit only field VALUES, never the
 * field-NAME labels ("Name", "Scientific name", ...) or the category headers ("Basic", "Taxonomy", ...) -- which are
 * the dialog's structure, not data. Headless: {@link NodeEditPanel} is a {@code JPanel}, so it (and the pure
 * {@code isPathEditable} logic) can be exercised without a display.
 */
public final class NodeEditPanelTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeEditPanel: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            final PhylogenyNode node = new PhylogenyNode();
            node.setName( "MyTip" );
            final Taxonomy tax = new Taxonomy();
            tax.setScientificName( "Homo sapiens" );
            node.getNodeData().addTaxonomy( tax );

            final NodeEditPanel panel = new NodeEditPanel( node, null, null );
            final JTree jt = panel.getJTreeForTest();
            final DefaultMutableTreeNode root = (DefaultMutableTreeNode) jt.getModel().getRoot();

            final DefaultMutableTreeNode basic = childByLabel( root, NodePanel.BASIC );
            if ( basic == null ) {
                return fail( "the 'Basic' category should be present" );
            }
            final DefaultMutableTreeNode name_label = childByLabel( basic, NodePanel.NODE_NAME );
            if ( ( name_label == null ) || ( name_label.getChildCount() != 1 ) ) {
                return fail( "the 'Name' field should have exactly one value child" );
            }
            final DefaultMutableTreeNode name_value = (DefaultMutableTreeNode) name_label.getChildAt( 0 );

            // the VALUE under a field label is editable ...
            if ( !editable( jt, name_value ) ) {
                return fail( "a field's VALUE must be editable" );
            }
            // ... but the field-NAME label ("Name"), the category header ("Basic") and the root are NOT (the bug)
            if ( editable( jt, name_label ) ) {
                return fail( "the field-NAME label ('Name') must NOT be editable" );
            }
            if ( editable( jt, basic ) ) {
                return fail( "a category header ('Basic') must NOT be editable" );
            }
            if ( editable( jt, root ) ) {
                return fail( "the root node must NOT be editable" );
            }

            // A non-empty field's value does NOT auto-open the editor (that path already works by clicking) ...
            if ( panel.wouldOpenEditorForEmptyFieldForTest( new TreePath( name_value.getPath() ) ) ) {
                return fail( "a NON-empty field value must not auto-open the editor" );
            }
            // ... nor does a field-name label.
            if ( panel.wouldOpenEditorForEmptyFieldForTest( new TreePath( name_label.getPath() ) ) ) {
                return fail( "a field-name label must not auto-open the editor" );
            }

            // an EMPTY field must stay editable AND be the one that auto-opens the inline editor (a blank tree row has
            // no clickable hit region to start editing on, so empty fields could not otherwise be filled in)
            final PhylogenyNode nameless = new PhylogenyNode(); // no name -> the "Name" field is empty
            final NodeEditPanel panel2 = new NodeEditPanel( nameless, null, null );
            final JTree jt2 = panel2.getJTreeForTest();
            final DefaultMutableTreeNode basic2 = childByLabel( (DefaultMutableTreeNode) jt2.getModel().getRoot(),
                                                                NodePanel.BASIC );
            final DefaultMutableTreeNode name_label2 = childByLabel( basic2, NodePanel.NODE_NAME );
            if ( ( name_label2 == null ) || ( name_label2.getChildCount() != 1 ) ) {
                return fail( "the nameless node's 'Name' field should have one value child" );
            }
            final DefaultMutableTreeNode empty_value = (DefaultMutableTreeNode) name_label2.getChildAt( 0 );
            if ( !String.valueOf( empty_value.getUserObject() ).isEmpty() ) {
                return fail( "a nameless node's 'Name' value should be empty" );
            }
            if ( !editable( jt2, empty_value ) ) {
                return fail( "an EMPTY field value must remain editable" );
            }
            if ( !panel2.wouldOpenEditorForEmptyFieldForTest( new TreePath( empty_value.getPath() ) ) ) {
                return fail( "an EMPTY field must auto-open the inline editor so it can be filled in" );
            }
            return true;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean editable( final JTree jt, final DefaultMutableTreeNode node ) {
        return jt.isPathEditable( new TreePath( node.getPath() ) );
    }

    private static DefaultMutableTreeNode childByLabel( final DefaultMutableTreeNode parent, final String label ) {
        if ( parent == null ) {
            return null;
        }
        for( int i = 0; i < parent.getChildCount(); ++i ) {
            final DefaultMutableTreeNode c = (DefaultMutableTreeNode) parent.getChildAt( i );
            if ( label.equals( String.valueOf( c.getUserObject() ) ) ) {
                return c;
            }
        }
        return null;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [NodeEditPanelTest] " + message );
        return false;
    }

    private NodeEditPanelTest() {
        // not instantiable
    }
}

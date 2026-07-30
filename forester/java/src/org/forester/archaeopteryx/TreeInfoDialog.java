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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;

import org.forester.phylogeny.Phylogeny;
import org.forester.util.ForesterUtil;

/**
 * A small modal editor for a tree's <b>Name</b> and <b>Description</b> -- the two human-readable,
 * figure-facing metadata fields on a {@link Phylogeny}. Before this dialog the description was
 * effectively write-only-by-features (e.g. the tree-mutation provenance sentence appended by editing
 * operations) and read-only-by-user (only visible buried in the "Basic Tree Information" stats dump),
 * and the name had no rename affordance at all.
 * <p>
 * On OK, if either field actually changed, {@link #apply} pushes an undo checkpoint (name and
 * description are carried by {@link Phylogeny#copy()}, so undo restores them), applies the change,
 * renames the current tab to match the name so the two cannot drift, marks the panel edited, and
 * repaints. Reachable from the View menu and from a tab right-click / double-click (see
 * {@link MainPanel}).
 */
final class TreeInfoDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    /** Opens the editor for the currently displayed tree of {@code mf}; a no-op if there is no non-empty tree. */
    static void showFor( final MainFrame mf ) {
        final TreePanel tp = mf.getCurrentTreePanel();
        if ( ( tp == null ) || ( tp.getPhylogeny() == null ) || tp.getPhylogeny().isEmpty() ) {
            return;
        }
        if ( tp.isCurrentTreeIsSubtree() ) {
            // The displayed tree is a transient sub-tree; editing its name/description would be discarded on
            // returning to the whole tree (and its undo history is cleared at the navigation boundary). Steer
            // the user back rather than silently losing the edit.
            JOptionPane.showMessageDialog( mf,
                                           "Return to the whole tree (\"Return to super-tree\") to edit its name "
                                                   + "and description.",
                                           "Sub-tree displayed", JOptionPane.INFORMATION_MESSAGE );
            return;
        }
        new TreeInfoDialog( mf, tp ).setVisible( true );
    }

    private TreeInfoDialog( final MainFrame mf, final TreePanel tp ) {
        super( mf, "Tree Info", true );
        final Phylogeny phy = tp.getPhylogeny();
        final JTextField name_field = new JTextField( phy.getName() == null ? "" : phy.getName(), 28 );
        final JTextArea desc_area = new JTextArea( phy.getDescription() == null ? "" : phy.getDescription(), 8, 28 );
        desc_area.setLineWrap( true );
        desc_area.setWrapStyleWord( true );

        final JPanel name_row = new JPanel( new BorderLayout( 6, 0 ) );
        name_row.add( new JLabel( "Name:" ), BorderLayout.WEST );
        name_row.add( name_field, BorderLayout.CENTER );

        final JPanel desc_head = new JPanel( new FlowLayout( FlowLayout.LEFT, 0, 0 ) );
        desc_head.add( new JLabel( "Description:" ) );
        final JPanel center = new JPanel( new BorderLayout( 0, 4 ) );
        center.add( desc_head, BorderLayout.NORTH );
        center.add( new JScrollPane( desc_area, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                                     ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER ), BorderLayout.CENTER );

        final JPanel form = new JPanel( new BorderLayout( 0, 8 ) );
        form.setBorder( BorderFactory.createEmptyBorder( 12, 14, 8, 14 ) );
        form.add( name_row, BorderLayout.NORTH );
        form.add( center, BorderLayout.CENTER );

        final JButton ok = new JButton( "OK" );
        ok.addActionListener( e -> {
            apply( mf, tp, name_field.getText(), desc_area.getText() );
            dispose();
        } );
        final JButton cancel = new JButton( "Cancel" );
        cancel.addActionListener( e -> dispose() );
        final JPanel south = new JPanel( new FlowLayout( FlowLayout.RIGHT ) );
        south.add( cancel );
        south.add( ok );

        setLayout( new BorderLayout() );
        add( form, BorderLayout.CENTER );
        add( south, BorderLayout.SOUTH );
        getRootPane().setDefaultButton( ok );
        pack();
        setMinimumSize( new Dimension( 400, 320 ) );
        setLocationRelativeTo( mf );
    }

    /**
     * Applies edited name/description to {@code tp}'s tree iff either actually changed. The name is a one-line
     * label (tab title, figure title), so its whitespace is cleaned up -- any run of whitespace (incl. tabs or
     * newlines from a paste) is collapsed to a single space and the ends trimmed; the (multi-line) description is
     * only trimmed, to preserve its formatting. Both sides are normalized before comparison, so re-confirming an
     * unedited field is a no-op. An empty name is not a meaningful rename (it would blank the tab label and the
     * on-disk tree name), so the existing name is kept in that case and only the description may change. On a real
     * change it checkpoints for undo before mutating, renames the current tab to match the name, marks the panel
     * edited, and repaints. Package-visible and static so it can be exercised in tests without driving the Swing
     * widgets. Returns {@code true} iff a change was applied.
     */
    static boolean apply( final MainFrame mf, final TreePanel tp, final String raw_name, final String raw_desc ) {
        final Phylogeny phy = tp.getPhylogeny();
        if ( ( phy == null ) || phy.isEmpty() || tp.isCurrentTreeIsSubtree() ) {
            return false; // no tree, or a transient sub-tree whose metadata edit would be discarded (see showFor)
        }
        // normalize BOTH sides so a stored value with stray/collapsible whitespace does not read as changed
        final String old_name = cleanName( phy.getName() );
        final String old_desc = phy.getDescription() == null ? "" : phy.getDescription().trim();
        final String typed_name = cleanName( raw_name );
        final String new_desc = raw_desc == null ? "" : raw_desc.trim();
        // an empty typed name is not a rename -- keep the current name (only the description may change)
        final String new_name = typed_name.isEmpty() ? old_name : typed_name;
        if ( new_name.equals( old_name ) && new_desc.equals( old_desc ) ) {
            return false; // nothing changed -- no undo checkpoint, no edited flag, no repaint
        }
        tp.pushUndoCheckpoint( "Edit Tree Info" ); // Phylogeny.copy() carries name + description, so undo restores them
        phy.setName( new_name );
        phy.setDescription( new_desc );
        if ( !new_name.equals( old_name ) && ( mf.getMainPanel() != null ) ) {
            // keep the tab label and the tree name in sync so they cannot drift (the tab title is what a
            // subsequent save writes back onto the tree name); undo/redo re-syncs it the other way
            mf.getMainPanel().setTitleOfSelectedTab( new_name );
        }
        tp.setEdited( true );
        tp.repaint();
        return true;
    }

    /**
     * Cleans up a tree name for use as a one-line label: collapses every run of whitespace (spaces, tabs, and
     * any newlines from a paste) into a single space and trims the ends. Null becomes "".
     */
    private static String cleanName( final String s ) {
        return s == null ? "" : ForesterUtil.collapseWhiteSpace( s ).trim();
    }
}

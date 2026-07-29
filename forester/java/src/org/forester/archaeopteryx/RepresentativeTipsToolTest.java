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
import java.awt.Dimension;
import java.awt.GraphicsEnvironment;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.NodeDataExporter;
import org.forester.archaeopteryx.tools.RepresentativeTipSelector;
import org.forester.archaeopteryx.tools.RepresentativeTipSelector.RepresentativePick;
import org.forester.archaeopteryx.tools.RepresentativeTipSelector.SelectionResult;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test for the "Select Representative Tips" tool: on a real {@link MainFrame}/{@link TreePanel}
 * it drives the two observable behaviors the (modal) menu handler relies on — highlighting the representatives
 * via the search-hit mechanism ({@code setFoundNodes0}), and extracting the representatives into a new tab
 * (copy → recompute → prune → {@code addPhylogenyInNewTab}) — and checks the menu item is present with a
 * tooltip. Guarded to a no-op on a headless box (needs FlatLaf via {@code createInstance}); run standalone or
 * as part of the non-headless suite.
 */
public final class RepresentativeTipsToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RepresentativeTipsTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { twoCherries() }, conf,
                                                                        "reps" ) );
            final boolean[] ok = { true };
            // gating: distance-cutoff is disabled (target preselected) when the tree has no branch lengths;
            // the "protect selected" checkbox is enabled+on only when tips are selected
            SwingUtilities.invokeAndWait( () -> {
                final MainFrameApplication.RepSelInputs no_bl = MainFrameApplication.buildRepSelInputs( false, true,
                                                                                                        0.05, 100, 0, 0 );
                if ( no_bl._cutoff_rb.isEnabled() ) {
                    System.out.println( "  [RepresentativeTipsToolTest] cutoff radio should be disabled without branch lengths" );
                    ok[ 0 ] = false;
                }
                if ( !no_bl._target_rb.isSelected() ) {
                    System.out.println( "  [RepresentativeTipsToolTest] target radio should be preselected without branch lengths" );
                    ok[ 0 ] = false;
                }
                if ( no_bl._protect_cb.isEnabled() ) {
                    System.out.println( "  [RepresentativeTipsToolTest] protect checkbox should be disabled with no selection" );
                    ok[ 0 ] = false;
                }
                final MainFrameApplication.RepSelInputs with_bl = MainFrameApplication.buildRepSelInputs( true, true,
                                                                                                          0.05, 100, 0, 0 );
                if ( !with_bl._cutoff_rb.isEnabled() || !with_bl._cutoff_rb.isSelected() ) {
                    System.out.println( "  [RepresentativeTipsToolTest] cutoff radio should be enabled+selected with branch lengths" );
                    ok[ 0 ] = false;
                }
                final MainFrameApplication.RepSelInputs with_sel = MainFrameApplication.buildRepSelInputs( true, true,
                                                                                                           0.05, 100, 0, 3 );
                if ( !with_sel._protect_cb.isEnabled() || !with_sel._protect_cb.isSelected() ) {
                    System.out.println( "  [RepresentativeTipsToolTest] protect checkbox should be enabled+on with a selection" );
                    ok[ 0 ] = false;
                }
                // The panel must be laid out compactly: a vertical BoxLayout so each row keeps its own height,
                // NOT a GridLayout(0,1) that inflates every row (labels/spacers included) to the tallest child
                // (the combo box) and made the dialog ~100px taller than needed -- clipping OK/Cancel off the
                // bottom on smaller/scaled displays. Assert the panel is shorter than that grid equivalent.
                final Component[] rows = with_sel._panel.getComponents();
                int row_sum = 0;
                int tallest = 0;
                for ( final Component row : rows ) {
                    final int h = row.getPreferredSize().height;
                    row_sum += h;
                    tallest = Math.max( tallest, h );
                }
                final int grid_equivalent = rows.length * tallest; // what GridLayout(0,1) would force
                if ( with_sel._panel.getPreferredSize().height >= grid_equivalent ) {
                    System.out.println( "  [RepresentativeTipsToolTest] input panel is not compact (rows all forced "
                            + "to the tallest height, like the old GridLayout)" );
                    ok[ 0 ] = false;
                }
                // fitToScreen: a form that already fits is returned unchanged (no scrollbar on normal displays).
                if ( MainFrameApplication.fitToScreen( with_sel._panel, 2000 ) != with_sel._panel ) {
                    System.out.println( "  [RepresentativeTipsToolTest] fitToScreen wrapped a form that already fits" );
                    ok[ 0 ] = false;
                }
                // ...but a form taller than the given max is scroll-capped to it, so the surrounding dialog's
                // OK/Cancel row can never be pushed off-screen. Use a deterministically over-tall stand-in.
                final javax.swing.JPanel tall = new javax.swing.JPanel();
                tall.setPreferredSize( new Dimension( 300, 1000 ) );
                final Component capped = MainFrameApplication.fitToScreen( tall, 500 );
                if ( capped == tall ) {
                    System.out.println( "  [RepresentativeTipsToolTest] fitToScreen did not cap an over-tall form" );
                    ok[ 0 ] = false;
                }
                else if ( !( capped instanceof JScrollPane ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] fitToScreen should wrap in a JScrollPane" );
                    ok[ 0 ] = false;
                }
                else if ( capped.getPreferredSize().height > 500 ) {
                    System.out.println( "  [RepresentativeTipsToolTest] fitToScreen did not cap the form to the given "
                            + "max height (OK/Cancel could still be clipped)" );
                    ok[ 0 ] = false;
                }
                // clampFullyOnScreen pulls a window that would spill its bottom (under the Dock) back inside the
                // usable area -- the real remaining cause of the clipped OK button -- without pushing it off top.
                final Rectangle usable = new Rectangle( 0, 25, 1440, 850 );
                final Point pulled_up = MainFrameApplication.clampFullyOnScreen( new Point( 200, 800 ),
                                                                                 new Dimension( 400, 300 ), usable );
                if ( ( ( pulled_up.y + 300 ) > ( usable.y + usable.height ) ) || ( pulled_up.y < usable.y )
                        || ( pulled_up.x < usable.x ) || ( ( pulled_up.x + 400 ) > ( usable.x + usable.width ) ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] clampFullyOnScreen left the window outside the "
                            + "usable area (OK/Cancel could still be clipped): " + pulled_up );
                    ok[ 0 ] = false;
                }
                final Point already_ok = MainFrameApplication.clampFullyOnScreen( new Point( 100, 100 ),
                                                                                   new Dimension( 400, 300 ), usable );
                if ( ( already_ok.x != 100 ) || ( already_ok.y != 100 ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] clampFullyOnScreen moved a window that already "
                            + "fits: " + already_ok );
                    ok[ 0 ] = false;
                }
                // extracted-tree name: "<parent>_<N>reps", singular "_1rep", "tree" fallback when unnamed
                if ( !"flaviviridae_247reps".equals( MainFrameApplication.representativeTreeName( "flaviviridae", 247 ) )
                        || !"flaviviridae_1rep".equals( MainFrameApplication.representativeTreeName( "flaviviridae", 1 ) )
                        || !"tree_3reps".equals( MainFrameApplication.representativeTreeName( "", 3 ) )
                        || !"tree_3reps".equals( MainFrameApplication.representativeTreeName( null, 3 ) ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] representativeTreeName format is wrong: "
                            + MainFrameApplication.representativeTreeName( "flaviviridae", 247 ) );
                    ok[ 0 ] = false;
                }
                // parent name comes from the DISPLAYED tab title ("mammals") when the internal name is empty, so
                // the base is preserved: mammals -> mammals_233reps; the "[N]" placeholder tab is ignored
                if ( !"mammals".equals( MainFrameApplication.resolveParentTreeName( "", "mammals" ) )                // tab title used when phy unnamed
                        || !"mammals".equals( MainFrameApplication.resolveParentTreeName( "mammals", "mammals" ) )   // both agree
                        || !"the_internal_name".equals( MainFrameApplication.resolveParentTreeName( "the_internal_name", null ) )
                        || ( MainFrameApplication.resolveParentTreeName( null, "[2]" ) != null )                     // placeholder tab ignored
                        || ( MainFrameApplication.resolveParentTreeName( "", "" ) != null ) ) {                      // nothing usable
                    System.out.println( "  [RepresentativeTipsToolTest] resolveParentTreeName picked the wrong parent name" );
                    ok[ 0 ] = false;
                }
                // end to end: unnamed phylogeny + "mammals" tab + 233 reps -> mammals_233reps
                if ( !"mammals_233reps".equals( MainFrameApplication.representativeTreeName(
                        MainFrameApplication.resolveParentTreeName( "", "mammals" ), 233 ) ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] parent-name + count did not compose to mammals_233reps" );
                    ok[ 0 ] = false;
                }
                // a file extension (dot + 1-5 chars) is stripped before "_Nreps"; longer/absent suffixes are kept
                if ( !"mammals".equals( MainFrameApplication.stripShortExtension( "mammals.xml" ) )                  // 3-char ext
                        || !"mammals".equals( MainFrameApplication.stripShortExtension( "mammals.nexus" ) )          // 5-char ext
                        || !"mammals".equals( MainFrameApplication.stripShortExtension( "mammals.h" ) )              // 1-char ext
                        || !"tree.v2".equals( MainFrameApplication.stripShortExtension( "tree.v2.xml" ) )            // only last suffix
                        || !"mammals".equals( MainFrameApplication.stripShortExtension( "mammals" ) )                // no suffix -> unchanged
                        || !"data.superlong".equals( MainFrameApplication.stripShortExtension( "data.superlong" ) ) // 9-char ext kept
                        || ( MainFrameApplication.stripShortExtension( null ) != null ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] stripShortExtension is wrong: "
                            + MainFrameApplication.stripShortExtension( "mammals.xml" ) );
                    ok[ 0 ] = false;
                }
                // and the whole name: "mammals.xml" -> "mammals_233reps"
                if ( !"mammals_233reps".equals( MainFrameApplication.representativeTreeName( "mammals.xml", 233 ) ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] representativeTreeName did not strip the .xml extension" );
                    ok[ 0 ] = false;
                }
                // provenance description written to the extracted tree: mode + value + pick + counts + parent name
                final String d_cut = MainFrameApplication.representativeTreeDescription( true, 0.05, 0,
                        RepresentativePick.MEDOID, 233, "mammals", 1000 );
                if ( !d_cut.equals( "Used the distance-cutoff (maximum distance 0.05, medoid representative) algorithm "
                        + "to select 233 representative tips from tree named \"mammals\" with 1000 tips." ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] cutoff-mode description is wrong: " + d_cut );
                    ok[ 0 ] = false;
                }
                final String d_tgt = MainFrameApplication.representativeTreeDescription( false, 0.0, 50,
                        RepresentativePick.LONGEST_BRANCH, 1, "mammals", 1 );
                if ( !d_tgt.equals( "Used the target-count (target 50, longest-branch representative) algorithm to "
                        + "select 1 representative tip from tree named \"mammals\" with 1 tip." ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] target-mode/singular description is wrong: " + d_tgt );
                    ok[ 0 ] = false;
                }
                // no parent name -> "tree"
                if ( !MainFrameApplication.representativeTreeDescription( true, 0.1, 0, RepresentativePick.MEDOID, 5,
                        null, 20 ).contains( "from tree named \"tree\" with 20 tips." ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] unnamed-parent description did not fall back to \"tree\"" );
                    ok[ 0 ] = false;
                }
            } );
            SwingUtilities.invokeAndWait( () -> {
                final MainPanel mp = mf[ 0 ].getMainPanel();
                final TreePanel tp = mp.getCurrentTreePanel();
                final Phylogeny phy = tp.getPhylogeny();

                // menu item present + tooltip
                final JMenuItem item = toolsItem( mf[ 0 ].getJMenuBar(), "Select Representative Tips" );
                if ( item == null ) {
                    System.out.println( "  [RepresentativeTipsToolTest] Tools menu item not found" );
                    ok[ 0 ] = false;
                }
                else if ( ( item.getToolTipText() == null )
                        || !item.getToolTipText().toLowerCase().contains( "representative" ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] menu item missing tooltip" );
                    ok[ 0 ] = false;
                }
                // pruning items are labelled "...Tips" (renamed from "...Nodes"); guard against a regression
                if ( ( toolsItem( mf[ 0 ].getJMenuBar(), "Delete Selected Tips" ) == null )
                        || ( toolsItem( mf[ 0 ].getJMenuBar(), "Retain Selected Tips" ) == null ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] \"Delete/Retain Selected Tips\" menu items not found" );
                    ok[ 0 ] = false;
                }
                if ( toolsItem( mf[ 0 ].getJMenuBar(), "Selected Nodes" ) != null ) {
                    System.out.println( "  [RepresentativeTipsToolTest] old \"...Selected Nodes\" menu label still present" );
                    ok[ 0 ] = false;
                }

                // protection via the real found-node path: select tip B (otherwise dropped in favor of A in
                // its cherry), resolve the selection to protected tips, and confirm B survives
                final PhylogenyNode b = tipByName( phy, "B" );
                final Set<Long> sel = new HashSet<>();
                sel.add( b.getId() );
                tp.setFoundNodes0( sel );
                final Set<Long> protected_ids = new HashSet<>();
                for( final PhylogenyNode t : NodeDataExporter.selectedExternalTips( phy,
                        tp.getFoundNodesAsListOfPhylogenyNodes() ) ) {
                    protected_ids.add( t.getId() );
                }
                final SelectionResult prot = RepresentativeTipSelector.selectByCutoff( phy, 0.02,
                                                                                       RepresentativePick.MEDOID,
                                                                                       protected_ids );
                if ( !prot.representativeIds().contains( b.getId() ) || ( prot.getProtectedKeptCount() != 1 ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] protected tip B was not kept via the found-node path" );
                    ok[ 0 ] = false;
                }
                tp.setFoundNodes0( null );

                // compute representatives (two cherries at cutoff 0.02 -> 2 groups)
                final SelectionResult res = RepresentativeTipSelector.selectByCutoff( phy, 0.02,
                                                                                      RepresentativePick.MEDOID );
                if ( res.getClusterCount() != 2 ) {
                    System.out.println( "  [RepresentativeTipsToolTest] expected 2 groups, got "
                            + res.getClusterCount() );
                    ok[ 0 ] = false;
                }

                // highlight via the search-hit mechanism, exactly as the handler does
                tp.setFoundNodes0( new HashSet<>( res.representativeIds() ) );
                final Set<Long> found = tp.getFoundNodes0();
                if ( ( found == null ) || !found.equals( res.representativeIds() ) || ( found.size() != 2 ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] representatives not highlighted" );
                    ok[ 0 ] = false;
                }

                // extract into a new tab: recompute on a copy, prune the non-representatives, add the tab
                final int tabs_before = mp.getTabbedPane().getTabCount();
                final Phylogeny copy = phy.copy();
                final SelectionResult on_copy = RepresentativeTipSelector.selectByCutoff( copy, 0.02,
                                                                                          RepresentativePick.MEDOID );
                final Set<Long> keep = on_copy.representativeIds();
                final Set<Long> to_delete = new HashSet<>();
                for( final PhylogenyNode ext : copy.getExternalNodes() ) {
                    if ( !keep.contains( ext.getId() ) ) {
                        to_delete.add( ext.getId() );
                    }
                }
                PhylogenyMethods.deleteExternalNodesNegativeSelection( to_delete, copy );
                mp.addPhylogenyInNewTab( copy, conf, "reps (representatives)", null );

                if ( mp.getTabbedPane().getTabCount() != ( tabs_before + 1 ) ) {
                    System.out.println( "  [RepresentativeTipsToolTest] a new tab was not added" );
                    ok[ 0 ] = false;
                }
                final Phylogeny shown = mp.getCurrentPhylogeny();
                if ( shown.getNumberOfExternalNodes() != 2 ) {
                    System.out.println( "  [RepresentativeTipsToolTest] extracted tree should have 2 tips, has "
                            + shown.getNumberOfExternalNodes() );
                    ok[ 0 ] = false;
                }
                // the original tree is untouched (non-destructive)
                if ( phy.getNumberOfExternalNodes() != 4 ) {
                    System.out.println( "  [RepresentativeTipsToolTest] original tree was mutated" );
                    ok[ 0 ] = false;
                }

                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static JMenuItem toolsItem( final JMenuBar bar, final String contains ) {
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m != null ) && "Tools".equals( m.getText() ) ) {
                for( int j = 0; j < m.getItemCount(); ++j ) {
                    final JMenuItem it = m.getItem( j );
                    if ( ( it != null ) && ( it.getText() != null ) && it.getText().contains( contains ) ) {
                        return it;
                    }
                }
            }
        }
        return null;
    }

    private static PhylogenyNode tipByName( final Phylogeny phy, final String name ) {
        for( final PhylogenyNode n : phy.getExternalNodes() ) {
            if ( name.equals( n.getName() ) ) {
                return n;
            }
        }
        return null;
    }

    private static Phylogeny twoCherries() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode x = child( root, null, 0.4 );
        final PhylogenyNode y = child( root, null, 0.4 );
        child( x, "A", 0.01 );
        child( x, "B", 0.01 );
        child( y, "C", 0.01 );
        child( y, "D", 0.01 );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode child( final PhylogenyNode parent, final String name, final double dist ) {
        final PhylogenyNode n = new PhylogenyNode();
        if ( name != null ) {
            n.setName( name );
        }
        parent.addAsChild( n );
        n.setDistanceToParent( dist );
        return n;
    }
}

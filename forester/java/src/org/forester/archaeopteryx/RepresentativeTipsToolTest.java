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
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
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

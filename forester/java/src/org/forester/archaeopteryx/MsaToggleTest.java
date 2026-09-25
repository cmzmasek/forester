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
import java.awt.image.BufferedImage;
import java.util.Set;

import javax.swing.JCheckBox;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;

/**
 * Tests the "Sequence Alignment" checkbox on the control panel -- the user-facing on/off for the alignment drawn
 * beside the tips, which before this existed only in the Settings dialog (so a tree that opened with its alignment
 * shown was awkward to quiet down).
 * <p>
 * What is checked: the checkbox is OFFERED exactly when the loaded tree carries an alignment the renderer could
 * actually draw (and stays hidden otherwise, including for the near misses that look like an alignment but are not);
 * a real click on it reaches the tab and changes what is drawn; each TAB keeps its own answer across a tab round
 * trip; a tree whose tips arrive aligned opens with the box already ticked; and Reset to Defaults puts both the tab
 * and the shared widget back to off.
 * <p>
 * The data-presence half needs no display and runs everywhere; the GUI half is a no-op on a headless box.
 */
public final class MsaToggleTest {

    private static final String[] ALIGN  = { "MKQLEDPFGH-WYVAST", "MKQIEDPFGY-WYVAST", "LRQMEDANGH-WFVCST" };
    /** A second, clearly different alignment, so the two tabs cannot pass by accident of looking alike. */
    private static final String[] ALIGN2 = { "CCCCWWWWYYYY-HHHH", "CCCCWWWWFFFF-HHHH", "CCCCWWWWYYYY-HHHR" };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MsaToggle: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };
        // ---- (1) WHEN the checkbox is offered: the data-presence scan behind the row's visibility ----------------
        // The pair the rule needs: the tree that must offer it, and the near misses that must not.
        if ( !presence( alignedTree( ALIGN ) ) ) {
            fail( ok, "a tree whose tips carry an aligned sequence must offer the Sequence Alignment checkbox" );
        }
        if ( presence( plainTree() ) ) {
            fail( ok, "a tree with no sequences at all must NOT offer the Sequence Alignment checkbox" );
        }
        // an UNALIGNED molecular sequence is not an alignment -- alignmentLength() ignores it, so nothing would draw
        if ( presence( sequenceTree( "MKQLEDPFGH", false ) ) ) {
            fail( ok, "an UNALIGNED molecular sequence must not offer the checkbox -- nothing would be drawn" );
        }
        // aligned, but empty: hasAlignedSequences/alignmentLength both refuse it, so the checkbox must too
        if ( presence( sequenceTree( "", true ) ) ) {
            fail( ok, "an EMPTY aligned sequence must not offer the checkbox" );
        }
        // aligned, but only on an INTERNAL node: the renderer builds the alignment from the EXTERNAL nodes, so this
        // draws nothing. This is the case a scan over all nodes (the scan IS pre-order over all nodes) gets wrong.
        if ( presence( internalOnlyAlignedTree() ) ) {
            fail( ok, "an aligned sequence on an INTERNAL node only must not offer the checkbox -- the renderer "
                    + "reads the external nodes, so there is nothing to draw" );
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return ok[ 0 ]; // the rest drives a real frame
        }
        try {
            guiChecks( ok );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    private static void guiChecks( final boolean[] ok ) throws Exception {
        final Phylogeny a = alignedTree( ALIGN );
        final Phylogeny b = alignedTree( ALIGN2 );
        final MainFrame[] mf = new MainFrame[ 1 ];
        // two tabs, both carrying an alignment; createInstance selects the LAST one
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { a, b }, new Configuration(), "msatoggle" ) );
        final MainPanel mp = mf[ 0 ].getMainPanel();
        final ControlPanel cp = mp.getControlPanel();
        final JTabbedPane tabs = mp.getTabbedPane();
        SwingUtilities.invokeAndWait( () -> {
            mf[ 0 ].getOptions().setShowOverview( false );
            mf[ 0 ].getOptions().setMsaColumnWidth( 10 );
        } );
        if ( tabs.getTabCount() != 2 ) {
            fail( ok, "expected two tabs, got " + tabs.getTabCount() );
            return;
        }
        final TreePanel tp0 = mp.getTreePanels().get( 0 );
        final TreePanel tp1 = mp.getTreePanels().get( 1 );

        // ---- (2) the checkbox exists on the panel, under the label the Settings dialog uses ----------------------
        final JCheckBox cb = cp.checkboxForTest( DisplayOption.SHOW_MSA );
        if ( cb == null ) {
            fail( ok, "there is no Sequence Alignment checkbox on the control panel" );
            return;
        }
        if ( !"Sequence Alignment".equals( cb.getText() ) ) {
            fail( ok, "the checkbox label must match the Settings dialog's wording, got \"" + cb.getText() + "\"" );
        }
        if ( !cp.isCheckboxRowVisibleForTest( DisplayOption.SHOW_MSA ) ) {
            fail( ok, "the row must be SHOWN for a tree that carries an alignment" );
        }

        // ---- (3) a tree that arrives aligned opens with the alignment already on --------------------------------
        if ( !cb.isSelected() || !tp1.isShowMsa() ) {
            fail( ok, "a tree whose tips carry an alignment must open with the box ticked and the tab showing it" );
        }

        // ---- (4) a REAL click switches the drawing off, and back on ---------------------------------------------
        // Drives the widget + the listener, not the setter, so an inert checkbox cannot pass.
        final BufferedImage[] img = new BufferedImage[ 3 ];
        SwingUtilities.invokeAndWait( () -> {
            tp1.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp1.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
            tp1.setSize( 760, 460 );
            img[ 0 ] = render( tp1 );
            cb.doClick(); // -> off
        } );
        if ( cb.isSelected() || tp1.isShowMsa() ) {
            fail( ok, "clicking the checkbox off must reach the tab -- the checkbox must not be inert" );
        }
        SwingUtilities.invokeAndWait( () -> img[ 1 ] = render( tp1 ) );
        if ( !differ( img[ 0 ], img[ 1 ] ) ) {
            fail( ok, "switching the alignment off must change what is drawn" );
        }
        SwingUtilities.invokeAndWait( () -> {
            cb.doClick(); // -> back on
            img[ 2 ] = render( tp1 );
        } );
        if ( !cb.isSelected() || !tp1.isShowMsa() ) {
            fail( ok, "clicking the checkbox back on must reach the tab" );
        }
        if ( !differ( img[ 1 ], img[ 2 ] ) ) {
            fail( ok, "switching the alignment back on must change what is drawn" );
        }

        // ---- (5) each TAB keeps its own answer, and a tab round trip re-seeds the shared widget -----------------
        // tab 1 off, tab 0 left on: neither may follow the other.
        SwingUtilities.invokeAndWait( () -> cb.doClick() ); // current tab (1) -> off
        if ( tp1.isShowMsa() ) {
            fail( ok, "tab 1 should be off at this point" );
        }
        if ( !tp0.isShowMsa() ) {
            fail( ok, "switching the alignment off in tab 1 must NOT switch it off in tab 0" );
        }
        SwingUtilities.invokeAndWait( () -> tabs.setSelectedIndex( 0 ) );
        if ( !cp.isShowMsa() ) {
            fail( ok, "arriving on a tab that shows its alignment must re-seed the checkbox to ON" );
        }
        SwingUtilities.invokeAndWait( () -> tabs.setSelectedIndex( 1 ) );
        if ( cp.isShowMsa() ) {
            fail( ok, "arriving back on the tab that hides its alignment must re-seed the checkbox to OFF" );
        }
        if ( tp1.isShowMsa() || !tp0.isShowMsa() ) {
            fail( ok, "a tab round trip must not change either tab's own answer" );
        }

        // ---- (6) the row is HIDDEN for a tree with no alignment, even with aligned trees open elsewhere ---------
        // File -> New opens an empty tab: the checkbox must collapse away rather than sit there doing nothing.
        if ( mf[ 0 ]._new_item == null ) {
            fail( ok, "could not reach File -> New" );
        }
        else {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ]._new_item.doClick() );
            SwingUtilities.invokeAndWait( () -> cp.updateDataCheckboxVisibility( true ) );
            if ( cp.isCheckboxRowVisibleForTest( DisplayOption.SHOW_MSA ) ) {
                fail( ok, "the row must be HIDDEN on a tab whose tree carries no alignment" );
            }
            SwingUtilities.invokeAndWait( () -> tabs.setSelectedIndex( 1 ) );
            SwingUtilities.invokeAndWait( () -> cp.updateDataCheckboxVisibility( true ) );
            if ( !cp.isCheckboxRowVisibleForTest( DisplayOption.SHOW_MSA ) ) {
                fail( ok, "...and come back when a tab that HAS an alignment is current again" );
            }
        }

        // ---- (7) Reset to Defaults: every tab back to off, and the shared widget must say so --------------------
        // A stale "on" in the widget is not cosmetic -- the next push writes the widget onto the tab.
        SwingUtilities.invokeAndWait( () -> {
            tp0.setShowMsa( true );
            tp1.setShowMsa( true );
            mf[ 0 ].resetToDefaults();
        } );
        if ( tp0.isShowMsa() || tp1.isShowMsa() ) {
            fail( ok, "Reset to Defaults must turn the alignment off in EVERY tab" );
        }
        if ( cp.isShowMsa() ) {
            fail( ok, "Reset to Defaults must re-seed the checkbox, not just the per-tab flags" );
        }
        SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
    }

    /** Whether the control panel would OFFER the alignment checkbox for this tree. */
    private static boolean presence( final Phylogeny phy ) {
        final Set<DisplayOption> present = AptxUtil.scanForDataPresence( phy );
        return present.contains( DisplayOption.SHOW_MSA );
    }

    private static BufferedImage render( final TreePanel tp ) {
        return AptxUtil.renderPhylogenyToImage( 760, 460, tp, tp.getOptions(), false, 1, false );
    }

    private static boolean differ( final BufferedImage a, final BufferedImage b ) {
        if ( ( a == null ) || ( b == null ) || ( a.getWidth() != b.getWidth() )
                || ( a.getHeight() != b.getHeight() ) ) {
            return true;
        }
        for( int x = 0; x < a.getWidth(); ++x ) {
            for( int y = 0; y < a.getHeight(); ++y ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    private static Phylogeny alignedTree( final String[] align ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < align.length; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( "seq_" + ( i + 1 ) );
            tip.setDistanceToParent( 0.1 + ( i * 0.05 ) );
            tip.getNodeData().addSequence( aligned( align[ i ], true ) );
            root.addAsChild( tip );
        }
        return wrap( root );
    }

    /** Tips with names and branch lengths but no sequence data of any kind. */
    private static Phylogeny plainTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( final String name : new String[] { "a", "b", "c" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.2 );
            root.addAsChild( tip );
        }
        return wrap( root );
    }

    /** Tips carrying a molecular sequence of the given text and aligned-ness. */
    private static Phylogeny sequenceTree( final String seq, final boolean is_aligned ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( final String name : new String[] { "a", "b", "c" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.2 );
            tip.getNodeData().addSequence( aligned( seq, is_aligned ) );
            root.addAsChild( tip );
        }
        return wrap( root );
    }

    /** An aligned sequence hung on the INTERNAL node, with the tips carrying none. */
    private static Phylogeny internalOnlyAlignedTree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode inner = new PhylogenyNode();
        inner.setName( "ancestor" );
        inner.setDistanceToParent( 0.1 );
        inner.getNodeData().addSequence( aligned( ALIGN[ 0 ], true ) );
        for( final String name : new String[] { "a", "b" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.2 );
            inner.addAsChild( tip );
        }
        root.addAsChild( inner );
        final PhylogenyNode out = new PhylogenyNode();
        out.setName( "out" );
        out.setDistanceToParent( 0.4 );
        root.addAsChild( out );
        return wrap( root );
    }

    private static Sequence aligned( final String seq, final boolean is_aligned ) {
        final Sequence s = new Sequence();
        s.setMolecularSequence( seq );
        s.setMolecularSequenceAligned( is_aligned );
        return s;
    }

    private static Phylogeny wrap( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [MsaToggleTest] " + message );
        ok[ 0 ] = false;
    }

    private MsaToggleTest() {
        // not instantiable
    }
}

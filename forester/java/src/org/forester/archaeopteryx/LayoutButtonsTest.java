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

import javax.swing.JFrame;
import javax.swing.JToggleButton;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.archaeopteryx.Options.TREE_ORIENTATION;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The control panel's layout row replaced two Settings-dialog dropdowns (tree style + orientation) with one
 * exclusive group of five buttons -- the five primary display types. This drives the REAL buttons and asserts:
 * <ul>
 * <li>each button lands the live tree on the right graphics type AND the right orientation, including the two
 * radial layouts, where orientation is a no-op and must be left alone;</li>
 * <li>the rectangular SUB-STYLE chosen in Settings (Euro / Rounded / Triangular) survives a round trip out to a
 * radial layout and back -- returning to "rectangular" must not silently reset the tree to Square;</li>
 * <li>a style change made from the Settings side re-seeds the row, so the two never disagree, and a style picked
 * while a RADIAL layout is showing is remembered rather than yanking the user back to rectangular;</li>
 * <li>"A" (aligned phylogram) is disabled in unrooted -- there is no column or ring to pin labels to there -- and
 * an active "A" is shown as "P" rather than sitting selected and inert, but comes BACK on the way out (the
 * override is on the button only; the tab's stored choice is the memory);</li>
 * <li>the row is re-seeded on Reset to Defaults, and on File -&gt; New, which forces a rectangular tree without
 * going through {@code typeChanged} (a stale always-visible control is exactly the trap the theme radios used to
 * fall into);</li>
 * <li>a tab switch cannot let the outgoing tab's P/A/C state overwrite the incoming tab's own stored choice.</li>
 * </ul>
 * Guarded to a no-op on a headless box (it needs a real {@code MainFrame}).
 */
public final class LayoutButtonsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LayoutButtons: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        return run( new Phylogeny[] {}, LayoutButtonsTest::noTreeStillSwitches )
                && run( new Phylogeny[] { tree() }, LayoutButtonsTest::exercise )
                && run( new Phylogeny[] { tree(), tree() }, LayoutButtonsTest::perTabDisplayTypeSurvivesTabSwitch )
                && run( new Phylogeny[] { tree() }, LayoutButtonsTest::newTreeResyncsLayoutRow );
    }

    /** Opens a frame on the given trees, runs {@code body} on the EDT, and always disposes the frame. */
    private static boolean run( final Phylogeny[] trees, final java.util.function.BiConsumer<MainFrame, boolean[]> body ) {
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( trees, new Configuration(),
                                                                        "layout-buttons" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    body.accept( mf[ 0 ], ok );
                }
                finally {
                    ( ( JFrame ) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * A tab switch has to re-seed the incoming tab's stored P/A/C choice BEFORE working out what that layout
     * allows, not after. An unrooted tab whose stored choice is "A" is the case that catches the wrong order:
     * re-seeding writes "A" back onto the buttons, so unless the "A is dead in unrooted" pass runs afterwards the
     * tab lands with "A" both DISABLED and SELECTED -- and since the paint reads the buttons, the tree would then
     * draw as an aligned phylogram in the one layout that cannot align. Also checks the plain per-tab property:
     * two tabs with different choices keep them across switches, and the stored "A" is still there to come back.
     */
    /**
     * The layout buttons must visibly work in an EMPTY window too. MainFrame.typeChanged guards its whole apply
     * chain -- including the zoom-row relabel -- behind a current-tree check, so without an explicit fallback a
     * user clicking "circular" before opening a tree saw the layout row light up and nothing else change:
     * X+/X- kept their zoom labels instead of becoming the rotate pair.
     */
    private static void noTreeStillSwitches( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 900, 600 );
        ( ( JFrame ) frame ).validate();
        final ControlPanel cp = frame.getMainPanel().getControlPanel();
        cp.getLayoutButton( LayoutIcon.Kind.CIRCULAR ).doClick();
        if ( cp.selectedLayoutKind() != LayoutIcon.Kind.CIRCULAR ) {
            fail( ok, "no tree: the layout row must still track the click, got " + cp.selectedLayoutKind() );
        }
        javax.swing.Icon icon = cp.getZoomInXButtonForTest().getIcon();
        if ( !( icon instanceof ControlButtonIcon )
                || ( ( ( ControlButtonIcon ) icon ).getKind() != ControlButtonIcon.Kind.ROTATE_CW ) ) {
            fail( ok, "no tree: choosing circular must still turn X+ into the rotate button" );
        }
        if ( cp.getFitButtonIconKind() == ControlButtonIcon.Kind.FIT_WIDTH ) {
            fail( ok, "no tree: choosing circular must switch the fit button to the label-direction flip" );
        }
        cp.getLayoutButton( LayoutIcon.Kind.ROOT_LEFT ).doClick();
        icon = cp.getZoomInXButtonForTest().getIcon();
        if ( icon != null ) {
            fail( ok, "no tree: returning to rectangular must restore the X+ text button" );
        }
        if ( cp.getFitButtonIconKind() != ControlButtonIcon.Kind.FIT_WIDTH ) {
            fail( ok, "no tree: returning to rectangular must restore the fit-width glyph, got "
                    + cp.getFitButtonIconKind() );
        }
    }

    private static void perTabDisplayTypeSurvivesTabSwitch( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 1000, 700 );
        ( ( JFrame ) frame ).validate();
        final ControlPanel cp = frame.getMainPanel().getControlPanel();
        final javax.swing.JTabbedPane tabs = frame.getMainPanel().getTabbedPane();
        if ( tabs.getTabCount() != 2 ) {
            fail( ok, "precondition: this exercise needs two tabs, got " + tabs.getTabCount() );
            return;
        }
        // tab 1: choose "A" while rectangular (so the tab's STORED type is ALIGNED), then send it unrooted --
        // the button drops to "P" while the store keeps ALIGNED, ready for the trip back out
        tabs.setSelectedIndex( 1 );
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        cp.getDisplayAsAlignedPhylogramRb().doClick();
        click( cp, LayoutIcon.Kind.UNROOTED );
        if ( cp.getDisplayAsAlignedPhylogramRb().isSelected() ) {
            fail( ok, "precondition: tab 1 should be showing 'P' while unrooted" );
        }
        // tab 0: a rectangular cladogram, so the two tabs disagree about everything
        tabs.setSelectedIndex( 0 );
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        cp.getDisplayAsCladogramRb().doClick();
        if ( cp.getTreeDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM ) {
            fail( ok, "precondition: tab 0 should be a cladogram, got " + cp.getTreeDisplayType() );
        }
        // Arriving back at the unrooted tab, re-seeding its stored ALIGNED must NOT leave "A" lit: it is disabled
        // in unrooted, and a disabled-but-selected button would also make the tree paint as an aligned phylogram
        // in a layout that cannot align. (This is what breaks if the re-seed runs AFTER the enable pass.)
        tabs.setSelectedIndex( 1 );
        if ( cp.getDisplayAsAlignedPhylogramRb().isEnabled() ) {
            fail( ok, "'A' must be disabled on returning to the unrooted tab" );
        }
        if ( cp.getDisplayAsAlignedPhylogramRb().isSelected() ) {
            fail( ok, "'A' must not be left SELECTED on an unrooted tab -- it is disabled there" );
        }
        if ( cp.getTreeDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM ) {
            fail( ok, "the unrooted tab must read back as a plain phylogram, got " + cp.getTreeDisplayType() );
        }
        // ... while the tab it came from keeps its own choice
        tabs.setSelectedIndex( 0 );
        if ( cp.getTreeDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM ) {
            fail( ok, "the rectangular tab's cladogram must survive the round trip, got " + cp.getTreeDisplayType() );
        }
        // and tab 1's stored "A" is still intact: bring it back to a rectangular layout and it returns
        tabs.setSelectedIndex( 1 );
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        if ( !cp.getDisplayAsAlignedPhylogramRb().isSelected() ) {
            fail( ok, "tab 1's stored 'A' must come back when it leaves unrooted, got " + cp.getTreeDisplayType() );
        }
    }

    /**
     * File -&gt; New forces the new tree rectangular directly (not through {@code typeChanged}), so it has to
     * re-seed the layout row by hand -- otherwise the row stays lit on whatever radial button was active while
     * the new tree is drawn rectangular.
     */
    private static void newTreeResyncsLayoutRow( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 1000, 700 );
        ( ( JFrame ) frame ).validate();
        final ControlPanel cp = frame.getMainPanel().getControlPanel();
        click( cp, LayoutIcon.Kind.CIRCULAR );
        if ( cp.selectedLayoutKind() != LayoutIcon.Kind.CIRCULAR ) {
            fail( ok, "precondition: the layout row should be on circular before File -> New" );
            return;
        }
        if ( frame._new_item == null ) {
            fail( ok, "could not reach the File -> New menu item" );
            return;
        }
        frame._new_item.doClick();
        final TreePanel fresh = frame.getCurrentTreePanel();
        if ( ( fresh == null ) || ( fresh.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) ) {
            fail( ok, "precondition: File -> New should open a rectangular tree, got "
                    + ( ( fresh == null ) ? "no tree" : fresh.getPhylogenyGraphicsType().toString() ) );
            return;
        }
        if ( cp.selectedLayoutKind() != LayoutIcon.Kind.ROOT_LEFT ) {
            fail( ok, "File -> New must leave the layout row on the rectangular button it actually drew, got "
                    + cp.selectedLayoutKind() );
        }
    }

    private static void exercise( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 1000, 700 );
        ( ( JFrame ) frame ).validate();
        final ControlPanel cp = frame.getMainPanel().getControlPanel();
        final TreePanel tp = frame.getCurrentTreePanel();
        if ( ( cp == null ) || ( tp == null ) ) {
            fail( ok, "could not reach the control panel / tree panel" );
            return;
        }
        // all five buttons exist and exactly one is lit
        for ( final LayoutIcon.Kind kind : LayoutIcon.Kind.values() ) {
            if ( cp.getLayoutButton( kind ) == null ) {
                fail( ok, "the control panel is missing the " + kind + " layout button" );
                return;
            }
        }
        if ( cp.selectedLayoutKind() == null ) {
            fail( ok, "one of the five layout buttons must be lit at startup" );
        }
        // 1. each button drives BOTH halves of the display type: the style and (for the rectangular three) the
        //    orientation. The radial two must leave the orientation exactly as they found it.
        click( cp, LayoutIcon.Kind.ROOT_TOP );
        check( ok, frame, tp, "root-top", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_TOP );
        click( cp, LayoutIcon.Kind.ROOT_BOTTOM );
        check( ok, frame, tp, "root-bottom", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_BOTTOM );
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        check( ok, frame, tp, "root-left", PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR, TREE_ORIENTATION.ROOT_LEFT );
        click( cp, LayoutIcon.Kind.CIRCULAR );
        check( ok, frame, tp, "circular", PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, TREE_ORIENTATION.ROOT_LEFT );
        click( cp, LayoutIcon.Kind.UNROOTED );
        check( ok, frame, tp, "unrooted", PHYLOGENY_GRAPHICS_TYPE.UNROOTED, TREE_ORIENTATION.ROOT_LEFT );
        // an orientation set BEFORE a radial trip must survive it untouched
        click( cp, LayoutIcon.Kind.ROOT_BOTTOM );
        click( cp, LayoutIcon.Kind.CIRCULAR );
        if ( frame.getOptions().getTreeOrientation() != TREE_ORIENTATION.ROOT_BOTTOM ) {
            fail( ok, "circular must not disturb the tree orientation, got "
                    + frame.getOptions().getTreeOrientation() );
        }
        // 2. a rectangular sub-style survives the round trip out to a radial layout and back
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        frame._rounded_type_cbmi.doClick(); // what the Settings "Rectangular style" dropdown does
        if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
            fail( ok, "picking Rounded in Settings should make the tree ROUNDED, got "
                    + tp.getPhylogenyGraphicsType() );
        }
        // 3. ... and that Settings-side change must re-seed the row rather than leaving it lit on the wrong button
        if ( cp.selectedLayoutKind() != LayoutIcon.Kind.ROOT_LEFT ) {
            fail( ok, "a Settings style change must leave the layout row on root-left, got "
                    + cp.selectedLayoutKind() );
        }
        click( cp, LayoutIcon.Kind.CIRCULAR );
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
            fail( ok, "returning to rectangular must restore the Rounded sub-style, got "
                    + tp.getPhylogenyGraphicsType() );
        }
        // 4. "A" (aligned phylogram) needs somewhere to pin the labels -- a right column (rectangular) or the
        //    outer ring (circular). UNROOTED has neither, so the button must go dead there rather than sit
        //    there doing nothing, and an active "A" must fall back to the plain phylogram unrooted really draws.
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        cp.getDisplayAsAlignedPhylogramRb().doClick();
        if ( !cp.getDisplayAsAlignedPhylogramRb().isSelected() ) {
            fail( ok, "precondition: 'A' should be selectable in a rectangular layout" );
        }
        click( cp, LayoutIcon.Kind.UNROOTED );
        if ( cp.getDisplayAsAlignedPhylogramRb().isEnabled() ) {
            fail( ok, "'A' (aligned phylogram) must be disabled in the unrooted layout -- nothing to align to" );
        }
        if ( cp.getDisplayAsAlignedPhylogramRb().isSelected() ) {
            fail( ok, "an active 'A' must fall back when entering unrooted, not stay selected and inert" );
        }
        if ( !cp.getDisplayAsUnalignedPhylogramRb().isSelected() ) {
            fail( ok, "the fallback from 'A' in unrooted must be 'P' (what unrooted actually draws)" );
        }
        // circular CAN align (labels pinned to the outer ring), so it must stay live there
        click( cp, LayoutIcon.Kind.CIRCULAR );
        if ( !cp.getDisplayAsAlignedPhylogramRb().isEnabled() ) {
            fail( ok, "'A' must stay enabled in circular -- labels align on the outer ring there" );
        }
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        if ( !cp.getDisplayAsAlignedPhylogramRb().isEnabled() ) {
            fail( ok, "'A' must come back on leaving unrooted" );
        }
        // 5. the rectangular sub-style can be set from Settings at ANY time. While a radial layout is showing it
        //    must NOT switch the layout out from under the user -- it is remembered and applies on the way back.
        click( cp, LayoutIcon.Kind.CIRCULAR );
        cp.setRectangularStyle( PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE );
        if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            fail( ok, "picking a rectangular style while circular must not leave circular, got "
                    + tp.getPhylogenyGraphicsType() );
        }
        if ( cp.getRectangularStyle() != PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
            fail( ok, "the pending rectangular style should be EURO_STYLE, got " + cp.getRectangularStyle() );
        }
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.EURO_STYLE ) {
            fail( ok, "the style chosen while circular must apply on returning to rectangular, got "
                    + tp.getPhylogenyGraphicsType() );
        }
        // ... and in a rectangular layout the same call applies immediately
        cp.setRectangularStyle( PHYLOGENY_GRAPHICS_TYPE.ROUNDED );
        if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.ROUNDED ) {
            fail( ok, "setting a rectangular style while rectangular must apply at once, got "
                    + tp.getPhylogenyGraphicsType() );
        }
        // 6. Leaving unrooted must bring "A" BACK. The fallback in (4) only overrides the BUTTON; the tab's
        //    stored choice stays ALIGNED, so a deliberate "A" is not silently converted into "P" for good by a
        //    detour through unrooted.
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        cp.getDisplayAsAlignedPhylogramRb().doClick();
        click( cp, LayoutIcon.Kind.UNROOTED );
        if ( cp.getDisplayAsAlignedPhylogramRb().isSelected() ) {
            fail( ok, "precondition: 'A' should be shown as 'P' while unrooted" );
        }
        click( cp, LayoutIcon.Kind.ROOT_LEFT );
        if ( !cp.getDisplayAsAlignedPhylogramRb().isSelected() ) {
            fail( ok, "'A' must be restored on leaving unrooted -- a detour through unrooted must not discard it" );
        }
        // 7. Reset to Defaults re-seeds the row (an always-visible control left stale is the classic trap here)
        click( cp, LayoutIcon.Kind.UNROOTED );
        frame.resetToDefaults();
        if ( cp.selectedLayoutKind() != LayoutIcon.Kind.ROOT_LEFT ) {
            fail( ok, "after Reset to Defaults the layout row must be back on root-left, got "
                    + cp.selectedLayoutKind() );
        }
        if ( tp.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR ) {
            fail( ok, "after Reset to Defaults the live tree must be RECTANGULAR, got "
                    + tp.getPhylogenyGraphicsType() );
        }
    }

    /** Fires the real button the way a user click does (setSelected alone would not run the listener). */
    private static void click( final ControlPanel cp, final LayoutIcon.Kind kind ) {
        final JToggleButton b = cp.getLayoutButton( kind );
        if ( b != null ) {
            b.doClick();
        }
    }

    private static void check( final boolean[] ok, final MainFrame frame, final TreePanel tp, final String what,
                               final PHYLOGENY_GRAPHICS_TYPE type, final TREE_ORIENTATION orientation ) {
        if ( tp.getPhylogenyGraphicsType() != type ) {
            fail( ok, what + " should set the tree style to " + type + ", got " + tp.getPhylogenyGraphicsType() );
        }
        if ( frame.getOptions().getTreeOrientation() != orientation ) {
            fail( ok, what + " should leave the orientation at " + orientation + ", got "
                    + frame.getOptions().getTreeOrientation() );
        }
        final ControlPanel cp = frame.getMainPanel().getControlPanel();
        if ( cp.selectedLayoutKind() != ControlPanel.layoutKindFor( type, orientation ) ) {
            fail( ok, what + " should light the " + ControlPanel.layoutKindFor( type, orientation )
                    + " button, got " + cp.selectedLayoutKind() );
        }
    }

    /** A small tree WITH branch lengths, so the P/A/C row is live and the radial layouts have something to draw. */
    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode inner = new PhylogenyNode();
        inner.setDistanceToParent( 0.2 );
        for ( final String name : new String[] { "a", "b" } ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.3 );
            inner.addAsChild( tip );
        }
        final PhylogenyNode c = new PhylogenyNode();
        c.setName( "c" );
        c.setDistanceToParent( 0.7 );
        root.addAsChild( inner );
        root.addAsChild( c );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        return phy;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [LayoutButtonsTest] " + msg );
        ok[ 0 ] = false;
    }

    private LayoutButtonsTest() {
    }
}

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
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.NODE_LABEL_DIRECTION;
import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The left-panel zoom cluster is relabeled and rewired for the RADIAL (circular/unrooted) layouts (where the two zoom
 * axes collapse to one -- both drive the single radial diameter since 0.11.9): Y+/Y- become a plain +/- zoom, the
 * now-redundant X-/X+ become rotate counter-clockwise / clockwise, "E" (vertical label spacing, a no-op in a fan) is
 * greyed out, and "W" (fit-width, redundant with "F" radially) becomes a node-label-direction flip. Switching back to a
 * rectangular layout reverts every button. Asserts the relabel (via the real MainFrame.typeChanged hook), that the
 * rotate glyphs actually render in the button font, that X-/X+ clicks rotate in the correct direction, that Y+ still
 * zooms (did not become a rotate), and that the "L" button flips the node-label direction (keeping the "Radial Labels"
 * checkbox in sync). Headful; a green no-op when headless.
 */
public final class RadialZoomControlTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RadialZoomControl: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return relabelAndRevertOk() && rotateButtonsOk() && zoomAndLabelOk();
    }

    /** For CIRCULAR and UNROOTED (reached through the real Type-menu path, MainFrame.typeChanged): Y+/Y- -> +/-,
     *  X+/X- -> the CW/CCW rotate glyphs (which must be renderable by the button font -- a tofu guard), E disabled,
     *  W -> "L". Then switching back to RECTANGULAR must restore Y+/Y-, X+/X- (native "+"/"-" or "X+"/"X-"), an
     *  enabled E, and the fit-width glyph. */
    private static boolean relabelAndRevertOk() {
        final boolean[] ok = { true };
        withFrame( "scale-axis.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                frame.typeChanged( gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ? frame._circular_type_cbmi
                        : frame._unrooted_type_cbmi );
                expect( ok, gt + " Y+ -> '+'", "+", cp.getZoomInYButtonForTest().getText() );
                expect( ok, gt + " Y- -> '-'", "-", cp.getZoomOutYButtonForTest().getText() );
                // the rotate buttons are DRAWN glyphs now (no font dependency, so no text fallback to mirror):
                // in a radial layout each must carry its own rotation icon and no leftover text
                final JButton x_in = cp.getZoomInXButtonForTest();
                final JButton x_out = cp.getZoomOutXButtonForTest();
                expectRotateIcon( ok, gt + " X+", x_in, ControlButtonIcon.Kind.ROTATE_CW );
                expectRotateIcon( ok, gt + " X-", x_out, ControlButtonIcon.Kind.ROTATE_CCW );
                // in radial the fit button becomes the label-direction flip and, like the theme toggle, shows
                // the state it will switch TO: labels are horizontal by default, so it offers LABELS_RADIAL
                if ( cp.getFitButtonIconKind() != ControlButtonIcon.Kind.LABELS_RADIAL ) {
                    fail( ok, gt + ": the label-flip button should offer LABELS_RADIAL, got "
                            + cp.getFitButtonIconKind() );
                }
                if ( ( cp.getFitWidthButtonForTest().getText() != null )
                        && !cp.getFitWidthButtonForTest().getText().isEmpty() ) {
                    fail( ok, gt + ": the label-flip button must not keep text beside its glyph" );
                }
                if ( cp.getExpandButtonForTest().isEnabled() ) {
                    fail( ok, "E must be greyed out (disabled) in " + gt + " (it does nothing radially)" );
                }
                // clockwise and counter-clockwise are mirror images -- exactly the pair a reader can most easily
                // be shown the wrong way round, so the two buttons must not carry the SAME glyph
                if ( ( x_in.getIcon() instanceof ControlButtonIcon ) && ( x_out.getIcon() instanceof ControlButtonIcon )
                        && ( ( ( ControlButtonIcon ) x_in.getIcon() ).getKind() == ( ( ControlButtonIcon ) x_out
                                .getIcon() ).getKind() ) ) {
                    fail( ok, "the two rotate buttons must not carry the same glyph" );
                }
            }
            // revert to rectangular (pin ROOT_LEFT so W reverts to "W" deterministically -- a persisted vertical
            // orientation would legitimately show "H", which the OrientationZoomControlTest already covers): every
            // radial re-label must be undone
            tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
            frame.typeChanged( frame._rectangular_type_cbmi );
            expect( ok, "revert Y+ -> 'Y+'", "Y+", cp.getZoomInYButtonForTest().getText() );
            expect( ok, "revert Y- -> 'Y-'", "Y-", cp.getZoomOutYButtonForTest().getText() );
            expect( ok, "revert X+", "X+", cp.getZoomInXButtonForTest().getText() );
            expect( ok, "revert X-", "X-", cp.getZoomOutXButtonForTest().getText() );
            // finding from review: the ICONS must be asserted cleared/restored too, not just the text -- a
            // regression dropping setIcon(null) would leave a stale rotate glyph beside the restored X+/X- text
            // and a text-only assertion would never notice
            if ( cp.getFitButtonIconKind() != ControlButtonIcon.Kind.FIT_WIDTH ) {
                fail( ok, "revert: the fit button must carry FIT_WIDTH again, got " + cp.getFitButtonIconKind() );
            }
            if ( ( cp.getZoomInXButtonForTest().getIcon() != null )
                    || ( cp.getZoomOutXButtonForTest().getIcon() != null ) ) {
                fail( ok, "revert: the rotate glyphs must be cleared when X+/X- get their text back" );
            }
            if ( !cp.getExpandButtonForTest().isEnabled() ) {
                fail( ok, "E must be re-enabled when reverting to a rectangular layout" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** In a radial layout, the X+/X- buttons rotate the fan: X+ CLOCKWISE (the starting angle increases -- screen y is
     *  down, so a larger angle sweeps clockwise), X- COUNTER-CLOCKWISE (decreases, wrapping into [0,2pi)). Driven by a
     *  real doClick gesture through actionPerformed, so it also proves the isRadialLayout dispatch. */
    private static boolean rotateButtonsOk() {
        final boolean[] ok = { true };
        withFrame( "scale-axis.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            // clockwise: from a mid angle (no wrap) the starting angle must INCREASE
            tp.setStartingAngle( Math.PI );
            cp.getZoomInXButtonForTest().doClick();
            if ( tp.getStartingAngle() <= Math.PI ) {
                fail( ok, "X+ must rotate clockwise (increase the starting angle) -- got " + tp.getStartingAngle()
                        + " from " + Math.PI );
            }
            // counter-clockwise: from the same mid angle the starting angle must DECREASE
            tp.setStartingAngle( Math.PI );
            cp.getZoomOutXButtonForTest().doClick();
            if ( tp.getStartingAngle() >= Math.PI ) {
                fail( ok, "X- must rotate counter-clockwise (decrease the starting angle) -- got "
                        + tp.getStartingAngle() + " from " + Math.PI );
            }
            // counter-clockwise past 0 must WRAP into [0,2pi), never go negative
            tp.setStartingAngle( 0.0 );
            tp.rotateRadial( false );
            if ( tp.getStartingAngle() <= 0.0 ) {
                fail( ok, "a counter-clockwise rotate past 0 must wrap positive -- got " + tp.getStartingAngle() );
            }
            // in a RECTANGULAR layout the X buttons must NOT rotate (they zoom): rotateRadial is a no-op there
            tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setStartingAngle( Math.PI );
            tp.rotateRadial( true );
            if ( tp.getStartingAngle() != Math.PI ) {
                fail( ok, "rotateRadial must be a no-op outside a radial layout -- got " + tp.getStartingAngle() );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Two radial re-wirings that must NOT be confused with rotation: the Y+ button still ZOOMS (grows the single
     *  radial diameter), and the "L" (former "W") button FLIPS the node-label direction, keeping the "Radial Labels"
     *  checkbox in sync. */
    private static boolean zoomAndLabelOk() {
        final boolean[] ok = { true };
        withFrame( "scale-axis.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            frame.typeChanged( frame._circular_type_cbmi );
            cp.showWhole();
            final int d0 = tp.radialDiameter();
            cp.getZoomInYButtonForTest().doClick();
            if ( tp.radialDiameter() <= d0 ) {
                fail( ok, "Y+ must still ZOOM IN (grow the radial diameter) in a radial layout -- " + d0 + " -> "
                        + tp.radialDiameter() );
            }
            // label-direction flip via the "L" button
            final NODE_LABEL_DIRECTION dir0 = o.getNodeLabelDirection();
            cp.getFitWidthButtonForTest().doClick();
            final NODE_LABEL_DIRECTION dir1 = o.getNodeLabelDirection();
            if ( dir1 == dir0 ) {
                fail( ok, "the 'L' button must flip the node-label direction (still " + dir0 + ")" );
            }
            if ( frame._label_direction_cbmi.isSelected() != ( dir1 == NODE_LABEL_DIRECTION.RADIAL ) ) {
                fail( ok, "the 'Radial Labels' checkbox must track the flipped direction (" + dir1 + " vs cbmi="
                        + frame._label_direction_cbmi.isSelected() + ")" );
            }
            // ... and the state GLYPH must follow the flip: it shows the direction it will switch TO, so after
            // flipping to radial labels it must offer LABELS_HORIZONTAL (and vice versa)
            final ControlButtonIcon.Kind offered = cp.getFitButtonIconKind();
            final ControlButtonIcon.Kind expected = ( dir1 == NODE_LABEL_DIRECTION.RADIAL )
                    ? ControlButtonIcon.Kind.LABELS_HORIZONTAL : ControlButtonIcon.Kind.LABELS_RADIAL;
            if ( offered != expected ) {
                fail( ok, "after the flip the label button must offer " + expected + ", got " + offered );
            }
            cp.getFitWidthButtonForTest().doClick();
            if ( o.getNodeLabelDirection() != dir0 ) {
                fail( ok, "a second 'L' click must flip the node-label direction back to " + dir0 + " (got "
                        + o.getNodeLabelDirection() + ")" );
            }
            // the Alt+W keyboard shortcut the "L" tooltip advertises must ALSO flip labels in a radial layout
            // (matching the button); before this it ran fit-width, contradicting the tooltip
            final NODE_LABEL_DIRECTION before_key = o.getNodeLabelDirection();
            dispatchAltW( tp );
            if ( o.getNodeLabelDirection() == before_key ) {
                fail( ok, "Alt+W must flip the node-label direction in a radial layout (still " + before_key + ")" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Dispatch Alt+W to the tree panel's key listeners (the app's hand-dispatched shortcut path). */
    private static void dispatchAltW( final TreePanel tp ) {
        final KeyEvent e = new KeyEvent( tp, KeyEvent.KEY_PRESSED, System.currentTimeMillis(),
                InputEvent.ALT_DOWN_MASK, KeyEvent.VK_W, 'W' );
        for ( final KeyListener kl : tp.getKeyListeners() ) {
            kl.keyPressed( e );
        }
    }

    /** The button must carry that rotation glyph, and NO text: a stale label beside an icon is the classic bug. */
    private static void expectRotateIcon( final boolean[] ok, final String what, final JButton b,
                                          final ControlButtonIcon.Kind kind ) {
        if ( !( b.getIcon() instanceof ControlButtonIcon )
                || ( ( ( ControlButtonIcon ) b.getIcon() ).getKind() != kind ) ) {
            fail( ok, what + " must show the " + kind + " glyph, got "
                    + ( ( b.getIcon() == null ) ? "no icon" : b.getIcon().getClass().getSimpleName() ) );
        }
        if ( ( b.getText() != null ) && !b.getText().isEmpty() ) {
            fail( ok, what + " must not keep its text label beside the glyph, got \"" + b.getText() + "\"" );
        }
    }

    private static void expect( final boolean[] ok, final String what, final String expected, final String actual ) {
        if ( !expected.equals( actual ) ) {
            fail( ok, what + ": expected '" + expected + "' but got '" + actual + "'" );
        }
    }

    private interface FrameBody {
        void run( MainFrame frame, TreePanel tp, Options o ) throws Exception;
    }

    private static void withFrame( final String demo, final FrameBody body, final boolean[] ok ) {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
            if ( !file.exists() ) {
                fail( ok, "demo tree missing: " + file.getAbsolutePath() );
                return;
            }
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "radialzoom" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    body.run( frame, frame.getMainPanel().getCurrentTreePanel(), frame.getOptions() );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [RadialZoomControlTest] " + msg );
        ok[ 0 ] = false;
    }

    private RadialZoomControlTest() {
    }
}

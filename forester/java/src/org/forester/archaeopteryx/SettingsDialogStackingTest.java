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
import java.awt.Dialog;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.AbstractButton;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.Timer;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * The Settings dialog never covers a dialog it opened (Christian, 2026-11-05: "The large Settings menu can cover a
 * smaller dialog box which has the focus and the Settings menu stays unresponsive"). The tree-font chooser it opened
 * had NO owner window, so macOS could stack the main window's family -- Settings with it -- above the chooser while
 * the chooser kept the focus. A dialog owned by Settings stays above it.
 * <p>
 * Settings is opened as the user opens it, from its menu. Then:
 * <ol>
 * <li>every button of Settings whose label ends in an ellipsis (the convention for "opens a dialog") is clicked,
 * the buttons iterated from the dialog itself: each opens a modal dialog (the drive took), OWNED by Settings and
 * lying within its bounds (centred over it, not over the main window: Settings is first moved aside, off the main
 * window's centre, where it opens), which is then closed unanswered -- a cancel,
 * so nothing is reset or set. "Clear Cache" is not clicked (it deletes the taxonomy cache when confirmed);</li>
 * <li>ONE Settings dialog: picking the menu again while it is open brings that dialog forward instead of stacking
 * a second copy over it; once closed -- by its Close button or its close box -- it is GONE (disposed, not hidden:
 * a hidden one stayed a window, re-skinned on every theme switch, one more per open/close cycle), and the menu
 * builds a NEW one (seeded from the current tab).</li>
 * </ol>
 * A dialog that never appears fails instead of hanging the suite: the watcher gives up after ~9 s and closes any modal.
 * Needs a display: a no-op when headless.
 */
public final class SettingsDialogStackingTest {

    private static final String ELLIPSIS = "…";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SettingsDialogStacking: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            final Configuration conf = new Configuration();
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { tree() }, conf, "settings stacking" ) );
            boolean ok = true;
            ok &= dialogsBelongToSettingsOk( mf[ 0 ] );
            ok &= oneSettingsDialogOk( mf[ 0 ] );
            return ok;
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
        finally {
            try {
                SwingUtilities.invokeAndWait( () -> {
                    for( final SettingsDialog d : settingsDialogs( mf[ 0 ], false ) ) {
                        d.dispose();
                    }
                    if ( mf[ 0 ] != null ) {
                        ( (JFrame) mf[ 0 ] ).dispose();
                    }
                } );
            }
            catch ( final Exception e ) {
                e.printStackTrace();
            }
        }
    }

    /** (1) Every ellipsis button of Settings opens a modal dialog owned by Settings, within its bounds. */
    private static boolean dialogsBelongToSettingsOk( final MainFrame mf ) throws Exception {
        pickSettingsMenu( mf );
        final List<SettingsDialog> open = settingsDialogs( mf, true );
        if ( !ck( open.size() == 1, "precondition: the Settings menu opened one dialog, found " + open.size() ) ) {
            return false;
        }
        final SettingsDialog settings = open.get( 0 );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            // Settings opens centred on the main window, where a dialog centred on the MAIN window lies too: move it
            // aside, so "over Settings" can only mean centred on Settings
            settings.setLocation( settings.getX() + 400, settings.getY() );
            final Rectangle s = settings.getBounds();
            final Rectangle m = mf.getBounds();
            ok[ 0 ] &= ck( Math.abs( s.getCenterX() - m.getCenterX() ) >= 300, "precondition: Settings sits aside the main "
                    + "window's centre: " + s + " / " + m );
            final List<AbstractButton> buttons = new ArrayList<>();
            collectButtons( settings.getContentPane(), buttons );
            final List<String> clicked = new ArrayList<>();
            for( final AbstractButton b : buttons ) {
                if ( !b.getText().contains( ELLIPSIS ) ) {
                    continue;
                }
                clicked.add( b.getText() );
                final Map<Window, Boolean> before = new IdentityHashMap<>();
                for( final Window w : Window.getWindows() ) {
                    before.put( w, Boolean.TRUE );
                }
                final Dialog[] seen = { null };
                final Rectangle[] seen_bounds = { null };
                final long start = System.currentTimeMillis();
                final Timer watcher = new Timer( 100, null );
                watcher.addActionListener( e -> {
                    for( final Window w : Window.getWindows() ) {
                        if ( ( w instanceof Dialog ) && ( w != settings ) && !before.containsKey( w ) && w.isShowing()
                                && ( (Dialog) w ).isModal() ) {
                            seen[ 0 ] = (Dialog) w;
                            seen_bounds[ 0 ] = w.getBounds();
                            watcher.stop();
                            w.dispose(); // closed unanswered: a cancel
                            return;
                        }
                    }
                    if ( ( System.currentTimeMillis() - start ) > 9000 ) {
                        watcher.stop();
                        closeAnyModal();
                    }
                } );
                watcher.start();
                b.doClick(); // a modal dialog runs its own event loop, in which the watcher fires
                watcher.stop();
                final String what = "\"" + b.getText() + "\"";
                if ( !ck( seen[ 0 ] != null, what + ": precondition: a dialog opened" ) ) {
                    ok[ 0 ] = false;
                    continue;
                }
                ok[ 0 ] &= ck( seen[ 0 ].getOwner() == settings, what + ": the dialog is owned by Settings, not by "
                        + describe( seen[ 0 ].getOwner(), mf ) );
                ok[ 0 ] &= ck( settings.getBounds().contains( seen_bounds[ 0 ] ), what
                        + ": the dialog lies over Settings " + settings.getBounds() + ", at " + seen_bounds[ 0 ] );
            }
            // the ones that exist today: the tree font, the two time-axis calibrations, Reset to Defaults
            ok[ 0 ] &= ck( clicked.size() >= 4, "precondition: the ellipsis buttons were reached, " + clicked );
            ok[ 0 ] &= ck( clicked.contains( "Choose" + ELLIPSIS ), "precondition: the tree font's \"Choose"
                    + ELLIPSIS + "\" was among them, " + clicked );
        } );
        SwingUtilities.invokeAndWait( () -> closeButton( settings ).doClick() );
        return ok[ 0 ];
    }

    /** (2) One Settings dialog: the menu brings an open one forward; a closed one is built anew. */
    private static boolean oneSettingsDialogOk( final MainFrame mf ) throws Exception {
        boolean ok = true;
        pickSettingsMenu( mf );
        final List<SettingsDialog> first = settingsDialogs( mf, true );
        if ( !ck( first.size() == 1, "precondition: the menu opened one Settings dialog, found " + first.size() ) ) {
            return false;
        }
        pickSettingsMenu( mf );
        final List<SettingsDialog> again = settingsDialogs( mf, true );
        ok &= ck( again.size() == 1, "the menu picked again while Settings is open shows ONE Settings dialog, found "
                + again.size() );
        ok &= ck( ( again.size() == 1 ) && ( again.get( 0 ) == first.get( 0 ) ),
                  "the menu picked again brings the open Settings dialog forward" );
        // closed by its Close button: gone, not hidden
        SwingUtilities.invokeAndWait( () -> closeButton( first.get( 0 ) ).doClick() );
        ok &= ck( settingsDialogs( mf, true ).isEmpty(), "precondition: Settings closed" );
        ok &= ck( !first.get( 0 ).isDisplayable(), "a Settings dialog closed by its Close button is disposed, not hidden" );
        pickSettingsMenu( mf );
        final List<SettingsDialog> reopened = settingsDialogs( mf, true );
        ok &= ck( ( reopened.size() == 1 ) && ( reopened.get( 0 ) != first.get( 0 ) ),
                  "a closed Settings dialog reopens as a NEW one (seeded from the current tab)" );
        // closed by the window's own close box: gone too
        if ( reopened.size() == 1 ) {
            final SettingsDialog d = reopened.get( 0 );
            SwingUtilities.invokeAndWait( () -> d.dispatchEvent( new WindowEvent( d, WindowEvent.WINDOW_CLOSING ) ) );
            ok &= ck( !d.isDisplayable(), "a Settings dialog closed by its close box is disposed, not hidden" );
        }
        // and no closed one lingers: after the cycles above, no Settings dialog is left as a window at all
        ok &= ck( settingsDialogs( mf, false ).isEmpty(), "no closed Settings dialog lingers as a window: "
                + settingsDialogs( mf, false ).size() + " left" );
        return ok;
    }

    private static JButton closeButton( final SettingsDialog d ) {
        final List<AbstractButton> buttons = new ArrayList<>();
        collectButtons( d.getContentPane(), buttons );
        for( final AbstractButton b : buttons ) {
            if ( "Close".equals( b.getText() ) ) {
                return (JButton) b;
            }
        }
        throw new IllegalStateException( "no Close button on the Settings dialog" );
    }

    /** Picks "Settings" on the menu bar as a click does (its menu listener opens the dialog later on the EDT). */
    private static void pickSettingsMenu( final MainFrame mf ) throws Exception {
        SwingUtilities.invokeAndWait( () -> {
            mf._settings_jmenu.setSelected( false );
            mf._settings_jmenu.setSelected( true );
        } );
        SwingUtilities.invokeAndWait( () -> {
            // let the listener's queued open run
        } );
        SwingUtilities.invokeAndWait( () -> mf._settings_jmenu.setSelected( false ) );
    }

    /** The Settings dialogs of {@code mf} that are windows: showing, or ({@code showing_only} false) merely
     *  displayable -- a disposed one is neither. */
    private static List<SettingsDialog> settingsDialogs( final MainFrame mf, final boolean showing_only ) {
        final List<SettingsDialog> found = new ArrayList<>();
        for( final Window w : Window.getWindows() ) {
            if ( ( w instanceof SettingsDialog ) && ( w.getOwner() == mf ) && w.isDisplayable()
                    && ( !showing_only || w.isShowing() ) ) {
                found.add( (SettingsDialog) w );
            }
        }
        return found;
    }

    private static void collectButtons( final Component c, final List<AbstractButton> out ) {
        if ( c instanceof JButton ) {
            out.add( (AbstractButton) c );
        }
        if ( c instanceof Container ) {
            for( final Component child : ( (Container) c ).getComponents() ) {
                collectButtons( child, out );
            }
        }
    }

    private static void closeAnyModal() {
        for( final Window w : Window.getWindows() ) {
            if ( ( w instanceof Dialog ) && ( (Dialog) w ).isModal() && w.isShowing() ) {
                w.dispose();
            }
        }
    }

    private static String describe( final Window w, final MainFrame mf ) {
        if ( w == null ) {
            return "nothing";
        }
        if ( w == mf ) {
            return "the main window";
        }
        return w.getClass().getName() + ( w.isDisplayable() ? "" : " (hidden)" );
    }

    private static Phylogeny tree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "A" );
        a.setDistanceToParent( 1.0 );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "B" );
        b.setDistanceToParent( 2.0 );
        root.addAsChild( a );
        root.addAsChild( b );
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            System.out.println( "  [SettingsDialogStackingTest] " + msg );
        }
        return cond;
    }
}

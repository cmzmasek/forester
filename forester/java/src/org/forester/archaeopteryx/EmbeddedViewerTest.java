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

import java.awt.Dimension;
import java.awt.GraphicsEnvironment;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.EnumMap;
import java.util.Map;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.LookAndFeel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.plaf.metal.MetalLookAndFeel;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * What a viewer from {@link MainFrameApplication#createEmbeddedInstance} promises its host beyond what
 * {@code org.forester.test.EmbeddingHostTest} sees through public members: the look and feel is the host's, the
 * standalone app's saved state is neither read nor written, {@code end()} closes the viewer's windows and silences
 * late work, and a radial tab that is not in front is fitted with the others.
 */
public final class EmbeddedViewerTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "EmbeddedViewer: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final String prev_dir = System.getProperty( GuiPreferences.DIR_PROPERTY );
        // the theme is a java.util.prefs key, not redirected by DIR_PROPERTY: snapshot it, restore it
        final Configuration.UI prev_ui = new Configuration().getUi();
        final LookAndFeel prev_laf = UIManager.getLookAndFeel();
        Path temp = null;
        try {
            temp = Files.createTempDirectory( "aptx-embedded" );
            System.setProperty( GuiPreferences.DIR_PROPERTY, temp.toString() );
            boolean ok = hostsLookAndFeelOk();
            ok = savedStateOk( temp ) && ok; // one failed check must not skip the next
            ok = endOk() && ok;
            ok = radialTabOk() && ok;
            return ok;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        finally {
            if ( prev_dir == null ) {
                System.clearProperty( GuiPreferences.DIR_PROPERTY );
            }
            else {
                System.setProperty( GuiPreferences.DIR_PROPERTY, prev_dir );
            }
            Configuration.saveUiPreference( prev_ui );
            try {
                SwingUtilities.invokeAndWait( () -> {
                    try {
                        if ( UIManager.getLookAndFeel() != prev_laf ) {
                            UIManager.setLookAndFeel( prev_laf );
                        }
                    }
                    catch ( final Exception e ) {
                        e.printStackTrace();
                    }
                } );
            }
            catch ( final Exception e ) {
                e.printStackTrace();
            }
            deleteQuietly( temp );
        }
    }

    /** The host's look and feel stands: at creation, on a theme switch and on Reset to Defaults; the saved theme is
     *  not written; the control panel offers no theme toggle. */
    private static boolean hostsLookAndFeelOk() throws Exception {
        final boolean[] ok = { true };
        Configuration.saveUiPreference( Configuration.UI.FLAT_LIGHT ); // so that a write of DARK would show
        final LookAndFeel[] hosts = { null };
        SwingUtilities.invokeAndWait( () -> {
            try {
                UIManager.setLookAndFeel( new MetalLookAndFeel() ); // a host that is not FlatLaf
            }
            catch ( final Exception e ) {
                ok[ 0 ] = fail( "could not install the host's look and feel: " + e );
            }
            hosts[ 0 ] = UIManager.getLookAndFeel();
        } );
        if ( !( hosts[ 0 ] instanceof MetalLookAndFeel ) ) {
            return fail( "fixture: the host's look and feel must be Metal, is " + hosts[ 0 ] );
        }
        final MainFrameApplication viewer = embedded();
        SwingUtilities.invokeAndWait( () -> {
            if ( UIManager.getLookAndFeel() != hosts[ 0 ] ) {
                ok[ 0 ] = fail( "creating an embedded viewer must leave the host's look and feel installed, found "
                        + UIManager.getLookAndFeel() );
            }
            if ( viewer.getMainPanel().getControlPanel().getThemeToggleButton() != null ) {
                ok[ 0 ] = fail( "an embedded viewer's control panel must offer no theme toggle" );
            }
            viewer.setDarkMode( true );
            if ( UIManager.getLookAndFeel() != hosts[ 0 ] ) {
                ok[ 0 ] = fail( "a theme switch in an embedded viewer must leave the host's look and feel, found "
                        + UIManager.getLookAndFeel() );
            }
            if ( viewer.getConfiguration().getUi() != Configuration.UI.FLAT_LIGHT ) {
                ok[ 0 ] = fail( "an embedded viewer's theme stays the light one, is "
                        + viewer.getConfiguration().getUi() );
            }
            viewer.resetToDefaults();
            if ( UIManager.getLookAndFeel() != hosts[ 0 ] ) {
                ok[ 0 ] = fail( "Reset to Defaults in an embedded viewer must leave the host's look and feel, found "
                        + UIManager.getLookAndFeel() );
            }
            viewer.end();
        } );
        if ( new Configuration().getUi() != Configuration.UI.FLAT_LIGHT ) {
            ok[ 0 ] = fail( "an embedded viewer must not write the saved theme, found " + new Configuration().getUi() );
        }
        return ok[ 0 ];
    }

    /** The standalone app's saved settings and directories are not read, Reset to Defaults does not delete its
     *  settings file, and there is no Open Recent. */
    private static boolean savedStateOk( final Path dir ) throws Exception {
        final boolean[] ok = { true };
        final Path settings = dir.resolve( GuiPreferences.SETTINGS_FILE );
        final boolean default_scale = Options.createInstance().isShowScale();
        final Options saved = Options.createInstance();
        saved.setShowScale( !default_scale );
        new GuiPreferences().saveFrom( saved );
        final File saved_dir = Files.createDirectory( dir.resolve( "last-used" ) ).toFile();
        final Map<DirectoryPreferences.Category, File> dirs = new EnumMap<>( DirectoryPreferences.Category.class );
        dirs.put( DirectoryPreferences.Category.OPEN, saved_dir );
        new DirectoryPreferences().saveFrom( dirs );
        // the drive took: a standalone viewer WOULD read both
        if ( MainFrameApplication.optionsWithSavedPreferences().isShowScale() == default_scale ) {
            return fail( "fixture: the saved settings must flip Show Scale for a standalone viewer" );
        }
        final Map<DirectoryPreferences.Category, File> read = new EnumMap<>( DirectoryPreferences.Category.class );
        new DirectoryPreferences().applyTo( read );
        if ( !saved_dir.equals( read.get( DirectoryPreferences.Category.OPEN ) ) ) {
            return fail( "fixture: the saved Open directory must be read back, got " + read );
        }
        final MainFrameApplication viewer = embedded();
        SwingUtilities.invokeAndWait( () -> {
            if ( viewer.getOptions().isShowScale() != default_scale ) {
                ok[ 0 ] = fail( "an embedded viewer must not read the saved display settings" );
            }
            if ( saved_dir.equals( viewer.getCurrentDir( DirectoryPreferences.Category.OPEN ) ) ) {
                ok[ 0 ] = fail( "an embedded viewer must not read the saved directories" );
            }
            final JMenu file = menu( viewer.getJMenuBar(), "File" );
            if ( file == null ) {
                ok[ 0 ] = fail( "fixture: no File menu" );
            }
            else if ( hasItem( file, "Open Recent" ) ) {
                ok[ 0 ] = fail( "an embedded viewer's File menu must have no Open Recent" );
            }
            viewer.resetToDefaults();
            viewer.end();
        } );
        if ( !Files.exists( settings ) ) {
            ok[ 0 ] = fail( "Reset to Defaults in an embedded viewer must not delete the saved settings file" );
        }
        return ok[ 0 ];
    }

    /** end() closes the windows the viewer opened and leaves the host's; a window shown afterwards is disposed, and a
     *  late process update adds no menu. */
    private static boolean endOk() throws Exception {
        final boolean[] ok = { true };
        final MainFrameApplication viewer = embedded();
        final JFrame[] host = { null };
        final SettingsDialog[] settings = { null };
        final JDialog[] plain = { null };
        try {
            SwingUtilities.invokeAndWait( () -> {
                host[ 0 ] = host( viewer, 900, 700 );
                settings[ 0 ] = new SettingsDialog( viewer );
                viewer.showViewerWindow( settings[ 0 ] );
                plain[ 0 ] = new JDialog( host[ 0 ], "a viewer window" );
                viewer.showViewerWindow( plain[ 0 ] );
            } );
            SwingUtilities.invokeAndWait( () -> {} );
            final int[] menus = { 0 };
            final long[] process = { 0 };
            SwingUtilities.invokeAndWait( () -> {
                if ( !settings[ 0 ].isShowing() || !plain[ 0 ].isShowing() ) {
                    ok[ 0 ] = fail( "fixture: the viewer's windows must be showing before end()" );
                }
                // the drive took: on a live viewer a running process DOES add a menu
                menus[ 0 ] = viewer.getJMenuBar().getMenuCount();
                process[ 0 ] = viewer.getProcessPool().addProcess( "probe" );
                viewer.updateProcessMenu();
            } );
            SwingUtilities.invokeAndWait( () -> {} );
            SwingUtilities.invokeAndWait( () -> {
                if ( viewer.getJMenuBar().getMenuCount() != ( menus[ 0 ] + 1 ) ) {
                    ok[ 0 ] = fail( "fixture: a running process must add a menu to a live viewer" );
                }
                viewer.getProcessPool().removeProcess( process[ 0 ] );
                viewer.updateProcessMenu();
            } );
            SwingUtilities.invokeAndWait( () -> {} );
            SwingUtilities.invokeAndWait( () -> {
                if ( viewer.getJMenuBar().getMenuCount() != menus[ 0 ] ) {
                    ok[ 0 ] = fail( "fixture: the process menu must go when the process ends" );
                }
                viewer.end();
            } );
            SwingUtilities.invokeAndWait( () -> {} );
            final JDialog[] late = { null };
            SwingUtilities.invokeAndWait( () -> {
                if ( !viewer.isClosed() ) {
                    ok[ 0 ] = fail( "end() must close the viewer" );
                }
                if ( settings[ 0 ].isDisplayable() ) {
                    ok[ 0 ] = fail( "end() must close the viewer's Settings dialog" );
                }
                if ( plain[ 0 ].isDisplayable() ) {
                    ok[ 0 ] = fail( "end() must close every window the viewer opened" );
                }
                if ( !host[ 0 ].isDisplayable() ) {
                    ok[ 0 ] = fail( "end() must leave the host's window open" );
                }
                late[ 0 ] = new JDialog( host[ 0 ], "late" );
                late[ 0 ].pack(); // displayable, as a window built by late work would be
                viewer.showViewerWindow( late[ 0 ] );
                if ( late[ 0 ].isVisible() || late[ 0 ].isDisplayable() ) {
                    ok[ 0 ] = fail( "a window shown after end() must be disposed, not shown" );
                }
                viewer.getProcessPool().addProcess( "late" );
                viewer.updateProcessMenu();
            } );
            SwingUtilities.invokeAndWait( () -> {} );
            SwingUtilities.invokeAndWait( () -> {
                if ( viewer.getJMenuBar().getMenuCount() != menus[ 0 ] ) {
                    ok[ 0 ] = fail( "a process update after end() must change nothing, the menu bar has "
                            + viewer.getJMenuBar().getMenuCount() + " menus, had " + menus[ 0 ] );
                }
                if ( viewer.isProcessAnimationRunningForTest() ) {
                    ok[ 0 ] = fail( "a process update after end() must not start the activity animation" );
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> {
                for( final java.awt.Window w : new java.awt.Window[] { settings[ 0 ], plain[ 0 ], host[ 0 ] } ) {
                    if ( w != null ) {
                        w.dispose();
                    }
                }
            } );
        }
        return ok[ 0 ];
    }

    /** A radial tab laid out before the host gave the viewer a size, and not in front when it does, is fitted to the
     *  viewport like the tab in front: it must not keep the fallback diameter. */
    private static boolean radialTabOk() throws Exception {
        final boolean[] ok = { true };
        final MainFrameApplication viewer = embedded();
        final JFrame[] host = { null };
        try {
            SwingUtilities.invokeAndWait( () -> {
                final ControlPanel cp = viewer.getMainPanel().getControlPanel();
                viewer.getMainPanel().getTabbedPane().setSelectedIndex( 0 );
                cp.getLayoutButton( LayoutIcon.Kind.CIRCULAR ).doClick();
                if ( !viewer.getMainPanel().getTreePanels().get( 0 ).isRadialLayout() ) {
                    ok[ 0 ] = fail( "fixture: tab 0 must be radial" );
                }
                viewer.getMainPanel().getTabbedPane().setSelectedIndex( 1 );
                if ( viewer.getMainPanel().getTreePanels().get( 1 ).isRadialLayout() ) {
                    ok[ 0 ] = fail( "fixture: tab 1 must stay rectangular" );
                }
                host[ 0 ] = host( viewer, 1100, 520 ); // a viewport whose short side is not the fallback diameter
            } );
            SwingUtilities.invokeAndWait( () -> {} );
            SwingUtilities.invokeAndWait( () -> {} );
            SwingUtilities.invokeAndWait( () -> {
                if ( viewer.getMainPanel().getTabbedPane().getSelectedIndex() != 1 ) {
                    ok[ 0 ] = fail( "fixture: tab 1 must be in front" );
                }
                final Dimension viewport = viewer.getMainPanel().getSizeOfViewport();
                final int expected = Math.min( viewport.width, viewport.height );
                if ( ( expected < 200 ) || ( expected == 600 ) ) {
                    ok[ 0 ] = fail( "fixture: the viewport's short side must be a usable diameter other than the "
                            + "fallback 600, is " + expected + " (" + viewport + ")" );
                }
                final int diameter = viewer.getMainPanel().getTreePanels().get( 0 ).radialDiameter();
                if ( diameter != expected ) {
                    ok[ 0 ] = fail( "a radial tab not in front must be fitted to the viewport: diameter " + diameter
                            + ", viewport " + viewport );
                }
                viewer.end();
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> {
                if ( host[ 0 ] != null ) {
                    host[ 0 ].dispose();
                }
            } );
        }
        return ok[ 0 ];
    }

    private static MainFrameApplication embedded() throws Exception {
        return MainFrameApplication.createEmbeddedInstance( new Phylogeny[] { tree( "((A,B),(C,D),(E,F),(G,H));" ),
                tree( "(((P,Q),R),((S,T),(U,(V,W))));" ) }, new Configuration(), "embedded" );
    }

    /** Hosts the viewer as another program would; call on the event thread. */
    private static JFrame host( final MainFrameApplication viewer, final int w, final int h ) {
        final JFrame host = new JFrame( "Host" );
        host.setJMenuBar( viewer.getJMenuBar() );
        host.getContentPane().add( viewer );
        host.setSize( w, h );
        host.setVisible( true );
        return host;
    }

    private static Phylogeny tree( final String newick ) throws Exception {
        return ParserBasedPhylogenyFactory.getInstance().create( newick, new NHXParser() )[ 0 ];
    }

    private static JMenu menu( final JMenuBar bar, final String title ) {
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m != null ) && title.equals( m.getText() ) ) {
                return m;
            }
        }
        return null;
    }

    private static boolean hasItem( final JMenu menu, final String text ) {
        for( int i = 0; i < menu.getItemCount(); ++i ) {
            if ( ( menu.getItem( i ) != null ) && text.equals( menu.getItem( i ).getText() ) ) {
                return true;
            }
        }
        return false;
    }

    private static void deleteQuietly( final Path dir ) {
        if ( dir == null ) {
            return;
        }
        try ( java.util.stream.Stream<Path> walk = Files.walk( dir ) ) {
            walk.sorted( java.util.Comparator.reverseOrder() ).forEach( p -> p.toFile().delete() );
        }
        catch ( final Exception e ) {
            // a temp directory left behind is harmless
        }
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [EmbeddedViewerTest] " + message );
        return false;
    }
}

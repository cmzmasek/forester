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
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Window;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.ButtonModel;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.LookAndFeel;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.plaf.ComponentUI;

import com.formdev.flatlaf.FlatLaf;

/**
 * The window opens wearing the theme it was asked for, all of it (Christian, 2026-11-05: launched in the dark theme,
 * the menus "File", "Edit", "Analysis", etc. "become grey on grey when mousing over"; a switch to light and back
 * cured it). A top-level menu's hover fill is DERIVED at paint time from the first opaque thing behind it -- FlatLaf's
 * dark {@code MenuBar.hoverBackground} is lighten(background, 10%, derived) -- and that was the frame's root pane,
 * which the JFrame constructor had built under the platform's look-and-feel before ours was installed: #D5D5D5
 * behind #DDDDDD text (measured; #54595C after a theme switch).
 * <p>
 * Each window is built by {@code createInstance}, as the application builds it, while ANOTHER look-and-feel is in
 * place: the platform's own (what a launch starts from; the reported case) and the other FlatLaf theme. Then:
 * <ol>
 * <li>no component of the frame, its menu bar and every menu's popup wears a UI delegate that is not FlatLaf's
 * (the whole tree, iterated -- not a list of the components that went wrong);</li>
 * <li>every top-level menu, painted hovered, shows its title at a contrast of at least 3:1 (the WCAG floor for UI
 * text) against the fill -- after the precondition that the hover drive took (the hovered paint is a fill, the plain
 * one is not);</li>
 * <li>launch paints what the theme switch paints: each hovered title is pixel-identical after a round trip through
 * the other theme -- the re-install and re-skin of every window that {@link MainFrame#setDarkMode} (the control
 * panel's toggle) performs, the path that always worked; done here by those two steps rather than through
 * {@code setDarkMode}, which would also write the developer's REAL theme preference (java.util.prefs).</li>
 * </ol>
 * The installed look-and-feel is restored. Needs a display: a no-op when headless.
 */
public final class LaunchThemeTest {

    /** The WCAG contrast floor for user-interface text, pinned as a literal. */
    private static final double MIN_CONTRAST = 3.0;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LaunchTheme: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final LookAndFeel installed = UIManager.getLookAndFeel();
        boolean ok = true;
        try {
            ok &= launchOk( Configuration.UI.FLAT_DARK, null );
            ok &= launchOk( Configuration.UI.FLAT_DARK, Configuration.UI.FLAT_LIGHT );
            ok &= launchOk( Configuration.UI.FLAT_LIGHT, Configuration.UI.FLAT_DARK );
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            ok = false;
        }
        finally {
            try {
                SwingUtilities.invokeAndWait( () -> {
                    if ( installed instanceof FlatLaf ) {
                        MainFrame.installLookAndFeel( ( (FlatLaf) installed ).isDark() ? Configuration.UI.FLAT_DARK
                                : Configuration.UI.FLAT_LIGHT );
                    }
                } );
            }
            catch ( final Exception e ) {
                e.printStackTrace();
                ok = false;
            }
        }
        return ok;
    }

    /**
     * Builds a window asking for {@code ui} while {@code before} is installed ({@code null}: the platform's own
     * look-and-feel, as at launch), and checks it.
     */
    private static boolean launchOk( final Configuration.UI ui, final Configuration.UI before ) throws Exception {
        final String when = ui + " built under " + ( ( before == null ) ? "the platform look-and-feel" : before );
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            try {
                if ( before == null ) {
                    UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName() );
                }
                else {
                    MainFrame.installLookAndFeel( before );
                }
            }
            catch ( final Exception e ) {
                ok[ 0 ] = fail( when + ": could not install the starting look-and-feel: " + e );
            }
            // precondition: the starting look-and-feel is NOT the one asked for (else nothing is exercised)
            final LookAndFeel start = UIManager.getLookAndFeel();
            final boolean start_is_target = ( start instanceof FlatLaf )
                    && ( ( (FlatLaf) start ).isDark() == ( ui == Configuration.UI.FLAT_DARK ) );
            ok[ 0 ] &= ck( !start_is_target, when + ": precondition: the window is built under another look-and-feel, was "
                    + start.getName() );
        } );
        if ( !ok[ 0 ] ) {
            return false;
        }
        final Configuration conf = new Configuration();
        conf.setUi( ui );
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication.createInstance( null, conf, "launch theme" ) );
        SwingUtilities.invokeAndWait( () -> {
            try {
                final MainFrame frame = mf[ 0 ];
                ok[ 0 ] &= everyUiIsFlatLaf( frame, when );
                // each check stands on its own: a failed contrast does not skip the comparison with the switch
                final Map<String, BufferedImage> at_launch = new HashMap<>();
                ok[ 0 ] &= hoveredTitlesOk( frame, when + ", at launch", at_launch );
                // the switch's re-install + re-skin of every window, there and back (no preference write)
                reskin( ( ui == Configuration.UI.FLAT_DARK ) ? Configuration.UI.FLAT_LIGHT : Configuration.UI.FLAT_DARK );
                reskin( ui );
                final Map<String, BufferedImage> after = new HashMap<>();
                ok[ 0 ] &= hoveredTitlesOk( frame, when + ", after a theme round trip", after );
                ok[ 0 ] &= ck( at_launch.keySet().equals( after.keySet() ), when + ": the same menus before and after: "
                        + at_launch.keySet() + " / " + after.keySet() );
                for( final String title : at_launch.keySet() ) {
                    final BufferedImage b = after.get( title );
                    ok[ 0 ] &= ck( ( b != null ) && samePixels( at_launch.get( title ), b ), when + ": \"" + title
                            + "\" hovered at launch paints as it does after a theme switch" );
                }
            }
            finally {
                ( (JFrame) mf[ 0 ] ).dispose();
            }
        } );
        return ok[ 0 ];
    }

    /** What {@link MainFrame#setDarkMode} does to the windows: install the theme, re-skin every window. */
    private static void reskin( final Configuration.UI ui ) {
        MainFrame.installLookAndFeel( ui );
        for( final Window w : Window.getWindows() ) {
            SwingUtilities.updateComponentTreeUI( w );
        }
    }

    /** (1) Every component of the frame, its menu bar and every menu popup has a FlatLaf UI delegate (or none). */
    private static boolean everyUiIsFlatLaf( final JFrame frame, final String when ) {
        final List<Component> all = new ArrayList<>();
        collect( frame, all );
        int with_ui = 0;
        boolean ok = true;
        for( final Component c : all ) {
            final ComponentUI ui = uiOf( c );
            if ( ui == null ) {
                continue;
            }
            ++with_ui;
            if ( !ui.getClass().getName().startsWith( "com.formdev.flatlaf." ) ) {
                ok = fail( when + ": " + c.getClass().getName() + " wears " + ui.getClass().getName() );
            }
        }
        // a sweep that found nothing to look at proves nothing: the frame has its root pane, panels, menus, items
        ok &= ck( with_ui >= 100, when + ": precondition: the sweep reached the window's components, found " + with_ui );
        return ok;
    }

    private static void collect( final Component c, final List<Component> all ) {
        all.add( c );
        if ( c instanceof JMenu ) {
            collect( ( (JMenu) c ).getPopupMenu(), all );
        }
        if ( c instanceof Container ) {
            for( final Component child : ( (Container) c ).getComponents() ) {
                collect( child, all );
            }
        }
    }

    /** A component's installed UI delegate, or null when it has none (a layered pane, an AWT component). */
    private static ComponentUI uiOf( final Component c ) {
        if ( !( c instanceof JComponent ) ) {
            return null;
        }
        try {
            final Object ui = c.getClass().getMethod( "getUI" ).invoke( c );
            return ( ui instanceof ComponentUI ) ? (ComponentUI) ui : null;
        }
        catch ( final NoSuchMethodException e ) {
            return null;
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    /**
     * (2) Every showing top-level menu, painted hovered: its title at >= 3:1 against the fill. The hovered paints go
     * into {@code paints} by title, whatever the verdict.
     */
    private static boolean hoveredTitlesOk( final JFrame frame, final String when,
                                            final Map<String, BufferedImage> paints ) {
        final JMenuBar bar = frame.getJMenuBar();
        if ( !ck( bar != null, when + ": the window has a menu bar" ) ) {
            return false;
        }
        boolean ok = true;
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final JMenu m = bar.getMenu( i );
            if ( ( m == null ) || !m.isShowing() || ( m.getWidth() <= 0 ) || ( m.getHeight() <= 0 ) ) {
                continue; // the found counter, the glue, a hidden process or error menu
            }
            final BufferedImage plain = paint( m, false );
            final BufferedImage hovered = paint( m, true );
            // precondition: the drive took -- hovered is a fill over most of the title, plain is not
            final double hovered_fill = opaqueFraction( hovered );
            final double plain_fill = opaqueFraction( plain );
            ok &= ck( ( hovered_fill >= 0.5 ) && ( plain_fill < hovered_fill ), when + ": precondition: \"" + m.getText()
                    + "\" hovered paints a fill (opaque " + pct( hovered_fill ) + ", plain " + pct( plain_fill ) + ")" );
            final int fill = modeOpaque( hovered );
            final int text = farthestFrom( hovered, fill );
            final double contrast = contrast( fill, text );
            ok &= ck( contrast >= MIN_CONTRAST, when + ": \"" + m.getText() + "\" hovered: text " + hex( text ) + " on "
                    + hex( fill ) + ", contrast " + String.format( "%.2f", contrast ) + " (needs >= " + MIN_CONTRAST + ")" );
            paints.put( m.getText(), hovered );
        }
        ok &= ck( paints.size() >= 6, when + ": precondition: the top-level menus were painted, " + paints.keySet() );
        return ok;
    }

    private static BufferedImage paint( final JMenu m, final boolean hover ) {
        final ButtonModel model = m.getModel();
        final boolean was = model.isRollover();
        model.setRollover( hover );
        try {
            final BufferedImage img = new BufferedImage( m.getWidth(), m.getHeight(), BufferedImage.TYPE_INT_ARGB );
            final Graphics2D g = img.createGraphics();
            try {
                m.paint( g );
            }
            finally {
                g.dispose();
            }
            return img;
        }
        finally {
            model.setRollover( was );
        }
    }

    private static double opaqueFraction( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                if ( ( img.getRGB( x, y ) >>> 24 ) == 255 ) {
                    ++n;
                }
            }
        }
        return n / (double) ( img.getWidth() * img.getHeight() );
    }

    /** The most frequent fully opaque colour: the fill. */
    private static int modeOpaque( final BufferedImage img ) {
        final Map<Integer, Integer> counts = new HashMap<>();
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int p = img.getRGB( x, y );
                if ( ( p >>> 24 ) == 255 ) {
                    counts.merge( p & 0xFFFFFF, 1, Integer::sum );
                }
            }
        }
        int best = 0;
        int best_n = -1;
        for( final Map.Entry<Integer, Integer> e : counts.entrySet() ) {
            if ( e.getValue() > best_n ) {
                best = e.getKey();
                best_n = e.getValue();
            }
        }
        return best;
    }

    /** The opaque pixel of the greatest contrast against the fill: the core of the title's glyphs. */
    private static int farthestFrom( final BufferedImage img, final int fill ) {
        int best = fill;
        double best_c = 1.0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int p = img.getRGB( x, y );
                if ( ( p >>> 24 ) == 255 ) {
                    final double c = contrast( fill, p & 0xFFFFFF );
                    if ( c > best_c ) {
                        best_c = c;
                        best = p & 0xFFFFFF;
                    }
                }
            }
        }
        return best;
    }

    /** WCAG 2 contrast ratio of two sRGB colours. */
    static double contrast( final int a, final int b ) {
        final double la = luminance( a );
        final double lb = luminance( b );
        return ( Math.max( la, lb ) + 0.05 ) / ( Math.min( la, lb ) + 0.05 );
    }

    private static double luminance( final int rgb ) {
        return ( 0.2126 * channel( ( rgb >> 16 ) & 0xFF ) ) + ( 0.7152 * channel( ( rgb >> 8 ) & 0xFF ) )
                + ( 0.0722 * channel( rgb & 0xFF ) );
    }

    private static double channel( final int v ) {
        final double c = v / 255.0;
        return ( c <= 0.03928 ) ? ( c / 12.92 ) : Math.pow( ( c + 0.055 ) / 1.055, 2.4 );
    }

    private static boolean samePixels( final BufferedImage a, final BufferedImage b ) {
        if ( ( a.getWidth() != b.getWidth() ) || ( a.getHeight() != b.getHeight() ) ) {
            return false;
        }
        for( int y = 0; y < a.getHeight(); ++y ) {
            for( int x = 0; x < a.getWidth(); ++x ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    return false;
                }
            }
        }
        return true;
    }

    private static String hex( final int rgb ) {
        return String.format( "#%06X", rgb & 0xFFFFFF );
    }

    private static String pct( final double f ) {
        return String.format( "%.0f%%", 100 * f );
    }

    private static boolean ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            System.out.println( "  [LaunchThemeTest] " + msg );
        }
        return cond;
    }

    private static boolean fail( final String msg ) {
        return ck( false, msg );
    }
}

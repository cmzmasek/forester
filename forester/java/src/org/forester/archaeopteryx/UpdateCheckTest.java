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
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;

import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Tests for the launch-time update check: version comparison and tag parsing (pure), a check against a local file
 * standing in for the GitHub API (newer / same / older / garbage), a check against a dead address that must end
 * silently -- nothing on stderr, nothing in the error log, no callback -- and the Help-menu line the main frame
 * adds when a newer release is known (headful part; skipped when headless).
 */
public final class UpdateCheckTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "UpdateCheck: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return versions() && parsing() && localFile() && failsSilently() && helpMenuLine();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean versions() {
        return check( "0.11.140 > 0.11.135", UpdateCheck.isNewer( "0.11.140", "0.11.135" ) )
                && check( "0.12 > 0.11.135", UpdateCheck.isNewer( "0.12", "0.11.135" ) )
                && check( "1.0.0 > 0.99.99", UpdateCheck.isNewer( "1.0.0", "0.99.99" ) )
                && check( "v-prefix accepted", UpdateCheck.isNewer( "v0.11.136", "0.11.135" ) )
                && check( "same is not newer", !UpdateCheck.isNewer( "0.11.135", "0.11.135" ) )
                && check( "same with trailing zero is not newer", !UpdateCheck.isNewer( "0.11.135.0", "0.11.135" ) )
                && check( "older is not newer", !UpdateCheck.isNewer( "0.11.134", "0.11.135" ) )
                && check( "older major is not newer", !UpdateCheck.isNewer( "0.9.200", "0.11.1" ) )
                && check( "numeric, not lexical", UpdateCheck.isNewer( "0.11.1000", "0.11.999" ) )
                && check( "garbage is never newer", !UpdateCheck.isNewer( "latest", "0.11.135" ) )
                && check( "partly garbage is never newer", !UpdateCheck.isNewer( "0.11.135-beta", "0.11.135" ) )
                && check( "empty is never newer", !UpdateCheck.isNewer( "", "0.11.135" ) )
                && check( "null is never newer", !UpdateCheck.isNewer( null, "0.11.135" ) )
                && check( "current garbage -> false", !UpdateCheck.isNewer( "0.11.136", "dev" ) )
                && eq( "normalize v", "0.11.140", UpdateCheck.normalize( " v0.11.140 " ) )
                && eq( "normalize blank", null, UpdateCheck.normalize( "  " ) )
                && eq( "current version is itself a plain dotted number (or the check is dead)", true,
                       UpdateCheck.isNewer( "999.0", AptxConstants.VERSION ) );
    }

    private static boolean parsing() {
        final String json = "{\"url\":\"https://api.github.com/x\",\"tag_name\":\"0.11.140\",\"name\":\"0.11.140\","
                + "\"prerelease\":false}";
        return eq( "tag_name", "0.11.140", UpdateCheck.latestTagFromJson( json ) )
                && eq( "spaces around the colon", "1.2.3", UpdateCheck.latestTagFromJson( "{ \"tag_name\" : \"1.2.3\" }" ) )
                && eq( "no tag -> null", null, UpdateCheck.latestTagFromJson( "{\"message\":\"Not Found\"}" ) )
                && eq( "empty -> null", null, UpdateCheck.latestTagFromJson( "" ) )
                && eq( "null -> null", null, UpdateCheck.latestTagFromJson( null ) );
    }

    /** A local file plays the API: the whole pipeline (fetch -> parse -> compare -> callback on the EDT) runs. */
    private static boolean localFile() throws Exception {
        final File f = File.createTempFile( "aptx-release", ".json" );
        f.deleteOnExit();
        final String url = f.toURI().toString();
        boolean ok = true;
        Files.write( f.toPath(), "{\"tag_name\":\"v0.11.999\",\"prerelease\":false}".getBytes( StandardCharsets.UTF_8 ) );
        ok = ok && eq( "newer from file", "0.11.999", UpdateCheck.newerVersionAt( url, "0.11.135" ) );
        ok = ok && check( "fetch returns the body", UpdateCheck.fetch( url ).contains( "0.11.999" ) );
        final AtomicReference<String> got = new AtomicReference<>();
        final CountDownLatch latch = new CountDownLatch( 1 );
        final Thread t = UpdateCheck.start( url, "0.11.135", 0, v -> {
            got.set( v + ( SwingUtilities.isEventDispatchThread() ? " on EDT" : " OFF EDT" ) );
            latch.countDown();
        } );
        ok = ok && check( "daemon thread (must never keep the JVM alive)", t.isDaemon() )
                && check( "callback arrives", latch.await( 10, TimeUnit.SECONDS ) )
                && eq( "callback value, on the EDT", "0.11.999 on EDT", got.get() );
        Files.write( f.toPath(), "{\"tag_name\":\"0.11.135\"}".getBytes( StandardCharsets.UTF_8 ) );
        ok = ok && eq( "same version -> nothing", null, UpdateCheck.newerVersionAt( url, "0.11.135" ) );
        Files.write( f.toPath(), "{\"tag_name\":\"0.11.100\"}".getBytes( StandardCharsets.UTF_8 ) );
        ok = ok && eq( "older -> nothing", null, UpdateCheck.newerVersionAt( url, "0.11.135" ) );
        Files.write( f.toPath(), "<html>rate limited</html>".getBytes( StandardCharsets.UTF_8 ) );
        ok = ok && eq( "garbage -> nothing", null, UpdateCheck.newerVersionAt( url, "0.11.135" ) );
        // the once-per-launch guard
        UpdateCheck.resetForTest();
        final int[] calls = { 0 };
        UpdateCheck.startAtLaunch( "0.11.135", v -> calls[ 0 ]++ ); // (real URL; may or may not answer -- not asserted)
        UpdateCheck.startAtLaunch( "0.11.135", v -> calls[ 0 ]++ );
        UpdateCheck.resetForTest();
        return ok;
    }

    /** No network / dead host / nonsense URL: null, no callback, and NOT ONE BYTE on stderr or in the error log. */
    private static boolean failsSilently() throws Exception {
        final PrintStream old_err = System.err;
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        System.setErr( new PrintStream( captured, true, "UTF-8" ) );
        boolean ok = true;
        try {
            // a closed port on localhost: connection refused, immediately
            ok = ok && eq( "refused connection -> null", null, UpdateCheck.newerVersionAt( "http://127.0.0.1:1/x", "0.1" ) );
            ok = ok && eq( "nonsense URL -> null", null, UpdateCheck.newerVersionAt( "not a url at all", "0.1" ) );
            ok = ok && eq( "missing file -> null", null,
                           UpdateCheck.newerVersionAt( new File( "/definitely/not/here.json" ).toURI().toString(), "0.1" ) );
            ok = ok && eq( "null url -> null", null, UpdateCheck.newerVersionAt( null, "0.1" ) );
            final CountDownLatch never = new CountDownLatch( 1 );
            final Thread t = UpdateCheck.start( "http://127.0.0.1:1/x", "0.1", 0, v -> never.countDown() );
            t.join( 10000 );
            ok = ok && check( "thread ended", !t.isAlive() ) && check( "no callback on failure", never.getCount() == 1 );
            // a callback that throws must not escape either
            final File f = File.createTempFile( "aptx-release", ".json" );
            f.deleteOnExit();
            Files.write( f.toPath(), "{\"tag_name\":\"9.9.9\"}".getBytes( StandardCharsets.UTF_8 ) );
            final Thread t2 = UpdateCheck.start( f.toURI().toString(), "0.1", 0, v -> {
                throw new IllegalStateException( "listener bug" );
            } );
            t2.join( 10000 );
            SwingUtilities.invokeAndWait( () -> {
            } ); // drain the EDT so the throwing callback has run
            ok = ok && check( "thread ended", !t2.isAlive() );
        }
        finally {
            System.setErr( old_err );
        }
        final String err = captured.toString( "UTF-8" );
        ok = ok && eq( "nothing on stderr", "", err );
        final ErrorLog log = ErrorLog.instance();
        ok = ok && check( "nothing in the error log", ( log == null ) || !log.hasEntries() );
        return ok;
    }

    private static boolean helpMenuLine() throws Exception {
        boolean ok = eq( "label", "New version available: 0.11.140", MainFrame.updateAvailableLabel( "0.11.140" ) );
        if ( GraphicsEnvironment.isHeadless() ) {
            return ok;
        }
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        root.setName( "r" );
        phy.setRoot( root );
        phy.setRooted( true );
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "update" ) );
        final boolean[] okk = { ok };
        SwingUtilities.invokeAndWait( () -> {
            final MainFrame f = mf[ 0 ];
            // OPT-IN: a fresh window must start with the check switched OFF
            if ( f.getOptions().isCheckForUpdatesAtLaunch() || f._check_for_updates_cbmi.isSelected() ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] the update check must be OFF by default (opt-in)" );
            }
            // How wide the Help menu is on its own, BEFORE the notice line is added.
            final javax.swing.JPopupMenu popup = f._help_jmenu.getPopupMenu();
            final int menu_w_without_notice = popup.getPreferredSize().width;
            final int before = f._help_jmenu.getItemCount();
            if ( f._update_available_item != null ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] no update line before a check says so" );
            }
            f.showUpdateAvailable( "0.11.140" );
            final JMenuItem item = f._update_available_item;
            if ( ( item == null ) || ( f._help_jmenu.getItem( 0 ) != item ) ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] the update line must be the FIRST Help item" );
            }
            else if ( !"New version available: 0.11.140".equals( item.getText() ) ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] wrong text: " + item.getText() );
            }
            if ( f._help_jmenu.getItemCount() != before + 2 ) { // the line + its separator
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] expected one line + one separator added" );
            }
            // quiet: the menu's own font and colour, no bold, no accent
            final JMenuItem plain = f._help_jmenu.getItem( 2 );
            if ( ( item != null ) && ( plain != null ) ) {
                if ( item.getFont().isBold() && !plain.getFont().isBold() ) {
                    okk[ 0 ] = false;
                    System.out.println( "  [UpdateCheckTest] the update line must not be bold" );
                }
                if ( !plain.getForeground().equals( item.getForeground() ) ) {
                    okk[ 0 ] = false;
                    System.out.println( "  [UpdateCheckTest] the update line must not be coloured: "
                            + item.getForeground() );
                }
            }
            // The reported problem was that the notice was CUT OFF the first time the menu was shown. It is
            // inserted into a menu that was built long before, so the safe property -- and the one that makes it
            // independent of when and how the popup re-measures itself -- is that the line does not need the
            // menu to get any wider than it already is. The old label ("Update available: Archaeopteryx x.y.z",
            // bold) was ~50% wider than the entire rest of the Help menu; this one fits inside it.
            if ( item.getPreferredSize().width > menu_w_without_notice ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] the update line (" + item.getPreferredSize().width
                        + " px) must fit the Help menu's own width (" + menu_w_without_notice
                        + " px), or it can be clipped: " + item.getText() );
            }
            f.showUpdateAvailable( "0.11.141" ); // a repeat refreshes, never duplicates
            if ( ( f._help_jmenu.getItemCount() != before + 2 )
                    || !"New version available: 0.11.141".equals( f._update_available_item.getText() ) ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] a second notice must refresh the same line" );
            }
            // the setting round-trips through the Settings checkbox item (UI -> Options -> UI)
            f._check_for_updates_cbmi.setSelected( false );
            f.updateOptions( f.getOptions() );
            if ( f.getOptions().isCheckForUpdatesAtLaunch() ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] unticking the setting must reach Options" );
            }
            f.getOptions().setCheckForUpdatesAtLaunch( true );
            f.applyOptionsToMenuStates( f.getOptions() );
            if ( !f._check_for_updates_cbmi.isSelected() ) {
                okk[ 0 ] = false;
                System.out.println( "  [UpdateCheckTest] Options -> checkbox must re-seed the control" );
            }
            ( (JFrame) f ).dispose();
        } );
        return okk[ 0 ];
    }

    private static boolean eq( final String what, final Object expected, final Object actual ) {
        if ( ( expected == null ) ? ( actual == null ) : expected.equals( actual ) ) {
            return true;
        }
        System.out.println( "  [UpdateCheckTest] " + what + ": expected <" + expected + "> but got <" + actual + ">" );
        return false;
    }

    private static boolean check( final String what, final boolean condition ) {
        if ( condition ) {
            return true;
        }
        System.out.println( "  [UpdateCheckTest] " + what );
        return false;
    }

    private UpdateCheckTest() {
        // not instantiable
    }
}

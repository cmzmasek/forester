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

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * The crash log an INSTALLED Archaeopteryx needs, because it has no console to print to.
 * <p>
 * The behaviour that matters is the de-duplication: the bug that motivated this threw inside {@code paintComponent},
 * so it recurred on every repaint. Writing that trace out each time would bury the log (and could fill a disk),
 * which is why a repeat only increments a counter.
 */
public final class ErrorLogTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ErrorLog: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return signature() && dedup() && errorVsChatter() && stderrTee() && sizeCap() && location()
                && tailAndConsole() && badPropertyCannotBlockStartup() && menu() && markerAfterStartupError();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ErrorLogTest] " + msg );
        return false;
    }

    private static Path tmp( final String name ) throws Exception {
        final Path dir = Files.createTempDirectory( "aptx-errorlog" );
        dir.toFile().deleteOnExit();
        return dir.resolve( name );
    }

    private static Throwable thrown( final String msg ) {
        try {
            throw new IllegalStateException( msg );
        }
        catch ( final IllegalStateException e ) {
            return e;
        }
    }

    private static boolean signature() {
        final Throwable a = thrown( "same" );
        final Throwable b = thrown( "same" );
        // same type + message + top frame -> the same failure, even though these are different objects
        if ( !ErrorLog.signature( a ).equals( ErrorLog.signature( b ) ) ) {
            return fail( "two throws of the same failure should share a signature" );
        }
        if ( ErrorLog.signature( a ).equals( ErrorLog.signature( thrown( "different" ) ) ) ) {
            return fail( "a different message is a different failure" );
        }
        if ( !"null".equals( ErrorLog.signature( null ) ) ) {
            return fail( "a null throwable must not blow up the logger" );
        }
        return true;
    }

    /** The point of the whole class: a repainting loop must not rewrite its trace thousands of times. */
    private static boolean dedup() {
        try {
            final Path f = tmp( "dedup.log" );
            final ErrorLog log = new ErrorLog( f, null );
            final Throwable boom = thrown( "paint blew up" );
            for( int i = 0; i < 500; ++i ) {
                log.logThrowable( "Uncaught in AWT-EventQueue-0", boom );
            }
            if ( log.pendingRepeatsForTest() != 499 ) {
                return fail( "499 repeats should be pending, got " + log.pendingRepeatsForTest() );
            }
            log.flushRepeats();
            final String text = Files.readString( f, StandardCharsets.UTF_8 );
            final int traces = text.split( "IllegalStateException", -1 ).length - 1;
            if ( traces > 3 ) {
                return fail( "the trace should be written ONCE, not per repaint -- found " + traces
                        + " mentions" );
            }
            if ( !text.contains( "repeated 499 more times" ) ) {
                return fail( "the repeats should be summarised in one line: " + text );
            }
            // a DIFFERENT failure flushes the pending count and writes its own trace
            log.logThrowable( "Uncaught in AWT-EventQueue-0", thrown( "something else" ) );
            final String text2 = Files.readString( f, StandardCharsets.UTF_8 );
            if ( !text2.contains( "something else" ) ) {
                return fail( "a different failure must still be written in full" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "exception: " + e );
        }
    }

    /** Only a real failure raises the "there is something to send" signal; stderr chatter does not. */
    private static boolean errorVsChatter() {
        try {
            final Path f = tmp( "notify.log" );
            final int[] notified = { 0 };
            final ErrorLog log = new ErrorLog( f, () -> notified[ 0 ]++ );
            final PrintStream tee = log.tee( new PrintStream( PrintStream.nullOutputStream() ) );
            tee.println( "WARNING: some library grumbling" );
            tee.flush();
            if ( notified[ 0 ] != 0 ) {
                return fail( "a stderr warning is captured but must NOT be announced as an error" );
            }
            if ( log.hasErrors() ) {
                return fail( "stderr chatter is not an error" );
            }
            log.logThrowable( "Uncaught in AWT-EventQueue-0", thrown( "real failure" ) );
            if ( notified[ 0 ] != 1 ) {
                return fail( "the first real failure should notify exactly once, got " + notified[ 0 ] );
            }
            log.logThrowable( "Uncaught in AWT-EventQueue-0", thrown( "another failure" ) );
            if ( notified[ 0 ] != 1 ) {
                return fail( "only the FIRST failure notifies, got " + notified[ 0 ] );
            }
            if ( !log.hasErrors() || !log.hasEntries() ) {
                return fail( "after a failure the log should report both" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "exception: " + e );
        }
    }

    /** The tee must reach BOTH the original stream and the file -- and a line at a time, not a file write per byte. */
    private static boolean stderrTee() {
        try {
            final Path f = tmp( "tee.log" );
            final ErrorLog log = new ErrorLog( f, null );
            final java.io.ByteArrayOutputStream original = new java.io.ByteArrayOutputStream();
            final PrintStream tee = log.tee( new PrintStream( original, true, StandardCharsets.UTF_8 ) );
            tee.println( "first line" );
            tee.println( "second line" );
            tee.flush();
            final String passed_through = original.toString( StandardCharsets.UTF_8 );
            if ( !passed_through.contains( "first line" ) || !passed_through.contains( "second line" ) ) {
                return fail( "the tee must still write to the original stderr: " + passed_through );
            }
            final String logged = Files.readString( f, StandardCharsets.UTF_8 );
            if ( !logged.contains( "first line" ) || !logged.contains( "second line" ) ) {
                return fail( "the tee must also write to the log: " + logged );
            }
            // an unterminated write is held back until its newline, rather than dribbling out a byte at a time
            tee.print( "no newline yet" );
            tee.flush();
            if ( Files.readString( f, StandardCharsets.UTF_8 ).contains( "no newline yet" ) ) {
                return fail( "a partial line should be buffered until its newline" );
            }
            tee.println();
            if ( !Files.readString( f, StandardCharsets.UTF_8 ).contains( "no newline yet" ) ) {
                return fail( "the buffered line should be written once its newline arrives" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "exception: " + e );
        }
    }

    /** A runaway loop must not fill the user's disk. */
    private static boolean sizeCap() {
        try {
            final Path f = tmp( "cap.log" );
            final ErrorLog log = new ErrorLog( f, null );
            final StringBuilder big = new StringBuilder();
            for( int i = 0; i < 4000; ++i ) {
                big.append( "0123456789012345678901234567890123456789012345678901234567890123\n" );
            }
            final PrintStream tee = log.tee( new PrintStream( PrintStream.nullOutputStream() ) );
            for( int i = 0; i < 12; ++i ) {
                tee.print( big ); // ~3 MB in total, comfortably past the 2 MB cap
            }
            tee.flush();
            final long size = Files.size( f );
            if ( size > ( 2 * ErrorLog.MAX_BYTES ) ) {
                return fail( "the log grew to " + size + " bytes, past twice the " + ErrorLog.MAX_BYTES
                        + "-byte cap" );
            }
            if ( !Files.readString( f, StandardCharsets.UTF_8 ).contains( "log restarted" ) ) {
                return fail( "restarting the log should say so" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "exception: " + e );
        }
    }

    /** The log lives beside the settings, so it follows the cache-dir override instead of the real home dir. */
    private static boolean location() {
        final String saved = System.getProperty( GuiPreferences.DIR_PROPERTY );
        try {
            System.setProperty( GuiPreferences.DIR_PROPERTY, "/tmp/aptx-log-location-probe" );
            final Path p = ErrorLog.defaultFile();
            if ( !p.toString().startsWith( "/tmp/aptx-log-location-probe" )
                    || !p.getFileName().toString().equals( ErrorLog.FILE_NAME ) ) {
                return fail( "the log should honour the cache-dir override, got " + p );
            }
            System.clearProperty( GuiPreferences.DIR_PROPERTY );
            if ( !ErrorLog.defaultFile().toString().contains( GuiPreferences.DEFAULT_DIR ) ) {
                return fail( "without an override the log belongs beside the settings, got "
                        + ErrorLog.defaultFile() );
            }
            return true;
        }
        finally {
            if ( saved != null ) {
                System.setProperty( GuiPreferences.DIR_PROPERTY, saved );
            }
            else {
                System.clearProperty( GuiPreferences.DIR_PROPERTY );
            }
        }
    }

    /**
     * The user-facing half: Help -> Show Error Log is always there, and the first logged error raises a quiet
     * menu-bar marker that can be dismissed.
     * <p>
     * Deliberately does NOT call {@link ErrorLog#install} -- that replaces {@code System.err}, sets the JVM-wide
     * uncaught handler and registers a shutdown hook, which would leak into every test that runs afterwards.
     */
    private static boolean menu() {
        if ( java.awt.GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            javax.swing.SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new org.forester.phylogeny.Phylogeny[ 0 ], new Configuration(), "errlog" ) );
            final boolean[] ok = { true };
            javax.swing.SwingUtilities.invokeAndWait( () -> {
                if ( helpItem( mf[ 0 ], "Show Error Log" ) == null ) {
                    ok[ 0 ] = fail( "Help should offer \"Show Error Log\" -- an installed user has no console, so "
                            + "this is the only way to reach a stack trace" );
                }
                final int before = mf[ 0 ].getJMenuBar().getMenuCount();
                mf[ 0 ].showErrorIndicator();
                final javax.swing.JMenu marker = menuByText( mf[ 0 ], "error logged" );
                if ( marker == null ) {
                    ok[ 0 ] = fail( "the first logged error should raise a menu-bar marker" );
                    return;
                }
                mf[ 0 ].showErrorIndicator(); // a repainting failure logs many times -- the marker must not stack up
                if ( mf[ 0 ].getJMenuBar().getMenuCount() != ( before + 1 ) ) {
                    ok[ 0 ] = fail( "the marker must appear ONCE however many errors are logged, menus went from "
                            + before + " to " + mf[ 0 ].getJMenuBar().getMenuCount() );
                }
                javax.swing.JMenuItem dismiss = null;
                for( int i = 0; i < marker.getItemCount(); ++i ) {
                    if ( ( marker.getItem( i ) != null ) && "Dismiss".equals( marker.getItem( i ).getText() ) ) {
                        dismiss = marker.getItem( i );
                    }
                }
                if ( dismiss == null ) {
                    ok[ 0 ] = fail( "the marker should be dismissable" );
                    return;
                }
                dismiss.doClick();
                if ( menuByText( mf[ 0 ], "error logged" ) != null ) {
                    ok[ 0 ] = fail( "dismissing the marker should remove it" );
                }
            } );
            javax.swing.SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Two things that would silently lose output: the tail of an unterminated stderr line, and the console.
     * <p>
     * Installing a default uncaught-exception handler SUPPRESSES the JVM's own "Exception in thread ..." print, so
     * without echoing, anyone running from a terminal (a developer, CI) would see nothing at all where they used to
     * see a stack trace.
     */
    private static boolean tailAndConsole() {
        try {
            final Path f = tmp( "tail.log" );
            final ErrorLog log = new ErrorLog( f, null );
            final java.io.ByteArrayOutputStream console = new java.io.ByteArrayOutputStream();
            final PrintStream original = new PrintStream( console, true, StandardCharsets.UTF_8 );
            final PrintStream tee = log.tee( original );
            tee.print( "a dying message with no newline" );
            tee.flush();
            log.close(); // what the shutdown hook runs
            if ( !Files.readString( f, StandardCharsets.UTF_8 ).contains( "a dying message with no newline" ) ) {
                return fail( "closing the log must write the last unterminated line, not drop it" );
            }
            // and a logged failure must still reach a real console
            final Path f2 = tmp( "console.log" );
            final ErrorLog log2 = new ErrorLog( f2, null );
            final java.io.ByteArrayOutputStream console2 = new java.io.ByteArrayOutputStream();
            log2.tee( new PrintStream( console2, true, StandardCharsets.UTF_8 ) ); // captures it as the console
            log2.logThrowable( "Uncaught in AWT-EventQueue-0", thrown( "visible please" ) );
            if ( !console2.toString( StandardCharsets.UTF_8 ).contains( "visible please" ) ) {
                return fail( "a logged failure must still be echoed to the console: "
                        + console2.toString( StandardCharsets.UTF_8 ) );
            }
            if ( !Files.readString( f2, StandardCharsets.UTF_8 ).contains( "visible please" ) ) {
                return fail( "...and still be written to the file" );
            }
            // exactly once in the file, though -- the echo must not come back through the tee
            final String logged = Files.readString( f2, StandardCharsets.UTF_8 );
            if ( ( logged.split( "visible please", -1 ).length - 1 ) > 2 ) {
                return fail( "the console echo must not be fed back through the tee and logged twice" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "exception: " + e );
        }
    }

    /**
     * A malformed {@code archaeopteryx.cache.dir} must not stop the application from starting.
     * <p>
     * {@code defaultFile()} builds a {@link Path} from that property, and {@code Paths.get} throws the UNCHECKED
     * {@code InvalidPathException} on a bad value. Thrown from {@code install()} -- which runs before anything else
     * in {@code main} -- it would propagate straight out and the user would get no window and no message at all: a
     * logger that prevents the program it is meant to diagnose from launching.
     */
    private static boolean badPropertyCannotBlockStartup() {
        final String saved = System.getProperty( GuiPreferences.DIR_PROPERTY );
        try {
            System.setProperty( GuiPreferences.DIR_PROPERTY, "bad\u0000path" ); // NUL is invalid on every platform
            ErrorLog.resetForTest();
            final ErrorLog log = ErrorLog.install( null ); // must NOT throw
            if ( log != null ) {
                return fail( "an unusable log directory should yield no logger, not a half-built one" );
            }
            if ( ErrorLog.instance() != null ) {
                return fail( "a failed install must not leave an instance behind" );
            }
            return true;
        }
        catch ( final Throwable t ) {
            return fail( "install() must never throw -- a bad cache dir would stop the app launching: " + t );
        }
        finally {
            ErrorLog.resetForTest();
            if ( saved != null ) {
                System.setProperty( GuiPreferences.DIR_PROPERTY, saved );
            }
            else {
                System.clearProperty( GuiPreferences.DIR_PROPERTY );
            }
        }
    }

    /**
     * The log is installed BEFORE any window exists, so an error thrown during start-up fires the one-shot
     * callback when there is no menu bar to mark. A frame built afterwards must still raise the marker, or that
     * error is never surfaced at all -- which would defeat the point of installing the log first.
     */
    private static boolean markerAfterStartupError() {
        if ( java.awt.GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Path f = tmp( "startup.log" );
            final ErrorLog log = new ErrorLog( f, null );
            log.logThrowable( "during start-up, before any window", thrown( "early boom" ) );
            ErrorLog.setInstanceForTest( log );
            final boolean[] ok = { true };
            final MainFrame[] mf = new MainFrame[ 1 ];
            try {
                javax.swing.SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                        .createInstance( new org.forester.phylogeny.Phylogeny[ 0 ], new Configuration(), "early" ) );
                javax.swing.SwingUtilities.invokeAndWait( () -> {
                    if ( menuByText( mf[ 0 ], "error logged" ) == null ) {
                        ok[ 0 ] = fail( "a frame opened AFTER a start-up error must still show the marker" );
                        return;
                    }
                    // and dismissing must not make it unrecoverable for the rest of the session
                    final javax.swing.JMenu marker = menuByText( mf[ 0 ], "error logged" );
                    for( int i = 0; i < marker.getItemCount(); ++i ) {
                        if ( ( marker.getItem( i ) != null ) && "Dismiss".equals( marker.getItem( i ).getText() ) ) {
                            marker.getItem( i ).doClick();
                        }
                    }
                    mf[ 0 ].showErrorIndicator();
                    if ( menuByText( mf[ 0 ], "error logged" ) == null ) {
                        ok[ 0 ] = fail( "after dismissing, a LATER error must be able to raise the marker again" );
                    }
                } );
                javax.swing.SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            }
            finally {
                ErrorLog.resetForTest();
            }
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            ErrorLog.resetForTest();
            e.printStackTrace();
            return false;
        }
    }

    private static javax.swing.JMenuItem helpItem( final MainFrame mf, final String text ) {
        final javax.swing.JMenuBar bar = mf.getJMenuBar();
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final javax.swing.JMenu m = bar.getMenu( i );
            if ( ( m == null ) || !"Help".equals( m.getText() ) ) {
                continue;
            }
            for( int j = 0; j < m.getItemCount(); ++j ) {
                if ( ( m.getItem( j ) != null ) && text.equals( m.getItem( j ).getText() ) ) {
                    return m.getItem( j );
                }
            }
        }
        return null;
    }

    private static javax.swing.JMenu menuByText( final MainFrame mf, final String contains ) {
        final javax.swing.JMenuBar bar = mf.getJMenuBar();
        for( int i = 0; i < bar.getMenuCount(); ++i ) {
            final javax.swing.JMenu m = bar.getMenu( i );
            if ( ( m != null ) && ( m.getText() != null ) && m.getText().contains( contains ) ) {
                return m;
            }
        }
        return null;
    }

    private ErrorLogTest() {
    }
}

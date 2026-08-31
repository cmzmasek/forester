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

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

/**
 * Where a crash goes when nobody is watching a terminal.
 * <p>
 * An INSTALLED Archaeopteryx has no console: the Windows launcher jpackage builds is a windowed one (no console at
 * all), a macOS {@code .app} hands its stdio to launchd, and a Linux desktop entry sends it to the session journal.
 * So a stack trace -- from an uncaught exception or from any of the {@code printStackTrace()} calls scattered
 * through the code -- simply vanished, and a user could only report "it stopped drawing". This writes both to
 * {@code archaeopteryx.log} beside the persisted settings, where a user can find it and attach it to a report.
 * <p>
 * Two things make it safe to leave switched on:
 * <ul>
 * <li><b>Repeats are counted, not rewritten.</b> An exception thrown inside {@code paintComponent} recurs on every
 * repaint -- thousands of identical traces within seconds. The first is written in full; the rest increment a
 * counter that is flushed as one "repeated N more times" line when something else happens or the app exits.</li>
 * <li><b>The file is capped.</b> Past {@link #MAX_BYTES} it is restarted with a note, so a runaway loop cannot fill
 * the user's disk.</li>
 * </ul>
 * The pure parts (signature, de-duplication, the cap) are headless-testable; installation touches global JVM state
 * and is done once from the application entry point.
 */
final class ErrorLog {

    static final String FILE_NAME = "archaeopteryx.log";
    /** Past this the log restarts (with a note). A runaway paint exception must not fill the disk. */
    static final long   MAX_BYTES = 2L * 1024 * 1024;

    private static final DateTimeFormatter STAMP = DateTimeFormatter.ofPattern( "yyyy-MM-dd HH:mm:ss" );
    private static ErrorLog                INSTANCE;

    private final Path             _file;
    private final Runnable         _on_first_error;
    /** The real stderr, captured before the tee replaced it. Traces are echoed here so a developer running from a
     *  terminal still SEES them: installing a default uncaught-exception handler suppresses the JVM's own
     *  "Exception in thread ..." output, which would otherwise make the console go quiet. */
    private volatile PrintStream   _console;
    private String                 _last_signature;
    private int                    _repeats;
    private boolean                _error_logged;
    private java.io.Writer         _out;
    private long                   _written_bytes;
    /** Line buffer for the System.err tee: a PrintStream can call through one BYTE at a time, and a file write per
     *  byte would be both absurdly slow and unreadable. Flushed on newline. */
    private final StringBuilder    _err_buffer = new StringBuilder();

    ErrorLog( final Path file, final Runnable on_first_error ) {
        _file = file;
        _on_first_error = on_first_error;
    }

    /** The log file: beside the persisted settings, so it moves with the {@code archaeopteryx.cache.dir} override
     *  (which is what keeps a test run out of the user's real directory). */
    static Path defaultFile() {
        final String override = System.getProperty( GuiPreferences.DIR_PROPERTY );
        final Path dir = ( ( override != null ) && !override.trim().isEmpty() )
                ? Paths.get( override.trim() )
                : Paths.get( System.getProperty( "user.home", "." ), GuiPreferences.DEFAULT_DIR );
        return dir.resolve( FILE_NAME );
    }

    /**
     * Installs the global handler and tees {@code System.err} into the log. Idempotent; never throws -- a logger
     * that brings down the application it is meant to diagnose would be worse than no logger.
     *
     * @param on_first_error run (once) the first time anything is logged, so the UI can offer the log to the user
     */
    static synchronized ErrorLog install( final Runnable on_first_error ) {
        if ( INSTANCE != null ) {
            return INSTANCE;
        }
        final ErrorLog log;
        try {
            // defaultFile() reads a system property and builds a Path: a malformed archaeopteryx.cache.dir throws
            // InvalidPathException, and that must not be what stops the application from starting.
            log = new ErrorLog( defaultFile(), on_first_error );
        }
        catch ( final Throwable t ) {
            return null;
        }
        INSTANCE = log;
        try {
            // Catches EDT throws too (verified): an uncaught exception on AWT-EventQueue-0 reaches the default
            // handler, which is exactly the paint-loop case this exists for.
            Thread.setDefaultUncaughtExceptionHandler( ( t, e ) -> log.logThrowable( "Uncaught in " + t.getName(), e ) );
            System.setErr( log.tee( System.err ) );
            // a repainting failure is still being counted when the user quits -- write that tail out
            Runtime.getRuntime().addShutdownHook( new Thread( log::close, "aptx-errorlog-close" ) );
        }
        catch ( final Throwable ignored ) {
            // a logger must never be the reason the app fails to start
        }
        return log;
    }

    /** The installed log, or null when it has not been installed (tests, the CLI). */
    static synchronized ErrorLog instance() {
        return INSTANCE;
    }

    /** For tests: forget the installed instance. */
    static synchronized void resetForTest() {
        INSTANCE = null;
    }

    /** For tests: stand in as the installed instance WITHOUT touching System.err, the default handler or the
     *  shutdown hook -- installing those for real would leak into every test that runs afterwards. */
    static synchronized void setInstanceForTest( final ErrorLog log ) {
        INSTANCE = log;
    }

    Path getFile() {
        return _file;
    }

    /** Whether anything has been written -- what the UI uses to decide there is something worth offering. */
    synchronized boolean hasEntries() {
        return _error_logged || Files.exists( _file );
    }

    /** Whether a real FAILURE (not just stderr chatter) has been logged this session. */
    synchronized boolean hasErrors() {
        return _error_logged;
    }

    /**
     * The de-duplication key for a throwable: its type, message and TOP STACK FRAME. Deliberately not the whole
     * trace -- the same paint bug reached through different repaint call stacks is still the same bug, and writing
     * it out again per repaint is what would drown the log.
     */
    static String signature( final Throwable t ) {
        if ( t == null ) {
            return "null";
        }
        final StackTraceElement[] st = t.getStackTrace();
        return t.getClass().getName() + "|" + t.getMessage() + "|" + ( ( st.length > 0 ) ? st[ 0 ].toString() : "" );
    }

    synchronized void logThrowable( final String context, final Throwable t ) {
        final String sig = signature( t );
        if ( sig.equals( _last_signature ) ) {
            ++_repeats; // same failure again (a repainting loop) -- count it, do not rewrite the trace
            return;
        }
        flushRepeats();
        _last_signature = sig;
        final StringBuilder sb = new StringBuilder();
        sb.append( LocalDateTime.now().format( STAMP ) ).append( "  " ).append( context ).append( '\n' );
        final java.io.StringWriter sw = new java.io.StringWriter();
        t.printStackTrace( new java.io.PrintWriter( sw ) );
        sb.append( sw ).append( '\n' );
        write( sb.toString(), true );
        // ...and to the real console, if there is one. Writing to System.err would come straight back through the
        // tee and be logged twice, so this goes to the stream captured before the tee was installed.
        final PrintStream console = _console;
        if ( console != null ) {
            console.print( sb );
        }
    }

    /** Writes the pending "repeated N more times" line, if any. */
    synchronized void flushRepeats() {
        if ( _repeats > 0 ) {
            final int n = _repeats;
            _repeats = 0;
            write( "    ... the same failure repeated " + n + " more time" + ( n == 1 ? "" : "s" ) + '\n', false );
        }
    }

    /**
     * @param is_error whether this is a real failure. Only a failure raises the "there is something to send"
     *                 signal -- plain stderr chatter (a JVM or library warning) is worth capturing but is not
     *                 worth telling the user their session went wrong.
     */
    private synchronized void write( final String text, final boolean is_error ) {
        try {
            // The writer is opened once and kept open, flushed per write. Opening and closing the file per line
            // (Files.writeString with APPEND) is ruinous when stderr is chatty -- measured at minutes for a few MB.
            if ( _out == null ) {
                if ( _file.getParent() != null ) {
                    Files.createDirectories( _file.getParent() );
                }
                _written_bytes = Files.exists( _file ) ? Files.size( _file ) : 0;
                _out = Files.newBufferedWriter( _file, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
                                                StandardOpenOption.APPEND );
            }
            if ( _written_bytes > MAX_BYTES ) {
                // restart rather than grow without bound: a repainting failure could otherwise fill the disk
                _out.close();
                _out = Files.newBufferedWriter( _file, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
                                                StandardOpenOption.TRUNCATE_EXISTING );
                final String note = LocalDateTime.now().format( STAMP ) + "  [log restarted: it had grown past "
                        + ( MAX_BYTES / ( 1024 * 1024 ) ) + " MB]\n";
                _out.write( note );
                _written_bytes = note.getBytes( StandardCharsets.UTF_8 ).length;
            }
            _out.write( text );
            _out.flush(); // a crash must not take the last entry with it
            // bytes, not chars: the cap is a byte budget and the seed came from Files.size(), so counting UTF-16
            // chars would let a log full of non-ASCII text overshoot it severalfold
            _written_bytes += text.getBytes( StandardCharsets.UTF_8 ).length;
        }
        catch ( final IOException ignored ) {
            _out = null; // try again next time; writing to System.err here would recurse straight back in
            return;
        }
        if ( is_error ) {
            final boolean first = !_error_logged;
            _error_logged = true;
            if ( first && ( _on_first_error != null ) ) {
                try {
                    _on_first_error.run();
                }
                catch ( final Throwable ignored ) {
                    // the notification must never take the application down
                }
            }
        }
    }

    /** Flushes any pending repeat count and closes the file. Called from a shutdown hook by {@link #install}. */
    synchronized void close() {
        flushErrBuffer(); // a final System.err.print with no trailing newline would otherwise be dropped
        flushRepeats();
        try {
            if ( _out != null ) {
                _out.close();
            }
        }
        catch ( final IOException ignored ) {
            // nothing useful to do while shutting down
        }
        _out = null;
    }

    /**
     * A {@link PrintStream} that writes to {@code original} AND to the log, so the {@code printStackTrace()} calls
     * already scattered through the code end up in the file too, without touching all 186 of them.
     */
    PrintStream tee( final PrintStream original ) {
        _console = original; // the real stderr: traces are echoed here so a terminal still shows them
        final OutputStream both = new OutputStream() {

            @Override
            public void write( final int b ) {
                original.write( b );
                writeRaw( new byte[] { (byte) b }, 0, 1 );
            }

            @Override
            public void write( final byte[] b, final int off, final int len ) {
                original.write( b, off, len );
                writeRaw( b, off, len );
            }

            @Override
            public void flush() {
                original.flush();
            }
        };
        return new PrintStream( both, true, StandardCharsets.UTF_8 );
    }

    /** Largest unterminated stderr run held in memory. A writer that never emits a newline must not be able to
     *  grow the buffer without bound (the file has a cap; this needs one too). */
    private static final int ERR_BUFFER_MAX = 64 * 1024;

    /** Buffers tee'd stderr and writes it a LINE at a time (see {@link #_err_buffer}). */
    private synchronized void writeRaw( final byte[] b, final int off, final int len ) {
        _err_buffer.append( new String( b, off, len, StandardCharsets.UTF_8 ) );
        int nl;
        while ( ( nl = _err_buffer.indexOf( "\n" ) ) >= 0 ) {
            final String line = _err_buffer.substring( 0, nl + 1 );
            _err_buffer.delete( 0, nl + 1 );
            write( line, false ); // stderr text is captured, but it does not by itself mean "an error happened"
        }
        if ( _err_buffer.length() > ERR_BUFFER_MAX ) {
            flushErrBuffer(); // no newline in sight -- write what we have rather than hoard it
        }
    }

    /** Writes whatever partial line is buffered. */
    private synchronized void flushErrBuffer() {
        if ( _err_buffer.length() > 0 ) {
            final String rest = _err_buffer.toString();
            _err_buffer.setLength( 0 );
            write( rest, false );
        }
    }

    /** For tests: how many repeats are pending (not yet summarised into a line). */
    synchronized int pendingRepeatsForTest() {
        return _repeats;
    }
}

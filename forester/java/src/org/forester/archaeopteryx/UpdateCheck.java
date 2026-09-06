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

import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URI;
import java.net.URLConnection;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.SwingUtilities;

/**
 * The launch-time "is there a newer Archaeopteryx?" check. One request to the GitHub API for the LATEST official
 * (non-prerelease) release of the Archaeopteryx home repository, on a daemon thread, a little after the window is
 * up; if its version is newer than the running one, {@link MainFrame} is told on the EDT and adds a line to its
 * Help menu. Everything else -- no network, GitHub down, a proxy in the way, an unexpected answer, a bug in here --
 * ends silently: this class never throws, never logs, never prints. The parsing and comparison are pure statics so
 * they can be tested without a network.
 */
final class UpdateCheck {

    /** The latest OFFICIAL release of the home repository (this endpoint ignores prereleases, which is what the
     *  forester repository publishes constantly). */
    static final String  LATEST_RELEASE_API = "https://api.github.com/repos/cmzmasek/archaeopteryx/releases/latest";
    /** Where the Help-menu line takes the user. */
    static final String  RELEASES_PAGE      = "https://github.com/cmzmasek/archaeopteryx/releases";
    static final int     TIMEOUT_MS         = 5000;
    /** Let the window come up and the first tree paint before the network is touched. */
    static final long    START_DELAY_MS     = 1500;
    private static final int     MAX_BODY_BYTES = 64 * 1024;
    private static final Pattern TAG_NAME       = Pattern.compile( "\"tag_name\"\\s*:\\s*\"([^\"]*)\"" );
    private static final AtomicBoolean STARTED  = new AtomicBoolean( false );

    private UpdateCheck() {
        // not instantiable
    }

    /**
     * Starts the launch-time check, at most ONCE per JVM (a second window does not check again). {@code on_newer}
     * is called on the EDT with the newer version's number, only if there is one.
     */
    static void startAtLaunch( final String current_version, final Consumer<String> on_newer ) {
        if ( !STARTED.compareAndSet( false, true ) ) {
            return;
        }
        start( LATEST_RELEASE_API, current_version, START_DELAY_MS, on_newer );
    }

    /**
     * The check itself, parameterised for tests: fetch {@code api_url}, take its {@code tag_name}, compare with
     * {@code current_version}, and call {@code on_newer} on the EDT if it is newer. Returns the (daemon) thread.
     * Swallows EVERYTHING.
     */
    static Thread start( final String api_url, final String current_version, final long delay_ms,
                         final Consumer<String> on_newer ) {
        final Thread t = new Thread( () -> {
            try {
                if ( delay_ms > 0 ) {
                    Thread.sleep( delay_ms );
                }
                final String newer = newerVersionAt( api_url, current_version );
                if ( newer != null ) {
                    SwingUtilities.invokeLater( () -> {
                        try {
                            on_newer.accept( newer );
                        }
                        catch ( final Throwable ignore ) {
                            // the caller's problem is not allowed to become the user's
                        }
                    } );
                }
            }
            catch ( final Throwable ignore ) {
                // by design: no network, no GitHub, no answer -> nothing happens, nothing is said
            }
        }, "aptx-update-check" );
        t.setDaemon( true );
        t.setPriority( Thread.MIN_PRIORITY );
        t.start();
        return t;
    }

    /** The version at {@code api_url} if it is newer than {@code current_version}, else null. Never throws. */
    static String newerVersionAt( final String api_url, final String current_version ) {
        try {
            final String body = fetch( api_url );
            if ( body == null ) {
                return null;
            }
            final String latest = normalize( latestTagFromJson( body ) );
            return ( ( latest != null ) && isNewer( latest, current_version ) ) ? latest : null;
        }
        catch ( final Throwable ignore ) {
            return null;
        }
    }

    /** The body of {@code url} as text, or null on ANY problem (never throws). Short timeouts, capped size. */
    static String fetch( final String url ) {
        try {
            final URLConnection c = URI.create( url ).toURL().openConnection();
            c.setConnectTimeout( TIMEOUT_MS );
            c.setReadTimeout( TIMEOUT_MS );
            c.setUseCaches( false );
            // GitHub's API refuses requests without a User-Agent
            c.setRequestProperty( "User-Agent", "Archaeopteryx/" + AptxConstants.VERSION );
            c.setRequestProperty( "Accept", "application/vnd.github+json" );
            if ( c instanceof HttpURLConnection ) {
                ( (HttpURLConnection) c ).setInstanceFollowRedirects( true );
            }
            try ( InputStream in = c.getInputStream() ) {
                final byte[] buf = new byte[ 4096 ];
                final java.io.ByteArrayOutputStream out = new java.io.ByteArrayOutputStream();
                int n;
                while ( ( n = in.read( buf ) ) > 0 ) {
                    out.write( buf, 0, n );
                    if ( out.size() > MAX_BODY_BYTES ) {
                        return null; // not the small JSON document this expects
                    }
                }
                return out.toString( StandardCharsets.UTF_8 );
            }
            finally {
                if ( c instanceof HttpURLConnection ) {
                    ( (HttpURLConnection) c ).disconnect();
                }
            }
        }
        catch ( final Throwable ignore ) {
            return null;
        }
    }

    /** The {@code "tag_name"} value in a GitHub release JSON document, or null. */
    static String latestTagFromJson( final String json ) {
        if ( json == null ) {
            return null;
        }
        final Matcher m = TAG_NAME.matcher( json );
        return m.find() ? m.group( 1 ) : null;
    }

    /** "v0.11.140" / " 0.11.140 " -> "0.11.140"; null/blank -> null. */
    static String normalize( final String tag ) {
        if ( tag == null ) {
            return null;
        }
        String s = tag.trim();
        if ( s.startsWith( "v" ) || s.startsWith( "V" ) ) {
            s = s.substring( 1 );
        }
        return s.isEmpty() ? null : s;
    }

    /**
     * Whether dotted version {@code candidate} is strictly newer than {@code current}: numeric, component by
     * component, missing components count as 0 ("0.12" == "0.12.0"). A candidate that is not a plain dotted number
     * is never "newer" (a typo'd tag must not nag every user).
     */
    static boolean isNewer( final String candidate, final String current ) {
        final int[] a = parts( normalize( candidate ) );
        final int[] b = parts( normalize( current ) );
        if ( ( a == null ) || ( b == null ) ) {
            return false;
        }
        final int n = Math.max( a.length, b.length );
        for( int i = 0; i < n; ++i ) {
            final int x = ( i < a.length ) ? a[ i ] : 0;
            final int y = ( i < b.length ) ? b[ i ] : 0;
            if ( x != y ) {
                return x > y;
            }
        }
        return false;
    }

    private static int[] parts( final String v ) {
        if ( ( v == null ) || v.isEmpty() ) {
            return null;
        }
        final String[] ss = v.split( "\\." );
        final int[] out = new int[ ss.length ];
        for( int i = 0; i < ss.length; ++i ) {
            if ( !ss[ i ].matches( "\\d{1,6}" ) ) {
                return null;
            }
            out[ i ] = Integer.parseInt( ss[ i ] );
        }
        return out;
    }

    /** For tests only: forget that the launch check ran. */
    static void resetForTest() {
        STARTED.set( false );
    }
}

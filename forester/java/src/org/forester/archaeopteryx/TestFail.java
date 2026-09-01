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

/**
 * Says WHERE a GUI test failed.
 * <p>
 * Many of the headful tests signalled failure with a bare {@code ok[0] = false;} and no message, so the suite
 * printed "failed." and nothing else. That is survivable for a deterministic failure -- run it again and debug it
 * -- but it made the suite's INTERMITTENT failure effectively unpinnable: it fired roughly once in ten runs,
 * always in a different test, and never said which assertion had tripped.
 * <p>
 * {@link #here()} returns {@code false} (so it drops straight into {@code ok[0] = ...}) and prints the file and
 * line of the caller, taken from the stack rather than hand-written, so it cannot drift as the file is edited.
 */
final class TestFail {

    /** Returns false, having printed the calling test's class, method and line. */
    static boolean here() {
        return here( null );
    }

    /** As {@link #here()}, with an extra note. */
    static boolean here( final String note ) {
        final StackTraceElement[] st = new Throwable().getStackTrace();
        // [0] is this method; walk out of TestFail to the first frame that is not us
        StackTraceElement caller = null;
        for( final StackTraceElement e : st ) {
            if ( !TestFail.class.getName().equals( e.getClassName() ) ) {
                caller = e;
                break;
            }
        }
        final String where = ( caller == null ) ? "unknown location"
                : ( simpleName( caller.getClassName() ) + "." + caller.getMethodName() + "("
                        + caller.getFileName() + ":" + caller.getLineNumber() + ")" );
        System.out.println( "  [FAILED at] " + where + ( ( note == null ) ? "" : ( " -- " + note ) ) );
        return false;
    }

    private static String simpleName( final String fqn ) {
        final int dot = fqn.lastIndexOf( '.' );
        return ( dot < 0 ) ? fqn : fqn.substring( dot + 1 );
    }

    private TestFail() {
    }
}

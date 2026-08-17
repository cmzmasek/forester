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

import java.util.List;
import java.util.Map;

/**
 * Headless tests for {@link TanglegramAssociation}: TSV and CSV parsing, quoted cells, blank/short-row skipping,
 * many:many mappings, CR/LF endings, and empty input.
 */
public final class TanglegramAssociationTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramAssociation: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return tsvOk() && csvAndQuotesOk() && skipsAndManyToManyOk() && emptyOk();
    }

    private static boolean tsvOk() {
        // a leading UTF-8 BOM (as Excel / Sheets prepend) must NOT contaminate the first key
        final TanglegramAssociation a = TanglegramAssociation.parse( "\uFEFF" + "hostA\tlouseA\nhostB\tlouseB\n" );
        if ( ( a.pairCount() != 2 ) || ( a.leftKeyCount() != 2 ) ) {
            return fail( "TSV should parse 2 pairs / 2 left keys, got " + a.pairCount() + "/" + a.leftKeyCount() );
        }
        final Map<String, List<String>> m = a.leftToRight();
        if ( ( m.get( "hostA" ) == null ) || !m.get( "hostA" ).equals( java.util.Arrays.asList( "louseA" ) )
                || !m.get( "hostB" ).equals( java.util.Arrays.asList( "louseB" ) ) ) {
            return fail( "TSV mapping wrong (BOM not stripped?): " + m );
        }
        return true;
    }

    private static boolean csvAndQuotesOk() {
        // comma-delimited, surrounding quotes, extra whitespace, a CR/LF ending
        final TanglegramAssociation a = TanglegramAssociation
                .parse( "\"host A\" , \"louse A\"\r\n host B ,louse B\r\n" );
        if ( a.pairCount() != 2 ) {
            return fail( "CSV should parse 2 pairs, got " + a.pairCount() );
        }
        final Map<String, List<String>> m = a.leftToRight();
        if ( ( m.get( "host A" ) == null ) || !m.get( "host A" ).get( 0 ).equals( "louse A" ) ) {
            return fail( "quoted/space CSV cell not cleaned: " + m );
        }
        if ( ( m.get( "host B" ) == null ) || !m.get( "host B" ).get( 0 ).equals( "louse B" ) ) {
            return fail( "trimmed CSV cell wrong: " + m );
        }
        return true;
    }

    private static boolean skipsAndManyToManyOk() {
        // blank line, a one-column short row (skipped), a blank cell (skipped), a duplicate exact pair (deduped),
        // and a genuine many:many (hostA -> louseA AND louseB)
        final TanglegramAssociation a = TanglegramAssociation.parse( "hostA\tlouseA\n\nhostA\tlouseB\nonly_one_column\n"
                + "hostB\t\nhostA\tlouseA\n" );
        if ( a.pairCount() != 2 ) {
            return fail( "expected 2 distinct pairs (hostA->louseA, hostA->louseB; the exact-dup, blank and short "
                    + "rows dropped), got " + a.pairCount() );
        }
        final List<String> rights = a.leftToRight().get( "hostA" );
        if ( ( rights.size() != 2 ) || !rights.contains( "louseA" ) || !rights.contains( "louseB" ) ) {
            return fail( "hostA should map to exactly {louseA, louseB}, got " + rights );
        }
        if ( a.leftToRight().containsKey( "hostB" ) ) {
            return fail( "a row with a blank right cell must be skipped (hostB should be absent)" );
        }
        return true;
    }

    private static boolean emptyOk() {
        if ( ( TanglegramAssociation.parse( "" ).pairCount() != 0 )
                || ( TanglegramAssociation.parse( null ).pairCount() != 0 )
                || ( TanglegramAssociation.parse( "\n \n\t\n" ).pairCount() != 0 ) ) {
            return fail( "empty / whitespace-only / null input should yield 0 pairs" );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramAssociation test failed: " + message );
        return false;
    }
}

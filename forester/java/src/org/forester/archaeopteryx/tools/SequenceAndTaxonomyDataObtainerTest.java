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

package org.forester.archaeopteryx.tools;

import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Unit test for the pure report-building logic of {@link SequenceAndTaxonomyDataObtainer}
 * ({@code hasIssues} / {@code buildCompletionMessage}). These need no network or GUI, so they
 * run in the headless suite; the actual web-service phases are exercised only against live
 * UniProt/EMBL-GenBank and are not unit-tested here.
 */
public final class SequenceAndTaxonomyDataObtainerTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SequenceAndTaxonomyDataObtainer: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // 1. nothing wrong: not an issue, and a plain success message (both null and both empty).
        if ( SequenceAndTaxonomyDataObtainer.hasIssues( null, null, null, null ) ) {
            return fail( "all-null must not be an issue" );
        }
        if ( SequenceAndTaxonomyDataObtainer.hasIssues( empty(), null, empty(), null ) ) {
            return fail( "empty not-found sets must not be an issue" );
        }
        final String success = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( empty(), null, empty(), null );
        if ( !success.equals( "Sequence and taxonomy data successfully obtained." ) ) {
            return fail( "unexpected success message: " + success );
        }

        // 2. a sequence error is an issue and is reported (taxonomy section absent when clean).
        if ( !SequenceAndTaxonomyDataObtainer.hasIssues( null, "boom", null, null ) ) {
            return fail( "a sequence error must be an issue" );
        }
        final String se = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( null, "no connection", empty(), null );
        if ( !se.contains( "Sequence data: error - no connection" ) ) {
            return fail( "sequence error not reported: " + se );
        }
        if ( se.contains( "Taxonomy data" ) ) {
            return fail( "clean taxonomy section must be absent: " + se );
        }

        // 3. a taxonomy error is an issue and is reported.
        if ( !SequenceAndTaxonomyDataObtainer.hasIssues( null, null, null, "tax boom" ) ) {
            return fail( "a taxonomy error must be an issue" );
        }
        final String te = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( empty(), null, null, "tax down" );
        if ( !te.contains( "Taxonomy data: error - tax down" ) ) {
            return fail( "taxonomy error not reported: " + te );
        }

        // 4. a single not-found node uses the singular wording and lists the node.
        final String one = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( set( "node_x" ), null, empty(), null );
        if ( !SequenceAndTaxonomyDataObtainer.hasIssues( set( "node_x" ), null, null, null ) ) {
            return fail( "a single not-found node must be an issue" );
        }
        if ( !one.contains( "following node:" ) || !one.contains( "node_x" ) ) {
            return fail( "single not-found wording/listing wrong: " + one );
        }
        if ( one.contains( "total:" ) ) {
            return fail( "single not-found must not print a total: " + one );
        }

        // 5. multiple not-found nodes print a total and list each.
        final SortedSet<String> three = set( "a", "b", "c" );
        final String multi = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( three, null, empty(), null );
        if ( !multi.contains( "total: 3" ) || !multi.contains( "a" ) || !multi.contains( "b" )
                || !multi.contains( "c" ) ) {
            return fail( "multi not-found wording/listing wrong: " + multi );
        }

        // 6. a long not-found list is truncated to 20 entries with an ellipsis, total still exact.
        final SortedSet<String> many = new TreeSet<>();
        for( int i = 0; i < 25; i++ ) {
            many.add( String.format( "n%02d", i ) );
        }
        final String trunc = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( empty(), null, many, null );
        if ( !trunc.contains( "total: 25" ) || !trunc.contains( "..." ) ) {
            return fail( "long list not truncated/totalled: " + trunc );
        }
        if ( trunc.contains( "n20" ) || trunc.contains( "n24" ) ) {
            return fail( "truncated list must not include the 21st+ entries: " + trunc );
        }
        if ( !trunc.contains( "n00" ) || !trunc.contains( "n19" ) ) {
            return fail( "truncated list must include the first 20 entries: " + trunc );
        }

        // 7. issues in BOTH phases produce BOTH sections in one message.
        final String both = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( set( "seqnode" ), null,
                                                                                    null, "tax error" );
        if ( !both.contains( "Sequence data" ) || !both.contains( "seqnode" )
                || !both.contains( "Taxonomy data: error - tax error" ) ) {
            return fail( "combined message missing a section: " + both );
        }
        return true;
    }

    private static SortedSet<String> empty() {
        return new TreeSet<>();
    }

    private static SortedSet<String> set( final String... values ) {
        final SortedSet<String> s = new TreeSet<>();
        for( final String v : values ) {
            s.add( v );
        }
        return s;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "\nSequenceAndTaxonomyDataObtainerTest failed: " + msg );
        return false;
    }

    private SequenceAndTaxonomyDataObtainerTest() {
    }
}

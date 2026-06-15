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

import org.forester.ws.seqdb.SequenceTaxonomyResolver;

/**
 * Unit test for the pure report/commit logic of {@link SequenceAndTaxonomyDataObtainer}
 * ({@code shouldCommit} / {@code hasIssues} / {@code buildCompletionMessage}) over the new
 * {@link SequenceTaxonomyResolver.Result}. No network or GUI, so it runs in the headless suite; the
 * actual web-service resolution is exercised only against live UniProt/NCBI.
 */
public final class SequenceAndTaxonomyDataObtainerTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "SequenceAndTaxonomyDataObtainer: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // 1. clean success: not an issue, plain success message, and it commits (something was written).
        final SequenceTaxonomyResolver.Result clean = result( empty(), 3, 2, false, null );
        if ( SequenceAndTaxonomyDataObtainer.hasIssues( clean ) ) {
            return fail( "a clean run must not be an issue" );
        }
        if ( !SequenceAndTaxonomyDataObtainer.shouldCommit( clean ) ) {
            return fail( "a run that wrote data must commit" );
        }
        final String ok_msg = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( clean );
        if ( !ok_msg.startsWith( "Sequence and taxonomy data successfully obtained" )
                || !ok_msg.contains( "3 sequence(s)" ) || !ok_msg.contains( "2 taxonomy" ) ) {
            return fail( "unexpected success message: " + ok_msg );
        }

        // 2. nothing written -> do not commit (leave the tree + its edited-state untouched).
        if ( SequenceAndTaxonomyDataObtainer.shouldCommit( result( empty(), 0, 0, false, null ) ) ) {
            return fail( "a run that wrote nothing must not commit" );
        }

        // 3. a transport error is an issue and is reported.
        final SequenceTaxonomyResolver.Result err = result( empty(), 1, 0, false, "connection reset" );
        if ( !SequenceAndTaxonomyDataObtainer.hasIssues( err ) ) {
            return fail( "an error must be an issue" );
        }
        if ( !SequenceAndTaxonomyDataObtainer.buildCompletionMessage( err ).contains( "Error: connection reset" ) ) {
            return fail( "error not reported" );
        }

        // 4. cancellation is an issue and is reported; a partial write still commits.
        final SequenceTaxonomyResolver.Result cancelled = result( empty(), 5, 5, true, null );
        if ( !SequenceAndTaxonomyDataObtainer.hasIssues( cancelled )
                || !SequenceAndTaxonomyDataObtainer.shouldCommit( cancelled ) ) {
            return fail( "a cancelled-but-partial run is an issue but must still commit" );
        }
        if ( !SequenceAndTaxonomyDataObtainer.buildCompletionMessage( cancelled ).contains( "Cancelled" ) ) {
            return fail( "cancellation not reported" );
        }

        // 5. a single not-found node: singular wording, listed, no total.
        final String one = SequenceAndTaxonomyDataObtainer
                .buildCompletionMessage( result( set( "node_x" ), 0, 1, false, null ) );
        if ( !one.contains( "following node:" ) || !one.contains( "node_x" ) || one.contains( "total:" ) ) {
            return fail( "single not-found wording/listing wrong: " + one );
        }

        // 6. multiple not-found nodes: a total and each listed.
        final String multi = SequenceAndTaxonomyDataObtainer
                .buildCompletionMessage( result( set( "a", "b", "c" ), 0, 0, false, null ) );
        if ( !multi.contains( "total: 3" ) || !multi.contains( "a" ) || !multi.contains( "b" )
                || !multi.contains( "c" ) ) {
            return fail( "multi not-found wording/listing wrong: " + multi );
        }

        // 7. a long not-found list is truncated to 20 with an ellipsis; total stays exact.
        final SortedSet<String> many = new TreeSet<>();
        for( int i = 0; i < 25; i++ ) {
            many.add( String.format( "n%02d", i ) );
        }
        final String trunc = SequenceAndTaxonomyDataObtainer.buildCompletionMessage( result( many, 0, 0, false, null ) );
        if ( !trunc.contains( "total: 25" ) || !trunc.contains( "..." ) ) {
            return fail( "long list not truncated/totalled: " + trunc );
        }
        if ( trunc.contains( "n20" ) || trunc.contains( "n24" ) ) {
            return fail( "truncated list must not include the 21st+ entries: " + trunc );
        }
        if ( !trunc.contains( "n00" ) || !trunc.contains( "n19" ) ) {
            return fail( "truncated list must include the first 20 entries: " + trunc );
        }
        return true;
    }

    private static SequenceTaxonomyResolver.Result result( final SortedSet<String> not_found,
                                                           final int seq,
                                                           final int tax,
                                                           final boolean cancelled,
                                                           final String error ) {
        return new SequenceTaxonomyResolver.Result( not_found, seq, tax, cancelled, error );
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

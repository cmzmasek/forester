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

package org.forester.ws.seqdb;

import java.util.Arrays;
import java.util.List;

/**
 * Unit tests for {@link TaxonLineage} (the unified Spine-A taxon representation) -- entirely local (no network):
 * the rank&rarr;name and rank&rarr;tax-id lookups (case-insensitive, {@code domain}&harr;{@code superkingdom}
 * alias both ways, absent-rank), the own-taxon-wins-a-same-rank-collision consistency (at/taxIdAt never disagree),
 * {@code ownRankOrNull} (Linnaean-or-null), the lenient {@code isEmpty}, and {@code lineageNames}/{@code namesByRank}.
 */
public final class TaxonLineageTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TaxonLineage: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // a Felis-like lineage: ancestors with tax-ids (incl. a leading "no rank"), own taxon genus/Felis/9682
        final List<TaxonLineage.Ancestor> anc = Arrays.asList(
                new TaxonLineage.Ancestor( "cellular organisms", "no rank", "131567" ),
                new TaxonLineage.Ancestor( "Eukaryota", "superkingdom", "2759" ),
                new TaxonLineage.Ancestor( "Mammalia", "class", "40674" ),
                new TaxonLineage.Ancestor( "Carnivora", "order", "33554" ) );
        final TaxonLineage t = new TaxonLineage( "9682", "genus", "Felis", "cats", anc );

        // at(rank): name, case-insensitive; the own rank/name is folded in
        if ( !"Carnivora".equals( t.at( "order" ) ) || !"Carnivora".equals( t.at( "ORDER" ) )
                || !"Felis".equals( t.at( "genus" ) ) || !"Mammalia".equals( t.at( "class" ) ) ) {
            return fail( "at(rank) name lookup wrong" );
        }
        // taxIdAt(rank): id, folded own id; and null for a rank present without an id would apply (none here)
        if ( !"33554".equals( t.taxIdAt( "order" ) ) || !"40674".equals( t.taxIdAt( "class" ) )
                || !"9682".equals( t.taxIdAt( "genus" ) ) ) {
            return fail( "taxIdAt(rank) wrong" );
        }
        // domain<->superkingdom alias, for BOTH name and id, BOTH directions
        if ( !"Eukaryota".equals( t.at( "domain" ) ) || !"2759".equals( t.taxIdAt( "domain" ) ) ) {
            return fail( "at/taxIdAt(\"domain\") must alias to the superkingdom entry" );
        }
        final TaxonLineage domain_only = new TaxonLineage( null, null, null, null,
                Arrays.asList( new TaxonLineage.Ancestor( "Bacteria", "domain", "2" ) ) );
        if ( !"Bacteria".equals( domain_only.at( "superkingdom" ) )
                || !"2".equals( domain_only.taxIdAt( "superkingdom" ) ) ) {
            return fail( "at/taxIdAt(\"superkingdom\") must alias to a domain entry when present" );
        }
        // an absent rank -> null for both accessors; a null rank query -> null
        if ( ( t.at( "species" ) != null ) || ( t.taxIdAt( "species" ) != null ) || ( t.at( null ) != null ) ) {
            return fail( "an absent/null rank must resolve to null" );
        }
        // the "no rank" ancestor is skipped from the rank map, but present in the flat lineage
        if ( t.namesByRank().containsValue( "cellular organisms" ) ) {
            return fail( "a 'no rank' entry must be skipped from namesByRank" );
        }
        final List<String> ln = t.lineageNames();
        if ( ( ln.size() != 5 ) || !"cellular organisms".equals( ln.get( 0 ) )
                || !"Felis".equals( ln.get( ln.size() - 1 ) ) || !ln.contains( "Carnivora" ) ) {
            return fail( "lineageNames must be the full name path (incl. no-rank) ending in the taxon: " + ln );
        }

        // ownRankOrNull: a Linnaean rank comes back raw; "no rank"/"clade"/absent -> null (the Fetch writer relies
        // on this to avoid tax.setRank("no rank"))
        if ( !"genus".equals( t.ownRankOrNull() )
                || ( new TaxonLineage( null, "no rank", "x", null, null ).ownRankOrNull() != null )
                || ( new TaxonLineage( null, "clade", "x", null, null ).ownRankOrNull() != null )
                || ( new TaxonLineage( null, null, "x", null, null ).ownRankOrNull() != null ) ) {
            return fail( "ownRankOrNull must return a Linnaean rank, else null" );
        }

        // lenient isEmpty (the deliberate Spine-A decision): a scientific-name-only, rank-less taxon is NOT empty,
        // so fetch() caches it (killing the old re-prompt/re-fetch loop); a nameless, ancestorless taxon IS empty.
        if ( new TaxonLineage( null, null, "Foo", null, null ).isEmpty() ) {
            return fail( "a scientific-name-only taxon must NOT be empty (lenient isEmpty)" );
        }
        if ( !TaxonLineage.EMPTY.isEmpty() || !new TaxonLineage( null, null, null, null, null ).isEmpty() ) {
            return fail( "a nameless, ancestorless taxon must be empty" );
        }
        if ( domain_only.isEmpty() ) {
            return fail( "an ancestor-carrying taxon must not be empty" );
        }

        // DESYNC REGRESSION: the own taxon, folded in LAST over a same-rank ancestor, must replace the WHOLE
        // {name, tax-id} entry -- so an own name with a NULL own id does NOT leave the ancestor's id behind
        // (at and taxIdAt would otherwise point at different taxa, poisoning Spine B's tax-id keying).
        final TaxonLineage collide = new TaxonLineage( null, "order", "OwnOrder", null,
                Arrays.asList( new TaxonLineage.Ancestor( "AncOrder", "order", "111" ) ) );
        if ( !"OwnOrder".equals( collide.at( "order" ) ) || ( collide.taxIdAt( "order" ) != null ) ) {
            return fail( "own taxon must win a same-rank collision as one unit; at=" + collide.at( "order" )
                    + " taxIdAt=" + collide.taxIdAt( "order" ) );
        }
        return true;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  TaxonLineage test failed: " + msg );
        return false;
    }

    private TaxonLineageTest() {
    }
}

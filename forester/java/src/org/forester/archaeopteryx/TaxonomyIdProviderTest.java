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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.sdi.GSDIR;

/**
 * A gene tree whose taxonomy ids came from UniProt must reconcile against a species tree whose ids came from NCBI.
 * <p>
 * They are the SAME identifiers -- the UniProt taxonomy is the NCBI taxonomy and uses its taxids -- but they arrive
 * labelled with different providers depending on where they were read from, and the reconciliation keyed on
 * value+provider. Every tip was therefore unmappable, every tip was stripped, and the user got
 * "species could not be mapped between gene tree and species tree (based on taxonomy id)" instantly, without any
 * database ever being consulted. Reported against a real 59-tip 16S tree in which 58 of the 59 ids matched exactly
 * once the provider label was ignored.
 */
public final class TaxonomyIdProviderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TaxonomyIdProvider: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return normalization() && reconcilesAcrossProviders();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [TaxonomyIdProviderTest] " + msg );
        return false;
    }

    private static boolean normalization() {
        // every label for the NCBI taxonomy namespace collapses to one, INCLUDING no label at all
        for( final String p : new String[] { "ncbi", "NCBI", "uniprot", "UniProt", "ncbitaxon", "ncbi_taxonomy",
                                             "uniprot.taxonomy", " ncbi ", "", null } ) {
            if ( !"ncbi".equals( Identifier.normalizedProvider( p ) ) ) {
                return fail( "provider '" + p + "' names the NCBI taxonomy namespace, got "
                        + Identifier.normalizedProvider( p ) );
            }
        }
        // ...but a genuinely different namespace stays distinct, or a GTDB id could be mistaken for a taxid
        if ( "ncbi".equals( Identifier.normalizedProvider( "gtdb" ) ) ) {
            return fail( "a GTDB id must not be treated as an NCBI taxid" );
        }
        if ( !new Identifier( "1313", "uniprot" ).getValuePlusNormalizedProvider()
                .equals( new Identifier( "1313", "ncbi" ).getValuePlusNormalizedProvider() ) ) {
            return fail( "the same taxid under two labels must produce the same comparison key" );
        }
        // a different id is still a different key
        if ( new Identifier( "1313", "uniprot" ).getValuePlusNormalizedProvider()
                .equals( new Identifier( "1314", "ncbi" ).getValuePlusNormalizedProvider() ) ) {
            return fail( "different taxids must not collide" );
        }
        // the raw key is left alone -- it is what phyloXML round-trips and what users see
        if ( !"1313uniprot".equals( new Identifier( "1313", "uniprot" ).getValuePlusProvider() ) ) {
            return fail( "getValuePlusProvider must keep reporting the provider as stored" );
        }
        return true;
    }

    /** The reported failure, in miniature: gene tree ids labelled uniprot, species tree ids labelled ncbi. */
    private static boolean reconcilesAcrossProviders() {
        try {
            final GSDIR gsdir = new GSDIR( geneTree( "uniprot" ), speciesTree( "ncbi" ), true, true, true );
            if ( gsdir.getMinDuplicationsSum() < 0 ) {
                return fail( "reconciliation produced no result" );
            }
        }
        catch ( final Exception e ) {
            return fail( "a UniProt-labelled gene tree must reconcile against an NCBI-labelled species tree, got: "
                    + e.getMessage() );
        }
        // and the SAME trees under one label must of course still work (no regression for the ordinary case)
        try {
            new GSDIR( geneTree( "ncbi" ), speciesTree( "ncbi" ), true, true, true );
        }
        catch ( final Exception e ) {
            return fail( "the same-provider case must keep working: " + e.getMessage() );
        }
        return true;
    }

    /** ((A,B),(C,D)) over four taxids -- binary, as GSDIR requires. */
    private static Phylogeny geneTree( final String provider ) {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( pair( tip( "1313", provider ), tip( "257758", provider ) ) );
        root.addAsChild( pair( tip( "68892", provider ), tip( "113107", provider ) ) );
        return rooted( root );
    }

    private static Phylogeny speciesTree( final String provider ) {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( pair( tip( "1313", provider ), tip( "257758", provider ) ) );
        root.addAsChild( pair( tip( "68892", provider ), tip( "113107", provider ) ) );
        return rooted( root );
    }

    private static PhylogenyNode pair( final PhylogenyNode a, final PhylogenyNode b ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.addAsChild( a );
        n.addAsChild( b );
        return n;
    }

    private static PhylogenyNode tip( final String id, final String provider ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "t" + id );
        final Taxonomy t = new Taxonomy();
        t.setIdentifier( new Identifier( id, provider ) );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static Phylogeny rooted( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private TaxonomyIdProviderTest() {
    }
}

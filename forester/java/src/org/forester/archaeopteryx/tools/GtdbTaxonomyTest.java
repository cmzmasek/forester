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

import java.util.Map;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;

/** Headless unit test for {@link GtdbTaxonomy} (GTDB lineage parsing + tip application). */
public final class GtdbTaxonomyTest {

    public static void main( final String[] args ) {
        System.out.println( "GtdbTaxonomy: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            final String ecoli = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;"
                    + "f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli";

            // --- parse: all seven ranks, ordered domain..species, prefixes mapped ---
            final Map<String, String> lin = GtdbTaxonomy.parse( ecoli );
            if ( lin.size() != 7 ) {
                return fail( "expected 7 ranks, got " + lin.size() + " -> " + lin );
            }
            if ( !"Bacteria".equals( lin.get( "domain" ) ) || !"Pseudomonadota".equals( lin.get( "phylum" ) )
                    || !"Gammaproteobacteria".equals( lin.get( "class" ) )
                    || !"Enterobacterales".equals( lin.get( "order" ) )
                    || !"Enterobacteriaceae".equals( lin.get( "family" ) )
                    || !"Escherichia".equals( lin.get( "genus" ) )
                    || !"Escherichia coli".equals( lin.get( "species" ) ) ) {
                return fail( "rank->name mapping wrong: " + lin );
            }
            // order must be domain-first .. species-last (the most-specific-last invariant applyToNode relies on)
            final String[] order = lin.keySet().toArray( new String[ 0 ] );
            if ( !"domain".equals( order[ 0 ] ) || !"species".equals( order[ 6 ] ) ) {
                return fail( "lineage order must be domain..species, got " + java.util.Arrays.toString( order ) );
            }

            // --- parse: empty ranks skipped (a bare s__ / f__ with no name), trailing ';', whitespace tolerated ---
            final Map<String, String> partial = GtdbTaxonomy.parse( " d__Bacteria ; p__Firmicutes ; f__ ; g__ ; s__ ; " );
            if ( ( partial.size() != 2 ) || !"Bacteria".equals( partial.get( "domain" ) )
                    || !"Firmicutes".equals( partial.get( "phylum" ) ) || partial.containsKey( "family" ) ) {
                return fail( "empty ranks must be skipped, whitespace tolerated: " + partial );
            }
            if ( !GtdbTaxonomy.parse( null ).isEmpty() || !GtdbTaxonomy.parse( "no ranks here" ).isEmpty() ) {
                return fail( "null / rank-less input must parse to an empty lineage" );
            }

            // --- rankForPrefix + looksLikeGtdb ---
            if ( !"domain".equals( GtdbTaxonomy.rankForPrefix( 'd' ) )
                    || !"family".equals( GtdbTaxonomy.rankForPrefix( 'f' ) )
                    || ( GtdbTaxonomy.rankForPrefix( 'x' ) != null ) ) {
                return fail( "rankForPrefix mapping wrong" );
            }
            if ( !GtdbTaxonomy.looksLikeGtdb( ecoli ) || !GtdbTaxonomy.looksLikeGtdb( "s__Bacillus subtilis" )
                    || GtdbTaxonomy.looksLikeGtdb( "Escherichia coli" ) || GtdbTaxonomy.looksLikeGtdb( null ) ) {
                return fail( "looksLikeGtdb detection wrong" );
            }
            // must NOT be fooled by a stray letter+"__" MID-token (a rank tag must START a ';'-delimited field);
            // "photo__album" contains "o__", "sample_g__bin" contains "g__", but neither is a classification
            if ( GtdbTaxonomy.looksLikeGtdb( "photo__album" ) || GtdbTaxonomy.looksLikeGtdb( "sample_g__bin" )
                    || GtdbTaxonomy.looksLikeGtdb( "genus_data" ) ) {
                return fail( "looksLikeGtdb false-positive on a mid-token __" );
            }

            // --- applyToNode: writes gtdb:<rank> properties + a species taxonomy ---
            final PhylogenyNode n = new PhylogenyNode();
            final int written = GtdbTaxonomy.applyToNode( n, ecoli );
            if ( written != 7 ) {
                return fail( "applyToNode should report 7 ranks written, got " + written );
            }
            if ( !"Escherichia coli".equals( prop( n, "gtdb:species" ) )
                    || !"Enterobacteriaceae".equals( prop( n, "gtdb:family" ) )
                    || !"Bacteria".equals( prop( n, "gtdb:domain" ) ) ) {
                return fail( "applyToNode did not write the expected gtdb:* properties" );
            }
            if ( !n.getNodeData().isHasTaxonomy()
                    || !"Escherichia coli".equals( n.getNodeData().getTaxonomy().getScientificName() )
                    || !"species".equals( n.getNodeData().getTaxonomy().getRank() ) ) {
                return fail( "applyToNode must set a species-rank taxonomy from the most specific rank" );
            }

            // --- applyToNode: most-specific taxonomy when species is absent (genus is the deepest) ---
            final PhylogenyNode g = new PhylogenyNode();
            GtdbTaxonomy.applyToNode( g, "d__Archaea;p__Thermoproteota;g__Nitrosopumilus" );
            if ( !"Nitrosopumilus".equals( g.getNodeData().getTaxonomy().getScientificName() )
                    || !"genus".equals( g.getNodeData().getTaxonomy().getRank() )
                    || ( prop( g, "gtdb:species" ) != null ) ) {
                return fail( "with no species, the taxonomy must be the deepest present rank (genus)" );
            }

            // --- applyToNode: idempotent re-import (no duplicate gtdb:* properties) ---
            final int before = g.getNodeData().getProperties().size();
            GtdbTaxonomy.applyToNode( g, "d__Archaea;p__Thermoproteota;g__Nitrosopumilus" );
            if ( g.getNodeData().getProperties().size() != before ) {
                return fail( "re-import must replace, not duplicate, gtdb:* properties (was " + before + ", now "
                        + g.getNodeData().getProperties().size() + ")" );
            }

            // --- applyToNode: re-import a SHORTER/different lineage fully REPLACES (no stale deeper ranks) ---
            final PhylogenyNode re = new PhylogenyNode();
            GtdbTaxonomy.applyToNode( re, ecoli ); // 7 ranks incl. gtdb:class..gtdb:species
            GtdbTaxonomy.applyToNode( re, "d__Bacteria;p__Firmicutes" ); // shorter: only domain + phylum
            if ( !"Firmicutes".equals( prop( re, "gtdb:phylum" ) ) || !"Bacteria".equals( prop( re, "gtdb:domain" ) )
                    || ( prop( re, "gtdb:class" ) != null ) || ( prop( re, "gtdb:genus" ) != null )
                    || ( prop( re, "gtdb:species" ) != null ) ) {
                return fail( "a re-import of a shorter lineage must CLEAR the stale deeper ranks (class/genus/species)" );
            }
            if ( !"Firmicutes".equals( re.getNodeData().getTaxonomy().getScientificName() )
                    || !"phylum".equals( re.getNodeData().getTaxonomy().getRank() ) ) {
                return fail( "a re-import's taxonomy must become the new most-specific rank (phylum Firmicutes)" );
            }

            // --- applyToNode: an unparseable classification writes nothing ---
            final PhylogenyNode empty = new PhylogenyNode();
            if ( GtdbTaxonomy.applyToNode( empty, "not a gtdb string" ) != 0 ) {
                return fail( "applyToNode on a non-GTDB string must write nothing" );
            }

            // --- applyByTipName: matches external tips by name; unmatched / unparseable tips are left alone ---
            final org.forester.phylogeny.Phylogeny phy = new org.forester.phylogeny.Phylogeny();
            final PhylogenyNode root = new PhylogenyNode();
            final PhylogenyNode t1 = leaf( "GCF_001" ), t2 = leaf( "GCF_002" ), t3 = leaf( "GCF_003" );
            root.addAsChild( t1 );
            root.addAsChild( t2 );
            root.addAsChild( t3 );
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            final Map<String, String> table = new java.util.HashMap<String, String>();
            table.put( "GCF_001", ecoli );
            table.put( "GCF_002", "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis" );
            // GCF_003 is NOT in the table -> left untouched; "GCF_999" in the table matches no tip -> ignored
            table.put( "GCF_999", ecoli );
            final int annotated = GtdbTaxonomy.applyByTipName( phy, table );
            if ( annotated != 2 ) {
                return fail( "applyByTipName should annotate the 2 matched tips, got " + annotated );
            }
            if ( !"Pseudomonadota".equals( prop( t1, "gtdb:phylum" ) )
                    || !"Bacillota".equals( prop( t2, "gtdb:phylum" ) )
                    || t3.getNodeData().isHasTaxonomy() ) {
                return fail( "applyByTipName mis-assigned or touched an unmatched tip" );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected exception: " + e.getMessage() );
        }
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static String prop( final PhylogenyNode n, final String ref ) {
        if ( ( n.getNodeData() == null ) || ( n.getNodeData().getProperties() == null ) ) {
            return null;
        }
        for( final Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( ref.equals( p.getRef() ) ) {
                return p.getValue();
            }
        }
        return null;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [GtdbTaxonomyTest] " + msg );
        return false;
    }

    private GtdbTaxonomyTest() {
    }
}

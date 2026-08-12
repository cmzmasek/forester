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

package org.forester.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.forester.analysis.AncestralTaxonomyInference.InferenceResult;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.ws.seqdb.TaxonLineage;

/**
 * Headless unit tests for the pure ancestral-taxonomy inference engine. Tip lineages are supplied as fully
 * ranked {@link TaxonLineage}s (as the taxonomy service would return them), so no network is involved.
 */
public final class AncestralTaxonomyInferenceTest {

    public static boolean test() {
        return testMultifurcationAndTaxId() && testPreserveVsOverwrite() && testUnresolvedDescendant()
                && testNameOnly() && testCopyUpAndEdges() && testCommonTaxon();
    }

    /** {@code commonTaxonOf} = the deepest taxon shared by ALL a clade's tips, matched by tax-id/name (NOT by
     *  position), so it is robust to lineages of different depths (a full NCBI lineage vs a shallow stored one); an
     *  unresolved tip makes it null. Used to label a collapsed clade with its common taxon. */
    private static boolean testCommonTaxon() {
        final Map<PhylogenyNode, TaxonLineage> map = new HashMap<PhylogenyNode, TaxonLineage>();
        final PhylogenyNode felis = leaf( "Felis" );
        map.put( felis, lin( "class", "Mammalia", "40674", "order", "Carnivora", "33554", "family", "Felidae", "9681",
                             "genus", "Felis", "9682" ) );
        final PhylogenyNode canis = leaf( "Canis" );
        // a SHALLOW lineage (own = order Carnivora, no ancestors) -- different depth than Felis; a position-aligned
        // intersection would spuriously find nothing, but the set-based deepest-shared must still return Carnivora
        map.put( canis, lin( "order", "Carnivora", "33554" ) );
        final PhylogenyNode clade = internal( felis, canis );
        treeOf( clade );
        final TaxonLineage.Ancestor common = AncestralTaxonomyInference.commonTaxonOf( clade, map );
        if ( ( common == null ) || !"Carnivora".equals( common.getName() ) ) {
            return fail( "commonTaxonOf must find the deepest shared taxon across MISALIGNED lineages (Carnivora), got "
                    + ( common == null ? "null" : common.getName() ) );
        }
        // an unresolved descendant tip (absent from the map) -> no shared taxon (conservative)
        final PhylogenyNode felis2 = leaf( "Felis2" );
        final PhylogenyNode mystery = leaf( "Mystery" );
        final Map<PhylogenyNode, TaxonLineage> map2 = new HashMap<PhylogenyNode, TaxonLineage>();
        map2.put( felis2, lin( "order", "Carnivora", "33554", "genus", "Felis", "9682" ) );
        final PhylogenyNode clade2 = internal( felis2, mystery );
        treeOf( clade2 );
        if ( AncestralTaxonomyInference.commonTaxonOf( clade2, map2 ) != null ) {
            return fail( "an unresolved tip must make the common taxon null (conservative)" );
        }
        return true;
    }

    /** A 3-way multifurcation must infer the deepest level ALL three tips share (not just the first two), and the
     *  inferred taxon carries the rank + NCBI tax-id. */
    private static boolean testMultifurcationAndTaxId() {
        final Map<PhylogenyNode, TaxonLineage> map = new HashMap<PhylogenyNode, TaxonLineage>();
        final PhylogenyNode felis = leaf( "Felis" );
        map.put( felis, lin( "order", "Carnivora", "33554", "family", "Felidae", "9681", "genus", "Felis", "9682" ) );
        final PhylogenyNode panthera = leaf( "Panthera" );
        map.put( panthera,
                 lin( "order", "Carnivora", "33554", "family", "Felidae", "9681", "genus", "Panthera", "9688" ) );
        final PhylogenyNode canis = leaf( "Canis" );
        map.put( canis, lin( "order", "Carnivora", "33554", "family", "Canidae", "9608", "genus", "Canis", "9611" ) );
        final PhylogenyNode node = internal( felis, panthera, canis );
        final Phylogeny phy = treeOf( node );
        final InferenceResult r = AncestralTaxonomyInference.inferInternalTaxonomies( phy, map, false );
        if ( !node.getNodeData().isHasTaxonomy() ) {
            return fail( "the multifurcation node should have been assigned" );
        }
        final Taxonomy t = node.getNodeData().getTaxonomy();
        // binary-only (Felis + Panthera) would wrongly infer family Felidae; all three share only order Carnivora
        if ( !"Carnivora".equals( t.getScientificName() ) ) {
            return fail( "3-way clade must infer the deepest level ALL three share (Carnivora), got "
                    + t.getScientificName() );
        }
        if ( !"order".equals( t.getRank() ) ) {
            return fail( "inferred rank should be 'order', got " + t.getRank() );
        }
        if ( ( t.getIdentifier() == null ) || !"33554".equals( t.getIdentifier().getValue() )
                || !"ncbi".equalsIgnoreCase( t.getIdentifier().getProvider() ) ) {
            return fail( "inferred taxon should carry the NCBI tax-id 33554, got " + t.getIdentifier() );
        }
        if ( r.getAssigned() != 1 ) {
            return fail( "expected 1 internal node assigned, got " + r.getAssigned() );
        }
        return true;
    }

    /** An existing internal taxonomy is preserved when overwrite is off and replaced when it is on. */
    private static boolean testPreserveVsOverwrite() {
        // preserve
        Map<PhylogenyNode, TaxonLineage> map = new HashMap<PhylogenyNode, TaxonLineage>();
        PhylogenyNode[] pre = felisPanthera( map );
        PhylogenyNode canis = leaf( "Canis" );
        map.put( canis, lin( "order", "Carnivora", "33554", "family", "Canidae", "9608", "genus", "Canis", "9611" ) );
        PhylogenyNode preNode = internal( pre[ 0 ], pre[ 1 ] );
        annotate( preNode, "MyClade" );
        PhylogenyNode root = internal( preNode, canis );
        InferenceResult r = AncestralTaxonomyInference.inferInternalTaxonomies( treeOf( root ), map, false );
        if ( !"MyClade".equals( preNode.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "an existing internal taxonomy must be preserved when overwrite is off" );
        }
        if ( r.getSkippedExisting() != 1 ) {
            return fail( "expected 1 skipped-existing, got " + r.getSkippedExisting() );
        }
        if ( !"Carnivora".equals( root.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "the un-annotated root should still be inferred (Carnivora)" );
        }
        // overwrite -- fresh tree (the previous run mutated the root in place)
        map = new HashMap<PhylogenyNode, TaxonLineage>();
        pre = felisPanthera( map );
        preNode = internal( pre[ 0 ], pre[ 1 ] );
        annotate( preNode, "MyClade" );
        root = internal( preNode, leaf( "Canis" ) );
        map.put( root.getChildNode( 1 ),
                 lin( "order", "Carnivora", "33554", "family", "Canidae", "9608", "genus", "Canis", "9611" ) );
        AncestralTaxonomyInference.inferInternalTaxonomies( treeOf( root ), map, true );
        if ( !"Felidae".equals( preNode.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "overwrite must replace the existing taxonomy with the inferred one (Felidae), got "
                    + preNode.getNodeData().getTaxonomy().getScientificName() );
        }
        return true;
    }

    /** A single unresolved descendant tip blanks the inference of every node above it (conservative). */
    private static boolean testUnresolvedDescendant() {
        final Map<PhylogenyNode, TaxonLineage> map = new HashMap<PhylogenyNode, TaxonLineage>();
        final PhylogenyNode felis = leaf( "Felis" );
        map.put( felis, lin( "order", "Carnivora", "33554", "genus", "Felis", "9682" ) );
        final PhylogenyNode unknown = leaf( "x_unknown" );
        map.put( unknown, TaxonLineage.EMPTY );
        final PhylogenyNode root = internal( felis, unknown );
        final InferenceResult r = AncestralTaxonomyInference.inferInternalTaxonomies( treeOf( root ), map, false );
        if ( root.getNodeData().isHasTaxonomy() ) {
            return fail( "a node with an unresolved descendant must infer nothing" );
        }
        if ( ( r.getAssigned() != 0 ) || ( r.getUnresolvedTips() != 1 ) || ( r.getResolvedTips() != 1 ) ) {
            return fail( "expected assigned=0, tipsWithout=1, tipsWith=1; got " + r.getAssigned() + "/"
                    + r.getUnresolvedTips() + "/" + r.getResolvedTips() );
        }
        return true;
    }

    /** Lineages whose ancestors carry no tax-id still infer by name (and no bogus identifier is written). */
    private static boolean testNameOnly() {
        final Map<PhylogenyNode, TaxonLineage> map = new HashMap<PhylogenyNode, TaxonLineage>();
        final PhylogenyNode felis = leaf( "Felis" );
        map.put( felis, lin( "order", "Carnivora", null, "genus", "Felis", null ) );
        final PhylogenyNode panthera = leaf( "Panthera" );
        map.put( panthera, lin( "order", "Carnivora", null, "genus", "Panthera", null ) );
        final PhylogenyNode root = internal( felis, panthera );
        AncestralTaxonomyInference.inferInternalTaxonomies( treeOf( root ), map, false );
        final Taxonomy t = root.getNodeData().getTaxonomy();
        if ( ( t == null ) || !"Carnivora".equals( t.getScientificName() ) || !"order".equals( t.getRank() ) ) {
            return fail( "name-only lineages should still infer Carnivora (order)" );
        }
        if ( t.getIdentifier() != null ) {
            return fail( "no NCBI tax-id was available, so the inferred taxon must carry no identifier" );
        }
        return true;
    }

    /** Unary copy-up, plus the null / empty / all-unresolved edge cases. */
    private static boolean testCopyUpAndEdges() {
        // unary: root -> mid -> [Felis, Panthera]; mid infers Felidae, and its single-child parent root copies it up
        final Map<PhylogenyNode, TaxonLineage> map = new HashMap<PhylogenyNode, TaxonLineage>();
        final PhylogenyNode[] pre = felisPanthera( map );
        final PhylogenyNode mid = internal( pre[ 0 ], pre[ 1 ] );
        final PhylogenyNode root = internal( mid );
        AncestralTaxonomyInference.inferInternalTaxonomies( treeOf( root ), map, false );
        if ( !"Felidae".equals( mid.getNodeData().getTaxonomy().getScientificName() )
                || !"Felidae".equals( root.getNodeData().getTaxonomy().getScientificName() ) ) {
            return fail( "a single-child parent should copy up the child's inferred taxon (Felidae)" );
        }
        // null tree
        if ( AncestralTaxonomyInference.inferInternalTaxonomies( null, map, false ).getAssigned() != 0 ) {
            return fail( "null tree must assign nothing" );
        }
        // all tips unresolved
        final Map<PhylogenyNode, TaxonLineage> empty = new HashMap<PhylogenyNode, TaxonLineage>();
        final PhylogenyNode a = leaf( "a" );
        final PhylogenyNode b = leaf( "b" );
        final InferenceResult r = AncestralTaxonomyInference.inferInternalTaxonomies( treeOf( internal( a, b ) ),
                                                                                      empty, false );
        if ( ( r.getAssigned() != 0 ) || ( r.getUnresolvedTips() != 2 ) ) {
            return fail( "all-unresolved tree: expected assigned=0, tipsWithout=2, got " + r.getAssigned() + "/"
                    + r.getUnresolvedTips() );
        }
        return true;
    }

    // ----- fixtures -----

    /** Two cat tips (Felis, Panthera) sharing family Felidae/order Carnivora; registers them in {@code map}. */
    private static PhylogenyNode[] felisPanthera( final Map<PhylogenyNode, TaxonLineage> map ) {
        final PhylogenyNode felis = leaf( "Felis" );
        map.put( felis, lin( "order", "Carnivora", "33554", "family", "Felidae", "9681", "genus", "Felis", "9682" ) );
        final PhylogenyNode panthera = leaf( "Panthera" );
        map.put( panthera,
                 lin( "order", "Carnivora", "33554", "family", "Felidae", "9681", "genus", "Panthera", "9688" ) );
        return new PhylogenyNode[] { felis, panthera };
    }

    /** Build a {@link TaxonLineage} from (rank, name, taxId) triples root&rarr;deepest; the last triple is the taxon
     *  itself, the earlier ones its ancestors. */
    private static TaxonLineage lin( final String... triples ) {
        final int groups = triples.length / 3;
        final List<TaxonLineage.Ancestor> anc = new ArrayList<TaxonLineage.Ancestor>();
        String own_rank = null;
        String own_name = null;
        String own_id = null;
        for( int g = 0; g < groups; g++ ) {
            final String rank = triples[ ( g * 3 ) ];
            final String name = triples[ ( g * 3 ) + 1 ];
            final String id = triples[ ( g * 3 ) + 2 ];
            if ( g == ( groups - 1 ) ) {
                own_rank = rank;
                own_name = name;
                own_id = id;
            }
            else {
                anc.add( new TaxonLineage.Ancestor( name, rank, id ) );
            }
        }
        return new TaxonLineage( own_id, own_rank, own_name, null, anc );
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode internal( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode c : children ) {
            n.addAsChild( c );
        }
        return n;
    }

    private static void annotate( final PhylogenyNode node, final String sci ) {
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        node.getNodeData().setTaxonomy( t );
    }

    private static Phylogeny treeOf( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "AncestralTaxonomyInference test failed: " + msg );
        return false;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }
}

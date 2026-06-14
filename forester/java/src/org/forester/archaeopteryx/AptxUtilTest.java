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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * Unit tests for {@link AptxUtil}. Lives in the {@code org.forester.archaeopteryx} package. Run
 * standalone via {@link #main(String[])}, or as part of the suite via {@link #test()}.
 */
public final class AptxUtilTest {

    public static void main( final String[] args ) {
        System.out.println( "AptxUtil: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        return testHasAtLeastOneNodeWithDomainArchitecture() && testGraphicsExportTypes() && testColorizableRanks()
                && testRankCounts() && testRankCoverageCounts() && testNodePruningOutcome()
                && testBranchesToCollapse();
    }

    /**
     * The candidate selection behind the "Collapse Branches" tools: only internal, non-root branches
     * below the threshold qualify (by max confidence, or by branch length); external nodes and the
     * root never do. Tree: root -&gt; (i1[conf 50, bl 0.1], i2[conf 90, bl 2.0]), each with 2 leaves.
     */
    private static boolean testBranchesToCollapse() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode i1 = internalBranch( 50.0, 0.1 );
        final PhylogenyNode i2 = internalBranch( 90.0, 2.0 );
        root.addAsChild( i1 );
        root.addAsChild( i2 );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        // confidence threshold 70 -> only i1 (50 < 70); i2 (90) is kept
        final List<PhylogenyNode> conf70 = AptxUtil.branchesToCollapseByConfidence( phy, 70.0 );
        if ( ( conf70.size() != 1 ) || ( conf70.get( 0 ) != i1 ) ) {
            return fail( "confidence < 70 must select only i1; got " + conf70.size() );
        }
        if ( AptxUtil.branchesToCollapseByConfidence( phy, 95.0 ).size() != 2 ) {
            return fail( "confidence < 95 must select both internal branches" );
        }
        if ( !AptxUtil.branchesToCollapseByConfidence( phy, 0.0 ).isEmpty() ) {
            return fail( "confidence < 0 must select nothing" );
        }
        // branch length threshold 1.0 -> only i1 (0.1 < 1.0); i2 (2.0) kept
        final List<PhylogenyNode> bl1 = AptxUtil.branchesToCollapseByBranchLength( phy, 1.0 );
        if ( ( bl1.size() != 1 ) || ( bl1.get( 0 ) != i1 ) ) {
            return fail( "branch length < 1.0 must select only i1; got " + bl1.size() );
        }
        if ( AptxUtil.branchesToCollapseByBranchLength( phy, 3.0 ).size() != 2 ) {
            return fail( "branch length < 3.0 must select both internal branches" );
        }
        // null / empty trees yield no candidates (no crash)
        if ( !AptxUtil.branchesToCollapseByConfidence( null, 50.0 ).isEmpty()
                || !AptxUtil.branchesToCollapseByBranchLength( new Phylogeny(), 1.0 ).isEmpty() ) {
            return fail( "null/empty tree must yield no candidates" );
        }
        return true;
    }

    /** An internal node (two leaf children) carrying a bootstrap confidence and a branch length. */
    private static PhylogenyNode internalBranch( final double confidence, final double branch_length ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setDistanceToParent( branch_length );
        n.getBranchData().addConfidence( new Confidence( confidence, "bootstrap" ) );
        n.addAsChild( new PhylogenyNode() );
        n.addAsChild( new PhylogenyNode() );
        return n;
    }

    /**
     * getRankCounts must count ranks on external nodes too: species (and other low ranks) live on
     * leaves, and the colorizer anchors on them, so excluding leaves undercounted them. Tree:
     * ((l1=species, l2=species)order, l3=species) -- species is on 3 leaves, order on 1 internal node.
     */
    private static boolean testRankCounts() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode order_node = new PhylogenyNode();
        order_node.getNodeData().setTaxonomy( taxonomy( "Diptera", "order" ) );
        order_node.addAsChild( speciesLeaf( "Aedes aegypti" ) );
        order_node.addAsChild( speciesLeaf( "Culex pipiens" ) );
        root.addAsChild( order_node );
        root.addAsChild( speciesLeaf( "Apis mellifera" ) );
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        final Map<String, Integer> counts = AptxUtil.getRankCounts( phy );
        if ( ( counts.get( "species" ) == null ) || ( counts.get( "species" ) != 3 ) ) {
            return fail( "species on 3 leaves must count as 3 (was the leaf-excluding bug); got "
                    + counts.get( "species" ) );
        }
        if ( ( counts.get( "order" ) == null ) || ( counts.get( "order" ) != 1 ) ) {
            return fail( "order on 1 internal node must count as 1; got " + counts.get( "order" ) );
        }
        return true;
    }

    private static PhylogenyNode speciesLeaf( final String sci_name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( sci_name );
        n.getNodeData().setTaxonomy( taxonomy( sci_name, "species" ) );
        return n;
    }

    /**
     * The decision behind "Delete/Retain Selected Nodes": no tree (or a single leaf) and an empty
     * selection are caught with warnings; deleting every selected node is refused; retaining can
     * never empty the tree; normal prunes proceed.
     */
    private static boolean testNodePruningOutcome() {
        // fewer than two external nodes -> nothing meaningful to prune
        if ( AptxUtil.nodePruningOutcome( 0, 0, true ) != AptxUtil.NodePruningOutcome.NO_TREE ) {
            return fail( "no tree (0 external nodes) must be NO_TREE" );
        }
        if ( AptxUtil.nodePruningOutcome( 1, 0, false ) != AptxUtil.NodePruningOutcome.NO_TREE ) {
            return fail( "a single-leaf tree must be NO_TREE" );
        }
        // a real tree but nothing selected (NO_TREE is checked first, so use a valid tree size)
        if ( AptxUtil.nodePruningOutcome( 10, 0, true ) != AptxUtil.NodePruningOutcome.NO_SELECTION ) {
            return fail( "no selected external nodes must be NO_SELECTION (delete)" );
        }
        if ( AptxUtil.nodePruningOutcome( 10, 0, false ) != AptxUtil.NodePruningOutcome.NO_SELECTION ) {
            return fail( "no selected external nodes must be NO_SELECTION (retain)" );
        }
        // deleting all selected when all are selected would empty the tree
        if ( AptxUtil.nodePruningOutcome( 10, 10, true ) != AptxUtil.NodePruningOutcome.WOULD_REMOVE_ALL ) {
            return fail( "deleting all 10 of 10 must be WOULD_REMOVE_ALL" );
        }
        // retaining never empties the tree, even with everything selected
        if ( AptxUtil.nodePruningOutcome( 10, 10, false ) != AptxUtil.NodePruningOutcome.OK ) {
            return fail( "retaining all selected nodes keeps them -> OK" );
        }
        // normal prunes proceed; deleting down to a single node is allowed
        if ( AptxUtil.nodePruningOutcome( 10, 3, true ) != AptxUtil.NodePruningOutcome.OK ) {
            return fail( "delete 3 of 10 -> OK" );
        }
        if ( AptxUtil.nodePruningOutcome( 10, 3, false ) != AptxUtil.NodePruningOutcome.OK ) {
            return fail( "retain 3 of 10 -> OK" );
        }
        if ( AptxUtil.nodePruningOutcome( 10, 9, true ) != AptxUtil.NodePruningOutcome.OK ) {
            return fail( "delete 9 of 10 leaves 1 node -> OK" );
        }
        return true;
    }

    /**
     * getColorizableRanks must offer only ranks that can actually drive a colorization -- present on
     * internal nodes with a count &gt; 0 -- labeled "rank (count) (P% coverage)" in canonical order,
     * and an empty array (which the caller turns into a warning) when nothing is colorizable.
     */
    private static boolean testColorizableRanks() {
        final Map<String, Integer> no_coverage = new LinkedHashMap<String, Integer>();
        // nothing present -> nothing colorizable (caller shows a warning and returns)
        if ( AptxUtil.getColorizableRanks( null, no_coverage, 10 ).length != 0 ) {
            return fail( "null rank counts must yield no colorizable ranks" );
        }
        if ( AptxUtil.getColorizableRanks( new LinkedHashMap<String, Integer>(), no_coverage, 10 ).length != 0 ) {
            return fail( "empty rank counts must yield no colorizable ranks" );
        }
        // a count of 0 is present but not colorizable -> excluded
        final Map<String, Integer> zero_only = new LinkedHashMap<String, Integer>();
        zero_only.put( "family", 0 );
        if ( AptxUtil.getColorizableRanks( zero_only, no_coverage, 10 ).length != 0 ) {
            return fail( "a rank with count 0 must not be offered" );
        }
        // only count > 0 ranks appear, labeled "rank (count) (P% coverage)", canonical order (order < genus)
        final Map<String, Integer> counts = new LinkedHashMap<String, Integer>();
        counts.put( "genus", 4 ); // inserted out of canonical order on purpose
        counts.put( "order", 2 );
        counts.put( "family", 0 ); // present but zero -> excluded
        final Map<String, Integer> coverage = new LinkedHashMap<String, Integer>();
        coverage.put( "order", 6 ); // 6 of 10 external nodes -> 60%
        coverage.put( "genus", 3 ); // 3 of 10 -> 30%
        final String[] ranks = AptxUtil.getColorizableRanks( counts, coverage, 10 );
        if ( ranks.length != 2 ) {
            return fail( "only ranks with count > 0 should appear; got " + Arrays.toString( ranks ) );
        }
        if ( !ranks[ 0 ].equals( "order (2) (60% coverage)" ) || !ranks[ 1 ].equals( "genus (4) (30% coverage)" ) ) {
            return fail( "labels must be \"rank (count) (P% coverage)\" in canonical order; got "
                    + Arrays.toString( ranks ) );
        }
        // missing coverage entry, or zero total, yields 0% rather than a crash
        final Map<String, Integer> just_order = new LinkedHashMap<String, Integer>();
        just_order.put( "order", 2 );
        if ( !AptxUtil.getColorizableRanks( just_order, null, 10 )[ 0 ].equals( "order (2) (0% coverage)" ) ) {
            return fail( "absent coverage map must render 0% coverage" );
        }
        if ( !AptxUtil.getColorizableRanks( just_order, coverage, 0 )[ 0 ].equals( "order (2) (0% coverage)" ) ) {
            return fail( "zero external total must render 0% coverage (no divide-by-zero)" );
        }
        // a nonzero-but-tiny coverage reads "<1%", never "0%" (a rank that colors something is not 0%)
        final Map<String, Integer> tiny = new LinkedHashMap<String, Integer>();
        tiny.put( "order", 1 );
        if ( !AptxUtil.getColorizableRanks( just_order, tiny, 300 )[ 0 ].equals( "order (2) (<1% coverage)" ) ) {
            return fail( "1 of 300 covered (0.33%) must render <1% coverage, not 0%" );
        }
        if ( !AptxUtil.coveragePercentLabel( 0, 100 ).equals( "0%" )
                || !AptxUtil.coveragePercentLabel( 1, 300 ).equals( "<1%" ) // 0.33% rounds to 0 -> <1%
                || !AptxUtil.coveragePercentLabel( 1, 200 ).equals( "1%" )  // 0.5% rounds half-up to 1%
                || !AptxUtil.coveragePercentLabel( 6, 10 ).equals( "60%" )
                || !AptxUtil.coveragePercentLabel( null, 10 ).equals( "0%" ) ) {
            return fail( "coveragePercentLabel: wrong percent label" );
        }
        return true;
    }

    /**
     * getRankCoverageCounts: an external node is covered by rank R when a non-root ancestor (or
     * itself) carries rank R. Tree: ((l1,l2)order=Diptera, l3) -- order covers l1,l2 (2 of 3 leaves).
     */
    private static boolean testRankCoverageCounts() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode order_node = new PhylogenyNode();
        order_node.getNodeData().setTaxonomy( taxonomy( "Diptera", "order" ) );
        order_node.addAsChild( leaf( "l1" ) );
        order_node.addAsChild( leaf( "l2" ) );
        root.addAsChild( order_node );
        root.addAsChild( leaf( "l3" ) ); // outside any ranked subtree
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        final Map<String, Integer> coverage = AptxUtil.getRankCoverageCounts( phy );
        if ( ( coverage.get( "order" ) == null ) || ( coverage.get( "order" ) != 2 ) ) {
            return fail( "rank 'order' must cover exactly the 2 leaves under it; got " + coverage.get( "order" ) );
        }
        if ( coverage.containsKey( "family" ) ) {
            return fail( "a rank not present must not appear in the coverage map" );
        }
        // empty / null trees are handled without error
        if ( !AptxUtil.getRankCoverageCounts( null ).isEmpty() || !AptxUtil.getRankCoverageCounts( new Phylogeny() ).isEmpty() ) {
            return fail( "null/empty tree must yield empty coverage" );
        }
        return true;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static org.forester.phylogeny.data.Taxonomy taxonomy( final String sci_name, final String rank ) {
        final org.forester.phylogeny.data.Taxonomy t = new org.forester.phylogeny.data.Taxonomy();
        t.setScientificName( sci_name );
        try {
            t.setRank( rank );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
        return t;
    }

    /**
     * Pins the set of supported graphics-export formats. The legacy raster formats BMP and GIF were
     * removed (both are strictly dominated by PNG for tree figures: BMP is uncompressed and ~10x
     * larger, GIF's 256-color palette dithers anti-aliased text and branch colors). This guard fails
     * if either is reintroduced, or if one of the intended formats is dropped by accident.
     */
    private static boolean testGraphicsExportTypes() {
        final Set<GraphicsExportType> expected = EnumSet.of( GraphicsExportType.JPG,
                                                             GraphicsExportType.PDF,
                                                             GraphicsExportType.PNG,
                                                             GraphicsExportType.TIFF,
                                                             GraphicsExportType.SVG,
                                                             GraphicsExportType.EPS );
        final Set<GraphicsExportType> actual = new HashSet<GraphicsExportType>();
        for( final GraphicsExportType t : GraphicsExportType.values() ) {
            actual.add( t );
            final String suffix = t.toString();
            if ( "bmp".equals( suffix ) || "gif".equals( suffix ) ) {
                return fail( "removed export format reintroduced: " + t.name() + " (" + suffix + ")" );
            }
        }
        if ( !actual.equals( expected ) ) {
            return fail( "supported export formats changed: expected " + expected + " but found " + actual );
        }
        return true;
    }

    private static boolean testHasAtLeastOneNodeWithDomainArchitecture() {
        // a tree whose leaves carry no sequence/domain data
        final Phylogeny no_domains = new Phylogeny();
        final PhylogenyNode root0 = new PhylogenyNode();
        root0.addAsChild( new PhylogenyNode() );
        root0.addAsChild( new PhylogenyNode() );
        no_domains.setRoot( root0 );
        no_domains.externalNodesHaveChanged();
        if ( AptxUtil.isHasAtLeastOneNodeWithDomainArchitecture( no_domains ) ) {
            return fail( "a tree without domain architectures must not be detected as having them" );
        }
        // give a single leaf a sequence with a domain architecture
        final Phylogeny with_domains = new Phylogeny();
        final PhylogenyNode root1 = new PhylogenyNode();
        final PhylogenyNode leaf = new PhylogenyNode();
        final ArrayList<PhylogenyData> domains = new ArrayList<PhylogenyData>();
        domains.add( new ProteinDomain( "d0", 10, 20 ) );
        domains.add( new ProteinDomain( "d1", 30, 40 ) );
        final Sequence seq = new Sequence();
        seq.setDomainArchitecture( new DomainArchitecture( domains, 100 ) );
        leaf.getNodeData().setSequence( seq );
        root1.addAsChild( leaf );
        root1.addAsChild( new PhylogenyNode() );
        with_domains.setRoot( root1 );
        with_domains.externalNodesHaveChanged();
        if ( !AptxUtil.isHasAtLeastOneNodeWithDomainArchitecture( with_domains ) ) {
            return fail( "a tree with a domain architecture must be detected" );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [AptxUtilTest] " + message );
        return false;
    }

    private AptxUtilTest() {
        // not instantiable
    }
}

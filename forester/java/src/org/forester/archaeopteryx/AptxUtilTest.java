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

import java.awt.Color;
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
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

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
        return testHasAtLeastOneNodeWithDomainArchitecture() && testGraphicsExportTypes() && testRankChoices()
                && testRankCounts() && testRankCoverageCounts() && testNodePruningOutcome()
                && testBranchesToCollapse() && testConfigFileOption() && testScanForDataPresence()
                && testAssignDistinctColors() && testShortenLabel()
                && testInternalNamesLookLikeConfidenceValues();
    }

    /**
     * The load-time "internal labels look like support values?" heuristic: numeric internal labels in
     * [0,100] (bare or bracketed) on >=2 nodes look like confidence; a real clade name, an out-of-range
     * number, fewer than two labels, or none do not. The root label is ignored.
     */
    private static boolean testInternalNamesLookLikeConfidenceValues() {
        try {
            // isSupportLikeNumber units (incl. the bracketed form)
            if ( !AptxUtil.isSupportLikeNumber( "95" ) || !AptxUtil.isSupportLikeNumber( "0.95" )
                    || !AptxUtil.isSupportLikeNumber( "[87]" ) || !AptxUtil.isSupportLikeNumber( "100" )
                    || !AptxUtil.isSupportLikeNumber( "0" ) ) {
                return fail( "support-like numbers (incl. bracketed) must be recognized" );
            }
            if ( AptxUtil.isSupportLikeNumber( "150" ) || AptxUtil.isSupportLikeNumber( "-1" )
                    || AptxUtil.isSupportLikeNumber( "Clade" ) || AptxUtil.isSupportLikeNumber( "" ) ) {
                return fail( "out-of-range / non-numeric must not be support-like" );
            }
            // two numeric internal labels in range -> looks like confidence (root label ignored)
            if ( !AptxUtil.internalNamesLookLikeConfidenceValues(
                    Phylogeny.createInstanceFromNhxString( "((A,B)95,(C,D)87)myTree;" ) ) ) {
                return fail( "two numeric internal labels should look like confidence" );
            }
            if ( !AptxUtil.internalNamesLookLikeConfidenceValues(
                    Phylogeny.createInstanceFromNhxString( "((A,B)0.95,(C,D)0.87);" ) ) ) {
                return fail( "posterior-probability internal labels should look like confidence" );
            }
            // a real clade name vetoes the offer
            if ( AptxUtil.internalNamesLookLikeConfidenceValues(
                    Phylogeny.createInstanceFromNhxString( "((A,B)95,(C,D)Mammalia);" ) ) ) {
                return fail( "a non-numeric clade name must veto the offer" );
            }
            // only one labeled internal node -> below the >=2 threshold
            if ( AptxUtil.internalNamesLookLikeConfidenceValues(
                    Phylogeny.createInstanceFromNhxString( "((A,B)95,(C,D));" ) ) ) {
                return fail( "a single numeric internal label is not enough" );
            }
            // out-of-range number vetoes
            if ( AptxUtil.internalNamesLookLikeConfidenceValues(
                    Phylogeny.createInstanceFromNhxString( "((A,B)150,(C,D)87);" ) ) ) {
                return fail( "an out-of-range number must veto the offer" );
            }
            // no internal labels
            if ( AptxUtil.internalNamesLookLikeConfidenceValues(
                    Phylogeny.createInstanceFromNhxString( "((A,B),(C,D));" ) ) ) {
                return fail( "no internal labels -> no offer" );
            }
            // end-to-end on the BRACKETED form (built by hand, since a parser may strip "[95]" as a
            // comment): detector accepts it, then strip-brackets + transfer yields confidence 95 / empty name.
            final Phylogeny bracketed = new Phylogeny();
            final PhylogenyNode root = new PhylogenyNode();
            final PhylogenyNode i1 = new PhylogenyNode();
            i1.setName( "[95]" );
            i1.addAsChild( new PhylogenyNode() );
            i1.addAsChild( new PhylogenyNode() );
            final PhylogenyNode i2 = new PhylogenyNode();
            i2.setName( "[87]" );
            i2.addAsChild( new PhylogenyNode() );
            i2.addAsChild( new PhylogenyNode() );
            root.addAsChild( i1 );
            root.addAsChild( i2 );
            bracketed.setRoot( root );
            bracketed.externalNodesHaveChanged();
            if ( !AptxUtil.internalNamesLookLikeConfidenceValues( bracketed ) ) {
                return fail( "bracketed numeric internal labels should look like confidence" );
            }
            AptxUtil.stripBracketsFromInternalNames( bracketed );
            PhylogenyMethods.transferInternalNodeNamesToConfidence( bracketed, "" );
            if ( !i1.getBranchData().isHasConfidences()
                    || ( i1.getBranchData().getConfidences().get( 0 ).getValue() != 95.0 )
                    || !i1.getName().isEmpty() ) {
                return fail( "strip-brackets + transfer must yield confidence 95 and clear the name" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "unexpected exception: " + e );
        }
    }

    /**
     * shortenLabel returns short/exact-length/null labels unchanged, and an over-long label as its head
     * (with trailing whitespace trimmed) plus a single-character ellipsis -- the display-only treatment
     * for pasted-in UniProt/NCBI FASTA-header names.
     */
    private static boolean testShortenLabel() {
        if ( !"abc".equals( AptxUtil.shortenLabel( "abc", 10 ) ) || ( AptxUtil.shortenLabel( null, 10 ) != null ) ) {
            return fail( "short or null labels must be returned unchanged" );
        }
        if ( !"abcde".equals( AptxUtil.shortenLabel( "abcde", 5 ) ) ) {
            return fail( "a label exactly at the limit must not be shortened" );
        }
        if ( !"abcd…".equals( AptxUtil.shortenLabel( "abcdef", 5 ) ) ) {
            return fail( "an over-long label must become head + ellipsis; got " + AptxUtil.shortenLabel( "abcdef", 5 ) );
        }
        if ( !"ab…".equals( AptxUtil.shortenLabel( "ab cdef", 4 ) ) ) {
            return fail( "whitespace before the cut must be trimmed; got " + AptxUtil.shortenLabel( "ab cdef", 4 ) );
        }
        final String header = "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical s-adenosyl methionine domain-containing protein OS=Fusarium phyllophilum";
        final String shortened = AptxUtil.shortenLabel( header, 60 );
        if ( ( shortened.length() > 60 ) || !shortened.endsWith( "…" )
                || !header.startsWith( shortened.substring( 0, shortened.length() - 1 ) ) ) {
            return fail( "a UniProt header must shorten to a <=60-char head ending in the ellipsis; got " + shortened );
        }
        return true;
    }

    /**
     * The rank colorizer's palette: each taxon gets a color, the colors are all distinct, the result
     * is deterministic (same input -&gt; same colors) and ordered like the (sorted) input, and an
     * empty/null input yields an empty map.
     */
    private static boolean testAssignDistinctColors() {
        if (!AptxUtil.assignDistinctColors(null).isEmpty()
                || !AptxUtil.assignDistinctColors(new java.util.TreeSet<String>()).isEmpty()) {
            return fail("null/empty taxa must yield an empty color map");
        }
        final java.util.SortedSet<String> taxa = new java.util.TreeSet<String>(
                java.util.Arrays.asList("Carnivora", "Rodentia", "Primates", "Chiroptera", "Cetacea"));
        final Map<String, Color> a = AptxUtil.assignDistinctColors(taxa);
        if (a.size() != taxa.size()) {
            return fail("every taxon must get a color");
        }
        // all colors distinct
        if (new java.util.HashSet<Color>(a.values()).size() != taxa.size()) {
            return fail("colors must be pairwise distinct; got " + a.values());
        }
        // ordered like the sorted input
        if (!new java.util.ArrayList<String>(a.keySet()).equals(new java.util.ArrayList<String>(taxa))) {
            return fail("color map must follow the input's (sorted) order");
        }
        // deterministic
        final Map<String, Color> b = AptxUtil.assignDistinctColors(taxa);
        if (!a.equals(b)) {
            return fail("assignDistinctColors must be deterministic for the same taxa");
        }
        // a single taxon still gets a (valid) color
        final java.util.SortedSet<String> one = new java.util.TreeSet<String>(java.util.Arrays.asList("Carnivora"));
        if ((AptxUtil.assignDistinctColors(one).size() != 1)
                || (AptxUtil.assignDistinctColors(one).get("Carnivora") == null)) {
            return fail("a single taxon must get exactly one color");
        }
        return true;
    }

    /**
     * The control panel uses {@link AptxUtil#scanForDataPresence} to show only the "Display Data"
     * checkboxes whose data is actually present. A barren tree (topology + node names only) must
     * report only node names; a richly annotated node must light up every corresponding flag; and
     * a null/empty tree yields the empty set.
     */
    private static boolean testScanForDataPresence() {
        try {
            if (!AptxUtil.scanForDataPresence(null).isEmpty()
                    || !AptxUtil.scanForDataPresence(new Phylogeny()).isEmpty()) {
                return false;
            }
            // barren tree: two named leaves, no branch lengths, no annotations
            final Phylogeny barren = new Phylogeny();
            final PhylogenyNode br = new PhylogenyNode();
            final PhylogenyNode b1 = new PhylogenyNode();
            b1.setName("a");
            final PhylogenyNode b2 = new PhylogenyNode();
            b2.setName("b");
            br.addAsChild(b1);
            br.addAsChild(b2);
            barren.setRoot(br);
            barren.externalNodesHaveChanged();
            final Set<Integer> barren_set = AptxUtil.scanForDataPresence(barren);
            if (!barren_set.contains(Configuration.show_node_names) || (barren_set.size() != 1)) {
                return false;
            }
            // rich node: one leaf carrying many kinds of data
            final Phylogeny rich = new Phylogeny();
            final PhylogenyNode rr = new PhylogenyNode();
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName("leaf");
            leaf.setDistanceToParent(0.5);
            leaf.getBranchData().addConfidence(new Confidence(0.9, "bootstrap"));
            leaf.getBranchData().setBranchWidth(new BranchWidth(3.0));
            leaf.getBranchData().setBranchColor(new BranchColor(java.awt.Color.RED));
            leaf.getNodeData().setEvent(Event.createSingleDuplicationEvent());
            final Taxonomy tax = new Taxonomy();
            tax.setScientificName("Homo sapiens");
            tax.setCommonName("human");
            tax.setTaxonomyCode("HUMAN");
            tax.setRank("species");
            leaf.getNodeData().setTaxonomy(tax);
            final Sequence seq = new Sequence();
            seq.setName("p53");
            seq.setGeneName("TP53");
            seq.setMolecularSequence("MEEPQSDPSV");
            seq.setAccession(new Accession("P04637", "uniprot"));
            leaf.getNodeData().setSequence(seq);
            final PropertiesList props = new PropertiesList();
            props.addProperty(new Property("aptx:test", "1", "", "xsd:string", Property.AppliesTo.NODE));
            leaf.getNodeData().setProperties(props);
            leaf.getNodeData().setVector(Arrays.asList(1.0, 2.0, 3.0));
            rr.addAsChild(leaf);
            rr.addAsChild(new PhylogenyNode());
            rich.setRoot(rr);
            rich.externalNodesHaveChanged();
            final Set<Integer> s = AptxUtil.scanForDataPresence(rich);
            final int[] expected = { Configuration.show_node_names, Configuration.write_branch_length_values,
                    Configuration.write_confidence_values, Configuration.width_branches, Configuration.use_style,
                    Configuration.write_events, Configuration.show_taxonomy_scientific_names,
                    Configuration.show_taxonomy_common_names, Configuration.show_tax_code, Configuration.show_tax_rank,
                    Configuration.show_seq_names, Configuration.show_gene_names, Configuration.show_sequence_acc,
                    Configuration.show_properties, Configuration.show_vector_data };
            for (final int which : expected) {
                if (!s.contains(which)) {
                    return false;
                }
            }
            // exact set: every present flag is an expected one too, so an over-report (a spurious flag
            // that would re-introduce an inert checkbox) is caught, not just under-reporting.
            if (s.size() != expected.length) {
                return false;
            }
            // spot-check a few that must be absent (subsumed by the size check, kept for intent).
            // show_mol_seqs in particular must NOT be reported even though the rich node HAS a molecular
            // sequence: there is no MSA viewer, so that checkbox is deliberately never built/scanned.
            if (s.contains(Configuration.show_domain_architectures) || s.contains(Configuration.show_seq_symbols)
                    || s.contains(Configuration.show_binary_characters) || s.contains(Configuration.show_mol_seqs)) {
                return false;
            }
            // regression (renderer match): a legitimate zero-length branch and a source-only accession
            // are both drawn, so both must be reported present.
            final Phylogeny edge = new Phylogeny();
            final PhylogenyNode er = new PhylogenyNode();
            final PhylogenyNode el = new PhylogenyNode();
            el.setName("e");
            el.setDistanceToParent(0.0);
            final Sequence es = new Sequence();
            es.setAccession(new Accession("", "UniProt"));
            el.getNodeData().setSequence(es);
            er.addAsChild(el);
            er.addAsChild(new PhylogenyNode());
            edge.setRoot(er);
            edge.externalNodesHaveChanged();
            final Set<Integer> edge_set = AptxUtil.scanForDataPresence(edge);
            if (!edge_set.contains(Configuration.write_branch_length_values)
                    || !edge_set.contains(Configuration.show_sequence_acc)) {
                return false;
            }
            // "Shorten Labels" is offered only when an external name is over-long (a pasted-in
            // UniProt/NCBI header); the short-named barren/rich trees above already assert it stays absent.
            final Phylogeny longlbl = new Phylogeny();
            final PhylogenyNode lr = new PhylogenyNode();
            final PhylogenyNode ll = new PhylogenyNode();
            ll.setName("tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical s-adenosyl methionine "
                    + "domain-containing protein OS=Fusarium phyllophilum OX=47803 GN=FPHYL_5758 PE=4 SV=1");
            lr.addAsChild(ll);
            lr.addAsChild(new PhylogenyNode());
            longlbl.setRoot(lr);
            longlbl.externalNodesHaveChanged();
            if (!AptxUtil.scanForDataPresence(longlbl).contains(Configuration.shorten_labels)) {
                return false;
            }
            return true;
        }
        catch (final Exception e) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * The launcher uses this to reject the removed {@code -c} configuration-file option: {@code -c} /
     * {@code -config} (any case, with or without an attached {@code =value}) match; tree files,
     * {@code -open}, empty, and {@code null} do not.
     */
    private static boolean testConfigFileOption() {
        for( final String yes : new String[] { "-c", "-config", "-C", "-Config", "-c=foo", "-config=my.txt" } ) {
            if ( !AptxUtil.isConfigFileOption( yes ) ) {
                return false;
            }
        }
        for( final String no : new String[] { null, "", "  ", "-open", "mytree.xml", "-cache", "c", "config" } ) {
            if ( AptxUtil.isConfigFileOption( no ) ) {
                return false;
            }
        }
        return true;
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
     * getRankChoices offers, in canonical order: every rank present on the tree (labeled
     * "rank (count) (P% coverage)") AND the major Linnaean ranks even when absent (labeled plainly) --
     * so a genus-only tree can still be colored by, e.g., order, which the taxonomy DB will resolve.
     */
    private static boolean testRankChoices() {
        // with no present ranks, the major ranks are still offered (plain), in canonical order
        final List<String> majors = Arrays.asList( AptxUtil.getRankChoices( new LinkedHashMap<String, Integer>(),
                                                                            null,
                                                                            0 ) );
        for( final String r : new String[] { "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus",
                "species" } ) {
            if ( !majors.contains( r ) ) {
                return fail( "major rank '" + r + "' must always be offered; got " + majors );
            }
        }
        if ( !( ( majors.indexOf( "order" ) < majors.indexOf( "family" ) )
                && ( majors.indexOf( "family" ) < majors.indexOf( "genus" ) )
                && ( majors.indexOf( "genus" ) < majors.indexOf( "species" ) ) ) ) {
            return fail( "ranks must be offered in canonical order; got " + majors );
        }
        // a present rank is labeled with count + coverage; absent majors stay plain
        final Map<String, Integer> counts = new LinkedHashMap<String, Integer>();
        counts.put( "order", 2 );
        final Map<String, Integer> coverage = new LinkedHashMap<String, Integer>();
        coverage.put( "order", 6 ); // 6 of 10 external nodes -> 60%
        final List<String> choices = Arrays.asList( AptxUtil.getRankChoices( counts, coverage, 10 ) );
        if ( !choices.contains( "order (2) (60% coverage)" ) ) {
            return fail( "a present rank must be labeled \"rank (count) (P% coverage)\"; got " + choices );
        }
        if ( !choices.contains( "family" ) || !choices.contains( "genus" ) ) {
            return fail( "absent major ranks must still be offered plainly; got " + choices );
        }
        // a present NON-major rank is offered too (labeled), even though it is not in the always-offer set
        final Map<String, Integer> sub = new LinkedHashMap<String, Integer>();
        sub.put( "subfamily", 3 );
        final List<String> sub_choices = Arrays.asList( AptxUtil.getRankChoices( sub, null, 10 ) );
        if ( !sub_choices.contains( "subfamily (3) (0% coverage)" ) ) {
            return fail( "a present non-major rank must be offered; got " + sub_choices );
        }
        // a present NON-major rank with count 0 must NOT be offered (it is neither colorizable nor a major)
        final Map<String, Integer> zero_sub = new LinkedHashMap<String, Integer>();
        zero_sub.put( "subfamily", 0 );
        final List<String> zero_choices = Arrays.asList( AptxUtil.getRankChoices( zero_sub, null, 10 ) );
        if ( zero_choices.contains( "subfamily (0) (0% coverage)" ) || zero_choices.contains( "subfamily" ) ) {
            return fail( "a present non-major rank with count 0 must not be offered; got " + zero_choices );
        }
        // two present ranks come out in canonical order (order before genus), regardless of map insertion order
        final Map<String, Integer> two = new LinkedHashMap<String, Integer>();
        two.put( "genus", 4 ); // inserted out of canonical order on purpose
        two.put( "order", 2 );
        final List<String> two_choices = Arrays.asList( AptxUtil.getRankChoices( two, null, 10 ) );
        if ( two_choices.indexOf( "order (2) (0% coverage)" ) >= two_choices.indexOf( "genus (4) (0% coverage)" ) ) {
            return fail( "present ranks must be emitted in canonical order (order before genus); got " + two_choices );
        }
        // coveragePercentLabel (still used to build the labels above)
        if ( !AptxUtil.coveragePercentLabel( 0, 100 ).equals( "0%" )
                || !AptxUtil.coveragePercentLabel( 1, 300 ).equals( "<1%" ) // 0.33% rounds to 0 -> <1%
                || !AptxUtil.coveragePercentLabel( 1, 200 ).equals( "1%" )  // 0.5% rounds half-up to 1%
                || !AptxUtil.coveragePercentLabel( 6, 10 ).equals( "60%" )
                || !AptxUtil.coveragePercentLabel( null, 10 ).equals( "0%" ) ) {
            return fail( "coveragePercentLabel: wrong percent label" );
        }
        // indexOfRank: pre-selecting the last-used rank must match a bare rank against count-suffixed choices
        final String[] mixed = { "phylum", "order (2) (60% coverage)", "family", "genus (4) (0% coverage)" };
        if ( ( AptxUtil.indexOfRank( mixed, "order" ) != 1 ) || ( AptxUtil.indexOfRank( mixed, "genus" ) != 3 )
                || ( AptxUtil.indexOfRank( mixed, "family" ) != 2 ) ) {
            return fail( "indexOfRank must find a bare rank whether the choice is plain or count-suffixed" );
        }
        if ( AptxUtil.indexOfRank( mixed, "ORDER" ) != 1 ) {
            return fail( "indexOfRank must be case-insensitive" );
        }
        if ( ( AptxUtil.indexOfRank( mixed, "species" ) != -1 ) || ( AptxUtil.indexOfRank( mixed, null ) != -1 )
                || ( AptxUtil.indexOfRank( mixed, "" ) != -1 ) || ( AptxUtil.indexOfRank( null, "order" ) != -1 ) ) {
            return fail( "indexOfRank must return -1 for an absent/empty rank or null choices (first-use case)" );
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

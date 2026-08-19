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
import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.io.parsers.nexus.NexusPhylogeniesParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Anti-rot guard for the {@code forester/demo/} gallery: every committed demo tree (see
 * {@link org.forester.demo.DemoTreeGenerator}) must still parse cleanly AND still carry the data its feature needs to
 * be demonstrable. Headless -- pure phyloXML parsing + model inspection, no GUI. If a parser or the property/taxonomy
 * model changes in a way that breaks a demo, this fails in the always-run suite instead of the demo silently rotting.
 *
 * <p>When a new demo tree is added (per the "one demo tree per new feature" convention), add a shape check here.
 */
public final class DemoTreesTest {

    private static final String DEMO_DIR = System.getProperty( "user.dir" ) + File.separator + "forester"
            + File.separator + "demo" + File.separator;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DemoTrees: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        boolean ok = true;
        // size-by: a numeric property to scale by
        ok &= hasNumericRef( "size-by-property.xml", "data:read_count" );
        // color-by: a categorical property (palette) AND a numeric one (gradient)
        ok &= hasCategoricalRef( "color-by-property.xml", "data:host" );
        ok &= hasNumericRef( "color-by-property.xml", "data:year" );
        // annotation columns: categorical + numeric + text properties
        ok &= hasCategoricalRef( "annotation-columns.xml", "data:host" );
        ok &= hasNumericRef( "annotation-columns.xml", "data:viral_load" );
        ok &= hasCategoricalRef( "annotation-columns.xml", "data:clade" );
        // colorize by rank / clade bands: each TIP carries its 'order' in-tree, so it colorizes OFFLINE (no prompt)
        // into 3 clades -- even though the Carnivora clade root is deliberately mis-annotated 'Rodentia' (a tip's
        // own identity wins over the wrong internal annotation)
        ok &= hasRank( "colorize-by-rank.xml", "order" );
        ok &= colorizesOfflineInto( "colorize-by-rank.xml", "order", 3 );
        // scale axis: a phylogram with real branch lengths (so the labeled distance axis has ticks)
        ok &= hasBranchLengths( "scale-axis.xml" );
        // node age bars: internal nodes carry a phyloXML <date> with a min/max interval + branch lengths (dated tree)
        ok &= hasBranchLengths( "node-hpd-bars.xml" );
        ok &= hasInternalDateInterval( "node-hpd-bars.xml" );
        ok &= isDetectedTimeTree( "node-hpd-bars.xml", AptxUtil.TIME_TREE_KIND.DATED ); // dated -> auto-labelled
        // zebra stripes: enough tips for alternating row bands to be meaningful + a categorical + numeric column to track
        ok &= hasAtLeastTips( "zebra-stripes.xml", 8 );
        ok &= hasCategoricalRef( "zebra-stripes.xml", "data:host" );
        ok &= hasNumericRef( "zebra-stripes.xml", "data:reads" );
        // reverse tip order: a ladder with several ordered tips so the tip-order reversal is unmistakable
        ok &= hasAtLeastTips( "reverse-tip-order.xml", 6 );
        // root-on-top orientation: a phylogram with several tips + branch lengths + internal support, so the vertical
        // dendrogram shows 45-degree tip labels and upright support/branch-length numbers
        ok &= hasAtLeastTips( "root-on-top.xml", 8 );
        ok &= hasBranchLengths( "root-on-top.xml" );
        ok &= hasInternalConfidence( "root-on-top.xml" );
        // search emphasis: enough tips + a searchable token shared by a subset (so a search highlights several) +
        // internal confidence values (so "Dim Non-Matches" fading the support numbers too is demonstrable)
        ok &= hasAtLeastTips( "domain-architectures.xml", 6 );
        ok &= hasBranchLengths( "domain-architectures.xml" );
        ok &= hasDomainArchitectures( "domain-architectures.xml", 6 );

        ok &= hasAtLeastTips( "heatmap-matrix.xml", 6 );
        ok &= hasNumericRef( "heatmap-matrix.xml", "data:s1" );
        ok &= hasNumericRef( "heatmap-matrix.xml", "data:s6" );

        // import annotations: a plain-named tip tree + a companion CSV whose "name" column must join onto EVERY tip
        ok &= hasAtLeastTips( "import-annotations.xml", 12 );
        ok &= csvJoinMatchesAllTips( "import-annotations.xml", "import-annotations.csv" );

        ok &= hasAtLeastTips( "search-emphasis.xml", 12 );
        ok &= tipsContaining( "search-emphasis.xml", "kinase", 4 );
        ok &= hasInternalConfidence( "search-emphasis.xml" );

        // node style: several tips carry a per-node NodeVisualData (font + node mark), for "Use Visual Styles" and
        // the Node Style editor (click-to / Tools)
        ok &= hasAtLeastTips( "node-visual-styles.xml", 5 );
        ok &= hasNodeVisualStyle( "node-visual-styles.xml", 3 );
        ok &= hasNodeShape( "node-visual-styles.xml", NodeVisualData.NodeShape.DIAMOND ); // round-trips the diamond

        // infer ancestor taxonomies: six real-species tips (rank 'species'), NO taxonomy on the internal nodes,
        // ready for Analysis > Infer Ancestor Taxonomies (which resolves the tips online and fills the internals)
        ok &= hasAtLeastTips( "infer-ancestor-taxonomies.xml", 6 );
        ok &= hasRank( "infer-ancestor-taxonomies.xml", "species" );
        ok &= internalNodesHaveNoTaxonomy( "infer-ancestor-taxonomies.xml" );

        // BEAST / BEAST X output: a NEXUS tree whose [&...] annotations parse into HPD date intervals (Node Age
        // Bars), posterior confidences (support), and a numeric beast:rate property (Color-by)
        ok &= beastAnnotationsOk( "beast-annotations.nex" );

        // ancestral-state pie charts: a discrete-trait posterior (beast:location_set + _set_prob) on the internal
        // nodes, so the pie feature has a trait with a parseable multi-state distribution
        ok &= hasAncestralStateTrait( "ancestral-pie-charts.xml", "location" );

        // create tanglegram: a linked PAIR of trees over the same eight taxa (different topologies) -- every tip of
        // one must link to a tip of the other on the scientific-name key (a fully-connected, crossing tanglegram)
        ok &= hasAtLeastTips( "tanglegram-tree-a.xml", 8 );
        ok &= hasAtLeastTips( "tanglegram-tree-b.xml", 8 );
        ok &= tanglegramPairLinks( "tanglegram-tree-a.xml", "tanglegram-tree-b.xml" );
        // a categorical 'clade' group per tip, so the tanglegram's connectors can be coloured by clade
        ok &= hasCategoricalRef( "tanglegram-tree-a.xml", "data:clade" );
        // cross-field tanglegram: a gene tree (accession names, species in the taxonomy) + a species tree (species
        // names) -- the same species in DIFFERENT fields, so they link only on per-tree fields (sci-name <-> node-name)
        ok &= hasAtLeastTips( "tanglegram-gene-tree.xml", 8 );
        ok &= hasAtLeastTips( "tanglegram-species-tree.xml", 8 );
        ok &= tanglegramCrossFieldLinks( "tanglegram-gene-tree.xml", "tanglegram-species-tree.xml" );
        // association tanglegram: a host tree + a parasite tree whose tip names share NOTHING, linked through a
        // two-column mapping file -- so linkByAssociation fully connects them while a value join links nothing
        ok &= hasAtLeastTips( "tanglegram-host-tree.xml", 6 );
        ok &= hasAtLeastTips( "tanglegram-parasite-tree.xml", 6 );
        ok &= tanglegramAssociationLinks( "tanglegram-host-tree.xml", "tanglegram-parasite-tree.xml",
                                          "tanglegram-association.tsv" );
        // geologic time axis: a dated dinosaur phylogram (branch lengths = time in Ma; internal nodes carry a
        // <date>), detected as a DATED time tree, with extant tips at age 0 -- so the geologic time axis maps ages
        // exactly and the Cretaceous/Jurassic/Triassic period bands render beneath it
        ok &= hasBranchLengths( "dinosaur-time-tree.xml" );
        ok &= hasAtLeastTips( "dinosaur-time-tree.xml", 7 );
        ok &= isDetectedTimeTree( "dinosaur-time-tree.xml", AptxUtil.TIME_TREE_KIND.DATED );
        return ok;
    }

    /** The host/parasite demo trees (whose tip names share nothing) fully link through the association file, and a
     *  plain value join links nothing -- exactly the parasite-vs-host case the association feature exists for. */
    private static boolean tanglegramAssociationLinks( final String host_file, final String parasite_file,
                                                       final String assoc_file ) {
        final Phylogeny hosts = load( host_file );
        final Phylogeny lice = load( parasite_file );
        if ( ( hosts == null ) || ( lice == null ) ) {
            return false;
        }
        final java.util.Map<String, java.util.List<String>> assoc;
        try {
            assoc = TanglegramAssociation
                    .parse( java.nio.file.Files.readString( new File( DEMO_DIR + assoc_file ).toPath() ) ).leftToRight();
        }
        catch ( final Exception e ) {
            return note( "could not read association file " + assoc_file + ": " + e.getMessage() );
        }
        final TanglegramLinker.Result r = TanglegramLinker.linkByAssociation( hosts, lice,
                                                                              TanglegramLinker.LinkField.NODE_NAME,
                                                                              TanglegramLinker.LinkField.NODE_NAME,
                                                                              assoc );
        if ( ( r.getLinks().size() < 6 ) || !r.getUnmatchedA().isEmpty() || !r.getUnmatchedB().isEmpty() ) {
            return note( host_file + "/" + parasite_file + " should fully link through " + assoc_file + "; found "
                    + r.getLinks().size() + " links, " + r.getUnmatchedA().size() + "/" + r.getUnmatchedB().size()
                    + " unmatched" );
        }
        if ( !TanglegramLinker.link( hosts, lice, TanglegramLinker.LinkField.NODE_NAME ).getLinks().isEmpty() ) {
            return note( "the host/parasite tip names should NOT value-join (the point of the association demo)" );
        }
        return true;
    }

    /** The gene/species demo pair links only when a field is chosen PER tree (gene tree's scientific name to species
     *  tree's node name): every tip matches with none unmatched, and NO single shared field links them. */
    private static boolean tanglegramCrossFieldLinks( final String gene_file, final String species_file ) {
        final Phylogeny gene = load( gene_file );
        final Phylogeny species = load( species_file );
        if ( ( gene == null ) || ( species == null ) ) {
            return false;
        }
        final TanglegramLinker.Result r = TanglegramLinker.link( gene, species,
                                                                 TanglegramLinker.LinkField.SCIENTIFIC_NAME,
                                                                 TanglegramLinker.LinkField.NODE_NAME );
        if ( ( r.getLinks().size() < 8 ) || !r.getUnmatchedA().isEmpty() || !r.getUnmatchedB().isEmpty() ) {
            return note( gene_file + " and " + species_file
                    + " should fully link on scientific-name <-> node-name; found " + r.getLinks().size() + " links, "
                    + r.getUnmatchedA().size() + "/" + r.getUnmatchedB().size() + " unmatched" );
        }
        // the demo only makes sense if a single shared field does NOT link them (else per-tree fields aren't needed)
        if ( !TanglegramLinker.link( gene, species, TanglegramLinker.LinkField.NODE_NAME ).getLinks().isEmpty()
                || !TanglegramLinker.link( gene, species, TanglegramLinker.LinkField.SCIENTIFIC_NAME ).getLinks()
                        .isEmpty() ) {
            return note( gene_file + "/" + species_file + " should NOT link on any single shared field (the point of "
                    + "the cross-field demo)" );
        }
        return true;
    }

    /** Every tip of one tanglegram demo tree links to a tip of the other on the scientific-name key (no tip left
     *  unmatched), so the Create Tanglegram demo shows a fully-connected, crossing tanglegram. */
    private static boolean tanglegramPairLinks( final String a_file, final String b_file ) {
        final Phylogeny a = load( a_file );
        final Phylogeny b = load( b_file );
        if ( ( a == null ) || ( b == null ) ) {
            return false;
        }
        final TanglegramLinker.Result r = TanglegramLinker.link( a, b, TanglegramLinker.LinkField.SCIENTIFIC_NAME );
        if ( r.getLinks().size() < 8 ) {
            return note( a_file + " and " + b_file + " should share >= 8 scientific-name links, found "
                    + r.getLinks().size() );
        }
        if ( !r.getUnmatchedA().isEmpty() || !r.getUnmatchedB().isEmpty() ) {
            return note( "tanglegram demo pair should leave no tip unmatched (a: " + r.getUnmatchedA().size()
                    + ", b: " + r.getUnmatchedB().size() + ")" );
        }
        return true;
    }

    /** The ancestral-pie demo offers {@code trait} (a beast:&lt;trait&gt;_set_prob property) AND at least one internal
     *  node has a real multi-state (&ge;2) distribution, so a pie with several wedges is demonstrable. */
    private static boolean hasAncestralStateTrait( final String file_name, final String trait ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( !TreePanelUtil.ancestralStateTraits( phy ).contains( trait ) ) {
            return note( file_name + " must offer ancestral-state trait '" + trait + "' (beast:" + trait
                    + "_set_prob on the internal nodes)" );
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && ( TreePanelUtil.stateDistribution( n, trait ).size() >= 2 ) ) {
                return true;
            }
        }
        return note( file_name + " must have an internal node with a >=2-state distribution for '" + trait + "'" );
    }

    private static boolean hasInternalConfidence( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && !n.getBranchData().getConfidences().isEmpty() ) {
                return true;
            }
        }
        return note( file_name + " must carry an internal-node confidence value (for the dim-the-numbers demo)" );
    }

    private static boolean tipsContaining( final String file_name, final String token, final int min_matches ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        int matches = 0;
        for( final PhylogenyNode leaf : phy.getExternalNodes() ) {
            if ( ( leaf.getName() != null ) && leaf.getName().contains( token ) ) {
                ++matches;
            }
        }
        if ( matches < min_matches ) {
            return note( file_name + " must have at least " + min_matches + " tips whose name contains '" + token
                    + "' (a searchable subset for the emphasis demo), found " + matches );
        }
        return true;
    }

    /** The demo tree is detected as the expected kind of time tree (auto-label vs offer vs none). */
    private static boolean isDetectedTimeTree( final String file_name, final AptxUtil.TIME_TREE_KIND expected ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        final AptxUtil.TIME_TREE_KIND got = AptxUtil.detectTimeTree( phy );
        if ( got != expected ) {
            return note( file_name + " must be detected as " + expected + " time tree, was " + got );
        }
        return true;
    }

    /** At least one tip carries the given node {@code shape} (round-trips the shape string through phyloXML). */
    private static boolean hasNodeShape( final String file_name, final NodeVisualData.NodeShape shape ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for ( final PhylogenyNode n : phy.getExternalNodes() ) {
            final NodeVisualData vis = n.getNodeData().getNodeVisualData();
            if ( ( vis != null ) && ( vis.getShape() == shape ) ) {
                return true;
            }
        }
        return note( file_name + " must have a tip with node shape " + shape );
    }

    /** At least {@code min} tips carry a non-empty per-node visual style ({@link NodeVisualData}). */
    private static boolean hasNodeVisualStyle( final String file_name, final int min ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        int styled = 0;
        for ( final PhylogenyNode n : phy.getExternalNodes() ) {
            final NodeVisualData vis = n.getNodeData().getNodeVisualData();
            if ( ( vis != null ) && !vis.isEmpty() ) {
                styled++;
            }
        }
        if ( styled < min ) {
            return note( file_name + " must have at least " + min + " tips with a visual style; found " + styled );
        }
        return true;
    }

    private static boolean hasAtLeastTips( final String file_name, final int min_tips ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( phy.getNumberOfExternalNodes() < min_tips ) {
            return note( file_name + " must have at least " + min_tips + " tips (for meaningful zebra row bands)" );
        }
        return true;
    }

    private static boolean hasInternalDateInterval( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                    && ( n.getNodeData().getDate().getMax() != null ) ) {
                return true;
            }
        }
        return note( file_name + " must carry an INTERNAL-node <date> with a min/max interval (for HPD age bars)" );
    }

    /** At least {@code min_tips} external tips carry a sequence with a domain architecture of &ge;1 domain. */
    private static boolean hasDomainArchitectures( final String file_name, final int min_tips ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        int n = 0;
        for( final PhylogenyNode tip : phy.getExternalNodes() ) {
            if ( tip.getNodeData().isHasSequence()
                    && ( tip.getNodeData().getSequence().getDomainArchitecture() != null )
                    && ( tip.getNodeData().getSequence().getDomainArchitecture().getDomains().size() > 0 ) ) {
                ++n;
            }
        }
        if ( n < min_tips ) {
            return note( file_name + " must carry domain architectures on >= " + min_tips + " tips, found " + n );
        }
        return true;
    }

    private static boolean hasBranchLengths( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( !( org.forester.phylogeny.PhylogenyMethods.calculateMaxDistanceToRoot( phy ) > 0.0 ) ) {
            return note( file_name + " must carry branch lengths (a positive max distance to root) for a scale axis" );
        }
        return true;
    }

    /** Parses a demo file, asserting it is present, error-free and non-empty; null on any failure (with a message). */
    /** The companion CSV parses and its key column joins onto every tip of the tree (no unmatched rows, no tip left out). */
    private static boolean csvJoinMatchesAllTips( final String tree_file, final String csv_file ) {
        final Phylogeny phy = load( tree_file );
        if ( phy == null ) {
            return false;
        }
        final File csv = new File( DEMO_DIR + csv_file );
        if ( !csv.exists() ) {
            return note( csv_file + " is missing from the demo gallery (" + csv.getAbsolutePath() + ")" );
        }
        try {
            final NodeDataImporter.Table table = NodeDataImporter.parseTable( java.nio.file.Files.readString( csv.toPath() ) );
            final NodeDataImporter.MatchReport rep = NodeDataImporter.dryRun( phy, table, table.defaultKeyColumn(),
                    NodeDataImporter.MatchBy.TIP_NAME );
            if ( ( rep.getTipsWithoutRow() != 0 ) || ( rep.getRowsUnmatched() != 0 ) ) {
                return note( csv_file + " should join onto every tip of " + tree_file + " (tips without a row: "
                        + rep.getTipsWithoutRow() + ", unmatched rows: " + rep.getRowsUnmatched() + ")" );
            }
            if ( !rep.getPropertyColumns().contains( "host" ) || !rep.getPropertyColumns().contains( "reads" ) ) {
                return note( csv_file + " should import the host + reads columns: " + rep.getPropertyColumns() );
            }
            return true;
        }
        catch ( final Exception e ) {
            return note( csv_file + " could not be read/joined: " + e.getMessage() );
        }
    }

    private static Phylogeny load( final String file_name ) {
        final File file = new File( DEMO_DIR + file_name );
        if ( !file.exists() ) {
            return fail( file_name + " is missing from the demo gallery (" + file.getAbsolutePath() + ")" );
        }
        try {
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance().create( file, parser );
            if ( parser.getErrorCount() > 0 ) {
                return fail( file_name + " parsed with errors: " + parser.getErrorMessages() );
            }
            if ( ( phys == null ) || ( phys.length != 1 ) || ( phys[ 0 ] == null ) || phys[ 0 ].isEmpty() ) {
                return fail( file_name + " did not yield exactly one non-empty tree" );
            }
            return phys[ 0 ];
        }
        catch ( final Exception e ) {
            return fail( file_name + " could not be read: " + e.getMessage() );
        }
    }

    private static boolean hasNumericRef( final String file_name, final String ref ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( !PropertyColorScheme.numericRefs( phy ).contains( ref ) ) {
            return note( file_name + " must offer numeric property '" + ref + "' for its feature demo" );
        }
        return true;
    }

    private static boolean hasCategoricalRef( final String file_name, final String ref ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        boolean present = false;
        for( final PhylogenyNode leaf : phy.getExternalNodes() ) {
            if ( PropertyColorScheme.valueFor( leaf, ref ) != null ) {
                present = true;
                break;
            }
        }
        if ( !present ) {
            return note( file_name + " must carry property '" + ref + "' on at least one tip" );
        }
        // 'categorical' means it is NOT treated as a numeric gradient
        if ( PropertyColorScheme.numericRefs( phy ).contains( ref ) ) {
            return note( file_name + " property '" + ref + "' should be categorical, not numeric" );
        }
        return true;
    }

    private static boolean hasRank( final String file_name, final String rank ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() && rank.equals( n.getNodeData().getTaxonomy().getRank() ) ) {
                return true;
            }
        }
        return note( file_name + " must carry an in-tree taxonomy at rank '" + rank + "' (offline colorize)" );
    }

    /** The BEAST NEXUS demo parses its [&...] annotations into structured data: internal-node HPD date intervals,
     *  posterior confidences, and a numeric beast:rate property on every node. */
    private static boolean beastAnnotationsOk( final String file_name ) {
        final Phylogeny phy = loadNexusBeast( file_name );
        if ( phy == null ) {
            return false;
        }
        int intervals = 0;
        int posteriors = 0;
        int rates = 0;
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getMin() != null )
                    && ( n.getNodeData().getDate().getMax() != null ) ) {
                intervals++;
            }
            if ( !n.isExternal() && n.getBranchData().isHasConfidences() ) {
                posteriors++;
            }
            if ( ( n.getNodeData().getProperties() != null )
                    && !n.getNodeData().getProperties().getProperties( "beast:rate" ).isEmpty() ) {
                rates++;
            }
        }
        if ( intervals < 3 ) {
            return note( file_name + " internal nodes must carry HPD date intervals (height_95%_HPD), got " + intervals );
        }
        if ( posteriors < 3 ) {
            return note( file_name + " internal nodes must carry posterior confidences, got " + posteriors );
        }
        if ( rates < 5 ) {
            return note( file_name + " nodes must carry a numeric beast:rate property, got " + rates );
        }
        return true;
    }

    /** Load a NEXUS demo file with BEAST-style annotation parsing on (the phyloXML {@link #load} helper is XML-only). */
    private static Phylogeny loadNexusBeast( final String file_name ) {
        final File file = new File( DEMO_DIR + file_name );
        if ( !file.exists() ) {
            return fail( file_name + " is missing from the demo gallery (" + file.getAbsolutePath() + ")" );
        }
        try {
            final NexusPhylogeniesParser parser = new NexusPhylogeniesParser();
            parser.setParseBeastStyleExtendedTags( true );
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance().create( file, parser );
            if ( ( phys == null ) || ( phys.length != 1 ) || ( phys[ 0 ] == null ) || phys[ 0 ].isEmpty() ) {
                return fail( file_name + " did not yield exactly one non-empty tree" );
            }
            return phys[ 0 ];
        }
        catch ( final Exception e ) {
            return fail( file_name + " could not be read: " + e.getMessage() );
        }
    }

    /** No internal node carries a taxonomy, so the tree is ready for "Infer Ancestor Taxonomies" to fill them. */
    private static boolean internalNodesHaveNoTaxonomy( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.isExternal() && n.getNodeData().isHasTaxonomy() ) {
                return note( file_name + " internal nodes must carry NO taxonomy (ready for ancestral inference)" );
            }
        }
        return true;
    }

    /** Anti-rot guard for the colorize-by-rank demo: it colorizes OFFLINE (null service) at {@code rank} into
     *  exactly {@code groups} maximal clades / {@code groups} distinct legend taxa, with NO tip needing an online
     *  fetch (every tip self-resolves the rank in-tree). This does NOT exercise the Spine B "tip identity wins over
     *  a wrong ancestor" fix -- that needs a DB, and the demo's tips self-resolve at step 1 (so pre-fix code would
     *  pass this too); the fix guard is {@code TreePanelUtilTest.testTipIdentityWins}. */
    private static boolean colorizesOfflineInto( final String file_name, final String rank, final int groups ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( !TreePanelUtil.unresolvedTipTaxa( phy, rank, null ).isEmpty() ) {
            return note( file_name + " must colorize at '" + rank + "' OFFLINE (no tip should need an online fetch)" );
        }
        final Map<String, Color> legend = new HashMap<String, Color>();
        final int colored = TreePanelUtil.colorPhylogenyAccordingToRanks( phy, rank, null, legend );
        // groups maximal clades AND groups distinct legend taxa: each tip carries its own order, so each order
        // forms exactly one clade. (An edit that dropped the tips' self-resolving order would change these counts.)
        if ( ( colored != groups ) || ( legend.size() != groups ) ) {
            return note( file_name + " must colorize into " + groups + " clades / " + groups + " distinct taxa at '"
                    + rank + "', got " + colored + " clades / " + legend.size() + " taxa " + legend.keySet() );
        }
        return true;
    }

    /** null-returning failure note for the parse path (so load(...) can bail with a message). */
    private static Phylogeny fail( final String msg ) {
        note( msg );
        return null;
    }

    /** false-returning failure note for the shape checks. */
    private static boolean note( final String msg ) {
        System.out.println( "  [DemoTreesTest] " + msg );
        return false;
    }

    private DemoTreesTest() {
    }
}

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
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
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
        return ok;
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

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

package org.forester.demo;

import java.io.File;
import java.io.IOException;
import java.util.Locale;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import java.util.ArrayList;
import java.util.List;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.NodeVisualData;
import org.forester.phylogeny.data.NodeVisualData.NodeFill;
import org.forester.phylogeny.data.NodeVisualData.NodeShape;
import java.awt.Color;
import java.awt.Font;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Regenerates the synthetic demo trees under {@code forester/demo/}. Each tree is named after the feature it
 * demonstrates and is crafted so a user can load it (File &gt; Open) and immediately try that feature -- see
 * {@code forester/demo/README.md}. Trees are small, synthetic (no data-licensing issues) and deliberately shaped to
 * show the feature at its best.
 *
 * <p>This is a developer tool: it is run once and its committed {@code .xml} output is what ships in the demo gallery
 * (so users need no build step). {@code org.forester.archaeopteryx.DemoTreesTest} smoke-loads every generated file so
 * the demos cannot silently rot. Run from the repo root: {@code java org.forester.demo.DemoTreeGenerator}.
 */
public final class DemoTreeGenerator {

    // exemplary phyloXML types: real numbers are xsd:decimal, categories / free text are xsd:string
    private static final String NUM = "xsd:decimal";
    private static final String TEXT = "xsd:string";

    public static void main( final String[] args ) throws IOException, PhyloXmlDataFormatException {
        final File dir = new File( System.getProperty( "user.dir" ), "forester/demo" );
        if ( !dir.exists() && !dir.mkdirs() ) {
            System.err.println( "could not create demo directory: " + dir.getAbsolutePath() );
            System.exit( 1 );
        }
        write( dir, "size-by-property.xml", sizeByPropertyTree() );
        write( dir, "color-by-property.xml", colorByPropertyTree() );
        write( dir, "annotation-columns.xml", annotationColumnsTree() );
        write( dir, "colorize-by-rank.xml", colorizeByRankTree() );
        write( dir, "infer-ancestor-taxonomies.xml", inferAncestorTaxonomiesTree() );
        write( dir, "scale-axis.xml", scaleAxisTree() );
        write( dir, "node-hpd-bars.xml", hpdBarsTree() );
        write( dir, "zebra-stripes.xml", zebraStripesTree() );
        write( dir, "reverse-tip-order.xml", reverseTipOrderTree() );
        write( dir, "root-on-top.xml", rootOnTopTree() );
        write( dir, "domain-architectures.xml", domainArchitecturesTree() );
        write( dir, "heatmap-matrix.xml", heatmapMatrixTree() );
        write( dir, "import-annotations.xml", importAnnotationsTree() );
        writeText( dir, "import-annotations.csv", importAnnotationsCsv() );
        write( dir, "search-emphasis.xml", searchEmphasisTree() );
        write( dir, "node-visual-styles.xml", nodeVisualStylesTree() );
        writeText( dir, "beast-annotations.nex", beastAnnotationsNexus() );
        write( dir, "ancestral-pie-charts.xml", ancestralPieChartsTree() );
        write( dir, "tanglegram-tree-a.xml", tanglegramTreeA() );
        write( dir, "tanglegram-tree-b.xml", tanglegramTreeB() );
        write( dir, "tanglegram-gene-tree.xml", tanglegramGeneTree() );
        write( dir, "tanglegram-species-tree.xml", tanglegramSpeciesTree() );
        write( dir, "tanglegram-host-tree.xml", tanglegramHostTree() );
        write( dir, "tanglegram-parasite-tree.xml", tanglegramParasiteTree() );
        writeText( dir, "tanglegram-association.tsv", tanglegramAssociationTsv() );
        System.out.println( "Wrote demo trees to " + dir.getAbsolutePath() );
    }

    private static void write( final File dir, final String file_name, final Phylogeny phy ) throws IOException {
        new PhylogenyWriter().toPhyloXML( phy, 0, new File( dir, file_name ) );
        System.out.println( "  " + file_name + " (" + phy.getNumberOfExternalNodes() + " tips)" );
    }

    // ----- "BEAST / BEAST X output": a small time-calibrated NEXUS tree with FigTree-style [&...] annotations
    //       (posterior support, node height + 95% HPD interval, per-branch rate). Open it to see the annotations
    //       parsed into HPD bars, posterior support, and a Color-by 'beast:rate' gradient. Hand-authored NEXUS
    //       (not phyloXML), written verbatim -- an ultrametric 5-tip tree (tip heights 0, root height 2.1).
    private static String beastAnnotationsNexus() {
        return "#NEXUS\n" + "BEGIN TREES;\n" + "\tTREE beast_demo = [&R] ("
                + "(isolate_A[&height=0.0,rate=0.0031]:1.2,isolate_B[&height=0.0,rate=0.0028]:1.2)"
                + "[&posterior=0.99,height=1.2,height_95%_HPD={0.95,1.5},rate=0.0030]:0.9,"
                + "(isolate_C[&height=0.0,rate=0.0026]:0.8,"
                + "(isolate_D[&height=0.0,rate=0.0035]:0.5,isolate_E[&height=0.0,rate=0.0033]:0.5)"
                + "[&posterior=0.92,height=0.5,height_95%_HPD={0.35,0.7},rate=0.0034]:0.3)"
                + "[&posterior=0.81,height=0.8,height_95%_HPD={0.6,1.1},rate=0.0029]:1.3)"
                + "[&posterior=1.0,height=2.1,height_95%_HPD={1.8,2.5},rate=0.0030];\n" + "END;\n";
    }

    // ----- "Ancestral-state pie charts": a BEAST-style phylogeography tree. Each TIP carries its single sampled
    //       location (beast:location); each INTERNAL node carries a posterior distribution over locations
    //       (beast:location_set + beast:location_set_prob), as a BEAST discrete-trait analysis produces. Load ->
    //       "Ancestral pie: location" to draw a per-node pie (wedges = state probabilities; tips = solid discs).
    private static Phylogeny ancestralPieChartsTree() {
        // three-tip regional clades; each internal node a posterior over {Africa, Americas, Asia, Europe}
        final PhylogenyNode africa = pieInternal( 0.05, "{Africa,Europe}", "{0.9,0.1}",
                locTip( "EPI_ISL_401", "Africa" ), locTip( "EPI_ISL_402", "Africa" ), locTip( "EPI_ISL_403", "Africa" ) );
        final PhylogenyNode europe = pieInternal( 0.05, "{Europe,Africa,Asia}", "{0.7,0.2,0.1}",
                locTip( "EPI_ISL_511", "Europe" ), locTip( "EPI_ISL_512", "Europe" ), locTip( "EPI_ISL_513", "Europe" ) );
        final PhylogenyNode asia = pieInternal( 0.05, "{Asia,Europe}", "{0.8,0.2}",
                locTip( "EPI_ISL_621", "Asia" ), locTip( "EPI_ISL_622", "Asia" ), locTip( "EPI_ISL_623", "Asia" ) );
        final PhylogenyNode americas = pieInternal( 0.05, "{Americas,Europe,Asia}", "{0.6,0.3,0.1}",
                locTip( "EPI_ISL_731", "Americas" ), locTip( "EPI_ISL_732", "Americas" ),
                locTip( "EPI_ISL_733", "Americas" ) );
        // deeper internal nodes: increasingly mixed toward the (African) root
        final PhylogenyNode eurasia = pieInternal( 0.04, "{Europe,Asia,Africa}", "{0.45,0.35,0.2}", europe, asia );
        final PhylogenyNode old_world = pieInternal( 0.04, "{Africa,Europe,Asia}", "{0.55,0.3,0.15}", africa, eurasia );
        final PhylogenyNode root = pieInternal( 0.0, "{Africa,Europe,Asia,Americas}", "{0.5,0.25,0.15,0.1}",
                old_world, americas );
        return tree( root, "Ancestral-state pie charts (demo)",
                     "Synthetic BEAST phylogeography tree. Each TIP carries its single sampled location "
                             + "(beast:location); each INTERNAL node carries a posterior distribution over locations "
                             + "(beast:location_set + beast:location_set_prob), as a BEAST discrete-trait / "
                             + "phylogeography analysis produces. Pick \"Ancestral pie: location\" to draw an "
                             + "ancestral-state pie at each node -- wedges sized by state probability, tips as solid "
                             + "single-state discs -- with a state->color legend." );
    }

    /** An internal node with a branch length and a BEAST discrete-trait posterior: a brace-wrapped state set + a
     *  matching probability set (as BEAST writes them), under beast:&lt;trait&gt;_set / beast:&lt;trait&gt;_set_prob. */
    private static PhylogenyNode pieInternal( final double bl, final String state_set, final String prob_set,
                                              final PhylogenyNode... kids ) {
        final PhylogenyNode n = clade( bl, kids );
        cat( n, "beast:location_set", state_set );
        cat( n, "beast:location_set_prob", prob_set );
        return n;
    }

    /** A tip carrying a single observed location (beast:location), rendered as a solid single-state pie disc. */
    private static PhylogenyNode locTip( final String name, final String location ) {
        final PhylogenyNode n = leaf( name );
        cat( n, "beast:location", location );
        return n;
    }

    /** Write a companion plain-text data file (e.g. a CSV to import onto a demo tree). */
    private static void writeText( final File dir, final String file_name, final String content ) throws IOException {
        java.nio.file.Files.writeString( new File( dir, file_name ).toPath(), content );
        System.out.println( "  " + file_name + " (" + content.lines().count() + " lines)" );
    }

    // ----- "Import annotations": a plain-named tip tree with NO per-tip data, paired with a companion CSV to
    //       import onto it (File > Import Annotations). Demonstrates the CSV/TSV keyed join.
    private static Phylogeny importAnnotationsTree() {
        final PhylogenyNode root = clade( 0.0,
                clade( 0.05, leaf( "isolate_01" ), leaf( "isolate_02" ), leaf( "isolate_03" ) ),
                clade( 0.06, leaf( "isolate_04" ), leaf( "isolate_05" ), leaf( "isolate_06" ) ),
                clade( 0.05, leaf( "isolate_07" ), leaf( "isolate_08" ),
                       clade( 0.03, leaf( "isolate_09" ), leaf( "isolate_10" ) ) ),
                clade( 0.06, leaf( "isolate_11" ), leaf( "isolate_12" ) ) );
        return tree( root, "Import annotations (demo)",
                     "A plain 12-tip tree (isolate_01..isolate_12) carrying NO per-tip data. Pair it with the "
                             + "companion import-annotations.csv: File > Import Annotations, pick the CSV, match the "
                             + "\"name\" column against the tip name; the host/country/reads columns are joined onto "
                             + "the tips -- then Color by host, or show reads as an annotation column." );
    }

    /** The companion table for {@link #importAnnotationsTree()} (one quoted field with an embedded comma, for CSV realism). */
    private static String importAnnotationsCsv() {
        return "name,host,country,reads\n"
                + "isolate_01,mosquito,USA,1200\n"
                + "isolate_02,mosquito,USA,980\n"
                + "isolate_03,bat,China,1543\n"
                + "isolate_04,bat,China,760\n"
                + "isolate_05,pig,Vietnam,410\n"
                + "isolate_06,pig,Vietnam,1330\n"
                + "isolate_07,bird,\"Congo, DR\",275\n"
                + "isolate_08,bird,Kenya,1890\n"
                + "isolate_09,mosquito,Brazil,640\n"
                + "isolate_10,mosquito,Brazil,1120\n"
                + "isolate_11,pig,Vietnam,505\n"
                + "isolate_12,bat,China,1450\n";
    }

    // ----- "Size by property": one numeric property (sequencing read count) spanning ~3 orders of magnitude, so the
    //       area-proportional tip dots vary strongly. Load -> "Size by: read_count".
    private static Phylogeny sizeByPropertyTree() {
        final PhylogenyNode a = sizedLeaf( "A/HongKong/1997", 120 );
        final PhylogenyNode b = sizedLeaf( "A/goose/Guangdong/1996", 5200 );
        final PhylogenyNode c = sizedLeaf( "A/Thailand/2004", 450 );
        final PhylogenyNode d = sizedLeaf( "A/Vietnam/2004", 1800 );
        final PhylogenyNode e = sizedLeaf( "A/duck/Hunan/2016", 9800 );
        final PhylogenyNode f = sizedLeaf( "A/chicken/Jiangsu/2016", 24000 );
        final PhylogenyNode g = sizedLeaf( "A/Shanghai/2013", 51000 );
        final PhylogenyNode h = sizedLeaf( "A/Anhui/2013", 88000 );
        final PhylogenyNode root = clade( 0, clade( 0.05, a, b ),
                                          clade( 0.04, clade( 0.03, c, d ),
                                                 clade( 0.03, clade( 0.02, e, f ), clade( 0.02, g, h ) ) ) );
        return tree( root, "Size by property (demo)",
                     "Synthetic avian-influenza tree. Each tip carries a numeric 'read_count' property. "
                             + "Try Size by: read_count -- the tip symbol AREA is proportional to the value." );
    }

    // ----- "Node Style": a few tips carry a per-node visual style (font style/size/colour + node shape/fill/size/
    //       colour). Turn on "Use Visual Styles" to see them, then click a node -> "Node Style" (or select/search
    //       nodes -> Tools -> "Node Style for Selected Nodes...") to change any node's style.
    private static Phylogeny nodeVisualStylesTree() {
        final PhylogenyNode a = styledLeaf( "Homo sapiens", Font.BOLD, 18, new Color( 0xD0, 0x00, 0x00 ),
                                            NodeShape.CIRCLE, NodeFill.SOLID, 14f, new Color( 0xD0, 0x00, 0x00 ) );
        final PhylogenyNode b = styledLeaf( "Pan troglodytes", Font.ITALIC, 14, new Color( 0x00, 0x66, 0xCC ),
                                            NodeShape.RECTANGLE, NodeFill.SOLID, 10f, new Color( 0x00, 0x66, 0xCC ) );
        final PhylogenyNode c = leaf( "Gorilla gorilla" ); // unstyled -- the default look, for comparison
        final PhylogenyNode d = styledLeaf( "Mus musculus", Font.PLAIN, 12, new Color( 0x00, 0x88, 0x00 ),
                                            NodeShape.CIRCLE, NodeFill.NONE, 9f, new Color( 0x00, 0x88, 0x00 ) );
        final PhylogenyNode e = leaf( "Rattus norvegicus" ); // unstyled
        final PhylogenyNode root = clade( 0, clade( 0.05, clade( 0.03, a, b ), c ), clade( 0.04, d, e ) );
        return tree( root, "Node Style (demo)",
                     "A few tips carry a per-node visual style (font + node mark). Turn on \"Use Visual Styles\", "
                             + "then click a node -> \"Node Style\", or select/search nodes and use Tools -> \"Node "
                             + "Style for Selected Nodes...\", to change the font (style/size/colour) and node mark "
                             + "(shape/fill/size/colour)." );
    }

    private static PhylogenyNode styledLeaf( final String name, final int font_style, final int font_size,
                                             final Color font_color, final NodeShape shape, final NodeFill fill,
                                             final float node_size, final Color node_color ) {
        final PhylogenyNode n = leaf( name );
        final NodeVisualData vis = new NodeVisualData();
        vis.setFontName( "Source Sans 3" );
        vis.setFontStyle( font_style );
        vis.setFontSize( font_size );
        vis.setFontColor( font_color );
        vis.setShape( shape );
        vis.setFillType( fill );
        vis.setSize( node_size );
        vis.setNodeColor( node_color );
        n.getNodeData().setNodeVisualData( vis );
        return n;
    }

    private static PhylogenyNode sizedLeaf( final String name, final int read_count ) {
        final PhylogenyNode n = leaf( name );
        num( n, "data:read_count", Integer.toString( read_count ) );
        return n;
    }

    // ----- "Color by property": a categorical property (host) for the palette AND a numeric one (year) for the
    //       blue->red gradient. Load -> "Color by: host" (categories) or "Color by: year" (gradient).
    private static Phylogeny colorByPropertyTree() {
        final PhylogenyNode[] tips = {
                hostYear( "A/California/2009", "Human", 2009 ),
                hostYear( "A/swine/Iowa/2010", "Swine", 2010 ),
                hostYear( "A/duck/Guangdong/2013", "Avian", 2013 ),
                hostYear( "A/Shanghai/2013", "Human", 2013 ),
                hostYear( "A/equine/Mongolia/2011", "Equine", 2011 ),
                hostYear( "A/canine/Florida/2015", "Canine", 2015 ),
                hostYear( "A/chicken/Jiangsu/2017", "Avian", 2017 ),
                hostYear( "A/swine/Shandong/2018", "Swine", 2018 ),
                hostYear( "A/Cambodia/2021", "Human", 2021 ),
                hostYear( "A/duck/Vietnam/2024", "Avian", 2024 ) };
        final PhylogenyNode root = clade( 0, clade( 0.06, tips[ 0 ], tips[ 1 ], tips[ 4 ] ),
                                          clade( 0.05, clade( 0.03, tips[ 2 ], tips[ 3 ], tips[ 6 ] ),
                                                 clade( 0.04, tips[ 5 ], tips[ 7 ], tips[ 8 ], tips[ 9 ] ) ) );
        return tree( root, "Color by property (demo)",
                     "Synthetic influenza-surveillance tree. Each tip has a categorical 'host' and a numeric 'year'. "
                             + "Try Color by: host (a distinct color per host) or Color by: year (a numeric gradient)." );
    }

    private static PhylogenyNode hostYear( final String name, final String host, final int year ) {
        final PhylogenyNode n = leaf( name );
        cat( n, "data:host", host );
        num( n, "data:year", Integer.toString( year ) );
        return n;
    }

    // ----- "Annotation Columns": several properties of the four supported kinds so the tool can render a color strip
    //       (categorical), a heat map / bar (numeric) and a text column side by side. Load -> Tools > Annotation Columns.
    private static Phylogeny annotationColumnsTree() {
        final PhylogenyNode[] tips = {
                annotated( "isolate_01", "Human", "2.3.4.4b", "HA", 7.8 ),
                annotated( "isolate_02", "Avian", "2.3.4.4b", "HA", 6.1 ),
                annotated( "isolate_03", "Avian", "2.3.2.1c", "NA", 3.4 ),
                annotated( "isolate_04", "Swine", "1A.3.3.2", "HA", 5.9 ),
                annotated( "isolate_05", "Human", "2.3.4.4b", "PB2", 8.6 ),
                annotated( "isolate_06", "Swine", "1A.3.3.2", "NA", 2.2 ),
                annotated( "isolate_07", "Avian", "2.3.2.1c", "HA", 4.7 ),
                annotated( "isolate_08", "Human", "2.3.4.4b", "PB2", 9.1 ),
                annotated( "isolate_09", "Equine", "Fc1", "HA", 1.5 ),
                annotated( "isolate_10", "Avian", "2.3.4.4b", "NA", 6.8 ) };
        final PhylogenyNode root = clade( 0, clade( 0.05, tips[ 0 ], tips[ 1 ], tips[ 4 ], tips[ 7 ] ),
                                          clade( 0.05, clade( 0.03, tips[ 2 ], tips[ 6 ], tips[ 9 ] ),
                                                 clade( 0.03, tips[ 3 ], tips[ 5 ], tips[ 8 ] ) ) );
        return tree( root, "Annotation columns (demo)",
                     "Synthetic tree with four annotation kinds per tip: 'host' and 'segment' (categorical), "
                             + "'viral_load' (numeric) and 'clade' (text). Try Tools > Annotation Columns to render "
                             + "them as tip-aligned color-strip, heat-map/bar and text columns." );
    }

    private static PhylogenyNode annotated( final String name, final String host, final String clade,
                                            final String segment, final double viral_load ) {
        final PhylogenyNode n = leaf( name );
        cat( n, "data:host", host );
        cat( n, "data:segment", segment );
        cat( n, "data:clade", clade ); // free-text label -> text column
        num( n, "data:viral_load", Double.toString( viral_load ) );
        return n;
    }

    // ----- "Colorize Subtrees / Clade Bands by Taxonomic Rank": a mammal tree where each ORDER's clade root carries an
    //       in-tree rank annotation, so colorizing by 'order' works fully OFFLINE (no NCBI lookup needed for the demo).
    private static Phylogeny colorizeByRankTree() throws PhyloXmlDataFormatException {
        // Each TIP carries its OWN 'order' taxonomy (self-resolvable), so the tool colorizes OFFLINE with no
        // "resolve online?" prompt -- the species name is kept as the node name for identification.
        final PhylogenyNode carnivora = clade( 0.08,
                                               orderTip( "Felis catus", "Carnivora" ),
                                               orderTip( "Panthera leo", "Carnivora" ),
                                               orderTip( "Canis lupus", "Carnivora" ) );
        // DELIBERATELY WRONG internal annotation on the cats' clade root: because each cat tip carries its OWN
        // order, that self-identity is authoritative and the mis-annotation is IGNORED (the cats stay Carnivora) --
        // an illustration that a tip's own taxonomy wins over a wrong/partial internal-node annotation.
        taxon( carnivora, "Rodentia", "order" );
        final PhylogenyNode rodentia = clade( 0.09,
                                              orderTip( "Mus musculus", "Rodentia" ),
                                              orderTip( "Rattus norvegicus", "Rodentia" ) );
        final PhylogenyNode primates = clade( 0.07,
                                              orderTip( "Homo sapiens", "Primates" ),
                                              orderTip( "Pan troglodytes", "Primates" ),
                                              orderTip( "Macaca mulatta", "Primates" ) );
        final PhylogenyNode root = clade( 0, clade( 0.04, carnivora, rodentia ), primates );
        return tree( root, "Colorize by taxonomic rank (demo)",
                     "Synthetic mammal tree. Each TIP carries its order as an in-tree rank annotation (the species "
                             + "is the node name), so Tools > Colorize Subtrees via Taxonomic Rank -- and Annotate "
                             + "Clades by Rank -- work at 'order' (Carnivora, Rodentia, Primates) with NO online "
                             + "lookup. The Carnivora clade root is deliberately MIS-annotated 'Rodentia' to show "
                             + "that a tip's own identity wins over a wrong/partial internal annotation." );
    }

    /** A tip named for its species but carrying its ORDER as an in-tree rank annotation, so it self-resolves the
     *  colorization rank offline. */
    private static PhylogenyNode orderTip( final String species_name, final String order )
            throws PhyloXmlDataFormatException {
        final PhylogenyNode n = leaf( species_name );
        taxon( n, order, "order" );
        return n;
    }

    // ----- "Infer Ancestor Taxonomies": real, well-known species at the tips, NO taxonomy on the internal nodes.
    //       Run Analysis > Infer Ancestor Taxonomies and accept the online resolve -> each internal node is filled
    //       with the deepest taxon its descendants share (rank + NCBI tax-id). Needs an internet connection for the
    //       resolve (a tip's lineage is not stored in phyloXML, so this cannot be demonstrated offline).
    private static Phylogeny inferAncestorTaxonomiesTree() throws PhyloXmlDataFormatException {
        final PhylogenyNode primates = clade( 0.06, speciesTip( "Homo sapiens" ), speciesTip( "Pan troglodytes" ) );
        final PhylogenyNode rodents = clade( 0.05, speciesTip( "Mus musculus" ), speciesTip( "Rattus norvegicus" ) );
        final PhylogenyNode carnivores = clade( 0.05, speciesTip( "Felis catus" ), speciesTip( "Canis lupus" ) );
        final PhylogenyNode root = clade( 0, primates, clade( 0.04, rodents, carnivores ) );
        return tree( root, "Infer ancestral taxonomies (demo)",
                     "Synthetic vertebrate tree whose six tips are real, well-known species (Homo sapiens, Pan "
                             + "troglodytes, Mus musculus, Rattus norvegicus, Felis catus, Canis lupus) and whose "
                             + "internal nodes carry NO taxonomy. Run Analysis > Infer Ancestor Taxonomies and accept "
                             + "the online resolve: each internal node is filled with the deepest taxon its "
                             + "descendants share (e.g. Murinae, Carnivora, Boreoeutheria), with ranks and NCBI "
                             + "tax-ids where defined, so "
                             + "you can then Colorize / Annotate Clades by Rank on the inferred internal taxa. Requires "
                             + "an internet connection for the resolve." );
    }

    /** A tip carrying a real species scientific name (rank species), so the taxonomy service can resolve its lineage
     *  online for ancestral-taxonomy inference. */
    private static PhylogenyNode speciesTip( final String scientific_name ) throws PhyloXmlDataFormatException {
        final PhylogenyNode n = leaf( scientific_name );
        taxon( n, scientific_name, "species" );
        return n;
    }

    // ----- "Create Tanglegram": a linked PAIR of trees over the SAME eight taxa but with DIFFERENT topologies, so
    //       Analysis > Create Tanglegram draws crossing tip-to-tip connectors. Each tip's name IS its scientific
    //       name, so the pair links on either 'Node Name' or 'Taxonomy: Scientific Name'; each tip also carries a
    //       categorical 'clade' group, so the tanglegram's connectors can be coloured by clade. Open BOTH files.
    private static Phylogeny tanglegramTreeA() throws PhyloXmlDataFormatException {
        final PhylogenyNode root = clade( 0,
                                          clade( 0.05,
                                                 clade( 0.04, cladeTip( "Homo sapiens", "Primates" ),
                                                        cladeTip( "Pan troglodytes", "Primates" ) ),
                                                 clade( 0.04, cladeTip( "Mus musculus", "Rodentia" ),
                                                        cladeTip( "Rattus norvegicus", "Rodentia" ) ) ),
                                          clade( 0.05,
                                                 clade( 0.04, cladeTip( "Felis catus", "Carnivora" ),
                                                        cladeTip( "Canis lupus", "Carnivora" ) ),
                                                 clade( 0.04, cladeTip( "Gallus gallus", "Aves" ),
                                                        cladeTip( "Danio rerio", "Actinopterygii" ) ) ) );
        return tree( root, "Tanglegram tree A",
                     "One of a linked PAIR of demo trees for Analysis > Create Tanglegram. Its eight tips are the "
                             + "same eight species as 'tanglegram-tree-b.xml' but arranged in a different topology, "
                             + "so the tanglegram's tip-to-tip connectors cross. Open BOTH tanglegram-tree files, "
                             + "then run Analysis > Create Tanglegram and link on Node Name or Taxonomy: Scientific "
                             + "Name. Each tip also carries a categorical 'clade' group, so in the tanglegram window "
                             + "you can set Colour: clade to colour the connectors by clade." );
    }

    private static Phylogeny tanglegramTreeB() throws PhyloXmlDataFormatException {
        final PhylogenyNode root = clade( 0,
                                          clade( 0.05,
                                                 clade( 0.04, cladeTip( "Homo sapiens", "Primates" ),
                                                        cladeTip( "Mus musculus", "Rodentia" ) ),
                                                 clade( 0.04, cladeTip( "Pan troglodytes", "Primates" ),
                                                        cladeTip( "Felis catus", "Carnivora" ) ) ),
                                          clade( 0.05,
                                                 clade( 0.04, cladeTip( "Rattus norvegicus", "Rodentia" ),
                                                        cladeTip( "Canis lupus", "Carnivora" ) ),
                                                 clade( 0.04, cladeTip( "Danio rerio", "Actinopterygii" ),
                                                        cladeTip( "Gallus gallus", "Aves" ) ) ) );
        return tree( root, "Tanglegram tree B",
                     "The companion of 'tanglegram-tree-a.xml' for Analysis > Create Tanglegram: the same eight "
                             + "species in a different topology, so the linked connectors cross. Open both files "
                             + "and link on Node Name or Taxonomy: Scientific Name; set Colour: clade in the "
                             + "tanglegram window to colour the connectors by clade." );
    }

    /** A tanglegram tip: a species (name + taxonomy) plus a categorical 'clade' group, so a tanglegram's connectors
     *  can be coloured by clade. Its terminal branch length VARIES per species (deterministically, from the name) so
     *  the pair is non-ultrametric -- the tanglegram window's "Branch lengths" (aligned-phylogram) toggle then draws
     *  the branches to scale with dotted tip leaders, instead of a flat cladogram. */
    private static PhylogenyNode cladeTip( final String scientific_name, final String clade )
            throws PhyloXmlDataFormatException {
        final PhylogenyNode n = speciesTip( scientific_name );
        n.setDistanceToParent( 0.01 + ( ( Math.abs( scientific_name.hashCode() ) % 7 ) * 0.01 ) ); // 0.01..0.07
        cat( n, "data:clade", clade );
        return n;
    }

    // ----- "Create Tanglegram" with a link field PER tree: a GENE tree whose tips are gene accessions (the species
    //       is in each tip's Taxonomy) paired with a SPECIES tree whose tips ARE the species names. The same species
    //       live in DIFFERENT fields, so the two trees are linked on DIFFERENT fields: link the gene tree on
    //       'Taxonomy: Scientific Name' and the species tree on 'Node Name'. A single shared field links nothing.
    private static Phylogeny tanglegramGeneTree() throws PhyloXmlDataFormatException {
        final PhylogenyNode root = clade( 0,
                                          clade( 0.05,
                                                 clade( 0.04, geneTip( "ENSG001", "Homo sapiens" ),
                                                        geneTip( "ENSG002", "Pan troglodytes" ) ),
                                                 clade( 0.04, geneTip( "ENSG003", "Mus musculus" ),
                                                        geneTip( "ENSG004", "Rattus norvegicus" ) ) ),
                                          clade( 0.05,
                                                 clade( 0.04, geneTip( "ENSG005", "Felis catus" ),
                                                        geneTip( "ENSG006", "Canis lupus" ) ),
                                                 clade( 0.04, geneTip( "ENSG007", "Gallus gallus" ),
                                                        geneTip( "ENSG008", "Danio rerio" ) ) ) );
        return tree( root, "Tanglegram gene tree",
                     "A GENE tree for Analysis > Create Tanglegram, paired with 'tanglegram-species-tree.xml'. Its "
                             + "tip NAMES are gene accessions; each tip's species is in its Taxonomy (Scientific "
                             + "Name). The companion stores the same eight species as plain tip NAMES. Open BOTH "
                             + "files, run Analysis > Create Tanglegram, and link the trees on DIFFERENT fields "
                             + "holding the same value: 'Link first tree by: Taxonomy: Scientific Name' and 'Link "
                             + "second tree by: Node Name'. (A single shared field links nothing here.)" );
    }

    private static Phylogeny tanglegramSpeciesTree() {
        final PhylogenyNode root = clade( 0,
                                          clade( 0.05,
                                                 clade( 0.04, leaf( "Homo sapiens" ), leaf( "Mus musculus" ) ),
                                                 clade( 0.04, leaf( "Pan troglodytes" ), leaf( "Felis catus" ) ) ),
                                          clade( 0.05,
                                                 clade( 0.04, leaf( "Rattus norvegicus" ), leaf( "Canis lupus" ) ),
                                                 clade( 0.04, leaf( "Danio rerio" ), leaf( "Gallus gallus" ) ) ) );
        return tree( root, "Tanglegram species tree",
                     "A SPECIES tree for Analysis > Create Tanglegram, paired with 'tanglegram-gene-tree.xml'. Its "
                             + "tip NAMES are the species names. Link it on 'Node Name' while linking the gene tree "
                             + "on 'Taxonomy: Scientific Name' -- the same species stored in different fields." );
    }

    /** A gene-tree tip for the cross-field tanglegram demo: an accession NODE NAME plus the species in the taxonomy
     *  (Scientific Name), so the tip's node name and its species differ. */
    private static PhylogenyNode geneTip( final String accession, final String scientific_name )
            throws PhyloXmlDataFormatException {
        final PhylogenyNode n = leaf( accession );
        taxon( n, scientific_name, "species" );
        return n;
    }

    // ----- "Create Tanglegram" with an ASSOCIATION FILE: a cophylogenetic HOST tree (pocket gophers) and PARASITE
    //       tree (their chewing lice) whose tip names share NOTHING, linked through 'tanglegram-association.tsv'
    //       (host<TAB>louse). No field value join can pair them -- this is the parasite-vs-host case. Open BOTH trees,
    //       run Analysis > Create Tanglegram, link BOTH on Node Name, tick "Link by an association file" and choose
    //       tanglegram-association.tsv. The two topologies differ (a host switch), so a few connectors cross.
    private static Phylogeny tanglegramHostTree() {
        final PhylogenyNode root = clade( 0,
                                          clade( 0.05, leaf( "Thomomys_bottae" ), leaf( "Thomomys_talpoides" ) ),
                                          clade( 0.05,
                                                 clade( 0.04, leaf( "Geomys_bursarius" ),
                                                        leaf( "Cratogeomys_merriami" ) ),
                                                 clade( 0.04, leaf( "Orthogeomys_hispidus" ),
                                                        leaf( "Zygogeomys_trichopus" ) ) ) );
        return tree( root, "Tanglegram host tree (pocket gophers)",
                     "A HOST tree (pocket gophers) for Analysis > Create Tanglegram, paired with "
                             + "'tanglegram-parasite-tree.xml' (their chewing lice). The host and parasite tip names "
                             + "share nothing, so they are linked through the two-column 'tanglegram-association.tsv' "
                             + "(host<TAB>louse). Open BOTH trees, run Analysis > Create Tanglegram, link BOTH on "
                             + "Node Name, tick \"Link by an association file\" and choose tanglegram-association.tsv. "
                             + "The two topologies differ (a host switch), so some connectors cross." );
    }

    private static Phylogeny tanglegramParasiteTree() {
        final PhylogenyNode root = clade( 0,
                                          clade( 0.05, leaf( "Geomydoecus_actuosi" ), leaf( "Geomydoecus_thomomyus" ) ),
                                          clade( 0.05,
                                                 clade( 0.04, leaf( "Geomydoecus_panamensis" ),
                                                        leaf( "Geomydoecus_trichopi" ) ),
                                                 clade( 0.04, leaf( "Geomydoecus_ewingi" ),
                                                        leaf( "Geomydoecus_perotensis" ) ) ) );
        return tree( root, "Tanglegram parasite tree (chewing lice)",
                     "A PARASITE tree (chewing lice) for Analysis > Create Tanglegram, the companion of "
                             + "'tanglegram-host-tree.xml'. Its tips are louse names, unrelated to the host names, so "
                             + "the trees are linked through 'tanglegram-association.tsv' (host<TAB>louse), not by "
                             + "value. See the host-tree description for the steps." );
    }

    /** The two-column (host<TAB>louse) association table linking 'tanglegram-host-tree.xml' to
     *  'tanglegram-parasite-tree.xml'. No header row: one association per line, tab-separated. */
    private static String tanglegramAssociationTsv() {
        return "Thomomys_bottae\tGeomydoecus_actuosi\n" + "Thomomys_talpoides\tGeomydoecus_thomomyus\n"
                + "Geomys_bursarius\tGeomydoecus_ewingi\n" + "Cratogeomys_merriami\tGeomydoecus_perotensis\n"
                + "Orthogeomys_hispidus\tGeomydoecus_panamensis\n" + "Zygogeomys_trichopus\tGeomydoecus_trichopi\n";
    }

    // ----- "Zebra Stripes": a wider (16-tip) tree with a categorical 'host' + numeric 'reads' per tip, so the faint
    //       alternating row bands help track a label across to its Annotation Columns. Load -> Settings > Display >
    //       Zebra Stripes (optionally Tools > Annotation Columns).
    private static Phylogeny zebraStripesTree() {
        final String[] hosts = { "Human", "Avian", "Swine", "Bat" };
        final PhylogenyNode[] clades = new PhylogenyNode[ 4 ];
        int n = 1;
        for( int c = 0; c < 4; c++ ) {
            final PhylogenyNode[] leaves = new PhylogenyNode[ 4 ];
            for( int i = 0; i < 4; i++, n++ ) {
                final PhylogenyNode leaf = leaf( String.format( Locale.ROOT, "isolate_%02d", n ) );
                cat( leaf, "data:host", hosts[ c ] );
                num( leaf, "data:reads", Integer.toString( 100 * n ) );
                leaves[ i ] = leaf;
            }
            clades[ c ] = clade( 0.05, leaves );
        }
        final PhylogenyNode root = clade( 0, clade( 0.04, clades[ 0 ], clades[ 1 ] ),
                                          clade( 0.04, clades[ 2 ], clades[ 3 ] ) );
        return tree( root, "Zebra stripes (demo)",
                "Synthetic 16-tip tree with a categorical 'host' and numeric 'reads' per tip. Turn on Settings > "
                        + "Display > Zebra Stripes -- faint alternating row bands make it easy to track a tip label "
                        + "across to its Annotation Columns (Tools > Annotation Columns)." );
    }

    // ----- "Reverse Tip Order": an 8-tip ladder (caterpillar) tree with sequentially-numbered tips, so the tip order
    //       (and the staircase) visibly inverts when reversed. Load -> Settings > Display > Reverse Tip Order.
    private static Phylogeny reverseTipOrderTree() {
        // build the caterpillar bottom-up: tip_08/tip_07 at the deepest fork, then prepend tip_06..tip_01
        PhylogenyNode n = clade( 0.07, blLeaf( "tip_07", 0.06 ), blLeaf( "tip_08", 0.06 ) );
        for( int i = 6; i >= 1; i-- ) {
            n = clade( 0.07, blLeaf( String.format( Locale.ROOT, "tip_%02d", i ), 0.06 ), n );
        }
        return tree( n, "Reverse tip order (demo)",
                "Synthetic 8-tip ladder tree with sequentially-numbered tips (tip_01 first). Turn on Settings > "
                        + "Display > Reverse Tip Order to reverse the tip order -- the staircase inverts and tip_08 "
                        + "moves to the other end. Display-only; the tree data is unchanged." );
    }

    // ----- "Search emphasis" (Bold Found Labels / Dim Non-Matches): 15 tips of which four are "*_kinase", scattered
    //       across the tree, so searching "kinase" highlights a subset. Search box a for "kinase", then turn on
    //       Settings > Display > Bold Found Labels and/or Dim Non-Matches.
    private static Phylogeny searchEmphasisTree() {
        final PhylogenyNode a = conf( clade( 0.05, leaf( "AKT1_kinase" ), leaf( "actin" ), leaf( "tubulin" ),
                                             leaf( "MAPK1_kinase" ), leaf( "myosin" ) ), 88 );
        final PhylogenyNode b = conf( clade( 0.05, leaf( "src_kinase" ), leaf( "collagen" ), leaf( "keratin" ),
                                             leaf( "histone_H3" ), leaf( "ubiquitin" ) ), 92 );
        final PhylogenyNode c = conf( clade( 0.05, leaf( "CDK2_kinase" ), leaf( "GAPDH" ), leaf( "beta_globin" ),
                                             leaf( "insulin" ), leaf( "albumin" ) ), 76 );
        final PhylogenyNode root = clade( 0, conf( clade( 0.04, a, b ), 95 ), c );
        return tree( root, "Search emphasis (demo)",
                "Synthetic 15-tip tree; four tips are '*_kinase', scattered across the tree, with branch lengths and "
                        + "bootstrap support on the internal nodes. Search (box a) for \"kinase\", then turn on "
                        + "Settings > Display > Bold Found Labels and/or Dim Non-Matches -- the four hits go bold "
                        + "while the rest (names AND their support / branch-length numbers) fade toward the background." );
    }

    /** Adds a bootstrap confidence to a node's incoming branch (so the demo can show support values). */
    private static PhylogenyNode conf( final PhylogenyNode node, final double value ) {
        node.getBranchData().addConfidence( new Confidence( value, "bootstrap" ) );
        return node;
    }

    // ----- "Root at top / bottom" orientation: a small balanced phylogram with branch lengths and bootstrap support,
    //       so the vertical orientation shows a dendrogram with 45-degree tip labels and UPRIGHT support / branch-length
    //       numbers. Load -> Settings > Display > Orientation: root at top. (Later: pairs with columns below the tips.)
    private static Phylogeny rootOnTopTree() {
        final PhylogenyNode mammals = conf( clade( 0.10,
                conf( clade( 0.08, blLeaf( "Human", 0.06 ), blLeaf( "Chimp", 0.05 ), blLeaf( "Gorilla", 0.07 ) ), 96 ),
                conf( clade( 0.09, blLeaf( "Mouse", 0.11 ), blLeaf( "Rat", 0.10 ) ), 88 ) ), 99 );
        final PhylogenyNode birds = conf( clade( 0.12,
                conf( clade( 0.07, blLeaf( "Chicken", 0.09 ), blLeaf( "Turkey", 0.08 ) ), 90 ),
                conf( clade( 0.10, blLeaf( "Finch", 0.12 ), blLeaf( "Sparrow", 0.11 ) ), 84 ) ), 97 );
        final PhylogenyNode fish = conf( clade( 0.14, blLeaf( "Zebrafish", 0.16 ), blLeaf( "Salmon", 0.15 ),
                                                blLeaf( "Trout", 0.14 ) ), 91 );
        final PhylogenyNode root = clade( 0, conf( clade( 0.06, mammals, birds ), 100 ), fish );
        final Phylogeny phy = tree( root, "Root at top (demo)",
                "Synthetic 12-tip vertebrate phylogram with branch lengths (substitutions/site) and bootstrap support "
                        + "on the internal nodes. Set Settings > Display > Orientation to \"Root at Top\" (or \"Root at "
                        + "Bottom\"): the tree becomes a vertical dendrogram with vertical tip labels, while the support "
                        + "and branch-length numbers stay upright." );
        phy.setDistanceUnit( "substitutions/site" );
        return phy;
    }

    // ----- "Scale Axis": a phylogram whose branch lengths (substitutions/site) span a useful range, so the labeled
    //       distance axis shows nice ticks. Load -> Settings > Display > Scale Axis.
    private static Phylogeny scaleAxisTree() {
        final PhylogenyNode root = clade( 0,
                clade( 0.12, blLeaf( "Homo_sapiens", 0.20 ), blLeaf( "Pan_troglodytes", 0.35 ),
                       clade( 0.18, blLeaf( "Mus_musculus", 0.30 ), blLeaf( "Rattus_norvegicus", 0.15 ) ) ),
                clade( 0.10, blLeaf( "Gallus_gallus", 0.45 ), blLeaf( "Xenopus_laevis", 0.25 ),
                       clade( 0.22, blLeaf( "Danio_rerio", 0.40 ), blLeaf( "Drosophila_melanogaster", 0.20 ) ) ) );
        final Phylogeny phy = tree( root, "Scale axis (demo)",
                "Synthetic gene-family phylogram with branch lengths in substitutions/site (max depth ~0.7). Turn on "
                        + "Settings > Display > Scale Axis to read a labeled distance axis with ticks along the bottom." );
        phy.setDistanceUnit( "substitutions/site" );
        return phy;
    }

    // ----- "Node Age Bars (HPD)": a dated (ultrametric) mammal phylogram, branch lengths = time (My). Each internal
    //       node carries a phyloXML <date> (value = node age, min/max = its 95% interval), the native age model.
    //       Load -> view as phylogram -> Settings > Display > Node Age Bars (HPD). Root age 90 My; every root-to-tip = 90.
    private static Phylogeny hpdBarsTree() {
        final PhylogenyNode homo_pan = hpdClade( 2, 7, 6, 9, blLeaf( "Homo_sapiens", 7 ), blLeaf( "Pan_troglodytes", 7 ) );
        final PhylogenyNode homininae = hpdClade( 11, 9, 7, 12, blLeaf( "Gorilla_gorilla", 9 ), homo_pan );
        final PhylogenyNode hominidae = hpdClade( 54, 20, 16, 25, blLeaf( "Pongo_abelii", 20 ), homininae );
        final PhylogenyNode primates = hpdClade( 16, 74, 66, 80, blLeaf( "Macaca_mulatta", 74 ), hominidae );
        final PhylogenyNode root = hpdClade( 0, 90, 82, 98, blLeaf( "Mus_musculus", 90 ), primates );
        final Phylogeny phy = tree( root, "Node age HPD bars (demo)",
                "Synthetic dated (ultrametric) mammal phylogram; branch lengths are time in millions of years. Each "
                        + "internal node carries a phyloXML <date> (age + 95% interval). View as a phylogram and turn "
                        + "on Settings > Display > Node Age Bars (HPD)." );
        phy.setDistanceUnit( "My" );
        return phy;
    }

    /** An internal clade node with a branch length (My) and a phyloXML date: point age + its 95% interval (My). */
    private static PhylogenyNode hpdClade( final double branch_length, final int age, final int hpd_min,
                                           final int hpd_max, final PhylogenyNode... kids ) {
        final PhylogenyNode n = clade( branch_length, kids );
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", java.math.BigDecimal.valueOf( age ),
                java.math.BigDecimal.valueOf( hpd_min ), java.math.BigDecimal.valueOf( hpd_max ), "mya" ) );
        return n;
    }

    // ----- "Domain Architectures": a small protein-family phylogram; each tip carries a sequence with a Pfam-style
    //       domain architecture (colored boxes along the protein). Load -> view as phylogram -> Display Data: Domain
    //       Architectures. Also works in a vertical orientation (root at top/bottom): each tip's architecture becomes a
    //       vertical bar hanging off the tip.
    private static Phylogeny domainArchitecturesTree() {
        final PhylogenyNode root = clade( 0,
                clade( 0.15,
                       domainLeaf( "kinase_a", 0.22, 500, dom( "SH3", 10, 60 ), dom( "SH2", 75, 165 ),
                                   dom( "Pkinase", 185, 445 ) ),
                       domainLeaf( "kinase_b", 0.30, 450, dom( "SH2", 20, 110 ), dom( "Pkinase", 130, 400 ) ),
                       domainLeaf( "kinase_c", 0.26, 470, dom( "PH", 15, 120 ), dom( "Pkinase", 145, 410 ) ) ),
                clade( 0.12,
                       domainLeaf( "adaptor_d", 0.18, 300, dom( "SH3", 5, 55 ), dom( "SH3", 70, 120 ),
                                   dom( "SH2", 140, 235 ) ),
                       domainLeaf( "kinase_e", 0.34, 430, dom( "Pkinase", 30, 300 ), dom( "SAM", 330, 400 ) ),
                       domainLeaf( "receptor_f", 0.28, 650, dom( "Ig", 20, 110 ), dom( "Ig", 130, 220 ),
                                   dom( "Pkinase", 360, 630 ) ) ) );
        return tree( root, "Domain architectures (demo)",
                     "Synthetic protein-family phylogram; each tip carries a sequence with a Pfam-style domain "
                             + "architecture. View as a phylogram and turn on Display Data: Domain Architectures. Also "
                             + "works in a vertical orientation (root at top/bottom), where each tip's architecture "
                             + "becomes a vertical bar." );
    }

    /** A leaf carrying a sequence with a domain architecture of {@code total_length} amino acids. */
    private static PhylogenyNode domainLeaf( final String name, final double bl, final int total_length,
                                             final ProteinDomain... domains ) {
        final PhylogenyNode n = blLeaf( name, bl );
        final List<PhylogenyData> ds = new ArrayList<>();
        for( final ProteinDomain d : domains ) {
            ds.add( d );
        }
        final Sequence seq = new Sequence();
        seq.setName( name );
        seq.setDomainArchitecture( new DomainArchitecture( ds, total_length ) );
        n.getNodeData().addSequence( seq );
        return n;
    }

    /** A protein domain with a strong (well below the default 1e-3 threshold) e-value so it is drawn. */
    private static ProteinDomain dom( final String name, final int from, final int to ) {
        return new ProteinDomain( name, from, to, 1e-6 );
    }

    // ----- "Heat map matrix": each tip carries a row of numeric values across six samples (s1..s6). Add s1..s6 as
    //       "Heat map (matrix)" annotation columns -> ONE shared color scale across the whole grid (a clustergram).
    private static Phylogeny heatmapMatrixTree() {
        final PhylogenyNode root = clade( 0.0,
                clade( 0.10, matrixLeaf( "gene_A", 0.20, 8, 2, 1, 0, 3, 1 ),
                       matrixLeaf( "gene_B", 0.25, 9, 3, 2, 1, 4, 0 ) ),
                clade( 0.12, matrixLeaf( "gene_C", 0.18, 1, 7, 8, 2, 0, 1 ),
                       clade( 0.10, matrixLeaf( "gene_D", 0.15, 0, 6, 9, 3, 1, 2 ),
                              matrixLeaf( "gene_E", 0.20, 2, 5, 7, 4, 0, 1 ) ) ),
                clade( 0.15, matrixLeaf( "gene_F", 0.22, 3, 1, 0, 8, 9, 6 ),
                       matrixLeaf( "gene_G", 0.19, 2, 0, 1, 7, 8, 9 ),
                       matrixLeaf( "gene_H", 0.21, 4, 2, 0, 6, 7, 8 ) ) );
        return tree( root, "Heat map matrix (demo)",
                     "Synthetic tree where each tip carries an abundance value across six samples (s1..s6). View as a "
                             + "phylogram, then Tools > Annotation Columns and add s1..s6 as \"Heat map (matrix)\": "
                             + "they render on ONE shared color scale (a clustergram grid), best in a vertical "
                             + "orientation with Tip Labels Below Columns." );
    }

    /** A leaf carrying numeric sample values s1..sN (as data:s1 .. data:sN properties). */
    private static PhylogenyNode matrixLeaf( final String name, final double bl, final int... samples ) {
        final PhylogenyNode n = blLeaf( name, bl );
        for( int i = 0; i < samples.length; ++i ) {
            num( n, "data:s" + ( i + 1 ), Integer.toString( samples[ i ] ) );
        }
        return n;
    }

    private static PhylogenyNode blLeaf( final String name, final double branch_length ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( branch_length );
        return n;
    }

    // ---- small builders --------------------------------------------------------------------------------------------

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( 0.02 );
        return n;
    }

    private static PhylogenyNode clade( final double branch_length, final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setDistanceToParent( branch_length );
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static void num( final PhylogenyNode n, final String ref, final String value ) {
        addProperty( n, ref, value, NUM );
    }

    private static void cat( final PhylogenyNode n, final String ref, final String value ) {
        addProperty( n, ref, value, TEXT );
    }

    private static void addProperty( final PhylogenyNode n, final String ref, final String value,
                                     final String datatype ) {
        PropertiesList pl = n.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            n.getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( ref, value, "", datatype, AppliesTo.NODE ) );
    }

    private static void taxon( final PhylogenyNode n, final String scientific_name, final String rank )
            throws PhyloXmlDataFormatException {
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        t.setRank( rank );
        n.getNodeData().setTaxonomy( t );
    }

    private static Phylogeny tree( final PhylogenyNode root, final String name, final String description ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setName( name );
        phy.setDescription( description );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private DemoTreeGenerator() {
    }
}

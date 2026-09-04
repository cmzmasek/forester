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
        write( dir, "label-properties.xml", labelPropertiesTree() );
        write( dir, "symbol-columns.xml", symbolColumnsTree() );
        write( dir, "stacked-bar-columns.xml", stackedBarColumnsTree() );
        write( dir, "colorize-by-rank.xml", colorizeByRankTree() );
        write( dir, "infer-ancestor-taxonomies.xml", inferAncestorTaxonomiesTree() );
        write( dir, "gene-duplication.xml", reconciliationGeneTree() );
        write( dir, "bat-phylogeny.xml", batSpeciesTree() );
        write( dir, "animal-tree-of-life.xml", animalTreeOfLife() );
        write( dir, "lagomorph-time-tree.xml", lagomorphTimeTree() );
        write( dir, "scale-axis.xml", scaleAxisTree() );
        write( dir, "long-branch-break.xml", longBranchBreakTree() );
        write( dir, "date-in-labels.xml", dateInLabelsTree() );
        write( dir, "node-hpd-bars.xml", hpdBarsTree() );
        write( dir, "dinosaur-time-tree.xml", dinosaurTimeTree() );
        write( dir, "fossil-range-bars.xml", fossilRangeTree() );
        write( dir, "ammonite-time-tree.xml", ammoniteTimeTree() );
        write( dir, "tree-of-life-deep-time.xml", deepTimeTree() );
        write( dir, "sars-cov-2-time-tree.xml", sarsCov2TimeTree() );
        write( dir, "zebra-stripes.xml", zebraStripesTree() );
        write( dir, "reverse-tip-order.xml", reverseTipOrderTree() );
        write( dir, "root-on-top.xml", rootOnTopTree() );
        write( dir, "domain-architectures.xml", domainArchitecturesTree() );
        write( dir, "heatmap-matrix.xml", heatmapMatrixTree() );
        write( dir, "import-annotations.xml", importAnnotationsTree() );
        write( dir, "alignment.xml", alignmentTree() );
        write( dir, "gtdb-genomes.xml", gtdbGenomeTree() );
        writeText( dir, "gtdb-classifications.tsv", gtdbClassificationsTsv() );
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

    /** The companion table for {@link #importAnnotationsTree()}: a categorical 'host', a quoted 'country' (with an
     *  embedded comma, for CSV realism), a numeric 'reads', and a binary 'resistant' (yes / no, with two blanks =
     *  untested -- so it renders as a filled / hollow / nothing Symbol column after import). 'resistant' is the LAST
     *  column so 'reads' stays column index 3 (the import tool test renames that column by index). */
    private static String importAnnotationsCsv() {
        return "name,host,country,reads,resistant\n"
                + "isolate_01,mosquito,USA,1200,yes\n"
                + "isolate_02,mosquito,USA,980,no\n"
                + "isolate_03,bat,China,1543,yes\n"
                + "isolate_04,bat,China,760,\n"
                + "isolate_05,pig,Vietnam,410,no\n"
                + "isolate_06,pig,Vietnam,1330,yes\n"
                + "isolate_07,bird,\"Congo, DR\",275,no\n"
                + "isolate_08,bird,Kenya,1890,yes\n"
                + "isolate_09,mosquito,Brazil,640,\n"
                + "isolate_10,mosquito,Brazil,1120,yes\n"
                + "isolate_11,pig,Vietnam,505,no\n"
                + "isolate_12,bat,China,1450,yes\n";
    }

    // ----- "Import GTDB Taxonomy": a tree of Bacteria + Archaea genomes named ONLY by their assembly accession
    //       (no taxonomy), paired with a GTDB-Tk-style classification table (user_genome -> d__..;p__..;s__..).
    //       File > Import GTDB Taxonomy, pick gtdb-classifications.tsv -> each rank becomes a gtdb:<rank> property
    //       (+ a taxonomy at the most specific rank) -> Color by gtdb:phylum / gtdb:domain, or add an Annotation
    //       Column for gtdb:family. Entirely offline.
    private static Phylogeny gtdbGenomeTree() {
        // Bacteria
        final PhylogenyNode pseudomonadota = clade( 0.06, clade( 0.03, leaf( "GCF_000005845.2" ), // E. coli
                                                                  leaf( "GCF_000006945.2" ) ),      // Salmonella enterica
                                                    leaf( "GCF_000006765.1" ) );                    // Pseudomonas aeruginosa
        final PhylogenyNode bacillota = clade( 0.06, clade( 0.03, leaf( "GCF_000009045.1" ),        // Bacillus subtilis
                                                            leaf( "GCF_000013425.1" ) ),            // Staphylococcus aureus
                                               leaf( "GCF_000009065.1" ) );                          // Clostridioides difficile
        final PhylogenyNode actino_bactero_cyano = clade( 0.05,
                                                          clade( 0.04, leaf( "GCF_000195955.2" ),    // Mycobacterium tuberculosis
                                                                 leaf( "GCF_000203835.1" ) ),       // Streptomyces coelicolor
                                                          clade( 0.05, leaf( "GCF_000025985.1" ),    // Bacteroides fragilis
                                                                 leaf( "GCF_000011465.1" ) ) );      // Prochlorococcus marinus
        final PhylogenyNode bacteria = clade( 0.08, pseudomonadota,
                                              clade( 0.05, bacillota, actino_bactero_cyano ) );
        // Archaea
        final PhylogenyNode archaea = clade( 0.09, leaf( "GCF_000091665.1" ),                        // Methanocaldococcus jannaschii
                                             clade( 0.05, leaf( "GCF_000006805.1" ),                 // Halobacterium salinarum
                                                    clade( 0.05, leaf( "GCF_000018465.1" ),          // Nitrosopumilus maritimus
                                                           leaf( "GCF_000007765.1" ) ) ) );          // Saccharolobus solfataricus
        final PhylogenyNode root = clade( 0.0, bacteria, archaea );
        return tree( root, "Import GTDB Taxonomy (demo)",
                     "A tree of 14 Bacteria + Archaea genomes named only by their assembly accession -- no taxonomy. "
                             + "Pair it with the companion gtdb-classifications.tsv (a GTDB-Tk-style table): "
                             + "File > Import GTDB Taxonomy, pick the TSV. Each GTDB rank becomes a gtdb:<rank> "
                             + "property (+ a taxonomy at the most specific rank), so you can Color by gtdb:phylum "
                             + "or gtdb:domain, add an Annotation Column for gtdb:family, and search gtdb:genus -- "
                             + "entirely offline." );
    }

    /** The companion GTDB-Tk-style classification table for {@link #gtdbGenomeTree()}: two columns
     *  (user_genome &lt;TAB&gt; classification), one genome per row, the classification a 7-rank d__..;p__..;s__.. string. */
    private static String gtdbClassificationsTsv() {
        return "user_genome\tclassification\n"
                + "GCF_000005845.2\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n"
                + "GCF_000006945.2\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica\n"
                + "GCF_000006765.1\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa\n"
                + "GCF_000009045.1\td__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis\n"
                + "GCF_000013425.1\td__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus\n"
                + "GCF_000009065.1\td__Bacteria;p__Bacillota_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__Clostridioides;s__Clostridioides difficile\n"
                + "GCF_000195955.2\td__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacterium tuberculosis\n"
                + "GCF_000203835.1\td__Bacteria;p__Actinomycetota;c__Actinomycetia;o__Streptomycetales;f__Streptomycetaceae;g__Streptomyces;s__Streptomyces coelicolor\n"
                + "GCF_000025985.1\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides fragilis\n"
                + "GCF_000011465.1\td__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__PCC-6307;f__Cyanobiaceae;g__Prochlorococcus;s__Prochlorococcus marinus\n"
                + "GCF_000091665.1\td__Archaea;p__Methanobacteriota;c__Methanococci;o__Methanococcales;f__Methanocaldococcaceae;g__Methanocaldococcus;s__Methanocaldococcus jannaschii\n"
                + "GCF_000006805.1\td__Archaea;p__Halobacteriota;c__Halobacteria;o__Halobacteriales;f__Halobacteriaceae;g__Halobacterium;s__Halobacterium salinarum\n"
                + "GCF_000018465.1\td__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrosopumilales;f__Nitrosopumilaceae;g__Nitrosopumilus;s__Nitrosopumilus maritimus\n"
                + "GCF_000007765.1\td__Archaea;p__Thermoproteota;c__Thermoprotei;o__Sulfolobales;f__Sulfolobaceae;g__Saccharolobus;s__Saccharolobus solfataricus\n";
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
        final PhylogenyNode e = styledLeaf( "Rattus norvegicus", Font.PLAIN, 12, new Color( 0xC8, 0x7A, 0x00 ),
                                            NodeShape.DIAMOND, NodeFill.SOLID, 12f, new Color( 0xC8, 0x7A, 0x00 ) );
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
                hostYear( "A/duck/Vietnam/2024", "Avian", 2024 ),
                // environmental samples with NO host annotation: the legend's dashed "no value" row counts them
                hostYear( "A/environment/Hunan/2016", null, 2016 ),
                hostYear( "A/environment/Fujian/2019", null, 2019 ) };
        final PhylogenyNode root = clade( 0, clade( 0.06, tips[ 0 ], tips[ 1 ], tips[ 4 ], tips[ 10 ] ),
                                          clade( 0.05, clade( 0.03, tips[ 2 ], tips[ 3 ], tips[ 6 ], tips[ 11 ] ),
                                                 clade( 0.04, tips[ 5 ], tips[ 7 ], tips[ 8 ], tips[ 9 ] ) ) );
        return tree( root, "Color by property (demo)",
                     "Synthetic influenza-surveillance tree. Each tip has a categorical 'host' and a numeric 'year'; "
                             + "two environmental samples carry no host. Try Color by: host (a distinct color per "
                             + "host; the legend's dashed 'no value' row counts the host-less tips) or Color by: "
                             + "year (a numeric gradient)." );
    }

    private static PhylogenyNode hostYear( final String name, final String host, final int year ) {
        final PhylogenyNode n = leaf( name );
        if ( host != null ) { // null = no host annotation at all (feeds the legend's "no value" row)
            cat( n, "data:host", host );
        }
        num( n, "data:year", Integer.toString( year ) );
        return n;
    }

    // ----- "Annotation Columns": several properties of the four supported kinds so the tool can render a color strip
    //       (categorical), a heat map / bar (numeric) and a text column side by side. Load -> Tools > Annotation Fields.
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
                             + "'viral_load' (numeric) and 'clade' (text). Try Tools > Annotation Fields to render "
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

    // ----- "Properties in labels": tips carrying SIX properties each -- the case that makes a tip label unreadable
    //       when every property is appended to it. The chooser gives each field ONE role: two stay in the label
    //       (the accession and the passage history, which belong beside the name) and four become tip-aligned
    //       columns. Load -> Tools > Annotation Fields.
    private static Phylogeny labelPropertiesTree() {
        final PhylogenyNode[] tips = {
                isolate( "A/Vietnam/1203/2004", "EPI1731", "Human", "Viet Nam", "2.3.4.4b", "E3", 2004 ),
                isolate( "A/HongKong/156/1997", "EPI0885", "Human", "China", "2.3.2.1c", "E2", 1997 ),
                isolate( "A/duck/Hunan/795/2002", "EPI1442", "Avian", "China", "2.3.2.1c", "E1", 2002 ),
                isolate( "A/chicken/Egypt/07/2010", "EPI2984", "Avian", "Egypt", "2.2.1", "E4", 2010 ),
                isolate( "A/swine/Iowa/15/1930", "EPI0112", "Swine", "USA", "1A.3.3.2", "MDCK1", 1930 ),
                isolate( "A/turkey/Ontario/6118/1968", "EPI0671", "Avian", "Canada", "2.2.1", "E2", 1968 ),
                isolate( "A/equine/Prague/1/1956", "EPI0334", "Equine", "Czechia", "Fc1", "E5", 1956 ),
                isolate( "A/mallard/Alberta/24/2001", "EPI1390", "Avian", "Canada", "2.3.4.4b", "E1", 2001 ) };
        final PhylogenyNode root = clade( 0, clade( 0.05, tips[ 0 ], tips[ 1 ], tips[ 2 ], tips[ 3 ] ),
                                          clade( 0.05, clade( 0.03, tips[ 4 ], tips[ 6 ] ),
                                                 clade( 0.03, tips[ 5 ], tips[ 7 ] ) ) );
        return tree( root, "Properties in labels (demo)",
                     "Synthetic tree whose tips each carry SIX properties: 'accession' (a different value on every "
                             + "tip), 'host', 'country', 'clade', 'passage' and 'year'. Showing all six in the tip "
                             + "label makes it unreadable, so give each field ONE role in Tools > Annotation Fields: "
                             + "keep 'accession' and 'passage' in the label, render 'host', 'country' and 'clade' as "
                             + "color-strip columns and 'year' as a heat map -- and use the up/down buttons to set "
                             + "the order the label reads in. Note that 'accession' is offered for the LABEL even "
                             + "though it has a different value on every tip, which makes it useless as a color." );
    }

    private static PhylogenyNode isolate( final String name, final String accession, final String host,
                                          final String country, final String clade, final String passage,
                                          final int year ) {
        final PhylogenyNode n = leaf( name );
        cat( n, "data:accession", accession );
        cat( n, "data:host", host );
        cat( n, "data:country", country );
        cat( n, "data:clade", clade );
        cat( n, "data:passage", passage );
        num( n, "data:year", Integer.toString( year ) );
        return n;
    }

    // ----- "Symbol columns": tips carry a BINARY flag (yes / no / untested) plus a categorical host, so a SYMBOL
    //       annotation column renders filled / hollow / nothing marks -- a presence/absence (binary) mark column.
    private static Phylogeny symbolColumnsTree() {
        final PhylogenyNode[] tips = {
                symbolTip( "isolate_01", "yes", "Human" ),
                symbolTip( "isolate_02", "no", "Avian" ),
                symbolTip( "isolate_03", "yes", "Avian" ),
                symbolTip( "isolate_04", null, "Swine" ),   // untested -> nothing is drawn for this tip
                symbolTip( "isolate_05", "yes", "Human" ),
                symbolTip( "isolate_06", "no", "Swine" ),
                symbolTip( "isolate_07", "no", "Avian" ),
                symbolTip( "isolate_08", "yes", "Equine" ) };
        final PhylogenyNode root = clade( 0, clade( 0.05, tips[ 0 ], tips[ 1 ], tips[ 2 ], tips[ 3 ] ),
                                          clade( 0.05, tips[ 4 ], tips[ 5 ], tips[ 6 ], tips[ 7 ] ) );
        return tree( root, "Symbol columns (demo)",
                     "Synthetic tree whose tips carry a binary 'resistant' flag (yes / no / untested) and a "
                             + "categorical 'host'. Try Tools > Annotation Fields and pick the 'Symbol' render "
                             + "type: a present value draws a filled mark, an explicit 'no' a hollow mark, and an "
                             + "untested tip nothing -- a presence/absence (binary) mark column. A many-valued "
                             + "categorical field ('host') becomes distinct colored marks." );
    }

    private static PhylogenyNode symbolTip( final String name, final String resistant, final String host ) {
        final PhylogenyNode n = leaf( name );
        if ( resistant != null ) {
            cat( n, "data:resistant", resistant );
        }
        cat( n, "data:host", host );
        return n;
    }

    // ----- "Stacked bar columns": metagenomic samples whose tips carry read counts for three bacterial phyla. Setting
    //       all three to "Stacked bar" MERGES them into one segmented bar per tip (segment length = read count, so the
    //       bar's total length shows sequencing depth AND its segments show composition); a normalized bar fills the
    //       width and shows pure proportion. Load -> Tools > Annotation Fields. (Compositional / proportional data.)
    private static Phylogeny stackedBarColumnsTree() {
        final PhylogenyNode[] tips = {
                sample( "gut_01", 820, 140, 40 ),
                sample( "gut_02", 610, 300, 90 ),
                sample( "gut_03", 450, 380, 170 ),
                sample( "oral_01", 300, 250, 1450 ), // a much larger total -> its absolute bar fills the column width
                sample( "oral_02", 260, 210, 980 ),
                sample( "skin_01", 900, 60, 40 ),
                sample( "skin_02", 780, 120, 100 ),
                sample( "soil_01", 200, 700, 300 ),
                sample( "soil_02", 150, 640, 210 ),
                sample( "water_01", 120, 400, 480 ) };
        final PhylogenyNode root = clade( 0, clade( 0.05, tips[ 0 ], tips[ 1 ], tips[ 2 ] ),
                                          clade( 0.05, clade( 0.03, tips[ 3 ], tips[ 4 ] ),
                                                 clade( 0.03, tips[ 5 ], tips[ 6 ] ) ),
                                          clade( 0.05, tips[ 7 ], tips[ 8 ], tips[ 9 ] ) );
        return tree( root, "Stacked bar columns (demo)",
                     "Synthetic metagenomic-sample tree: each tip carries read counts for three bacterial phyla "
                             + "('firmicutes', 'bacteroidetes', 'proteobacteria'). Try Tools > Annotation Fields, set "
                             + "all three to 'Stacked bar' -- they merge into one segmented bar per tip (segment length "
                             + "= the read count, so the bar's total length shows sequencing depth AND its segments show "
                             + "composition). Tick 'Normalize stacked bars to 100%' for pure proportional composition." );
    }

    private static PhylogenyNode sample( final String name, final int firmicutes, final int bacteroidetes,
                                         final int proteobacteria ) {
        final PhylogenyNode n = leaf( name );
        num( n, "data:firmicutes", Integer.toString( firmicutes ) );
        num( n, "data:bacteroidetes", Integer.toString( bacteroidetes ) );
        num( n, "data:proteobacteria", Integer.toString( proteobacteria ) );
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

    // ----- "Gene duplications & speciations (using NCBI taxonomy)": a small gene-family tree that DOES NOT come with a
    //       species tree. Every tip is a gene (its node NAME is the gene label) tagged with the real species it comes
    //       from (its <taxonomy> scientific name). The family was duplicated in the common ancestor of the three
    //       mammals into two paralogs (BCL2-like clade + BCLX-like clade), each then following the species tree
    //       (human+chimp)+mouse. Run Analysis > Infer Duplications & Speciations (using NCBI taxonomy): a species tree
    //       is built AUTOMATICALLY from the tips' NCBI taxonomy and GSDIR marks the ancestral gene duplication (red) and
    //       the speciations (green). Needs an internet connection for the resolve (a tip's lineage is not stored in
    //       phyloXML, so this cannot be demonstrated offline).
    private static Phylogeny reconciliationGeneTree() throws PhyloXmlDataFormatException {
        final PhylogenyNode bcl2 = clade( 0.04,
                clade( 0.06, geneTip( "BCL2_human", "Homo sapiens" ), geneTip( "BCL2_chimp", "Pan troglodytes" ) ),
                geneTip( "BCL2_mouse", "Mus musculus" ) );
        final PhylogenyNode bclx = clade( 0.04,
                clade( 0.06, geneTip( "BCLX_human", "Homo sapiens" ), geneTip( "BCLX_chimp", "Pan troglodytes" ) ),
                geneTip( "BCLX_mouse", "Mus musculus" ) );
        final PhylogenyNode root = clade( 0, bcl2, bclx );
        return tree( root, "Gene duplications & speciations (demo)",
                     "Synthetic gene-family tree with NO accompanying species tree. Each of the six tips is a gene "
                             + "(node name = gene label, e.g. BCL2_human) tagged with the real species it comes from "
                             + "(taxonomy scientific name: Homo sapiens, Pan troglodytes, Mus musculus). The family was "
                             + "duplicated in the mammalian common ancestor into two paralogs (BCL2-like + BCLX-like), "
                             + "each following the species tree (human+chimp)+mouse. Run Analysis > Infer Duplications & "
                             + "Speciations (using NCBI taxonomy): a species tree is built AUTOMATICALLY from the tips' "
                             + "NCBI taxonomy and GSDIR marks the ancestral gene DUPLICATION and the SPECIATIONS. "
                             + "Requires an internet connection for the resolve." );
    }

    // ----- "Bat phylogeny": a larger (34-tip) species tree of the bats (order Chiroptera), spanning both suborders and
    //       eight families. Every tip carries a full <taxonomy> -- common name, scientific (Latin) name, rank, and a
    //       taxonomic SYNONYM (mostly the original Vespertilio/Pteropus/Phyllostoma combination the species was first
    //       described in). Every INTERNAL node is annotated with its clade name + rank (order / suborder / superfamily /
    //       family / subfamily / genus), so Colorize / Annotate Clades by Rank works OFFLINE at any of those ranks and the
    //       Internal Taxonomy Key summarises the backbone. Topology follows the accepted two-suborder classification
    //       (Yinpterochiroptera + Yangochiroptera); it is schematic (illustrative branch lengths), not a rigorous estimate.
    private static Phylogeny batSpeciesTree() throws PhyloXmlDataFormatException {
        // --- Yinpterochiroptera ---
        final PhylogenyNode pteropodidae = namedClade( 0.05, "Pteropodidae", "family",
                batTip( "Large flying fox", "Pteropus vampyrus", "Vespertilio vampyrus" ),
                batTip( "Egyptian fruit bat", "Rousettus aegyptiacus", "Pteropus aegyptiacus" ),
                batTip( "Straw-coloured fruit bat", "Eidolon helvum", "Vespertilio vampyrus helvus" ) );
        final PhylogenyNode rhinolophidae = namedClade( 0.05, "Rhinolophidae", "family",
                batTip( "Greater horseshoe bat", "Rhinolophus ferrumequinum", "Vespertilio ferrumequinum" ),
                batTip( "Lesser horseshoe bat", "Rhinolophus hipposideros", "Vespertilio hipposideros" ) );
        final PhylogenyNode megadermatidae = namedClade( 0.06, "Megadermatidae", "family",
                batTip( "Greater false vampire bat", "Lyroderma lyra", "Megaderma lyra" ) );
        final PhylogenyNode rhinolophoidea = namedClade( 0.04, "Rhinolophoidea", "superfamily",
                rhinolophidae, megadermatidae );
        final PhylogenyNode yinptero = namedClade( 0.04, "Yinpterochiroptera", "suborder",
                pteropodidae, rhinolophoidea );

        // --- Yangochiroptera : Noctilionoidea ---
        final PhylogenyNode noctilionidae = namedClade( 0.06, "Noctilionidae", "family",
                batTip( "Greater bulldog bat", "Noctilio leporinus", "Vespertilio leporinus" ) );
        final PhylogenyNode phyllostomidae = namedClade( 0.05, "Phyllostomidae", "family",
                batTip( "Common vampire bat", "Desmodus rotundus", "Phyllostoma rotundum" ),
                batTip( "Seba's short-tailed bat", "Carollia perspicillata", "Vespertilio perspicillatus" ),
                batTip( "Pallas's long-tongued bat", "Glossophaga soricina", "Vespertilio soricinus" ) );
        final PhylogenyNode noctilionoidea = namedClade( 0.04, "Noctilionoidea", "superfamily",
                noctilionidae, phyllostomidae );

        // --- Yangochiroptera : Vespertilionoidea ---
        final PhylogenyNode molossidae = namedClade( 0.05, "Molossidae", "family",
                batTip( "Velvety free-tailed bat", "Molossus molossus", "Vespertilio molossus" ),
                batTip( "European free-tailed bat", "Tadarida teniotis", "Cephalotes teniotis" ) );
        final PhylogenyNode miniopteridae = namedClade( 0.06, "Miniopteridae", "family",
                batTip( "Common bent-wing bat", "Miniopterus schreibersii", "Vespertilio schreibersii" ) );
        // Vespertilionidae (vesper bats -- the largest family), split into subfamilies with genus clades
        final PhylogenyNode myotis = namedClade( 0.04, "Myotis", "genus",
                batTip( "Greater mouse-eared bat", "Myotis myotis", "Vespertilio myotis" ),
                batTip( "Little brown bat", "Myotis lucifugus", "Vespertilio lucifugus" ),
                batTip( "Daubenton's bat", "Myotis daubentonii", "Vespertilio daubentonii" ),
                batTip( "Natterer's bat", "Myotis nattereri", "Vespertilio nattereri" ) );
        final PhylogenyNode myotinae = namedClade( 0.04, "Myotinae", "subfamily", myotis );
        final PhylogenyNode pipistrellus = namedClade( 0.03, "Pipistrellus", "genus",
                batTip( "Common pipistrelle", "Pipistrellus pipistrellus", "Vespertilio pipistrellus" ),
                batTip( "Nathusius's pipistrelle", "Pipistrellus nathusii", "Vespertilio nathusii" ),
                batTip( "Kuhl's pipistrelle", "Pipistrellus kuhlii", "Vespertilio kuhlii" ) );
        final PhylogenyNode nyctalus = namedClade( 0.03, "Nyctalus", "genus",
                batTip( "Common noctule", "Nyctalus noctula", "Vespertilio noctula" ),
                batTip( "Leisler's bat", "Nyctalus leisleri", "Vespertilio leisleri" ) );
        final PhylogenyNode eptesicus = namedClade( 0.03, "Eptesicus", "genus",
                batTip( "Serotine bat", "Eptesicus serotinus", "Vespertilio serotinus" ),
                batTip( "Big brown bat", "Eptesicus fuscus", "Vespertilio fuscus" ) );
        final PhylogenyNode lasiurus = namedClade( 0.03, "Lasiurus", "genus",
                batTip( "Eastern red bat", "Lasiurus borealis", "Vespertilio borealis" ),
                batTip( "Hoary bat", "Lasiurus cinereus", "Vespertilio cinereus" ) );
        final PhylogenyNode vespertilioninae = namedClade( 0.04, "Vespertilioninae", "subfamily",
                pipistrellus, nyctalus, eptesicus, lasiurus,
                batTip( "Parti-coloured bat", "Vespertilio murinus", "Vespertilio discolor" ),
                batTip( "Savi's pipistrelle", "Hypsugo savii", "Vespertilio savii" ),
                batTip( "Tricolored bat", "Perimyotis subflavus", "Vespertilio subflavus" ),
                batTip( "Evening bat", "Nycticeius humeralis", "Vespertilio humeralis" ),
                batTip( "Pallid bat", "Antrozous pallidus", "Vespertilio pallidus" ) );
        final PhylogenyNode plecotinae = namedClade( 0.04, "Plecotinae", "subfamily",
                batTip( "Brown long-eared bat", "Plecotus auritus", "Vespertilio auritus" ),
                batTip( "Western barbastelle", "Barbastella barbastellus", "Vespertilio barbastellus" ),
                batTip( "Townsend's big-eared bat", "Corynorhinus townsendii", "Plecotus townsendii" ) );
        final PhylogenyNode vespertilionidae = namedClade( 0.04, "Vespertilionidae", "family",
                myotinae, vespertilioninae, plecotinae );
        final PhylogenyNode vespertilionoidea = namedClade( 0.03, "Vespertilionoidea", "superfamily",
                molossidae, miniopteridae, vespertilionidae );

        final PhylogenyNode yangochiro = namedClade( 0.04, "Yangochiroptera", "suborder",
                noctilionoidea, vespertilionoidea );
        final PhylogenyNode root = namedClade( 0, "Chiroptera", "order", yinptero, yangochiro );
        return tree( root, "Bat phylogeny (Chiroptera, demo)",
                     "A 34-species tree of the bats (order Chiroptera) across both suborders and eight families. Every "
                             + "TIP carries a full taxonomy -- common name, scientific name, and a taxonomic synonym (the "
                             + "original combination, e.g. Vespertilio pipistrellus) -- so turn on Display Data: Taxonomy "
                             + "to see the italic Latin names, or search the synonyms. Every INTERNAL node is annotated "
                             + "with its clade + rank (order / suborder / superfamily / family / subfamily / genus), so "
                             + "Tools > Colorize Subtrees via Taxonomic Rank -- and Annotate Clades by Rank -- work OFFLINE "
                             + "at 'family' (Pteropodidae, Rhinolophidae, Phyllostomidae, Vespertilionidae, ...) or 'genus', "
                             + "and Settings > Display > Internal Taxonomy Key lists the backbone. Schematic (illustrative "
                             + "branch lengths), not a rigorous phylogeny." );
    }

    /** A bat (or any) species tip: node name = common name; taxonomy carries the scientific (Latin) name, the common
     *  name, rank 'species', and a taxonomic synonym (the original combination). A null synonym is omitted. */
    private static PhylogenyNode batTip( final String common_name, final String scientific_name, final String synonym )
            throws PhyloXmlDataFormatException {
        final PhylogenyNode n = leaf( common_name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        t.setCommonName( common_name );
        t.setRank( "species" );
        if ( synonym != null ) {
            t.getSynonyms().add( synonym );
        }
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    /** An internal clade node with a branch length AND a clade-name taxonomy at the given rank (order / family / ...),
     *  so the rank-based colorize / annotate tools and the Internal Taxonomy Key read it directly. */
    private static PhylogenyNode namedClade( final double branch_length, final String scientific_name, final String rank,
                                             final PhylogenyNode... kids ) throws PhyloXmlDataFormatException {
        final PhylogenyNode n = clade( branch_length, kids );
        taxon( n, scientific_name, rank );
        return n;
    }

    // ----- "Lagomorph time tree": a time-calibrated species tree of the rabbits, hares and pikas (order Lagomorpha),
    //       ages in Ma, all tips extant (0 Ma) so the axis reaches the present. Every tip carries a common + scientific
    //       name; every internal node is dated AND annotated with its clade + rank, so it auto-detects as a time tree and
    //       Colorize / Annotate Clades by Rank works offline. Turn on Settings > Overlays > Time axis: Geologic for the
    //       ICS geologic scale (the Lagomorpha crown is Eocene). Schematic divergence times, not a rigorous chronogram.
    private static Phylogeny lagomorphTimeTree() throws PhyloXmlDataFormatException {
        // --- Ochotonidae (pikas) : the Ochotona radiation (crown ~15 Ma, mid-Miocene). Only the genus Ochotona is
        //     extant, so the family crown IS the genus crown -- one node, labelled at family rank (the Leporidae vs
        //     Ochotonidae split is the headline for the family-rank colorize demo). ---
        final PhylogenyNode ochotonidae = datedTaxonNode( "Ochotonidae", "pikas", "family", 15,
                datedNamedTip( "American pika", "Ochotona princeps", 0 ),
                datedNamedTip( "Collared pika", "Ochotona collaris", 0 ),
                datedNamedTip( "Alpine pika", "Ochotona alpina", 0 ),
                datedNamedTip( "Plateau pika", "Ochotona curzoniae", 0 ),
                datedNamedTip( "Afghan pika", "Ochotona rufescens", 0 ) );

        // --- Leporidae (rabbits & hares), crown ~25 Ma (late Oligocene) ---
        final PhylogenyNode lepus = datedTaxonNode( "Lepus", "hares", "genus", 5, // hare crown radiation, Pliocene
                datedNamedTip( "European hare", "Lepus europaeus", 0 ),
                datedNamedTip( "Mountain hare", "Lepus timidus", 0 ),
                datedNamedTip( "Arctic hare", "Lepus arcticus", 0 ),
                datedNamedTip( "Snowshoe hare", "Lepus americanus", 0 ),
                datedNamedTip( "Black-tailed jackrabbit", "Lepus californicus", 0 ) );
        final PhylogenyNode sylvilagus = datedTaxonNode( "Sylvilagus", "cottontails", "genus", 7, // late Miocene
                datedNamedTip( "Eastern cottontail", "Sylvilagus floridanus", 0 ),
                datedNamedTip( "Desert cottontail", "Sylvilagus audubonii", 0 ) );
        // core Leporinae: hares + cottontails + the European rabbit + the pygmy rabbit
        final PhylogenyNode oryctolagus = datedNamedTip( "European rabbit", "Oryctolagus cuniculus", 0 );
        final PhylogenyNode brachylagus = datedNamedTip( "Pygmy rabbit", "Brachylagus idahoensis", 0 );
        final PhylogenyNode rabbits = datedTaxonNode( "", "", "", 10, oryctolagus, brachylagus );
        final PhylogenyNode hares_cottontails = datedTaxonNode( "", "", "", 12, lepus, sylvilagus );
        final PhylogenyNode crown_leporids = datedTaxonNode( "", "", "", 18, hares_cottontails, rabbits );
        // relict / basal leporid genera branching earlier
        final PhylogenyNode relicts = datedTaxonNode( "", "", "", 20,
                datedNamedTip( "Volcano rabbit", "Romerolagus diazi", 0 ),
                datedNamedTip( "Amami rabbit", "Pentalagus furnessi", 0 ),
                datedNamedTip( "Sumatran striped rabbit", "Nesolagus netscheri", 0 ),
                datedNamedTip( "Hispid hare", "Caprolagus hispidus", 0 ) );
        final PhylogenyNode leporidae = datedTaxonNode( "Leporidae", "rabbits & hares", "family", 25, crown_leporids,
                relicts );

        final PhylogenyNode root = datedTaxonNode( "Lagomorpha", "rabbits, hares & pikas", "order", 50, ochotonidae,
                leporidae );
        setTimeBranchLengths( root, 50 );
        final Phylogeny phy = tree( root, "Lagomorph time tree (demo)",
                                    "A time-calibrated species tree of the rabbits, hares and pikas (order Lagomorpha), "
                                            + "ages in Ma; all 18 tips are extant (0 Ma) so the axis reaches the present. "
                                            + "It auto-detects as a time tree; turn on Settings > Overlays > Time axis: "
                                            + "Geologic to draw the ICS geologic scale (the pika/rabbit split is Eocene, "
                                            + "the Lepus radiation Pliocene). Every tip carries a common + scientific name "
                                            + "and every internal node is dated and annotated with its clade + rank "
                                            + "(order Lagomorpha; families Leporidae, Ochotonidae; genera Lepus, "
                                            + "Sylvilagus, Ochotona), so Colorize / Annotate Clades by Rank works offline. "
                                            + "Schematic divergence times, not a rigorous chronogram." );
        phy.setDistanceUnit( "My" );
        return phy;
    }

    /** A dated tip with a common + scientific name: node name = common name, a phyloXML <date> (age in Ma), and a
     *  species-rank taxonomy. */
    private static PhylogenyNode datedNamedTip( final String common_name, final String scientific_name,
                                                final double age_ma ) throws PhyloXmlDataFormatException {
        final PhylogenyNode n = datedTip( common_name, age_ma );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        t.setCommonName( common_name );
        t.setRank( "species" );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    /** A dated internal node carrying a clade-name TAXONOMY at the given rank, plus an optional common name (blank
     *  sci-name/rank = an unlabelled dated internal node). The clade name lives ONLY in the taxonomy (not the node name
     *  too), so the label isn't drawn twice; taxonomy display shows the clade (+ common name) + rank, and the
     *  rank-colorize / annotate tools read it. */
    private static PhylogenyNode datedTaxonNode( final String scientific_name, final String common_name,
                                                 final String rank, final double age_ma,
                                                 final PhylogenyNode... kids ) throws PhyloXmlDataFormatException {
        final PhylogenyNode n = new PhylogenyNode();
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", java.math.BigDecimal.valueOf( age_ma ),
                null, null, "mya" ) );
        for ( final PhylogenyNode k : kids ) {
            n.addAsChild( k );
        }
        if ( ( scientific_name != null ) && ( scientific_name.length() > 0 ) ) {
            taxon( n, scientific_name, rank );
            if ( ( common_name != null ) && ( common_name.length() > 0 ) ) {
                n.getNodeData().getTaxonomy().setCommonName( common_name );
                // display the common name via the node NAME (like the tips) -- it differs from the italic scientific
                // name, so the clade shows e.g. "Leporidae rabbits & hares" with no doubling
                n.setName( common_name );
            }
        }
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
    //       alternating row bands help track a label across to its annotation columns. Load -> Settings > Display >
    //       Zebra Stripes (optionally Tools > Annotation Fields).
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
                        + "across to its annotation columns (Tools > Annotation Fields)." );
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

    // ----- "Break Long Branches": an ingroup with normal (~0.1-0.3) branches plus one distant, fast-evolving outgroup
    //       on a ~4.0 branch that squashes the rest. Turn on Settings > Layout > Break Long Branches.
    private static Phylogeny longBranchBreakTree() {
        final PhylogenyNode primates = conf( clade( 0.20, blLeaf( "Homo_sapiens", 0.18 ),
                                                    blLeaf( "Pan_troglodytes", 0.15 ) ), 98 );
        final PhylogenyNode rodents = conf( clade( 0.22, blLeaf( "Mus_musculus", 0.25 ),
                                                   blLeaf( "Rattus_norvegicus", 0.20 ) ), 95 );
        final PhylogenyNode mammals = conf( clade( 0.15, primates, rodents ), 99 );
        final PhylogenyNode others = conf( clade( 0.18, blLeaf( "Gallus_gallus", 0.30 ),
                                                  blLeaf( "Xenopus_laevis", 0.28 ) ), 88 );
        final PhylogenyNode root = clade( 0, mammals, others, blLeaf( "Fast_evolving_outgroup", 4.0 ) );
        final Phylogeny phy = tree( root, "Break long branches (demo)",
                "Synthetic phylogram whose ingroup branches are ~0.1-0.3 substitutions/site (with bootstrap support on "
                        + "the internal nodes) while one distant, fast-evolving outgroup sits on a ~4.0 branch that "
                        + "squashes the rest. Turn on Settings > Layout > Break Long Branches: the outgroup is drawn "
                        + "shortened with a break mark (its true length is still shown as the label) and the ingroup "
                        + "reclaims the freed width; the support values stay clear of the break mark." );
        phy.setDistanceUnit( "substitutions/site" );
        return phy;
    }

    // ----- "Extract Dates from Labels": tips carry the sampling DATE in the label (mixed formats), no structured <date>
    private static Phylogeny dateInLabelsTree() {
        final PhylogenyNode root = clade( 0,
                clade( 0.2, blLeaf( "hCoV-19/China/WH-01/2019-12-30", 0.1 ),
                       blLeaf( "hCoV-19/Japan/TY/2020-01", 0.2 ),
                       blLeaf( "hCoV-19/Italy/LOM-2020-02-20", 0.3 ) ),
                clade( 0.5, blLeaf( "hCoV-19/USA/CA-2020-03-15", 0.4 ), blLeaf( "hCoV-19/UK/ENG-2020-09-20", 0.9 ),
                       clade( 1.0, blLeaf( "hCoV-19/India/delta/2021-04-10", 0.5 ),
                              blLeaf( "hCoV-19/SouthAfrica/omicron/2021-11-25", 1.2 ),
                              blLeaf( "hCoV-19/USA/BA2/2022-02-14", 1.4 ) ) ),
                clade( 0.1, blLeaf( "A/swine/2019", 0.1 ), blLeaf( "GISAID/01-Dec-2021", 2.0 ) ) );
        return tree( root, "Dates in tip labels (demo)",
                "Synthetic SARS-CoV-2-like tree whose tips carry the SAMPLING DATE in the label (mixed formats: ISO "
                        + "2020-02-20, ISO year-month 2020-01, bare year 2019, month-name 01-Dec-2021) but NO "
                        + "structured <date>. Run Tools > Extract Dates from Labels... to pull each date into a <date> "
                        + "+ a Color-by-able data:date property; the tree then gets the Calendar axis." );
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

    // ----- "Geologic time axis": a schematic time-calibrated archosaur tree (ages in Ma), including Archaeopteryx in
    //       the Late Jurassic. Extant tips (Crocodylia, modern birds) sit at 0 Ma so the axis reaches the present;
    //       branch lengths are age differences (a proper time tree). Turn on Settings > Overlays > Time axis: Geologic.
    private static Phylogeny dinosaurTimeTree() {
        final PhylogenyNode croc = datedTip( "Crocodylia", 0 );
        final PhylogenyNode birds = datedTip( "Aves (modern birds)", 0 );
        final PhylogenyNode archaeopteryx = datedTip( "Archaeopteryx", 150 );
        final PhylogenyNode trex = datedTip( "Tyrannosaurus", 66 );
        final PhylogenyNode diplodocus = datedTip( "Diplodocus", 150 );
        final PhylogenyNode stegosaurus = datedTip( "Stegosaurus", 150 );
        final PhylogenyNode triceratops = datedTip( "Triceratops", 66 );
        final PhylogenyNode paraves = datedNode( "Paraves", 165, archaeopteryx, birds );
        final PhylogenyNode theropoda = datedNode( "Theropoda", 231, trex, paraves );
        final PhylogenyNode saurischia = datedNode( "Saurischia", 232, diplodocus, theropoda );
        final PhylogenyNode ornithischia = datedNode( "Ornithischia", 235, stegosaurus, triceratops );
        final PhylogenyNode dinosauria = datedNode( "Dinosauria", 240, ornithischia, saurischia );
        final PhylogenyNode root = datedNode( "Archosauria", 250, croc, dinosauria );
        setTimeBranchLengths( root, 250 ); // branch length = parent age - node age -> a proper time phylogram
        final Phylogeny phy = tree( root, "Dinosaur time tree (demo)",
                                    "A schematic, time-calibrated archosaur tree (ages in Ma), with Archaeopteryx in "
                                            + "the Late Jurassic. It is auto-detected as a time tree; turn on "
                                            + "Settings > Overlays > Time axis: Geologic to draw the coloured ICS "
                                            + "geologic axis (Period over Epoch). Extant tips (Crocodylia, modern "
                                            + "birds) are at 0 Ma so the axis reaches the present; the extinct taxa "
                                            + "sit at their ages (Late Jurassic / Late Cretaceous). Schematic, not a "
                                            + "rigorous phylogeny." );
        phy.setDistanceUnit( "My" );
        return phy;
    }

    /** A schematic horse (equid) phylogeny where each tip is a genus carrying a stratigraphic RANGE (First/Last
     *  Appearance Datum) as a phyloXML date value+min+max, so the "Fossil Range Bars (FAD/LAD)" overlay draws each
     *  taxon's known duration as a capped bar. Each genus is placed at its LAD (last appearance), so the range bar runs
     *  back over the terminal branch toward the older FAD and the tip label sits clear to the right; the extant Equus
     *  reaches 0 Ma, so the tree is a proper time phylogram and the geologic axis reaches the present. */
    private static Phylogeny fossilRangeTree() {
        // fossilTip( name, placement age = LAD, FAD (oldest), LAD (youngest) ) -- all in Ma
        final PhylogenyNode eohippus = fossilTip( "Eohippus", 48, 55, 48 );
        final PhylogenyNode mesohippus = fossilTip( "Mesohippus", 30, 38, 30 );
        final PhylogenyNode merychippus = fossilTip( "Merychippus", 11, 17, 11 );
        final PhylogenyNode pliohippus = fossilTip( "Pliohippus", 6, 12, 6 );
        final PhylogenyNode equus = fossilTip( "Equus", 0, 4.5, 0 ); // extant -> placed at the present, range back to 4.5 Ma
        final PhylogenyNode n4 = datedNode( "", 13, pliohippus, equus );
        final PhylogenyNode n3 = datedNode( "", 20, merychippus, n4 );
        final PhylogenyNode n2 = datedNode( "", 42, mesohippus, n3 );
        final PhylogenyNode root = datedNode( "Equidae", 58, eohippus, n2 );
        setTimeBranchLengths( root, 58 );
        final Phylogeny phy = tree( root, "Fossil ranges (equid evolution, demo)",
                                    "A schematic horse (Equidae) evolutionary tree (ages in Ma) where each genus carries "
                                            + "a stratigraphic range -- its First Appearance Datum (oldest) and Last "
                                            + "Appearance Datum (youngest) -- as a phyloXML date min/max. It is "
                                            + "auto-detected as a time tree and \"Fossil Range Bars (FAD/LAD)\" turns on, "
                                            + "drawing each taxon's known duration as a capped bar. Turn on Settings > "
                                            + "Overlays > Time axis: Geologic to read the ranges against the coloured ICS "
                                            + "geologic axis; the extant Equus reaches 0 Ma so the axis reaches the "
                                            + "present. Schematic, not a rigorous phylogeny or literal ranges." );
        phy.setDistanceUnit( "My" );
        return phy;
    }

    /** A schematic, FOSSIL-ONLY ammonite tree: every tip is an extinct Mesozoic genus and the youngest tips die at the
     *  end-Cretaceous (K-Pg, 66 Ma), so NO tip reaches the present (age 0). This exercises the fossil-only geologic-axis
     *  alignment (maxDistanceToRoot < rootAge): the bands must line up with the branches over [66, 250] Ma rather than
     *  wrongly pinning age 0 to the youngest tip. Each genus also carries a FAD/LAD range for the fossil range bars. */
    private static Phylogeny ammoniteTimeTree() {
        // fossilTip( name, placement age = LAD, FAD (oldest), LAD (youngest) ) -- all in Ma
        final PhylogenyNode ceratites = fossilTip( "Ceratites", 235, 245, 235 );       // Middle Triassic
        final PhylogenyNode psiloceras = fossilTip( "Psiloceras", 199, 201, 199 );     // earliest Jurassic
        final PhylogenyNode perisphinctes = fossilTip( "Perisphinctes", 152, 163, 152 ); // Late Jurassic
        final PhylogenyNode baculites = fossilTip( "Baculites", 66, 83, 66 );          // Late Cretaceous -> K-Pg
        final PhylogenyNode scaphites = fossilTip( "Scaphites", 66, 90, 66 );          // Late Cretaceous -> K-Pg
        final PhylogenyNode n4 = datedNode( "", 95, baculites, scaphites );
        final PhylogenyNode n3 = datedNode( "", 170, perisphinctes, n4 );
        final PhylogenyNode n2 = datedNode( "", 210, psiloceras, n3 );
        final PhylogenyNode root = datedNode( "Ammonoidea", 250, ceratites, n2 );
        setTimeBranchLengths( root, 250 );
        final Phylogeny phy = tree( root, "Ammonite time tree (fossil-only, demo)",
                                    "A schematic ammonite (Ammonoidea) tree (ages in Ma) in which EVERY genus is extinct "
                                            + "-- the youngest tips (Baculites, Scaphites) die at the end-Cretaceous "
                                            + "(K-Pg, 66 Ma), so NO tip reaches the present. Show it as a phylogram (P) "
                                            + "and turn on Settings > Overlays > Time axis: Geologic to read the tree "
                                            + "against the ICS intervals (Triassic / Jurassic / Cretaceous): the coloured "
                                            + "bands line up with the branches over 250-66 Ma, and the youngest tips sit "
                                            + "at the K-Pg boundary rather than at the present -- the fossil-only case. "
                                            + "Each genus also carries a FAD/LAD range (Fossil Range Bars). Schematic, "
                                            + "not a rigorous phylogeny or literal ranges." );
        phy.setDistanceUnit( "My" );
        return phy;
    }

    /** A schematic, time-calibrated tree of the three domains of life reaching back to ~3.8 Ga (LUCA, deep in the
     *  Archean) -- so the geologic axis picks the DEEP-TIME (Eon over Era) band pair and the Precambrian is not blank. */
    private static Phylogeny deepTimeTree() {
        final PhylogenyNode cyano = datedTip( "Cyanobacteria", 0 );
        final PhylogenyNode proteo = datedTip( "Proteobacteria", 0 );
        final PhylogenyNode firmi = datedTip( "Firmicutes", 0 );
        final PhylogenyNode eury = datedTip( "Euryarchaeota", 0 );
        final PhylogenyNode cren = datedTip( "Crenarchaeota", 0 );
        final PhylogenyNode plants = datedTip( "Plantae", 0 );
        final PhylogenyNode animals = datedTip( "Animalia", 0 );
        final PhylogenyNode fungi = datedTip( "Fungi", 0 );
        final PhylogenyNode gracilicutes = datedNode( "", 3000, proteo, firmi );
        final PhylogenyNode bacteria = datedNode( "Bacteria", 3400, cyano, gracilicutes );
        final PhylogenyNode archaea = datedNode( "Archaea", 3200, eury, cren );
        final PhylogenyNode opistho = datedNode( "Opisthokonta", 1100, fungi, animals );
        final PhylogenyNode eukaryota = datedNode( "Eukaryota", 1800, plants, opistho );
        final PhylogenyNode neomura = datedNode( "", 3500, archaea, eukaryota );
        final PhylogenyNode luca = datedNode( "LUCA", 3800, bacteria, neomura );
        setTimeBranchLengths( luca, 3800 );
        final Phylogeny phy = tree( luca, "Tree of life (deep time, demo)",
                                    "A schematic, time-calibrated tree of the three domains of life (ages in Ma), with "
                                            + "LUCA near 3.8 Ga in the Archean and the crown eukaryotes in the "
                                            + "Proterozoic. Turn on Settings > Overlays > Time axis: Geologic to draw "
                                            + "the coloured ICS geologic axis; because the tree reaches deep time it "
                                            + "adapts to the coarse band pair (Eon over Era), so the Precambrian "
                                            + "(Archean, Proterozoic) is banded rather than blank. Extant tips at 0 Ma "
                                            + "reach the present. Schematic (rough molecular-clock ages), not a "
                                            + "rigorous phylogeny." );
        phy.setDistanceUnit( "My" );
        return phy;
    }

    /** A schematic, time-calibrated SARS-CoV-2 phylodynamic tree with TIP DATES as calendar years (2019.98 .. 2022.5),
     *  so the CALENDAR (absolute-date) time axis maps the tree to calendar time. Non-ultrametric (tips sampled at
     *  different dates) -- the tip-dated case the calendar axis targets. */
    private static Phylogeny sarsCov2TimeTree() {
        final PhylogenyNode wuhan = calTip( "Wuhan-Hu-1", 2019.98 );
        final PhylogenyNode early = calTip( "B.1 (early 2020)", 2020.2 );
        final PhylogenyNode alpha = calTip( "Alpha (B.1.1.7)", 2020.95 );
        final PhylogenyNode delta = calTip( "Delta (B.1.617.2)", 2021.5 );
        final PhylogenyNode ba1 = calTip( "Omicron BA.1", 2021.95 );
        final PhylogenyNode ba2 = calTip( "Omicron BA.2", 2022.3 );
        final PhylogenyNode ba5 = calTip( "Omicron BA.5", 2022.5 );
        final PhylogenyNode ba2ba5 = calNode( "", 2022.0, ba2, ba5 );
        final PhylogenyNode omicron = calNode( "Omicron", 2021.0, ba1, ba2ba5 );
        final PhylogenyNode delta_omi = calNode( "", 2020.5, delta, omicron );
        final PhylogenyNode variants = calNode( "", 2020.3, alpha, delta_omi );
        final PhylogenyNode post_wuhan = calNode( "", 2020.0, early, variants );
        final PhylogenyNode root = calNode( "SARS-CoV-2 MRCA", 2019.9, wuhan, post_wuhan );
        setCalendarBranchLengths( root, 2019.9 );
        final Phylogeny phy = tree( root, "SARS-CoV-2 time tree (demo)",
                                    "A schematic, time-calibrated SARS-CoV-2 phylodynamic tree: the tips are named "
                                            + "lineages/variants with their sampling dates as calendar years (Wuhan-Hu-1 "
                                            + "late 2019 through the Omicron sub-lineages in 2022), so the tree is "
                                            + "TIP-DATED (non-ultrametric). Turn on Settings > Overlays > Time axis: "
                                            + "Calendar to draw the labelled calendar-year axis (the most-recent tip = "
                                            + "the present is derived automatically). Works in every layout, including "
                                            + "concentric year rings in the circular view. Schematic, not a rigorous "
                                            + "phylogeny." );
        phy.setDistanceUnit( "years" );
        return phy;
    }

    private static PhylogenyNode calTip( final String name, final double year ) {
        final PhylogenyNode n = leaf( name );
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", java.math.BigDecimal.valueOf( year ), null,
                null, "year" ) );
        return n;
    }

    private static PhylogenyNode calNode( final String name, final double year, final PhylogenyNode... kids ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", java.math.BigDecimal.valueOf( year ), null,
                null, "year" ) );
        for ( final PhylogenyNode k : kids ) {
            n.addAsChild( k );
        }
        return n;
    }

    /** Sets each node's branch length to (node calendar year - parent calendar year), so distance-from-root == elapsed
     *  calendar time (a proper tip-dated time phylogram). */
    private static void setCalendarBranchLengths( final PhylogenyNode node, final double parent_year ) {
        final double year = node.getNodeData().getDate().getValue().doubleValue();
        node.setDistanceToParent( year - parent_year );
        for ( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            setCalendarBranchLengths( node.getChildNode( i ), year );
        }
    }

    private static PhylogenyNode datedTip( final String name, final double age_ma ) {
        final PhylogenyNode n = leaf( name );
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", java.math.BigDecimal.valueOf( age_ma ),
                null, null, "mya" ) );
        return n;
    }

    private static PhylogenyNode datedNode( final String name, final double age_ma, final PhylogenyNode... kids ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", java.math.BigDecimal.valueOf( age_ma ),
                null, null, "mya" ) );
        for ( final PhylogenyNode k : kids ) {
            n.addAsChild( k );
        }
        return n;
    }

    /** A fossil tip placed at {@code placement_ma} (its point age on the time axis) carrying a stratigraphic RANGE as a
     *  phyloXML date value + min (LAD, youngest) + max (FAD, oldest) -- what the Fossil Range Bar spans. */
    private static PhylogenyNode fossilTip( final String name, final double placement_ma, final double fad_ma,
                                            final double lad_ma ) {
        final PhylogenyNode n = leaf( name );
        n.getNodeData().setDate( new org.forester.phylogeny.data.Date( "", java.math.BigDecimal.valueOf( placement_ma ),
                java.math.BigDecimal.valueOf( lad_ma ), java.math.BigDecimal.valueOf( fad_ma ), "mya" ) );
        return n;
    }

    /** Sets each node's branch length to (parent age - node age), so distance-from-root == elapsed time. */
    private static void setTimeBranchLengths( final PhylogenyNode node, final double parent_age_ma ) {
        final double age = node.getNodeData().getDate().getValue().doubleValue();
        node.setDistanceToParent( parent_age_ma - age );
        for ( int i = 0; i < node.getNumberOfDescendants(); ++i ) {
            setTimeBranchLengths( node.getChildNode( i ), age );
        }
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
                             + "phylogram, then Tools > Annotation Fields and add s1..s6 as \"Heat map (matrix)\": "
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

    // ---- sequence alignment (embedded aligned mol_seqs; shown with Settings -> Overlays -> Sequence Alignment) ------

    private static Phylogeny alignmentTree() {
        final String base = "MKTAYIAKQR-QISFVKSHFSRQLEERLGLIEVQ"; // 34 columns, one aligned gap (protein)
        final Phylogeny phy = new Phylogeny();
        phy.setName( "Sequence Alignment" );
        phy.setDescription( "Synthetic 34-column protein alignment (one aligned gap) beside a six-tip tree. Turn it on "
                + "with Settings > Overlays > Sequence Alignment. HOVER A RESIDUE for its alignment column, its number "
                + "within that sequence's own ungapped residues, its full name and physico-chemical class, and its "
                + "Kyte-Doolittle hydropathy -- the column number is a property of the alignment, the residue number "
                + "is the coordinate that maps back onto the real protein. Under the alignment, the CONSERVATION "
                + "TRACK shows a bar per column and the consensus residue beneath it, scored over the tips currently "
                + "displayed: collapse a clade and the profile re-scores for what is left. Switch between consensus "
                + "identity and information content under Settings > Overlays > Conservation measure." );
        phy.setRoot( clade( 0.0,
                clade( 0.05,
                        alignedLeaf( "Human", 0.03, base ),
                        alignedLeaf( "Chimp", 0.03, sub( base, 24, 'D' ) ) ),
                clade( 0.05,
                        alignedLeaf( "Gorilla", 0.04, sub( base, 12, 'V' ) ),
                        clade( 0.04,
                                alignedLeaf( "Mouse", 0.06, sub( sub( base, 20, 'K' ), 24, 'D' ) ),
                                clade( 0.03,
                                        alignedLeaf( "Chicken", 0.08, sub( sub( base, 3, 'G' ), 25, 'K' ) ),
                                        alignedLeaf( "Frog", 0.09, sub( sub( base, 2, 'S' ), 30, 'V' ) ) ) ) ) ) );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode alignedLeaf( final String name, final double branch_length, final String aligned_seq ) {
        final PhylogenyNode n = blLeaf( name, branch_length );
        final Sequence s = new Sequence();
        s.setMolecularSequence( aligned_seq );
        s.setMolecularSequenceAligned( true );
        n.getNodeData().addSequence( s );
        return n;
    }

    /** {@code s} with the character at {@code i} replaced by {@code c} (same length -- keeps the alignment rectangular). */
    private static String sub( final String s, final int i, final char c ) {
        final char[] a = s.toCharArray();
        a[ i ] = c;
        return new String( a );
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


    // ----- "Animal tree of life (sponges to mammals)": a 25-tip backbone of the Metazoa, built specifically for
    //       NESTED clade annotation. Every internal clade is annotated at PHYLUM, CLASS or ORDER, which is exactly the
    //       three-rank ladder "Annotate Clades by Rank" can draw at once (genus/family/order is the equivalent inside a
    //       single order -- see the bat demo). All offline: no database lookup is needed at any of the three ranks.
    //
    //       TOPOLOGY: the backbone follows the standard Bilateria split (Protostomia = Ecdysozoa + Spiralia,
    //       Deuterostomia = Echinodermata + Chordata). The one genuinely CONTESTED node is the animal root: this tree
    //       places Ctenophora as sister to all other animals, following the chromosome-scale synteny evidence of
    //       Schultz et al. (2023, Nature 618:110-117), but the competing Porifera-sister hypothesis is very much alive
    //       -- later modelling work continues to favour it -- so that node should be read as unresolved, and the tree's
    //       description says so. Branch lengths are illustrative, not estimated.
    private static Phylogeny animalTreeOfLife() throws PhyloXmlDataFormatException {
        // --- Ctenophora: comb jellies (placed sister to the rest here; see the note above) ---
        final PhylogenyNode ctenophora = namedClade( 0.10, "Ctenophora", "phylum",
                namedClade( 0.06, "Tentaculata", "class",
                        namedClade( 0.05, "Cydippida", "order",
                                animalTip( "Sea gooseberry", "Pleurobrachia bachei" ) ),
                        namedClade( 0.05, "Lobata", "order",
                                animalTip( "Warty comb jelly", "Mnemiopsis leidyi" ) ) ) );

        // --- Porifera: sponges ---
        final PhylogenyNode porifera = namedClade( 0.09, "Porifera", "phylum",
                namedClade( 0.06, "Demospongiae", "class",
                        namedClade( 0.05, "Haplosclerida", "order",
                                animalTip( "Great Barrier Reef sponge", "Amphimedon queenslandica" ) ) ),
                namedClade( 0.06, "Calcarea", "class",
                        namedClade( 0.05, "Leucosolenida", "order",
                                animalTip( "Calcareous sponge", "Sycon ciliatum" ) ) ) );

        // --- Cnidaria: corals, hydras, jellyfish ---
        final PhylogenyNode cnidaria = namedClade( 0.07, "Cnidaria", "phylum",
                namedClade( 0.05, "Anthozoa", "class",
                        namedClade( 0.04, "Scleractinia", "order",
                                animalTip( "Staghorn coral", "Acropora digitifera" ) ) ),
                namedClade( 0.05, "Hydrozoa", "class",
                        namedClade( 0.04, "Anthoathecata", "order",
                                animalTip( "Freshwater hydra", "Hydra vulgaris" ) ) ) );

        // --- Ecdysozoa: the moulting animals ---
        final PhylogenyNode arthropoda = namedClade( 0.06, "Arthropoda", "phylum",
                namedClade( 0.04, "Insecta", "class",
                        namedClade( 0.04, "Diptera", "order",
                                animalTip( "Fruit fly", "Drosophila melanogaster" ) ),
                        namedClade( 0.04, "Hymenoptera", "order",
                                animalTip( "Western honey bee", "Apis mellifera" ) ) ),
                namedClade( 0.05, "Arachnida", "class",
                        namedClade( 0.04, "Araneae", "order",
                                animalTip( "Common house spider", "Parasteatoda tepidariorum" ) ) ) );
        final PhylogenyNode nematoda = namedClade( 0.09, "Nematoda", "phylum",
                namedClade( 0.05, "Chromadorea", "class",
                        namedClade( 0.05, "Rhabditida", "order",
                                animalTip( "Roundworm", "Caenorhabditis elegans" ),
                                animalTip( "Nematode", "Pristionchus pacificus" ) ) ) );
        final PhylogenyNode ecdysozoa = clade( 0.04, arthropoda, nematoda );

        // --- Spiralia (Lophotrochozoa): molluscs and annelids ---
        final PhylogenyNode mollusca = namedClade( 0.06, "Mollusca", "phylum",
                namedClade( 0.05, "Gastropoda", "class",
                        namedClade( 0.04, "Patellogastropoda", "order",
                                animalTip( "Owl limpet", "Lottia gigantea" ) ) ),
                namedClade( 0.05, "Bivalvia", "class",
                        namedClade( 0.04, "Ostreida", "order",
                                animalTip( "Pacific oyster", "Magallana gigas" ) ) ) );
        final PhylogenyNode annelida = namedClade( 0.06, "Annelida", "phylum",
                namedClade( 0.05, "Polychaeta", "class",
                        namedClade( 0.04, "Capitellida", "order",
                                animalTip( "Bristle worm", "Capitella teleta" ) ) ),
                namedClade( 0.05, "Clitellata", "class",
                        namedClade( 0.04, "Rhynchobdellida", "order",
                                animalTip( "Freshwater leech", "Helobdella robusta" ) ) ) );
        final PhylogenyNode spiralia = clade( 0.04, mollusca, annelida );
        final PhylogenyNode protostomia = clade( 0.04, ecdysozoa, spiralia );

        // --- Deuterostomia: echinoderms and chordates ---
        final PhylogenyNode echinodermata = namedClade( 0.07, "Echinodermata", "phylum",
                namedClade( 0.05, "Echinoidea", "class",
                        namedClade( 0.04, "Camarodonta", "order",
                                animalTip( "Purple sea urchin", "Strongylocentrotus purpuratus" ) ) ),
                namedClade( 0.05, "Asteroidea", "class",
                        namedClade( 0.04, "Forcipulatida", "order",
                                animalTip( "Common starfish", "Asterias rubens" ) ) ) );
        final PhylogenyNode mammalia = namedClade( 0.03, "Mammalia", "class",
                namedClade( 0.03, "Monotremata", "order",
                        animalTip( "Platypus", "Ornithorhynchus anatinus" ) ),
                clade( 0.01,
                        namedClade( 0.03, "Rodentia", "order",
                                animalTip( "House mouse", "Mus musculus" ) ),
                        clade( 0.01,
                                namedClade( 0.03, "Primates", "order",
                                        animalTip( "Human", "Homo sapiens" ) ),
                                namedClade( 0.03, "Carnivora", "order",
                                        animalTip( "Dog", "Canis lupus familiaris" ) ) ) ) );
        final PhylogenyNode amniota = clade( 0.02, mammalia,
                namedClade( 0.04, "Aves", "class",
                        namedClade( 0.03, "Galliformes", "order",
                                animalTip( "Red junglefowl", "Gallus gallus" ) ) ) );
        final PhylogenyNode tetrapoda = clade( 0.02, amniota,
                namedClade( 0.05, "Amphibia", "class",
                        namedClade( 0.04, "Anura", "order",
                                animalTip( "Western clawed frog", "Xenopus tropicalis" ) ) ) );
        final PhylogenyNode vertebrata = clade( 0.02, tetrapoda,
                namedClade( 0.05, "Actinopterygii", "class",
                        namedClade( 0.04, "Cypriniformes", "order",
                                animalTip( "Zebrafish", "Danio rerio" ) ) ) );
        final PhylogenyNode chordata = namedClade( 0.04, "Chordata", "phylum", vertebrata,
                namedClade( 0.07, "Ascidiacea", "class",
                        namedClade( 0.05, "Phlebobranchia", "order",
                                animalTip( "Sea squirt", "Ciona intestinalis" ) ) ) );
        final PhylogenyNode deuterostomia = clade( 0.04, echinodermata, chordata );

        final PhylogenyNode bilateria = clade( 0.04, protostomia, deuterostomia );
        final PhylogenyNode parahoxozoa = clade( 0.03, cnidaria, bilateria );
        final PhylogenyNode rest_of_animals = clade( 0.03, porifera, parahoxozoa );
        final PhylogenyNode root = namedClade( 0, "Metazoa", "kingdom", ctenophora, rest_of_animals );
        return tree( root, "Animal tree of life (sponges to mammals, demo)",
                     "A 25-species backbone of the animals (Metazoa), from sponges and comb jellies out to the mammals. "
                             + "Built for NESTED clade annotation: every internal clade is annotated at PHYLUM, CLASS or "
                             + "ORDER, so Tools > Annotate Clades by Rank can draw all three at once as nested bars or "
                             + "brackets -- pick 'order' plus 'class' plus 'phylum' and each phylum's classes take shades "
                             + "of that phylum's colour. Also works for a single rank, and for Colorize Subtrees via "
                             + "Taxonomic Rank. Everything resolves OFFLINE; no database lookup is needed. "
                             + "CAVEAT on the deepest node: this tree places Ctenophora (comb jellies) as sister to all "
                             + "other animals, following the chromosome-synteny evidence of Schultz et al. 2023 (Nature "
                             + "618:110-117), but the competing Porifera-sister (sponges-first) hypothesis remains well "
                             + "supported and the animal root is not settled -- read that split as unresolved. Branch "
                             + "lengths are illustrative, not estimated." );
    }

    /** A tip of the animal tree: common name shown, scientific name + species rank for the taxonomy machinery. */
    private static PhylogenyNode animalTip( final String common_name, final String scientific_name )
            throws PhyloXmlDataFormatException {
        final PhylogenyNode n = leaf( common_name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        t.setCommonName( common_name );
        t.setRank( "species" );
        n.getNodeData().setTaxonomy( t );
        return n;
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

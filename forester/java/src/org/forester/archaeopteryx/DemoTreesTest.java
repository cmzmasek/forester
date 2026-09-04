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
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterConstants;
import org.forester.util.ForesterUtil;

/**
 * Anti-rot guard for the {@code forester/demo/} gallery: every committed demo tree (see
 * {@link org.forester.demo.DemoTreeGenerator}) must still parse cleanly -- the phyloXML demos additionally VALIDATE
 * against the phyloXML schema (see {@link #load}) -- AND still carry the data its feature needs to be demonstrable.
 * Headless -- parsing + model inspection, no GUI. If a parser or the property/taxonomy model changes in a way that
 * breaks a demo, this fails in the always-run suite instead of the demo silently rotting.
 *
 * <p>When a new demo tree is added (per the "one demo tree per new feature" convention), add a shape check here.
 */
public final class DemoTreesTest {

    private static final String DEMO_DIR = System.getProperty( "user.dir" ) + File.separator + "forester"
            + File.separator + "demo" + File.separator;

    /** Whether the phyloXML XSD is on the classpath (staged by Ant {@code copy_resources}). When present, {@link #load}
     *  VALIDATES each phyloXML demo against the schema; when absent (a raw IDE compile that skipped that step), it falls
     *  back to the lenient parser + a one-time note. A missing XSD is an ENVIRONMENT gap, NOT demo rot, so it must not
     *  hard-fail every demo -- mirrors {@link DemoTreesGalleryTest}, which SKIPS its bundle check rather than failing.
     *  Decided ONCE here, distinct from a genuine per-file parse/validation failure. */
    private static final boolean XSD_ON_CLASSPATH = xsdOnClasspath();

    private static boolean xsdOnClasspath() {
        final boolean present = PhyloXmlParser.class.getClassLoader()
                .getResource( ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE ) != null;
        if ( !present ) {
            System.out.println( "  [DemoTreesTest] NOTE: the phyloXML XSD (" + ForesterConstants.LOCAL_PHYLOXML_XSD_RESOURCE
                    + ") is not on the classpath (Ant copy_resources not staged?) -- checking the demos WITHOUT schema "
                    + "validation" );
        }
        return present;
    }

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
        ok &= hasHostlessTip( "color-by-property.xml" ); // feeds the legend's dashed "no value" row
        // annotation columns: categorical + numeric + text properties
        ok &= hasCategoricalRef( "annotation-columns.xml", "data:host" );
        ok &= hasNumericRef( "annotation-columns.xml", "data:viral_load" );
        ok &= hasCategoricalRef( "annotation-columns.xml", "data:clade" );
        // symbol columns: a binary flag + a categorical host, both rendered as shape (SYMBOL) marks
        ok &= hasAtLeastTips( "label-properties.xml", 8 );
        // the point of this demo is a tip carrying MANY fields, one of them different on every tip (so it is
        // offered for the label but not as a colour) -- a regression that thinned it out would gut the demo
        ok &= hasCategoricalRef( "label-properties.xml", "data:host" );
        ok &= hasCategoricalRef( "label-properties.xml", "data:passage" );
        ok &= hasNumericRef( "label-properties.xml", "data:year" );
        ok &= hasSixFieldsWithAUniqueOne( "label-properties.xml" );
        ok &= hasAtLeastTips( "symbol-columns.xml", 8 );
        ok &= hasCategoricalRef( "symbol-columns.xml", "data:resistant" );
        ok &= hasCategoricalRef( "symbol-columns.xml", "data:host" );
        // stacked bar columns: three numeric read-count fields that MERGE into one segmented (stacked) bar per tip
        ok &= hasAtLeastTips( "stacked-bar-columns.xml", 8 );
        ok &= hasNumericRef( "stacked-bar-columns.xml", "data:firmicutes" );
        ok &= hasNumericRef( "stacked-bar-columns.xml", "data:bacteroidetes" );
        ok &= hasNumericRef( "stacked-bar-columns.xml", "data:proteobacteria" );
        // colorize by rank / clade bands: each TIP carries its 'order' in-tree, so it colorizes OFFLINE (no prompt)
        // into 3 clades -- even though the Carnivora clade root is deliberately mis-annotated 'Rodentia' (a tip's
        // own identity wins over the wrong internal annotation)
        ok &= hasRank( "colorize-by-rank.xml", "order" );
        ok &= colorizesOfflineInto( "colorize-by-rank.xml", "order", 3 );
        // scale axis: a phylogram with real branch lengths (so the labeled distance axis has ticks)
        ok &= hasBranchLengths( "scale-axis.xml" );
        // break long branches: a phylogram with branch lengths AND a genuine outlier branch (> the auto cap), so the
        // "Break Long Branches" option actually has something to break, plus internal bootstrap support (so the demo
        // shows the support values stay clear of the break mark)
        ok &= hasBranchLengths( "long-branch-break.xml" );
        ok &= hasOutlierBranch( "long-branch-break.xml" );
        ok &= hasInternalConfidence( "long-branch-break.xml" );
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

        // sequence alignment: tips carry equal-length aligned molecular sequences (the alignment shown beside the tree)
        ok &= alignmentDemoOk( "alignment.xml" );

        // bat phylogeny: a large taxonomy tree (common + scientific names + synonyms at the tips, ranks on the clades)
        ok &= batTreeOk( "bat-phylogeny.xml" );

        // animal tree of life: the NESTED clade-levels demo -- three ranks annotated at once (phylum/class/order)
        ok &= animalTreeOk( "animal-tree-of-life.xml" );

        // lagomorph time tree: a dated species tree with clade + rank annotations, for the geologic axis + rank colorize
        ok &= lagomorphTreeOk( "lagomorph-time-tree.xml" );
        ok &= hasBranchLengths( "lagomorph-time-tree.xml" );
        ok &= isDetectedTimeTree( "lagomorph-time-tree.xml", AptxUtil.TIME_TREE_KIND.DATED );

        // filoviridae: a REAL filovirus phylogeny (Ebola/Marburg reps) with repseq:* metadata -> color by species
        ok &= hasAtLeastTips( "filoviridae-tree.xml", 13 );
        ok &= hasBranchLengths( "filoviridae-tree.xml" );
        ok &= hasCategoricalRef( "filoviridae-tree.xml", "repseq:species" ); // the demo colors by this
        ok &= hasCategoricalRef( "filoviridae-tree.xml", "repseq:host" );    // real host metadata present

        ok &= hasAtLeastTips( "heatmap-matrix.xml", 6 );
        ok &= hasNumericRef( "heatmap-matrix.xml", "data:s1" );
        ok &= hasNumericRef( "heatmap-matrix.xml", "data:s6" );

        // import annotations: a plain-named tip tree + a companion CSV whose "name" column must join onto EVERY tip
        ok &= hasAtLeastTips( "import-annotations.xml", 12 );
        ok &= csvJoinMatchesAllTips( "import-annotations.xml", "import-annotations.csv" );

        // import GTDB taxonomy: accession-named genome tree + a GTDB-Tk table whose classifications must apply to
        // EVERY tip, giving gtdb:<rank> properties + a species taxonomy
        ok &= hasAtLeastTips( "gtdb-genomes.xml", 14 );
        ok &= gtdbDemoOk( "gtdb-genomes.xml", "gtdb-classifications.tsv" );

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

        // gene duplication / speciation demo: a gene-family tree whose tips carry species-rank taxonomy, with at least
        // one species appearing on >1 tip (a paralog) -- the shape that makes "Infer Duplications & Speciations (using
        // NCBI taxonomy)" recover a gene duplication
        ok &= hasAtLeastTips( "gene-duplication.xml", 6 );
        ok &= hasRank( "gene-duplication.xml", "species" );
        ok &= hasParalogSpecies( "gene-duplication.xml" );

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
        // fossil range bars: a dated equid tree whose TIPS carry a <date> min/max (FAD/LAD) stratigraphic range --
        // the "Fossil Range Bars (FAD/LAD)" overlay draws each tip's known duration
        ok &= hasBranchLengths( "fossil-range-bars.xml" );
        ok &= hasAtLeastTips( "fossil-range-bars.xml", 5 );
        ok &= isDetectedTimeTree( "fossil-range-bars.xml", AptxUtil.TIME_TREE_KIND.DATED );
        ok &= hasExternalDateInterval( "fossil-range-bars.xml" );
        // fossil-only geologic alignment: an all-extinct ammonite tree whose youngest tips die at the K-Pg (66 Ma),
        // so NO tip reaches the present -- the geologic axis must align to the branches (maxDist < rootAge)
        ok &= hasBranchLengths( "ammonite-time-tree.xml" );
        ok &= hasAtLeastTips( "ammonite-time-tree.xml", 5 );
        ok &= isDetectedTimeTree( "ammonite-time-tree.xml", AptxUtil.TIME_TREE_KIND.DATED );
        ok &= hasExternalDateInterval( "ammonite-time-tree.xml" );
        ok &= allExternalDatesAtLeast( "ammonite-time-tree.xml", 60.0 ); // fossil-only: every tip older than 60 Ma
        // deep-time geologic axis: a dated tree of life reaching into the Archean (oldest <date> > 2500 Ma), so the
        // geologic axis adapts to the coarse Eon/Era band pair -- the Precambrian is banded, not blank
        ok &= hasBranchLengths( "tree-of-life-deep-time.xml" );
        ok &= hasAtLeastTips( "tree-of-life-deep-time.xml", 8 );
        ok &= isDetectedTimeTree( "tree-of-life-deep-time.xml", AptxUtil.TIME_TREE_KIND.DATED );
        ok &= oldestDateAtLeast( "tree-of-life-deep-time.xml", 2500.0 );
        // calendar time axis: a tip-dated SARS-CoV-2 tree whose <date> values are CALENDAR YEARS (2019.98..2022.5) --
        // detected DATED, and the max date (the most-recent tip = the present) is a plausible calendar year
        ok &= hasBranchLengths( "sars-cov-2-time-tree.xml" );
        ok &= hasAtLeastTips( "sars-cov-2-time-tree.xml", 7 );
        ok &= isDetectedTimeTree( "sars-cov-2-time-tree.xml", AptxUtil.TIME_TREE_KIND.DATED );
        ok &= oldestDateAtLeast( "sars-cov-2-time-tree.xml", 2000.0 );
        // Auspice / Nextstrain v2 JSON: a synthetic nCoV dataset -- parses via the AuspiceJsonParser, is a dated
        // (calendar) tree, and carries a discrete 'region' trait with per-node confidence (-> ancestral-state pies)
        ok &= auspiceDemoOk( "nextstrain-ncov.json" );
        // Extract Dates from Labels: the tips carry the date in the LABEL (no structured <date>), and a strict
        // majority parse via TipDateExtractor -- so Tools > Extract Dates from Labels has work to do
        ok &= dateInLabelsDemoOk( "date-in-labels.xml" );
        // GUARD: no demo tree carries an aptx: property on a NODE -- Archaeopteryx app state (import profile,
        // time-axis config, ...) must ride the PHYLOGENY level (a <property applies_to="phylogeny">), never node data
        ok &= noDemoCarriesNodeLevelAptxProperty();
        return ok;
    }

    /** No committed demo {@code .xml} carries an {@code aptx:}-prefixed property on any NODE: app state belongs at the
     *  PHYLOGENY level (a {@code <property applies_to="phylogeny">}, auto-hidden from displays), so a node-level
     *  {@code aptx:} property would be an accidental-write regression. */
    private static boolean noDemoCarriesNodeLevelAptxProperty() {
        final File[] xmls = new File( DEMO_DIR ).listFiles( ( d, name ) -> name.endsWith( ".xml" ) );
        if ( xmls == null ) {
            return note( "could not list the demo directory " + DEMO_DIR );
        }
        boolean ok = true;
        for( final File f : xmls ) {
            final Phylogeny phy = load( f.getName() );
            if ( phy == null ) {
                ok = false;
                continue;
            }
            for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                final PropertiesList pl = it.next().getNodeData().getProperties();
                if ( pl != null ) {
                    for( final Property p : pl.getProperties() ) {
                        if ( ( p.getRef() != null ) && p.getRef().startsWith( "aptx:" ) ) {
                            ok = note( f.getName() + " carries a NODE-level aptx: property (" + p.getRef()
                                    + ") -- app state must be a phylogeny-level <property applies_to=\"phylogeny\">" );
                        }
                    }
                }
            }
        }
        return ok;
    }

    /** The date-in-labels demo has &ge;10 tips whose LABELS carry a parseable date but which have NO structured
     *  {@code <date>} yet -- exactly what the Extract-Dates-from-Labels tool consumes. */
    private static boolean dateInLabelsDemoOk( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return note( file_name + " is missing / unparseable" );
        }
        if ( phy.getNumberOfExternalNodes() < 10 ) {
            return note( file_name + " must have >= 10 tips, has " + phy.getNumberOfExternalNodes() );
        }
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( ext.getNodeData().isHasDate() ) {
                return note( file_name + " tips must NOT carry a structured <date> (the tool extracts it from the label)" );
            }
        }
        if ( !org.forester.archaeopteryx.tools.TipDateExtractor.mostLabelsHaveDates( phy ) ) {
            return note( file_name + " tip labels must carry parseable dates" );
        }
        return true;
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

    /** The alignment demo: enough tips, they carry aligned molecular sequences, and all aligned rows are the same
     *  length (a real alignment) with a reasonable number of columns. */
    private static boolean alignmentDemoOk( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( phy.getNumberOfExternalNodes() < 6 ) {
            return note( file_name + " must have at least 6 tips" );
        }
        if ( !AptxUtil.hasAlignedSequences( phy ) ) {
            return note( file_name + " must carry aligned molecular sequences (the alignment)" );
        }
        int len = -1;
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( ext.getNodeData().isHasSequence() && ext.getNodeData().getSequence().isMolecularSequenceAligned() ) {
                final String mol = ext.getNodeData().getSequence().getMolecularSequence();
                if ( len < 0 ) {
                    len = mol.length();
                }
                else if ( mol.length() != len ) {
                    return note( file_name + " aligned rows must all be the same length (it is an alignment)" );
                }
            }
        }
        if ( len < 10 ) {
            return note( file_name + " alignment must have at least 10 columns" );
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

    /** The tree has a genuine outlier branch: at least one branch exceeds the auto break cap
     *  ({@code LONG_BRANCH_BREAK_MULTIPLIER} * median positive branch length). Delegates to the PRODUCTION helpers so
     *  the demo guard tracks the shipping threshold. */
    private static boolean hasOutlierBranch( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        final double cap = TreePanelUtil.longBranchBreakCap( phy, TreePanel.LONG_BRANCH_BREAK_MULTIPLIER );
        if ( cap <= 0 ) {
            return note( file_name + " must have positive branch lengths (a break cap)" );
        }
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            if ( it.next().getDistanceToParent() > cap ) {
                return true;
            }
        }
        return note( file_name + " must carry an outlier branch longer than the auto break cap " + cap );
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

    /** At least one EXTERNAL tip carries a {@code <date>} with a min/max interval (a fossil FAD/LAD range) -- what the
     *  Fossil Range Bars overlay draws. Delegates to the PRODUCTION predicate that drives the auto-enable, so the demo
     *  guard and the shipping code can't drift. */
    private static boolean hasExternalDateInterval( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( AptxUtil.isHasAtLeastOneExternalNodeWithDateInterval( phy ) ) {
            return true;
        }
        return note( file_name + " must carry an EXTERNAL-tip <date> with a min/max interval (for fossil range bars)" );
    }

    /** Every external tip's {@code <date>} value is at least {@code min_ma} Ma -- i.e. the tree is FOSSIL-ONLY (no extant
     *  tip at/near the present), the case the fossil-only geologic alignment targets. */
    private static boolean allExternalDatesAtLeast( final String file_name, final double min_ma ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( !n.getNodeData().isHasDate() || ( n.getNodeData().getDate().getValue() == null )
                    || ( n.getNodeData().getDate().getValue().doubleValue() < min_ma ) ) {
                return note( file_name + " must be fossil-only: every tip's <date> value >= " + min_ma + " Ma" );
            }
        }
        return true;
    }

    /** The oldest node {@code <date>} value is at least {@code min_ma} Ma (so the geologic axis picks the deep-time
     *  band pair -- e.g. Eon/Era for an Archean tree). */
    private static boolean oldestDateAtLeast( final String file_name, final double min_ma ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        double oldest = 0;
        for( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasDate() && ( n.getNodeData().getDate().getValue() != null ) ) {
                oldest = Math.max( oldest, n.getNodeData().getDate().getValue().doubleValue() );
            }
        }
        if ( oldest < min_ma ) {
            return note( file_name + " oldest <date> must be >= " + min_ma + " Ma (deep time), was " + oldest );
        }
        return true;
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

    /**
     * The animal backbone exists to demo THREE ranks annotated at once, so the thing to guard is that all three are
     * really there and really nest: every tip taxonomised, and phylum / class / order all present on the internal
     * clades (offline, no lookup). Without any one of them the demo silently degrades to fewer columns than the
     * gallery entry promises.
     */
    private static boolean animalTreeOk( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( phy.getNumberOfExternalNodes() < 20 ) {
            return note( file_name + " must span at least 20 animals, has " + phy.getNumberOfExternalNodes() );
        }
        for( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( !ext.getNodeData().isHasTaxonomy()
                    || ForesterUtil.isEmpty( ext.getNodeData().getTaxonomy().getScientificName() ) ) {
                return note( file_name + ": every tip needs a scientific name (missing on " + ext.getName() + ")" );
            }
        }
        for( final String rank : new String[] { "phylum", "class", "order" } ) {
            if ( !hasRank( file_name, rank ) ) {
                return false; // hasRank reports which rank is missing
            }
        }
        // sponges AND mammals: the span the demo is named for
        return hasClade( file_name, "Porifera" ) && hasClade( file_name, "Mammalia" ) && hasClade( file_name,
                                                                                                  "Chordata" );
    }

    /** Whether some node carries {@code clade} as its scientific name (an internal clade annotation). */
    private static boolean hasClade( final String file_name, final String clade ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final PhylogenyNode n : phy.getNodesViaScientificName( clade ) ) {
            if ( n != null ) {
                return true;
            }
        }
        return note( file_name + " must contain the clade '" + clade + "'" );
    }

    /** The bat phylogeny: &ge;30 species; every TIP carries a taxonomy with a scientific name, a common name AND a
     *  synonym; and the INTERNAL nodes carry rank annotations at 'family' and 'genus' (so Colorize / Annotate Clades by
     *  Rank works offline). Guards the taxonomy richness the demo advertises. */
    private static boolean batTreeOk( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( phy.getNumberOfExternalNodes() < 30 ) {
            return note( file_name + " must have at least 30 species, has " + phy.getNumberOfExternalNodes() );
        }
        for( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( !ext.getNodeData().isHasTaxonomy() ) {
                return note( file_name + ": every tip must carry a taxonomy (missing on " + ext.getName() + ")" );
            }
            final org.forester.phylogeny.data.Taxonomy t = ext.getNodeData().getTaxonomy();
            if ( ForesterUtil.isEmpty( t.getScientificName() ) || ForesterUtil.isEmpty( t.getCommonName() )
                    || t.getSynonyms().isEmpty() ) {
                return note( file_name + ": every tip must have a scientific name, common name AND a synonym (tip '"
                        + ext.getName() + "')" );
            }
        }
        return hasRank( file_name, "family" ) && hasRank( file_name, "genus" ); // clade annotations for rank colorize
    }

    /** The lagomorph time tree: &ge;15 species, a dated time tree with branch lengths, and internal clade annotations at
     *  'family' and 'genus' (so the geologic axis + offline rank colorize both work). */
    private static boolean lagomorphTreeOk( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        if ( phy.getNumberOfExternalNodes() < 15 ) {
            return note( file_name + " must have at least 15 species, has " + phy.getNumberOfExternalNodes() );
        }
        for( final PhylogenyNode ext : phy.getExternalNodes() ) {
            final org.forester.phylogeny.data.Taxonomy t = ext.getNodeData().isHasTaxonomy()
                    ? ext.getNodeData().getTaxonomy() : null;
            if ( ( t == null ) || ForesterUtil.isEmpty( t.getScientificName() )
                    || ForesterUtil.isEmpty( t.getCommonName() ) ) {
                return note( file_name + ": every tip must carry a scientific + common name (tip '" + ext.getName()
                        + "')" );
            }
        }
        return hasRank( file_name, "family" ) && hasRank( file_name, "genus" ); // clade annotations for rank colorize
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

    /** The GTDB-Tk-style companion table parses, its classification column joins onto EVERY tip by name, and applying
     *  it (via {@link org.forester.archaeopteryx.tools.GtdbTaxonomy}) writes a gtdb:phylum property + a taxonomy onto
     *  the tips -- so the "Import GTDB Taxonomy" demo lights up offline. */
    private static boolean gtdbDemoOk( final String tree_file, final String tsv_file ) {
        final Phylogeny phy = load( tree_file );
        if ( phy == null ) {
            return false;
        }
        // the genome tree must start WITHOUT taxonomy (the whole point of importing it)
        if ( phy.getFirstExternalNode().getNodeData().isHasTaxonomy() ) {
            return note( tree_file + " tips should carry no taxonomy before the GTDB import" );
        }
        final File tsv = new File( DEMO_DIR + tsv_file );
        if ( !tsv.exists() ) {
            return note( tsv_file + " is missing from the demo gallery (" + tsv.getAbsolutePath() + ")" );
        }
        try {
            final NodeDataImporter.Table table = NodeDataImporter.parseTable( java.nio.file.Files.readString( tsv.toPath() ) );
            // find the GTDB classification column and pair it with the other (key) column
            int class_col = -1;
            for ( int c = 0; ( class_col < 0 ) && ( c < table.getColumnCount() ); ++c ) {
                for ( int r = 0; r < table.getRowCount(); ++r ) {
                    if ( org.forester.archaeopteryx.tools.GtdbTaxonomy.looksLikeGtdb( table.getCell( r, c ) ) ) {
                        class_col = c;
                        break;
                    }
                }
            }
            if ( class_col < 0 ) {
                return note( tsv_file + " has no GTDB classification column (a d__..;p__..;s__.. value)" );
            }
            final int key_col = ( class_col == 0 ) ? 1 : 0;
            final java.util.Map<String, String> map = new java.util.HashMap<String, String>();
            for ( int r = 0; r < table.getRowCount(); ++r ) {
                map.put( table.getCell( r, key_col ), table.getCell( r, class_col ) );
            }
            final int annotated = org.forester.archaeopteryx.tools.GtdbTaxonomy.applyByTipName( phy, map );
            if ( annotated != phy.getNumberOfExternalNodes() ) {
                return note( tsv_file + " should annotate all " + phy.getNumberOfExternalNodes() + " tips of " + tree_file
                        + ", annotated " + annotated );
            }
            // every tip must now carry a gtdb:phylum property + a taxonomy
            for ( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
                final org.forester.phylogeny.PhylogenyNode tip = it.next();
                if ( !tip.getNodeData().isHasTaxonomy() ) {
                    return note( tree_file + " tip " + tip.getName() + " has no taxonomy after the GTDB import" );
                }
                if ( propValue( tip, "gtdb:phylum" ) == null ) {
                    return note( tree_file + " tip " + tip.getName() + " has no gtdb:phylum property after the GTDB import" );
                }
            }
            return true;
        }
        catch ( final Exception e ) {
            return note( tsv_file + " could not be read/applied: " + e.getMessage() );
        }
    }

    private static String propValue( final org.forester.phylogeny.PhylogenyNode n, final String ref ) {
        if ( ( n.getNodeData() == null ) || ( n.getNodeData().getProperties() == null ) ) {
            return null;
        }
        for ( final org.forester.phylogeny.data.Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( ref.equals( p.getRef() ) ) {
                return p.getValue();
            }
        }
        return null;
    }

    /** The Auspice/Nextstrain v2 JSON demo parses via {@link org.forester.io.parsers.json.AuspiceJsonParser}, is a
     *  dated calendar tree with &ge;5 tips, and carries a discrete 'region' trait with confidence (-> pies). */
    private static boolean auspiceDemoOk( final String file_name ) {
        final File file = new File( DEMO_DIR + file_name );
        if ( !file.exists() ) {
            return note( file_name + " is missing from the demo gallery (" + file.getAbsolutePath() + ")" );
        }
        try {
            final org.forester.io.parsers.json.AuspiceJsonParser parser = new org.forester.io.parsers.json.AuspiceJsonParser();
            parser.setSource( file );
            final Phylogeny[] phys = parser.parse();
            if ( ( phys == null ) || ( phys.length != 1 ) || phys[ 0 ].isEmpty() ) {
                return note( file_name + " did not yield one non-empty tree" );
            }
            final Phylogeny phy = phys[ 0 ];
            if ( phy.getNumberOfExternalNodes() < 15 ) {
                return note( file_name + " must have >= 15 tips, has " + phy.getNumberOfExternalNodes() );
            }
            if ( AptxUtil.detectTimeTree( phy ) != AptxUtil.TIME_TREE_KIND.DATED ) {
                return note( file_name + " must be detected as a dated time tree" );
            }
            if ( !TreePanelUtil.ancestralStateTraits( phy ).contains( "region" ) ) {
                return note( file_name + " must carry a 'region' ancestral-state trait (for pies)" );
            }
            // num_date is exposed as a numeric property so it can drive the "Color by" date gradient -- it must show
            // up in the colour-by refs AND be typed numeric (a gradient), which is what makes "Color by date" work
            if ( !hasTipProperty( phy, "nextstrain:num_date" )
                    || !PropertyColorScheme.colorableRefs( phy ).contains( "nextstrain:num_date" )
                    || !PropertyColorScheme.numericRefs( phy ).contains( "nextstrain:num_date" ) ) {
                return note( file_name + " must offer nextstrain:num_date as a NUMERIC 'Color by' ref (date gradient)" );
            }
            // host metadata with real variety (mostly human + at least one non-human -> a useful 'Color by host')
            final java.util.Set<String> hosts = tipPropertyValues( phy, "nextstrain:host" );
            if ( hosts.size() < 2 || !hosts.contains( "Homo sapiens" ) ) {
                return note( file_name + " must carry a 'host' trait with >= 2 distinct values incl. Homo sapiens, got "
                        + hosts );
            }
            return true;
        }
        catch ( final Exception e ) {
            return note( file_name + " could not be read: " + e.getMessage() );
        }
    }

    /** Parse a demo phyloXML file, asserting it is present, error-free and a single non-empty tree; null on any failure
     *  (with a message). Uses the XSD-VALIDATING parser when the schema is on the classpath ({@link #XSD_ON_CLASSPATH} --
     *  the Ant/CI path), so a schema-invalid demo (a mis-ordered/misplaced element, or app state written on a node
     *  instead of the phylogeny) fails the always-run suite instead of surfacing only as a silent runtime quirk; falls
     *  back to the lenient parser when the XSD isn't staged. */
    private static Phylogeny load( final String file_name ) {
        final File file = new File( DEMO_DIR + file_name );
        if ( !file.exists() ) {
            return fail( file_name + " is missing from the demo gallery (" + file.getAbsolutePath() + ")" );
        }
        try {
            final PhyloXmlParser parser = XSD_ON_CLASSPATH ? PhyloXmlParser.createPhyloXmlParserXsdValidating()
                                                           : PhyloXmlParser.createPhyloXmlParser();
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

    /**
     * The "properties in labels" demo only makes its point if every tip really does carry SIX fields, one of them
     * ('data:accession') with a DIFFERENT value on every tip. That unique field is precisely what
     * {@code colorableRefs} drops and what the Annotation Fields chooser must still offer for the label, so a
     * regression that thinned the demo down -- or that made the field constant -- would leave the gallery green
     * while the demo no longer demonstrated anything.
     */
    private static boolean hasSixFieldsWithAUniqueOne( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        final java.util.List<String> refs = TreePanelUtil.userVisiblePropertyRefs( phy );
        if ( refs.size() < 6 ) {
            return note( file_name + " should offer at least 6 annotation fields, has " + refs.size() );
        }
        final java.util.Set<String> unique_values = new java.util.HashSet<String>();
        for( final PhylogenyNode leaf : phy.getExternalNodes() ) {
            final String acc = PropertyColorScheme.valueFor( leaf, "data:accession" );
            if ( acc == null ) {
                return note( file_name + " every tip needs a 'data:accession' value" );
            }
            unique_values.add( acc );
        }
        if ( unique_values.size() != phy.getExternalNodes().size() ) {
            return note( file_name + " 'data:accession' must differ on every tip (that is what makes it a "
                    + "label-only field), got " + unique_values.size() + " distinct values for "
                    + phy.getExternalNodes().size() + " tips" );
        }
        if ( PropertyColorScheme.colorableRefs( phy ).contains( "data:accession" ) ) {
            return note( file_name + " 'data:accession' should NOT be colorable -- the demo's point is that the "
                    + "chooser offers it for the LABEL anyway" );
        }
        return true;
    }

    /** The color-by demo must keep at least one tip WITHOUT a host, so the legend's "no value" row shows. */
    private static boolean hasHostlessTip( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        for( final PhylogenyNode leaf : phy.getExternalNodes() ) {
            if ( PropertyColorScheme.valueFor( leaf, "data:host" ) == null ) {
                return true;
            }
        }
        return note( file_name + " must keep a host-less tip (demos the legend's 'no value' row)" );
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

    /** True iff at least one species scientific name appears on more than one EXTERNAL tip (a paralog) -- the property
     *  the taxonomy reconciliation needs to be able to recover a gene duplication. */
    private static boolean hasParalogSpecies( final String file_name ) {
        final Phylogeny phy = load( file_name );
        if ( phy == null ) {
            return false;
        }
        final HashMap<String, Integer> counts = new HashMap<>();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.getNodeData().isHasTaxonomy() ) {
                final String sci = n.getNodeData().getTaxonomy().getScientificName();
                if ( ( sci != null ) && !sci.isEmpty() ) {
                    counts.merge( sci, 1, Integer::sum );
                }
            }
        }
        for( final Integer c : counts.values() ) {
            if ( c >= 2 ) {
                return true;
            }
        }
        return note( file_name + " must have at least one species on >1 tip (a paralog) for reconciliation" );
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

    private static boolean hasTipProperty( final Phylogeny phy, final String ref ) {
        return !tipPropertyValues( phy, ref ).isEmpty();
    }

    /** Distinct values of a given node property across the external nodes (tips). */
    private static java.util.Set<String> tipPropertyValues( final Phylogeny phy, final String ref ) {
        final java.util.Set<String> vals = new java.util.HashSet<>();
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( ext.getNodeData().getProperties() != null ) {
                for ( final org.forester.phylogeny.data.Property p : ext.getNodeData().getProperties()
                        .getProperties() ) {
                    if ( ref.equals( p.getRef() ) ) {
                        vals.add( p.getValue() );
                    }
                }
            }
        }
        return vals;
    }

    /** false-returning failure note for the shape checks. */
    private static boolean note( final String msg ) {
        System.out.println( "  [DemoTreesTest] " + msg );
        return false;
    }

    private DemoTreesTest() {
    }
}

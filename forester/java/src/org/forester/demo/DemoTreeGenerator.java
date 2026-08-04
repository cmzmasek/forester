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
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
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
        write( dir, "scale-axis.xml", scaleAxisTree() );
        write( dir, "node-hpd-bars.xml", hpdBarsTree() );
        write( dir, "zebra-stripes.xml", zebraStripesTree() );
        write( dir, "flip-vertically.xml", flipVerticallyTree() );
        System.out.println( "Wrote demo trees to " + dir.getAbsolutePath() );
    }

    private static void write( final File dir, final String file_name, final Phylogeny phy ) throws IOException {
        new PhylogenyWriter().toPhyloXML( phy, 0, new File( dir, file_name ) );
        System.out.println( "  " + file_name + " (" + phy.getNumberOfExternalNodes() + " tips)" );
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
        final PhylogenyNode carnivora = orderClade( "Carnivora", 0.08,
                                                     species( "Felis catus" ), species( "Panthera leo" ),
                                                     species( "Canis lupus" ) );
        final PhylogenyNode rodentia = orderClade( "Rodentia", 0.09,
                                                   species( "Mus musculus" ), species( "Rattus norvegicus" ) );
        final PhylogenyNode primates = orderClade( "Primates", 0.07,
                                                   species( "Homo sapiens" ), species( "Pan troglodytes" ),
                                                   species( "Macaca mulatta" ) );
        final PhylogenyNode root = clade( 0, clade( 0.04, carnivora, rodentia ), primates );
        return tree( root, "Colorize by taxonomic rank (demo)",
                     "Synthetic mammal tree. Tips carry species taxonomy; each order's clade root is annotated with "
                             + "rank 'order' (Carnivora, Rodentia, Primates), so Tools > Colorize Subtrees via "
                             + "Taxonomic Rank -- and Annotate Clades by Rank -- work at 'order' with no online lookup." );
    }

    private static PhylogenyNode orderClade( final String order, final double bl, final PhylogenyNode... kids )
            throws PhyloXmlDataFormatException {
        final PhylogenyNode n = clade( bl, kids );
        taxon( n, order, "order" );
        return n;
    }

    private static PhylogenyNode species( final String scientific_name ) throws PhyloXmlDataFormatException {
        final PhylogenyNode n = leaf( scientific_name );
        taxon( n, scientific_name, "species" );
        return n;
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

    // ----- "Flip Vertically": an 8-tip ladder (caterpillar) tree with sequentially-numbered tips, so the tip order
    //       (and the staircase) visibly inverts top<->bottom when flipped. Load -> Settings > Display > Flip Vertically.
    private static Phylogeny flipVerticallyTree() {
        // build the caterpillar bottom-up: tip_08/tip_07 at the deepest fork, then prepend tip_06..tip_01
        PhylogenyNode n = clade( 0.07, blLeaf( "tip_07", 0.06 ), blLeaf( "tip_08", 0.06 ) );
        for( int i = 6; i >= 1; i-- ) {
            n = clade( 0.07, blLeaf( String.format( Locale.ROOT, "tip_%02d", i ), 0.06 ), n );
        }
        return tree( n, "Flip vertically (demo)",
                "Synthetic 8-tip ladder tree with sequentially-numbered tips (tip_01 at the top). Turn on Settings > "
                        + "Display > Flip Vertically to reverse the tip order top-to-bottom -- the staircase inverts "
                        + "and tip_08 moves to the top. Display-only; the tree data is unchanged." );
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

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

import java.io.File;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.forester.archaeopteryx.TreeFacts.Group;
import org.forester.archaeopteryx.TreeFacts.Histogram;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Point;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

/**
 * Headless tests for {@link TreeFacts}: every group's numbers on a synthetic tree (structure, branch lengths
 * incl. zero/negative/ultrametric, one support group per type, annotation coverage down to property counts, the
 * time-axis group), the histogram binning, the file facts (format sniffing, size text), and the text helpers.
 */
public final class TreeFactsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeFacts: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return structure() && branchLengths() && support() && coverage() && timeAxis() && histogram()
                    && fileFacts() && helpers() && edgeCases();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // ---------------------------------------------------------------- fixture
    /**
     * root -> ( X(0.1, boot 90 / prob 0.9) -> ( A(0.2), B(0.3), C(0) ), Y(0.4, boot 70) -> ( D(-1 = no length),
     * E(0.5) ) ). A: taxonomy with id + sequence with mol seq + domains + property p + q; B: taxonomy without id
     * + sequence without mol seq + property p; C: date + distribution; D: nothing; E: property q. X: named, dated,
     * duplication event; Y: unnamed, speciation.
     */
    static Phylogeny fixture() {
        final PhylogenyNode root = new PhylogenyNode();
        root.setName( "root" );
        final PhylogenyNode x = new PhylogenyNode();
        x.setName( "X" );
        x.setDistanceToParent( 0.1 );
        x.getBranchData().addConfidence( new Confidence( 90, "bootstrap" ) );
        x.getBranchData().addConfidence( new Confidence( 0.9, "probability" ) );
        x.getNodeData().setDate( new Date( "x", new BigDecimal( "50" ), null, null, "mya" ) );
        x.getNodeData().setEvent( new Event( 2, 0, 1 ) );
        final PhylogenyNode y = new PhylogenyNode();
        y.setDistanceToParent( 0.4 );
        y.getBranchData().addConfidence( new Confidence( 70, "bootstrap" ) );
        y.getNodeData().setEvent( new Event( 0, 1, 0 ) );
        root.addAsChild( x );
        root.addAsChild( y );
        final PhylogenyNode a = tip( "A", 0.2 );
        final Taxonomy ta = new Taxonomy();
        ta.setScientificName( "Homo sapiens" );
        ta.setIdentifier( new Identifier( "9606", "ncbi" ) );
        a.getNodeData().addTaxonomy( ta );
        final Sequence sa = new Sequence();
        sa.setName( "A_seq" );
        sa.setAccession( new Accession( "P1", "UniProt" ) );
        sa.setMolecularSequence( "MKT" );
        final DomainArchitecture da = new DomainArchitecture();
        da.addDomain( new ProteinDomain( "kinase", 1, 10, 0.001 ) );
        sa.setDomainArchitecture( da );
        a.getNodeData().addSequence( sa );
        prop( a, "p", "1" );
        prop( a, "q", "2" );
        prop( a, "q", "3" ); // the SAME ref twice on one node counts that node once
        final PhylogenyNode b = tip( "B", 0.3 );
        final Taxonomy tb = new Taxonomy();
        tb.setScientificName( "Mus musculus" );
        b.getNodeData().addTaxonomy( tb );
        final Sequence sb = new Sequence();
        sb.setName( "B_seq" );
        b.getNodeData().addSequence( sb );
        prop( b, "p", "4" );
        final PhylogenyNode c = tip( "C", 0 );
        c.getNodeData().setDate( new Date( "c", new BigDecimal( "0" ), null, null, "mya" ) );
        final List<Point> pts = new ArrayList<>();
        pts.add( new Point( "WGS84", new BigDecimal( "1" ), new BigDecimal( "2" ), null, null ) );
        c.getNodeData().setDistribution( new Distribution( "somewhere", pts ) );
        x.addAsChild( a );
        x.addAsChild( b );
        x.addAsChild( c );
        final PhylogenyNode d = new PhylogenyNode();
        d.setName( "D" ); // NO branch length
        final PhylogenyNode e = tip( "E", 0.5 );
        prop( e, "q", "5" );
        y.addAsChild( d );
        y.addAsChild( e );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode tip( final String name, final double bl ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( bl );
        return n;
    }

    private static void prop( final PhylogenyNode n, final String ref, final String value ) {
        PropertiesList pl = n.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            n.getNodeData().setProperties( pl );
        }
        pl.addProperty( new Property( "t:" + ref, value, "", "xsd:string", AppliesTo.NODE ) );
    }

    private static Group group( final List<Group> groups, final String title ) {
        for( final Group g : groups ) {
            if ( g.title.equals( title ) ) {
                return g;
            }
        }
        return null;
    }

    private static boolean is( final Group g, final String key, final String expected ) {
        final String v = ( g == null ) ? null : g.value( key );
        if ( !expected.equals( v ) ) {
            System.out.println( "  [TreeFactsTest] " + ( ( g == null ) ? "?" : g.title ) + " / " + key + ": expected \""
                    + expected + "\", got \"" + v + "\"" );
            return false;
        }
        return true;
    }

    // ---------------------------------------------------------------- tests
    private static boolean structure() {
        final List<Group> gs = TreeFacts.compute( fixture(), null, false, null );
        // group order: File, Structure, Branch lengths, Support (bootstrap), Support (probability), Coverage
        final StringBuilder titles = new StringBuilder();
        for( final Group g : gs ) {
            titles.append( g.title ).append( '|' );
        }
        if ( !titles.toString().equals( "File|Structure|Branch lengths|Support values (bootstrap)|"
                + "Support values (probability)|Annotation coverage|" ) ) {
            return TestFail.here( titles.toString() );
        }
        final Group s = group( gs, TreeFacts.STRUCTURE );
        boolean ok = true;
        ok &= is( s, "Tips", "5" );
        ok &= is( s, "Internal nodes", "3" );
        ok &= is( s, "Total nodes", "8" );
        ok &= is( s, "Branches", "7" );
        ok &= is( s, "Rooted", "yes" );
        ok &= is( s, "Rerootable", "yes" );
        ok &= is( s, "Branching", "1 polytomy" );
        ok &= is( s, "Depth", "2 (root to deepest tip)" );
        ok &= is( s, "Height", "0.9 (longest root-to-tip path)" );
        if ( !"5 tips".equals( s.detail ) ) {
            return TestFail.here( s.detail );
        }
        if ( s.value( "Collapsed clades" ) != null ) {
            return TestFail.here( "no collapsed clades" );
        }
        // a collapsed clade shows up
        final Phylogeny phy = fixture();
        phy.getRoot().getChildNode( 0 ).setCollapse( true );
        if ( !is( group( TreeFacts.compute( phy, null, false, null ), TreeFacts.STRUCTURE ), "Collapsed clades", "1" ) ) {
            return false;
        }
        return ok || TestFail.here();
    }

    private static boolean branchLengths() {
        final Group b = group( TreeFacts.compute( fixture(), null, false, null ), TreeFacts.BRANCH_LENGTHS );
        boolean ok = true;
        // 6 of 7 branches have a length (D has none): 0.1, 0.4, 0.2, 0.3, 0, 0.5
        ok &= is( b, "Branches with lengths", "6 of 7 (86%)" );
        ok &= is( b, "Median", "0.25" );
        ok &= is( b, "Minimum", "0" );
        ok &= is( b, "Maximum", "0.5" );
        ok &= is( b, "Sum", "1.5 (total tree length)" );
        ok &= is( b, "Zero-length", "1" );
        ok &= is( b, "Ultrametric", "no" );
        if ( b.value( "Negative" ) != null ) {
            return TestFail.here( "no negative lengths in the fixture" );
        }
        if ( !b.value( "Mean" ).startsWith( "0.25 ± " ) || !b.value( "Mean" ).endsWith( " (sd)" ) ) {
            return TestFail.here( b.value( "Mean" ) );
        }
        if ( !"median 0.25".equals( b.detail ) ) {
            return TestFail.here( b.detail );
        }
        if ( ( b.histogram == null ) || ( b.histogram.total() != 6 ) || ( b.histogram.min != 0 )
                || ( b.histogram.max != 0.5 ) ) {
            return TestFail.here( String.valueOf( b.histogram ) );
        }
        // a negative length is counted separately and noted; an ultrametric tree says so
        final Phylogeny neg = fixture();
        neg.getRoot().getChildNode( 1 ).setDistanceToParent( -0.2 );
        final Group bn = group( TreeFacts.compute( neg, null, false, null ), TreeFacts.BRANCH_LENGTHS );
        ok &= is( bn, "Negative", "1 (not drawn to scale)" );
        ok &= is( bn, "Branches with lengths", "5 of 7 (71%)" );
        final Phylogeny ultra = ultrametric();
        final Group bu = group( TreeFacts.compute( ultra, null, false, null ), TreeFacts.BRANCH_LENGTHS );
        ok &= is( bu, "Ultrametric", "yes (all tips equidistant from the root)" );
        // a tree without any branch lengths has no branch-length group (and no Height)
        final Phylogeny none = fixture();
        for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = none.iteratorPreorder(); it.hasNext(); ) {
            it.next().setDistanceToParent( org.forester.phylogeny.data.PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
        }
        final List<Group> gn = TreeFacts.compute( none, null, false, null );
        if ( ( group( gn, TreeFacts.BRANCH_LENGTHS ) != null ) || ( group( gn, TreeFacts.STRUCTURE ).value( "Height" ) != null ) ) {
            return TestFail.here();
        }
        return ok || TestFail.here();
    }

    private static Phylogeny ultrametric() {
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode x = new PhylogenyNode();
        x.setDistanceToParent( 1 );
        x.addAsChild( tip( "A", 2 ) );
        x.addAsChild( tip( "B", 2 ) );
        root.addAsChild( x );
        root.addAsChild( tip( "C", 3 ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean support() {
        final List<Group> gs = TreeFacts.compute( fixture(), null, false, null );
        final Group boot = group( gs, "Support values (bootstrap)" );
        final Group prob = group( gs, "Support values (probability)" );
        boolean ok = true;
        ok &= is( boot, "Type", "bootstrap" );
        ok &= is( boot, "Branches with support", "2" );
        ok &= is( boot, "Median", "80" );
        ok &= is( boot, "Minimum", "70" );
        ok &= is( boot, "Maximum", "90" );
        ok &= is( prob, "Branches with support", "1" );
        ok &= is( prob, "Median", "0.9" );
        if ( !prob.value( "Mean" ).equals( "0.9" ) ) { // a single value: no sd
            return TestFail.here( prob.value( "Mean" ) );
        }
        if ( ( boot.histogram != null ) || ( prob.histogram != null ) ) {
            return TestFail.here( "fewer than 4 values: no histogram" );
        }
        // untyped support only: the plain title, no Type line; a negative value is "unset" and skipped
        final Phylogeny phy = fixture();
        for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            it.next().getBranchData().getConfidences().clear();
        }
        phy.getRoot().getChildNode( 0 ).getBranchData().addConfidence( new Confidence( 55, "" ) );
        phy.getRoot().getChildNode( 1 ).getBranchData().addConfidence( new Confidence( -1, "" ) );
        final Group untyped = group( TreeFacts.compute( phy, null, false, null ), TreeFacts.SUPPORT );
        if ( ( untyped == null ) || ( untyped.value( "Type" ) != null ) ) {
            return TestFail.here();
        }
        ok &= is( untyped, "Branches with support", "1" );
        // no support at all: no support group
        final Phylogeny bare = ultrametric();
        for( final Group g : TreeFacts.compute( bare, null, false, null ) ) {
            if ( g.title.startsWith( TreeFacts.SUPPORT ) ) {
                return TestFail.here( g.title );
            }
        }
        return ok || TestFail.here();
    }

    private static boolean coverage() {
        final Group c = group( TreeFacts.compute( fixture(), null, false, null ), TreeFacts.COVERAGE );
        boolean ok = true;
        ok &= is( c, "Taxonomy", "2 of 5 (40%) tips" );
        ok &= is( c, "Taxonomy identifiers", "1 of 5 (20%) tips" );
        ok &= is( c, "Distinct taxonomies", "2" );
        ok &= is( c, "Sequences", "2 of 5 (40%) tips" );
        ok &= is( c, "Molecular sequences", "1 of 5 (20%) tips" );
        ok &= is( c, "Domain architectures", "1 of 5 (20%) tips" );
        ok &= is( c, "Dates", "1 of 5 (20%) tips, 1 internal nodes" );
        ok &= is( c, "Distributions", "1 of 5 (20%) tips" );
        ok &= is( c, "Named internal nodes", "2 of 3 (67%)" ); // root + X named, Y not
        ok &= is( c, "Events", "2 duplications, 1 speciation, 1 gene loss" );
        ok &= is( c, "Property t:p", "2 nodes" );
        ok &= is( c, "Property t:q", "2 nodes" ); // A (twice, once) + E
        if ( c.value( "Properties" ) != null ) {
            return TestFail.here( "'Properties: none' only without properties" );
        }
        if ( !"6 kinds".equals( c.detail ) ) {
            return TestFail.here( c.detail );
        }
        // a bare tree: zeros, 'none', 'no annotations'
        final Group bare = group( TreeFacts.compute( ultrametric(), null, false, null ), TreeFacts.COVERAGE );
        ok &= is( bare, "Taxonomy", "0 of 3 (0%) tips" );
        ok &= is( bare, "Sequences", "0 of 3 (0%) tips" );
        ok &= is( bare, "Properties", "none" );
        if ( ( bare.value( "Taxonomy identifiers" ) != null ) || ( bare.value( "Events" ) != null )
                || !"no annotations".equals( bare.detail ) ) {
            return TestFail.here( bare.detail );
        }
        return ok || TestFail.here();
    }

    private static boolean timeAxis() {
        // no panel info: no group
        if ( group( TreeFacts.compute( fixture(), null, false, null ), TreeFacts.TIME_AXIS ) != null ) {
            return TestFail.here();
        }
        // axis off and no dates anywhere: no group
        final TreeFacts.TimeAxis off = new TreeFacts.TimeAxis( Options.TIME_AXIS_TYPE.NONE, false, null, 0, 0 );
        if ( group( TreeFacts.compute( ultrametric(), null, false, off ), TreeFacts.TIME_AXIS ) != null ) {
            return TestFail.here();
        }
        // axis off but the tree carries dates: the group says so
        final Group g0 = group( TreeFacts.compute( fixture(), null, false, off ), TreeFacts.TIME_AXIS );
        boolean ok = true;
        ok &= is( g0, "Axis", "off" );
        ok &= is( g0, "Dated nodes", "2" );
        // geologic with a root age and a unit
        final TreeFacts.TimeAxis geo = new TreeFacts.TimeAxis( Options.TIME_AXIS_TYPE.GEOLOGIC, true, "mya", 50, 0 );
        final Group gg = group( TreeFacts.compute( fixture(), null, false, geo ), TreeFacts.TIME_AXIS );
        ok &= is( gg, "Axis", "Geologic (ICS)" );
        ok &= is( gg, "Dated nodes", "2 (a dated tree)" );
        ok &= is( gg, "Unit", "mya" );
        ok &= is( gg, "Root age", "50 Ma" );
        if ( !"Geologic (ICS)".equals( gg.detail ) || ( gg.value( "Most recent date" ) != null ) ) {
            return TestFail.here( gg.detail );
        }
        // calendar: the present date, not a root age
        final TreeFacts.TimeAxis cal = new TreeFacts.TimeAxis( Options.TIME_AXIS_TYPE.CALENDAR, true, "year", 0,
                                                              2024.5 );
        final Group gc = group( TreeFacts.compute( fixture(), null, false, cal ), TreeFacts.TIME_AXIS );
        ok &= is( gc, "Most recent date", "2024.5" );
        if ( gc.value( "Root age" ) != null ) {
            return TestFail.here();
        }
        return ok || TestFail.here();
    }

    private static boolean histogram() {
        final DescriptiveStatistics st = new BasicDescriptiveStatistics();
        for( final double v : new double[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } ) {
            st.addValue( v );
        }
        final Histogram h = TreeFacts.histogram( st );
        if ( ( h == null ) || ( h.counts.length != TreeFacts.HISTOGRAM_BINS ) || ( h.min != 0 ) || ( h.max != 12 ) ) {
            return TestFail.here();
        }
        if ( h.total() != 13 ) {
            return TestFail.here( "every value lands in a bin: " + h.total() );
        }
        // 12 bins of width 1 over [0,12]: 0..10 -> bins 0..10 one each; 11 -> bin 11; 12 (the max) -> the LAST bin
        for( int i = 0; i < 11; ++i ) {
            if ( h.counts[ i ] != 1 ) {
                return TestFail.here( "bin " + i + " = " + h.counts[ i ] );
            }
        }
        if ( ( h.counts[ 11 ] != 2 ) || ( h.maxCount() != 2 ) ) {
            return TestFail.here( "last bin closed: " + h.counts[ 11 ] );
        }
        if ( ( h.edge( 0 ) != 0 ) || ( h.edge( 12 ) != 12 ) || ( Math.abs( h.edge( 6 ) - 6 ) > 1e-9 ) ) {
            return TestFail.here();
        }
        // too few values, or all (nearly) the same: no histogram
        final DescriptiveStatistics few = new BasicDescriptiveStatistics();
        few.addValue( 1 );
        few.addValue( 2 );
        few.addValue( 3 );
        if ( TreeFacts.histogram( few ) != null ) {
            return TestFail.here();
        }
        final DescriptiveStatistics same = new BasicDescriptiveStatistics();
        for( int i = 0; i < 10; ++i ) {
            same.addValue( 0.5 + ( i * 0.000001 ) );
        }
        if ( ( TreeFacts.histogram( same ) != null ) || ( TreeFacts.histogram( null ) != null ) ) {
            return TestFail.here();
        }
        return true;
    }

    private static boolean fileFacts() throws Exception {
        // no file
        final Group none = TreeFacts.fileGroup( null, true );
        boolean ok = true;
        ok &= is( none, "Location", "not saved to a file yet" );
        ok &= is( none, "Unsaved changes", "yes" );
        if ( !none.detail.isEmpty() || ( none.value( "Path" ) != null ) ) {
            return TestFail.here();
        }
        // a real file: path, sniffed format, size, modified, detail = the file name
        final File xml = File.createTempFile( "treefacts", ".xml" );
        xml.deleteOnExit();
        Files.write( xml.toPath(), "\n<?xml version=\"1.0\"?>\n<phyloxml/>\n".getBytes( StandardCharsets.UTF_8 ) );
        final Group g = TreeFacts.fileGroup( xml, false );
        ok &= is( g, "Path", xml.getCanonicalPath() );
        ok &= is( g, "Format", "phyloXML" );
        ok &= is( g, "Size", xml.length() + " B" );
        ok &= is( g, "Unsaved changes", "no" );
        if ( !xml.getName().equals( g.detail ) || ( g.value( "Modified" ) == null )
                || !g.value( "Modified" ).matches( "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}" ) ) {
            return TestFail.here( g.detail + " / " + g.value( "Modified" ) );
        }
        // the format sniff, by content first, then by suffix
        ok &= sniff( ".nex", "#NEXUS\nBEGIN TREES;", "Nexus" );
        ok &= sniff( ".txt", "  begin taxa;", "Nexus" );
        ok &= sniff( ".tre", "((A:1,B:2):0.5,C:3);", "Newick" );
        ok &= sniff( ".txt", "((A:1,B:2)[&&NHX:S=x]:0.5,C:3);", "NHX" );
        ok &= sniff( ".json", "{\"tree\": {}}", "JSON" );
        ok &= sniff( ".nwk", "", "empty file" );
        ok &= sniff( ".nhx", "", "empty file" );
        ok &= sniff( ".weird", "garbage", "unknown" );
        ok &= sniff( ".nexus", "garbage", "Nexus" ); // unrecognizable content: the suffix decides
        ok &= sniff( ".phyloxml", "garbage", "phyloXML" );
        ok &= sniff( ".nh", "garbage", "Newick" );
        if ( !"unknown".equals( TreeFacts.fileFormatLabel( null ) )
                || !"unknown".equals( TreeFacts.fileFormatLabel( new File( "/no/such/file.xml" ) ) ) ) {
            return TestFail.here();
        }
        final File zip = File.createTempFile( "treefacts", ".xml.zip" );
        zip.deleteOnExit();
        Files.write( zip.toPath(), new byte[] { 0x50, 0x4b } );
        if ( !"zipped (phyloXML)".equals( TreeFacts.fileFormatLabel( zip ) ) ) {
            return TestFail.here( TreeFacts.fileFormatLabel( zip ) );
        }
        return ok || TestFail.here();
    }

    private static boolean sniff( final String suffix, final String content, final String expected ) throws Exception {
        final File f = File.createTempFile( "treefacts", suffix );
        f.deleteOnExit();
        Files.write( f.toPath(), content.getBytes( StandardCharsets.UTF_8 ) );
        final String got = TreeFacts.fileFormatLabel( f );
        if ( !expected.equals( got ) ) {
            System.out.println( "  [TreeFactsTest] sniff " + suffix + " \"" + content.replace( "\n", "\\n" )
                    + "\": expected " + expected + ", got " + got );
            return false;
        }
        return true;
    }

    private static boolean helpers() {
        if ( !"12 of 40 (30%)".equals( TreeFacts.ofText( 12, 40 ) ) || !"0 of 3 (0%)".equals( TreeFacts.ofText( 0, 3 ) )
                || !"3 of 3 (100%)".equals( TreeFacts.ofText( 3, 3 ) ) || !"1,234 of 5,000 (25%)".equals( TreeFacts.ofText( 1234, 5000 ) )
                || !"7".equals( TreeFacts.ofText( 7, 0 ) ) ) {
            return TestFail.here( TreeFacts.ofText( 12, 40 ) );
        }
        if ( !"512 B".equals( TreeFacts.humanSize( 512 ) ) || !"1 KB".equals( TreeFacts.humanSize( 1024 ) )
                || !"2.5 KB".equals( TreeFacts.humanSize( 2560 ) ) || !"3.2 MB".equals( TreeFacts.humanSize( 3355443 ) )
                || !"1 GB".equals( TreeFacts.humanSize( 1L << 30 ) ) ) {
            return TestFail.here( TreeFacts.humanSize( 3355443 ) );
        }
        if ( !"12,345".equals( TreeFacts.count( 12345 ) ) || !"0.123457".equals( TreeFacts.number( 0.1234567 ) )
                || !"5".equals( TreeFacts.number( 5.0 ) ) || !"0.02".equals( TreeFacts.number( 0.02 ) ) ) {
            return TestFail.here( TreeFacts.number( 0.1234567 ) );
        }
        final Group g = new Group( "T", null, new ArrayList<TreeFacts.Fact>(), null );
        if ( !"".equals( g.detail ) || ( g.value( "x" ) != null ) ) {
            return TestFail.here();
        }
        return true;
    }

    private static boolean edgeCases() {
        // a null or empty tree: only the file group
        final List<Group> gs = TreeFacts.compute( null, null, false, null );
        if ( ( gs.size() != 1 ) || !TreeFacts.FILE.equals( gs.get( 0 ).title ) ) {
            return TestFail.here();
        }
        if ( TreeFacts.compute( new Phylogeny(), null, false, null ).size() != 1 ) {
            return TestFail.here();
        }
        // a single node
        final Phylogeny one = new Phylogeny();
        final PhylogenyNode only = new PhylogenyNode();
        only.setName( "only" );
        one.setRoot( only );
        one.setRooted( true );
        final List<Group> g1 = TreeFacts.compute( one, null, false, null );
        final Group s1 = group( g1, TreeFacts.STRUCTURE );
        if ( !is( s1, "Tips", "1" ) || !is( s1, "Branches", "0" ) || !"1 tip".equals( s1.detail ) ) {
            return TestFail.here();
        }
        if ( group( g1, TreeFacts.BRANCH_LENGTHS ) != null ) {
            return TestFail.here();
        }
        final Group c1 = group( g1, TreeFacts.COVERAGE );
        if ( ( c1.value( "Named internal nodes" ) != null ) || !is( c1, "Taxonomy", "0 of 1 (0%) tips" ) ) {
            return TestFail.here();
        }
        // an unrooted tree
        final Phylogeny un = fixture();
        un.setRooted( false );
        if ( !is( group( TreeFacts.compute( un, null, false, null ), TreeFacts.STRUCTURE ), "Rooted", "no" ) ) {
            return false;
        }
        return true;
    }

    private TreeFactsTest() {
    }
}

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

import java.math.BigDecimal;
import java.net.URI;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import org.forester.archaeopteryx.NodeDataDraft.ConfidenceDraft;
import org.forester.archaeopteryx.NodeDataDraft.Problem;
import org.forester.archaeopteryx.NodeDataDraft.PropertyDraft;
import org.forester.archaeopteryx.NodeDataDraft.SequenceDraft;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Accession;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.data.BranchWidth;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Distribution;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Point;
import org.forester.phylogeny.data.Polygon;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Reference;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.data.Uri;

/**
 * Headless tests for {@link NodeDataDraft}: reading a node into a draft (without touching the node), validation of
 * every field rule, and the DIFF write-back -- only changed fields are applied, untouched precision and untouched
 * sub-elements survive, and sections emptied by the user are removed from the node.
 */
public final class NodeDataDraftTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeDataDraft: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return readRichNode() && openDoesNotMutate() && noChangeIsNoOp() && diffWriteKeepsPrecision()
                    && basicsWrite() && taxonomyWrite() && sequencesWrite() && eventsWrite() && dateWrite()
                    && distributionWrite() && referenceWrite() && propertiesWrite() && validation()
                    && invalidWriteRefused() && equalityAndProvenance() && helpers();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // ---------------------------------------------------------------- fixtures
    /** An internal node carrying one of everything the editor knows about. */
    static PhylogenyNode richNode() throws Exception {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "BRCA1 clade" );
        n.setDistanceToParent( 0.0123 );
        n.getBranchData().setBranchWidth( new BranchWidth( 2 ) );
        n.getBranchData().addConfidence( new Confidence( 95, "bootstrap" ) );
        n.getBranchData().addConfidence( new Confidence( 0.98, "probability", 0.01 ) );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( "Homo sapiens" );
        t.setTaxonomyCode( "HUMAN" );
        t.setRank( "species" );
        t.setIdentifier( new Identifier( "9606", "ncbi" ) );
        t.setCommonName( "human" );
        t.setAuthority( "Linnaeus, 1758" );
        t.getSynonyms().add( "man" );
        t.getSynonyms().add( "" ); // an empty synonym must not surface as a blank line
        t.addUri( new Uri( new URI( "https://example.org/tax/9606" ), "taxonomy page", "text/html" ) );
        n.getNodeData().addTaxonomy( t );
        final Sequence s1 = new Sequence();
        s1.setName( "BRCA1_HUMAN" );
        s1.setSymbol( "BRCA1" );
        s1.setGeneName( "BRCA1" );
        s1.setType( "protein" );
        s1.setAccession( new Accession( "P38398", "UniProt" ) );
        s1.setLocation( "chr17" );
        s1.setMolecularSequence( "MDLSALRV" );
        s1.setMolecularSequenceAligned( true );
        s1.addUri( new Uri( new URI( "https://example.org/seq/P38398" ) ) );
        s1.addAnnotation( new Annotation( "GO:0005634" ) ); // an EXTRA the form does not edit
        n.getNodeData().addSequence( s1 );
        final Sequence s2 = new Sequence();
        s2.setName( "BRCA1 mRNA" );
        s2.setType( "rna" );
        n.getNodeData().addSequence( s2 );
        n.getNodeData().setEvent( new Event( 1, 0, 2 ) );
        n.getNodeData().setDate( new Date( "split", new BigDecimal( "6.5" ), new BigDecimal( "5.0" ),
                                           new BigDecimal( "8.0" ), "mya" ) );
        final List<Point> pts = new ArrayList<>();
        pts.add( new Point( "WGS84", new BigDecimal( "37.77" ), new BigDecimal( "-122.42" ), new BigDecimal( "16" ),
                            "m" ) );
        n.getNodeData().setDistribution( new Distribution( "San Francisco", pts ) );
        n.getNodeData().setReference( new Reference( "Miki et al. 1994", "10.1126/science.7545954" ) );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "data:depth", "120", "METRIC:m", "xsd:decimal", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:habitat", "coastal", "", "xsd:string", AppliesTo.CLADE, "ref1" ) );
        n.getNodeData().setProperties( pl );
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "tipA" );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "tipB" );
        n.addAsChild( a );
        n.addAsChild( b );
        return n;
    }

    // ---------------------------------------------------------------- tests
    private static boolean readRichNode() throws Exception {
        final NodeDataDraft d = NodeDataDraft.from( richNode() );
        return eq( "name", "BRCA1 clade", d.name ) && eq( "branch length", "0.0123", d.branchLength )
                && eq( "width", "2", d.branchWidth ) && eq( "confidences", 2, d.confidences.size() )
                && eq( "conf 0", new ConfidenceDraft( "95", "bootstrap", "" ), d.confidences.get( 0 ) )
                && eq( "conf 1", new ConfidenceDraft( "0.98", "probability", "0.01" ), d.confidences.get( 1 ) )
                && eq( "sci name", "Homo sapiens", d.taxSciName ) && eq( "code", "HUMAN", d.taxCode )
                && eq( "rank", "species", d.taxRank ) && eq( "id", "9606", d.taxId )
                && eq( "provider", "ncbi", d.taxProvider ) && eq( "common", "human", d.taxCommonName )
                && eq( "authority", "Linnaeus, 1758", d.taxAuthority ) && eq( "synonyms", "man", d.taxSynonyms )
                && eq( "tax uris", "https://example.org/tax/9606", d.taxUris )
                && eq( "sequences", 2, d.sequences.size() ) && eq( "seq name", "BRCA1_HUMAN", d.sequences.get( 0 ).name )
                && eq( "symbol", "BRCA1", d.sequences.get( 0 ).symbol )
                && eq( "gene", "BRCA1", d.sequences.get( 0 ).geneName )
                && eq( "type", "protein", d.sequences.get( 0 ).type )
                && eq( "acc", "P38398", d.sequences.get( 0 ).accession )
                && eq( "source", "UniProt", d.sequences.get( 0 ).source )
                && eq( "location", "chr17", d.sequences.get( 0 ).location )
                && eq( "mol seq", "MDLSALRV", d.sequences.get( 0 ).molSeq )
                && eq( "aligned", true, d.sequences.get( 0 ).aligned )
                && eq( "seq uris", "https://example.org/seq/P38398", d.sequences.get( 0 ).uris )
                && eq( "seq 2 type", "rna", d.sequences.get( 1 ).type )
                && eq( "dups", "1", d.duplications ) && eq( "specs", "0", d.speciations )
                && eq( "losses", "2", d.geneLosses ) && eq( "date desc", "split", d.dateDesc )
                && eq( "date value", "6.5", d.dateValue ) && eq( "date min", "5.0", d.dateMin )
                && eq( "date max", "8.0", d.dateMax ) && eq( "date unit", "mya", d.dateUnit )
                && eq( "dist desc", "San Francisco", d.distDesc ) && eq( "datum", "WGS84", d.distDatum )
                && eq( "lat", "37.77", d.distLat ) && eq( "long", "-122.42", d.distLong )
                && eq( "alt", "16", d.distAlt ) && eq( "alt unit", "m", d.distAltUnit )
                && eq( "ref desc", "Miki et al. 1994", d.refDesc ) && eq( "doi", "10.1126/science.7545954", d.refDoi )
                && eq( "properties", 2, d.properties.size() ) && eq( "prop ref", "data:depth", d.properties.get( 0 ).ref )
                && eq( "prop unit", "METRIC:m", d.properties.get( 0 ).unit )
                && eq( "prop applies", AppliesTo.CLADE, d.properties.get( 1 ).appliesTo )
                && eq( "prop id ref", "ref1", d.properties.get( 1 ).idRef )
                && check( "no problems on a rich node", d.validate( true ).isEmpty() );
    }

    /** The old editor added empty Taxonomy/Sequence/Distribution/Reference objects to the node just to show rows. */
    private static boolean openDoesNotMutate() {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "bare" );
        final NodeDataDraft d = NodeDataDraft.from( n );
        return check( "bare draft is empty", !d.hasTaxonomy() && d.sequences.isEmpty() && !d.hasDate()
                && !d.hasDistribution() && !d.hasReference() && d.properties.isEmpty() && d.confidences.isEmpty() )
                && check( "reading must not add a taxonomy", !n.getNodeData().isHasTaxonomy() )
                && check( "reading must not add a sequence", !n.getNodeData().isHasSequence() )
                && check( "reading must not add a distribution", !n.getNodeData().isHasDistribution() )
                && check( "reading must not add a reference", !n.getNodeData().isHasReference() )
                && check( "reading must not add a date", !n.getNodeData().isHasDate() )
                && check( "reading must not add properties", !n.getNodeData().isHasProperties() )
                && check( "reading must not add an event", !n.getNodeData().isHasEvent() )
                && eq( "a root's branch length shows empty", "", d.branchLength );
    }

    private static boolean noChangeIsNoOp() throws Exception {
        final PhylogenyNode n = richNode();
        final NodeDataDraft base = NodeDataDraft.from( n );
        final Sequence s1 = n.getNodeData().getSequence( 0 );
        final Set<String> changed = base.copy().writeTo( n, base );
        return check( "no field changed -> nothing written", changed.isEmpty() )
                && check( "no-op write leaves the sequence object alone", n.getNodeData().getSequence( 0 ) == s1 )
                && eq( "re-read equals the baseline", base, NodeDataDraft.from( n ) );
    }

    /** Only CHANGED fields are written, so a branch length the user never touched keeps every digit. */
    private static boolean diffWriteKeepsPrecision() {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "x" );
        n.setDistanceToParent( 0.123456789012345 );
        final PhylogenyNode parent = new PhylogenyNode();
        parent.addAsChild( n );
        final NodeDataDraft base = NodeDataDraft.from( n );
        final NodeDataDraft d = base.copy();
        d.name = "renamed";
        final Set<String> changed = d.writeTo( n, base );
        return eq( "9 decimals shown", "0.123456789", base.branchLength )
                && eq( "only Basic changed", "[Basic]", changed.toString() )
                && eq( "name written", "renamed", n.getName() )
                && check( "untouched branch length keeps full precision",
                          n.getDistanceToParent() == 0.123456789012345 );
    }

    private static boolean basicsWrite() throws Exception {
        final PhylogenyNode n = richNode();
        final NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        d.branchLength = " 0.5 ";
        d.branchWidth = "1"; // the default -> the element goes away
        d.confidences.get( 0 ).value = "88";
        d.confidences.get( 0 ).type = ""; // -> "unknown"
        d.confidences.add( new ConfidenceDraft( "0.5", "posterior", "0.1" ) );
        d.confidences.add( new ConfidenceDraft() ); // an untouched "+ Add" row is ignored
        d.writeTo( n, base );
        boolean ok = check( "branch length", n.getDistanceToParent() == 0.5 )
                && check( "default width removes the element", n.getBranchData().getBranchWidth() == null )
                && eq( "confidences", 3, n.getBranchData().getConfidences().size() )
                && check( "conf 0 value", n.getBranchData().getConfidence( 0 ).getValue() == 88 )
                && eq( "empty type -> unknown", "unknown", n.getBranchData().getConfidence( 0 ).getType() )
                && check( "sd written", n.getBranchData().getConfidence( 2 ).getStandardDeviation() == 0.1 );
        // clearing the branch length restores the "no branch length" default
        final NodeDataDraft base2 = NodeDataDraft.from( n );
        d = base2.copy();
        d.branchLength = "";
        d.branchWidth = "3.5";
        d.confidences.clear();
        d.writeTo( n, base2 );
        return ok && check( "empty -> default branch length",
                            n.getDistanceToParent() == org.forester.phylogeny.data.PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT )
                && check( "width 3.5", n.getBranchData().getBranchWidth().getValue() == 3.5 )
                && check( "no confidences left", !n.getBranchData().isHasConfidences() );
    }

    private static boolean taxonomyWrite() throws Exception {
        // create on a node without one
        final PhylogenyNode n = new PhylogenyNode();
        NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        d.taxSciName = "Mus musculus";
        d.taxRank = "Species"; // any case in, phyloXML case out
        d.taxCode = "MOUSE";
        d.taxId = "10090";
        d.taxSynonyms = "mouse\n\nhouse mouse\n";
        d.taxUris = "https://example.org/a\nhttps://example.org/b";
        boolean ok = eq( "section", "[Taxonomy]", d.writeTo( n, base ).toString() )
                && check( "taxonomy created", n.getNodeData().isHasTaxonomy() )
                && eq( "sci name", "Mus musculus", n.getNodeData().getTaxonomy().getScientificName() )
                && eq( "rank lower-cased", "species", n.getNodeData().getTaxonomy().getRank() )
                && eq( "id value", "10090", n.getNodeData().getTaxonomy().getIdentifier().getValue() )
                && eq( "synonyms (blank lines dropped)", Arrays.asList( "mouse", "house mouse" ),
                       n.getNodeData().getTaxonomy().getSynonyms() )
                && eq( "two uris", 2, n.getNodeData().getTaxonomy().getUris().size() );
        // update in place keeps the extras the form does not edit (here: a URI's description, the lineage)
        final PhylogenyNode r = richNode();
        final Taxonomy t = r.getNodeData().getTaxonomy();
        t.setLineage( Arrays.asList( "Eukaryota", "Metazoa" ) );
        final Uri kept = t.getUris().get( 0 );
        base = NodeDataDraft.from( r );
        d = base.copy();
        d.taxCommonName = "person";
        d.writeTo( r, base );
        ok = ok && check( "same taxonomy object", r.getNodeData().getTaxonomy() == t )
                && eq( "common name", "person", t.getCommonName() )
                && check( "URI object reused (its description survives)", t.getUris().get( 0 ) == kept )
                && eq( "lineage kept", 2, t.getLineage().size() );
        // clearing every field removes the element -- unless a lineage still hangs on it
        base = NodeDataDraft.from( r );
        d = base.copy();
        d.taxId = d.taxProvider = d.taxCode = d.taxSciName = d.taxAuthority = d.taxCommonName = "";
        d.taxSynonyms = d.taxRank = d.taxUris = "";
        d.writeTo( r, base );
        ok = ok && check( "lineage keeps the (emptied) taxonomy", r.getNodeData().isHasTaxonomy() )
                && eq( "fields cleared", "", NodeDataDraft.nn( t.getScientificName() ) );
        t.setLineage( null );
        base = NodeDataDraft.from( r );
        d = base.copy();
        d.taxRank = "genus";
        d.writeTo( r, base );
        base = NodeDataDraft.from( r );
        d = base.copy();
        d.taxRank = "";
        d.writeTo( r, base );
        return ok && check( "all empty and no lineage -> taxonomy removed", !r.getNodeData().isHasTaxonomy() );
    }

    private static boolean sequencesWrite() throws Exception {
        final PhylogenyNode n = richNode();
        final Sequence s1 = n.getNodeData().getSequence( 0 );
        final Sequence s2 = n.getNodeData().getSequence( 1 );
        NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        // edit the first in place, add a third, leave a blank "+ Add" card
        d.sequences.get( 0 ).symbol = "BRCA1x";
        d.sequences.get( 0 ).molSeq = ">hdr\n10 MDLS ALRV\n20 EEVQ*"; // pasted junk gets cleaned
        d.sequences.get( 0 ).aligned = false;
        final SequenceDraft third = new SequenceDraft();
        third.name = "third";
        third.type = "DNA";
        third.accession = "X1";
        d.sequences.add( third );
        d.sequences.add( new SequenceDraft() );
        boolean ok = eq( "section", "[Sequences]", d.writeTo( n, base ).toString() )
                && eq( "three sequences (blank card ignored)", 3, n.getNodeData().getSequences().size() )
                && check( "first mutated IN PLACE", n.getNodeData().getSequence( 0 ) == s1 )
                && eq( "symbol", "BRCA1x", s1.getSymbol() )
                && eq( "mol seq cleaned", "hdrMDLSALRVEEVQ*", s1.getMolecularSequence() )
                && check( "aligned flag", !s1.isMolecularSequenceAligned() )
                && eq( "annotation survives an in-place edit", 1, s1.getAnnotations().size() )
                && check( "second untouched", n.getNodeData().getSequence( 1 ) == s2 )
                && eq( "third created", "third", n.getNodeData().getSequence( 2 ).getName() )
                && eq( "type lower-cased", "dna", n.getNodeData().getSequence( 2 ).getType() )
                && eq( "accession", "X1", n.getNodeData().getSequence( 2 ).getAccession().getValue() )
                && check( "draft now bound to the created object",
                          third.origin == n.getNodeData().getSequence( 2 ) );
        // remove the first; empty out the third (nothing else on it -> gone); the second survives as itself
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.sequences.remove( 0 );
        final SequenceDraft t = d.sequences.get( 1 );
        t.name = t.type = t.accession = "";
        d.writeTo( n, base );
        ok = ok && eq( "one sequence left", 1, n.getNodeData().getSequences().size() )
                && check( "the survivor is the former second object", n.getNodeData().getSequence( 0 ) == s2 );
        // an emptied sequence that still carries an annotation is KEPT (the form cannot show what it would lose)
        s2.addAnnotation( new Annotation( "GO:1" ) );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.sequences.get( 0 ).name = "";
        d.sequences.get( 0 ).type = "";
        d.writeTo( n, base );
        ok = ok && eq( "annotated sequence kept", 1, n.getNodeData().getSequences().size() );
        // removing the last one leaves the node with NO sequence element (not an empty list)
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.sequences.clear();
        d.writeTo( n, base );
        return ok && check( "no sequences -> isHasSequence false", !n.getNodeData().isHasSequence() )
                && check( "getSequences() is null (no empty shell)", n.getNodeData().getSequences() == null );
    }

    private static boolean eventsWrite() throws Exception {
        final PhylogenyNode n = richNode();
        NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        d.speciations = "3";
        boolean ok = eq( "section", "[Events]", d.writeTo( n, base ).toString() )
                && eq( "speciations", 3, n.getNodeData().getEvent().getNumberOfSpeciations() )
                && eq( "duplications kept", 1, n.getNodeData().getEvent().getNumberOfDuplications() );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.duplications = d.speciations = d.geneLosses = "";
        d.writeTo( n, base );
        ok = ok && check( "all cleared on a counts-only event -> removed", !n.getNodeData().isHasEvent() );
        // a TYPED event (here: a transfer) keeps its type when its counts are cleared, instead of vanishing --
        // and instead of being stamped "mixed" by the count setters
        n.getNodeData().setEvent( new Event( 0, 1, 0, "transfer" ) );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.duplications = d.speciations = d.geneLosses = "";
        d.writeTo( n, base );
        ok = ok && check( "typed event kept", n.getNodeData().isHasEvent() )
                && eq( "type kept", Event.EventType.transfer, n.getNodeData().getEvent().getEventType() )
                && eq( "counts gone", "", NodeDataDraft.from( n ).speciations );
        // an external node never gets events, even if a draft carries them
        final PhylogenyNode tip = new PhylogenyNode();
        tip.setName( "tip" );
        base = NodeDataDraft.from( tip );
        d = base.copy();
        d.duplications = "2";
        final Set<String> changed = d.writeTo( tip, base );
        return ok && check( "no events on a tip", !tip.getNodeData().isHasEvent() && changed.isEmpty() )
                && check( "tip draft validates without events", d.validate( false ).isEmpty() );
    }

    private static boolean dateWrite() throws Exception {
        final PhylogenyNode n = new PhylogenyNode();
        NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        d.dateValue = "12.5";
        d.dateUnit = "mya";
        boolean ok = eq( "section", "[Date]", d.writeTo( n, base ).toString() )
                && check( "date created", n.getNodeData().isHasDate() )
                && eq( "value", new BigDecimal( "12.5" ), n.getNodeData().getDate().getValue() )
                && check( "min null", n.getNodeData().getDate().getMin() == null );
        final Date date = n.getNodeData().getDate();
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.dateMin = "10";
        d.writeTo( n, base );
        ok = ok && check( "same date object", n.getNodeData().getDate() == date )
                && eq( "min", new BigDecimal( "10" ), date.getMin() );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.dateValue = d.dateMin = d.dateUnit = "";
        d.writeTo( n, base );
        return ok && check( "all cleared -> removed", !n.getNodeData().isHasDate() );
    }

    private static boolean distributionWrite() throws Exception {
        final PhylogenyNode n = new PhylogenyNode();
        NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        d.distLat = "10";
        d.distLong = "20";
        boolean ok = eq( "section", "[Distribution]", d.writeTo( n, base ).toString() )
                && check( "created", n.getNodeData().isHasDistribution() )
                && eq( "empty datum -> ?", Point.UNKNOWN_GEODETIC_DATUM,
                       n.getNodeData().getDistribution().getPoints().get( 0 ).getGeodeticDatum() )
                && eq( "datum reads back as empty", "", NodeDataDraft.from( n ).distDatum );
        // polygons the form does not edit survive a point edit
        final List<Point> poly = new ArrayList<>();
        poly.add( new Point( "WGS84", BigDecimal.ONE, BigDecimal.ONE ) );
        poly.add( new Point( "WGS84", BigDecimal.TEN, BigDecimal.ONE ) );
        poly.add( new Point( "WGS84", BigDecimal.TEN, BigDecimal.TEN ) );
        final Distribution old = n.getNodeData().getDistribution();
        final List<Polygon> polygons = new ArrayList<>();
        polygons.add( new Polygon( poly ) );
        n.getNodeData().setDistribution( new Distribution( old.getDesc(), old.getPoints(), polygons ) );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.distDesc = "somewhere";
        d.distAlt = "5";
        d.distAltUnit = "m";
        d.writeTo( n, base );
        ok = ok && eq( "desc", "somewhere", n.getNodeData().getDistribution().getDesc() )
                && eq( "alt", new BigDecimal( "5" ), n.getNodeData().getDistribution().getPoints().get( 0 ).getAltitude() )
                && eq( "polygon kept", 1, n.getNodeData().getDistribution().getPolygons().size() );
        // Distribution(desc, points) used to DROP the description -- fixed alongside this editor
        ok = ok && eq( "2-arg Distribution keeps its desc", "d", new Distribution( "d", new ArrayList<>() ).getDesc() );
        // clearing everything (with no polygons) removes the element
        n.getNodeData().setDistribution( new Distribution( "x", old.getPoints() ) );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.distDesc = d.distDatum = d.distLat = d.distLong = d.distAlt = d.distAltUnit = "";
        d.writeTo( n, base );
        return ok && check( "all cleared -> removed", !n.getNodeData().isHasDistribution() );
    }

    private static boolean referenceWrite() throws Exception {
        final PhylogenyNode n = new PhylogenyNode();
        NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        d.refDoi = "10.1000/xyz123";
        boolean ok = eq( "section", "[Reference]", d.writeTo( n, base ).toString() )
                && eq( "doi", "10.1000/xyz123", n.getNodeData().getReference().getDoi() );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.refDoi = "";
        d.writeTo( n, base );
        return ok && check( "cleared -> removed", !n.getNodeData().isHasReference() );
    }

    private static boolean propertiesWrite() throws Exception {
        final PhylogenyNode n = richNode();
        NodeDataDraft base = NodeDataDraft.from( n );
        NodeDataDraft d = base.copy();
        d.properties.get( 0 ).value = "130";
        d.properties.add( new PropertyDraft( "data:year", "2019", "", "xsd:integer", AppliesTo.NODE ) );
        d.properties.add( new PropertyDraft() ); // blank row ignored
        boolean ok = eq( "section", "[Properties]", d.writeTo( n, base ).toString() )
                && eq( "three properties", 3, n.getNodeData().getProperties().size() )
                && eq( "value", "130", n.getNodeData().getProperties().getProperties( "data:depth" ).get( 0 ).getValue() )
                && eq( "id_ref survives the rebuild", "ref1",
                       n.getNodeData().getProperties().getProperties( "data:habitat" ).get( 0 ).getIdRef() )
                && eq( "new one", "2019", n.getNodeData().getProperties().getProperties( "data:year" ).get( 0 ).getValue() );
        base = NodeDataDraft.from( n );
        d = base.copy();
        d.properties.clear();
        d.writeTo( n, base );
        return ok && check( "none left -> removed", !n.getNodeData().isHasProperties() );
    }

    private static boolean validation() throws Exception {
        final NodeDataDraft d = NodeDataDraft.from( richNode() );
        d.branchLength = "abc";
        d.branchWidth = "-1";
        d.confidences.get( 0 ).value = "";
        d.confidences.get( 1 ).sd = "x";
        d.taxCode = "human";
        d.taxRank = "kind";
        d.taxUris = "https://ok.example\nwww.no-scheme.example";
        d.sequences.get( 0 ).symbol = "has space";
        d.sequences.get( 1 ).type = "peptide";
        d.sequences.get( 1 ).uris = "not a url";
        d.duplications = "-2";
        d.speciations = "two";
        d.dateValue = "1e";
        d.distLat = "95";
        d.distLong = "-200";
        d.distAlt = "10";
        d.distAltUnit = "";
        d.refDoi = "no doi";
        d.properties.get( 0 ).value = "deep"; // xsd:decimal
        d.properties.get( 1 ).ref = "noprefix";
        d.properties.add( new PropertyDraft( "a:b", "v", "nounit", "xsd", AppliesTo.NODE ) );
        final List<Problem> ps = d.validate( true );
        final List<String> keys = new ArrayList<>();
        for( final Problem p : ps ) {
            keys.add( p.key );
        }
        final List<String> expected = Arrays.asList( NodeDataDraft.BRANCH_LENGTH,
                                                     NodeDataDraft.BRANCH_WIDTH,
                                                     NodeDataDraft.confidenceKey( 0, NodeDataDraft.CONF_VALUE ),
                                                     NodeDataDraft.confidenceKey( 1, NodeDataDraft.CONF_SD ),
                                                     NodeDataDraft.TAX_CODE,
                                                     NodeDataDraft.TAX_RANK,
                                                     NodeDataDraft.TAX_URIS,
                                                     NodeDataDraft.sequenceKey( 0, NodeDataDraft.SEQ_SYMBOL ),
                                                     NodeDataDraft.sequenceKey( 1, NodeDataDraft.SEQ_TYPE ),
                                                     NodeDataDraft.sequenceKey( 1, NodeDataDraft.SEQ_URIS ),
                                                     NodeDataDraft.EV_DUPLICATIONS,
                                                     NodeDataDraft.EV_SPECIATIONS,
                                                     NodeDataDraft.DATE_VALUE,
                                                     NodeDataDraft.DIST_LAT,
                                                     NodeDataDraft.DIST_LONG,
                                                     NodeDataDraft.DIST_ALT_UNIT,
                                                     NodeDataDraft.REF_DOI,
                                                     NodeDataDraft.propertyKey( 0, NodeDataDraft.PROP_VALUE ),
                                                     NodeDataDraft.propertyKey( 1, NodeDataDraft.PROP_REF ),
                                                     NodeDataDraft.propertyKey( 2, NodeDataDraft.PROP_UNIT ),
                                                     NodeDataDraft.propertyKey( 2, NodeDataDraft.PROP_TYPE ) );
        boolean ok = eq( "every rule fires exactly once, in field order", expected, keys );
        for( final Problem p : ps ) {
            ok = ok && check( "a message for " + p.key, !p.message.isEmpty() );
        }
        // events are NOT validated for an external node
        ok = ok && check( "tip: no event problems", !keysOf( d.validate( false ) ).contains( NodeDataDraft.EV_DUPLICATIONS ) );
        // the things that must PASS: a rank in any case, a symbol at the 20-char limit, blank rows
        final NodeDataDraft good = NodeDataDraft.from( richNode() );
        good.taxRank = "SPECIES";
        good.sequences.get( 0 ).symbol = "12345678901234567890";
        good.confidences.add( new ConfidenceDraft() );
        good.properties.add( new PropertyDraft() );
        good.distAlt = "";
        return ok && eq( "good draft", "[]", good.validate( true ).toString() );
    }

    private static boolean invalidWriteRefused() throws Exception {
        final PhylogenyNode n = richNode();
        final NodeDataDraft base = NodeDataDraft.from( n );
        final NodeDataDraft d = base.copy();
        d.name = "changed";
        d.branchLength = "bad";
        try {
            d.writeTo( n, base );
            return check( "an invalid draft must not write", false );
        }
        catch ( final IllegalStateException e ) {
            return eq( "node untouched", "BRCA1 clade", n.getName() );
        }
    }

    private static boolean equalityAndProvenance() throws Exception {
        final NodeDataDraft a = NodeDataDraft.from( richNode() );
        final NodeDataDraft b = a.copy();
        boolean ok = eq( "copy equals", a, b ) && eq( "hash", a.hashCode(), b.hashCode() );
        b.sequences.get( 1 ).location = "elsewhere";
        ok = ok && check( "a nested change breaks equality", !a.equals( b ) )
                && eq( "changed sections", "[Sequences]", b.changedSections( a, true ).toString() );
        final SequenceDraft x = new SequenceDraft();
        final SequenceDraft y = SequenceDraft.from( new Sequence() );
        ok = ok && eq( "sequence drafts compare by FIELDS, not origin", x, y );
        b.duplications = "9";
        b.name = "n";
        ok = ok && eq( "sections in page order", "[Basic, Sequences, Events]",
                       b.changedSections( a, true ).toString() )
                && eq( "a tip ignores the events diff", "[Basic, Sequences]", b.changedSections( a, false ).toString() );
        return ok && eq( "provenance", "Manually edited the node data (Basic, Sequences) of node \"BRCA1 clade\".",
                         NodeDataDraft.provenance( "BRCA1 clade", Arrays.asList( "Basic", "Sequences" ) ) );
    }

    private static boolean helpers() {
        return eq( "number: no trailing zeros", "0.5", NodeDataDraft.formatNumber( 0.5 ) )
                && eq( "number: leading zero", "0.000001", NodeDataDraft.formatNumber( 0.000001 ) )
                && eq( "number: integer", "95", NodeDataDraft.formatNumber( 95 ) )
                && eq( "number: 9 decimals max", "0.123456789", NodeDataDraft.formatNumber( 0.1234567891234 ) )
                && eq( "lines", Arrays.asList( "a", "b" ), NodeDataDraft.lines( " a \r\n\n b\n" ) )
                && eq( "lines of null", 0, NodeDataDraft.lines( null ).size() )
                && eq( "mol seq: keeps gaps and stops", "AC-G.T*?", NodeDataDraft.sanitizeMolSeq( " 1 AC-G .T*?\n" ) )
                && eq( "nn", "", NodeDataDraft.nn( null ) );
    }

    // ---------------------------------------------------------------- plumbing
    private static List<String> keysOf( final List<Problem> ps ) {
        final List<String> keys = new ArrayList<>();
        for( final Problem p : ps ) {
            keys.add( p.key );
        }
        return keys;
    }

    private static boolean eq( final String what, final Object expected, final Object actual ) {
        if ( ( expected == null ) ? ( actual == null ) : expected.equals( actual ) ) {
            return true;
        }
        System.out.println( "  [NodeDataDraftTest] " + what + ": expected <" + expected + "> but got <" + actual + ">" );
        return false;
    }

    private static boolean check( final String what, final boolean condition ) {
        if ( condition ) {
            return true;
        }
        System.out.println( "  [NodeDataDraftTest] " + what );
        return false;
    }

    private NodeDataDraftTest() {
        // not instantiable
    }
}

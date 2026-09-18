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

package org.forester.phylogeny.data;

import java.io.File;
import java.io.FileWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.Iterator;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The ordering contract of {@link PropertiesList}: a node's properties keep INSERTION order, which for every parser
 * and importer is the order the SOURCE stated them in.
 * <p>
 * This list used to re-sort itself by ref on every add. The assertions here pin the replacement from both ends --
 * that insertion order survives, and that nothing re-sorts -- because the old behaviour was invisible until a tree
 * carried many properties: it alphabetized a gene presence/absence matrix's columns on read and wrote them back out
 * alphabetized, destroying the author's grouping on a plain open-and-save, and it disagreed with Archaeopteryx.js,
 * whose reader and writer both preserve document order.
 * <p>
 * Also pinned deliberately: a property REPLACED by the remove-then-re-add path (the annotation importer, the GTDB
 * and tip-date tools) moves to the END of the list. The old sort silently returned it to its alphabetical slot; now
 * the relocation is real, so it is stated here rather than left to be discovered.
 */
public final class PropertiesListTest {

    private static final String STR = "xsd:string";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "PropertiesListTest: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return testInsertionOrder() && testNothingReSorts() && testDuplicateRefs() && testLookupsKeepOrder()
                && testCopyAndTextKeepOrder() && testPhyloXmlWritesInListOrder() && testDocumentOrderSurvivesParse()
                && testReplacedPropertyMovesToTail() && testEdges();
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [PropertiesListTest] " + msg );
        return false;
    }

    private static Property p( final String ref, final String value ) {
        return new Property( ref, value, "", STR, AppliesTo.NODE );
    }

    /** The refs of a list, in order, space-joined -- so a failure prints the order that was actually built. */
    private static String refs( final PropertiesList props ) {
        return refs( props.getProperties() );
    }

    private static String refs( final List<Property> list ) {
        final StringBuilder sb = new StringBuilder();
        for( final Property prop : list ) {
            if ( sb.length() > 0 ) {
                sb.append( " " );
            }
            sb.append( prop.getRef() );
        }
        return sb.toString();
    }

    // ---- insertion order is kept, even when it is the opposite of alphabetical ------------------------------------
    private static boolean testInsertionOrder() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( p( "data:zebra", "z" ) );
        props.addProperty( p( "data:aardvark", "a" ) );
        if ( !"data:zebra data:aardvark".equals( refs( props ) ) ) {
            return fail( "insertion order must survive (this list must NOT sort by ref): got [" + refs( props ) + "]" );
        }
        // the shape this exists for: a gene matrix grouped core / resistance / mobile, which alphabetizing scrambles
        final PropertiesList genes = new PropertiesList();
        final String[] in_order = { "meta:rpoB", "meta:gyrA", "meta:blaTEM", "meta:vanA", "meta:IS26", "meta:tnpA" };
        for( final String ref : in_order ) {
            genes.addProperty( p( ref, "4" ) );
        }
        final StringBuilder expected = new StringBuilder();
        for( final String ref : in_order ) {
            if ( expected.length() > 0 ) {
                expected.append( " " );
            }
            expected.append( ref );
        }
        if ( !expected.toString().equals( refs( genes ) ) ) {
            return fail( "a grouped gene matrix must keep its column order: expected [" + expected + "] got ["
                    + refs( genes ) + "]" );
        }
        if ( genes.size() != in_order.length ) {
            return fail( "size() should be " + in_order.length + ", got " + genes.size() );
        }
        return true;
    }

    /**
     * The deliberate NON-behaviour, asserted from the other side: adding a property that would sort first must not
     * disturb what is already there. A test that only checks the final order would still pass if the list sorted by
     * something other than ref (insertion time, say), so this pins that EXISTING entries never move.
     */
    private static boolean testNothingReSorts() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( p( "m:mid", "1" ) );
        props.addProperty( p( "z:last", "2" ) );
        final Property first_before = props.getProperties().get( 0 );
        final Property second_before = props.getProperties().get( 1 );
        props.addProperty( p( "a:first_alphabetically", "3" ) );
        if ( props.getProperties().get( 0 ) != first_before ) {
            return fail( "adding an alphabetically-first property moved an existing entry: [" + refs( props ) + "]" );
        }
        if ( props.getProperties().get( 1 ) != second_before ) {
            return fail( "adding a property disturbed the second existing entry: [" + refs( props ) + "]" );
        }
        if ( !"m:mid z:last a:first_alphabetically".equals( refs( props ) ) ) {
            return fail( "a new property must be APPENDED, got [" + refs( props ) + "]" );
        }
        return true;
    }

    // ---- several properties may share one ref, and they keep the order they were added in -------------------------
    private static boolean testDuplicateRefs() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( p( "meta:zeta", "Z1" ) );
        props.addProperty( p( "meta:alpha", "A1" ) );
        props.addProperty( p( "meta:zeta", "Z2" ) );
        if ( props.size() != 3 ) {
            return fail( "a repeated ref must NOT be dropped or merged: size " + props.size() + " [" + refs( props )
                    + "]" );
        }
        if ( !"meta:zeta meta:alpha meta:zeta".equals( refs( props ) ) ) {
            return fail( "repeated refs keep their positions: [" + refs( props ) + "]" );
        }
        final List<Property> zetas = props.getProperties( "meta:zeta" );
        if ( ( zetas.size() != 2 ) || !"Z1".equals( zetas.get( 0 ).getValue() )
                || !"Z2".equals( zetas.get( 1 ).getValue() ) ) {
            return fail( "both values of a repeated ref must come back, in insertion order" );
        }
        return true;
    }

    // ---- the by-ref lookups return their matches in list order -----------------------------------------------------
    private static boolean testLookupsKeepOrder() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( p( "ns:zulu", "1" ) );
        props.addProperty( p( "other:x", "2" ) );
        props.addProperty( p( "ns:alpha", "3" ) );
        if ( !"ns:zulu ns:alpha".equals( refs( props.getPropertiesWithGivenReferencePrefix( "ns:" ) ) ) ) {
            return fail( "a prefix lookup must keep list order, got ["
                    + refs( props.getPropertiesWithGivenReferencePrefix( "ns:" ) ) + "]" );
        }
        if ( props.getPropertiesWithGivenRef( "ns:alpha" ).size() != 1 ) {
            return fail( "getPropertiesWithGivenRef should find exactly one ns:alpha" );
        }
        if ( !props.getProperties( "other:x" ).get( 0 ).getValue().equals( "2" ) ) {
            return fail( "getProperties(ref) should find other:x" );
        }
        if ( !props.getProperties( "nope:none" ).isEmpty() ) {
            return fail( "an absent ref should yield an empty list, not a gap" );
        }
        return true;
    }

    // ---- copy() and the text forms follow the list, so nothing quietly re-sorts on the way out ---------------------
    private static boolean testCopyAndTextKeepOrder() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( p( "data:zebra", "z" ) );
        props.addProperty( p( "data:aardvark", "a" ) );
        final PropertiesList copy = ( PropertiesList ) props.copy();
        if ( !"data:zebra data:aardvark".equals( refs( copy ) ) ) {
            return fail( "copy() must preserve order, got [" + refs( copy ) + "]" );
        }
        final String text = props.asSimpleText().toString();
        if ( text.indexOf( "zebra" ) > text.indexOf( "aardvark" ) ) {
            return fail( "asSimpleText() must follow list order: [" + text + "]" );
        }
        return true;
    }

    // ---- phyloXML output is written in list order ------------------------------------------------------------------
    private static boolean testPhyloXmlWritesInListOrder() {
        try {
            final PropertiesList props = new PropertiesList();
            props.addProperty( p( "data:zebra", "z" ) );
            props.addProperty( p( "data:aardvark", "a" ) );
            final Writer w = new StringWriter();
            props.toPhyloXML( w, 0, "  " );
            final String xml = w.toString();
            if ( xml.indexOf( "data:zebra" ) > xml.indexOf( "data:aardvark" ) ) {
                return fail( "phyloXML must be written in list order, got: " + xml );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "writing phyloXML threw " + e );
        }
    }

    /**
     * The end-to-end claim: properties written in a phyloXML document come back in DOCUMENT order. This is the one
     * a user meets -- it is what makes an open-and-save round trip keep an author's column grouping.
     */
    private static boolean testDocumentOrderSurvivesParse() {
        File tmp = null;
        try {
            final String[] doc_order = { "meta:rpoB", "meta:gyrA", "meta:blaTEM", "meta:aardvark" };
            final StringBuilder xml = new StringBuilder();
            xml.append( "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" );
            xml.append( "<phyloxml xmlns=\"http://www.phyloxml.org\">\n<phylogeny rooted=\"true\">\n" );
            xml.append( "<clade><clade><name>t1</name>\n" );
            for( final String ref : doc_order ) {
                xml.append( "<property ref=\"" ).append( ref )
                        .append( "\" datatype=\"xsd:string\" applies_to=\"node\">v</property>\n" );
            }
            xml.append( "</clade></clade>\n</phylogeny>\n</phyloxml>\n" );
            tmp = File.createTempFile( "aptx_proplist_", ".xml" );
            tmp.deleteOnExit();
            final Writer fw = new FileWriter( tmp );
            fw.write( xml.toString() );
            fw.close();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( tmp, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final PhylogenyNode tip = phy.getExternalNodes().get( 0 );
            if ( ( tip.getNodeData() == null ) || ( tip.getNodeData().getProperties() == null ) ) {
                return fail( "the parsed tip carries no properties -- the fixture did not reach the parser" );
            }
            final String got = refs( tip.getNodeData().getProperties() );
            final StringBuilder expected = new StringBuilder();
            for( final String ref : doc_order ) {
                if ( expected.length() > 0 ) {
                    expected.append( " " );
                }
                expected.append( ref );
            }
            // "meta:aardvark" is last in the document and first alphabetically, so this cannot pass by coincidence
            if ( !expected.toString().equals( got ) ) {
                return fail( "phyloXML document order must survive parsing: expected [" + expected + "] got [" + got
                        + "]" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "the parse round trip threw " + e );
        }
        finally {
            if ( tmp != null ) {
                tmp.delete();
            }
        }
    }

    /**
     * The decided consequence of dropping the sort: the remove-then-re-add path used by the annotation importer and
     * the GTDB / tip-date tools appends, so a REPLACED property lands at the end rather than back in its old slot.
     * Pinned so the relocation stays a decision instead of a surprise.
     */
    private static boolean testReplacedPropertyMovesToTail() {
        final PropertiesList props = new PropertiesList();
        props.addProperty( p( "meta:a", "1" ) );
        props.addProperty( p( "meta:b", "2" ) );
        props.addProperty( p( "meta:c", "3" ) );
        for( final Iterator<Property> it = props.getProperties().iterator(); it.hasNext(); ) {
            if ( "meta:b".equals( it.next().getRef() ) ) {
                it.remove();
            }
        }
        props.addProperty( p( "meta:b", "22" ) );
        if ( !"meta:a meta:c meta:b".equals( refs( props ) ) ) {
            return fail( "a replaced property should move to the tail, got [" + refs( props ) + "]" );
        }
        if ( !"22".equals( props.getProperties( "meta:b" ).get( 0 ).getValue() ) ) {
            return fail( "the replacement value should win" );
        }
        return true;
    }

    // ---- edges -----------------------------------------------------------------------------------------------------
    private static boolean testEdges() {
        final PropertiesList empty = new PropertiesList();
        if ( ( empty.size() != 0 ) || !empty.getProperties().isEmpty() ) {
            return fail( "a new list is empty" );
        }
        if ( empty.asSimpleText().length() != 0 ) {
            return fail( "an empty list has no text" );
        }
        if ( !empty.getPropertiesWithGivenReferencePrefix( "x:" ).isEmpty() ) {
            return fail( "a prefix lookup on an empty list finds nothing" );
        }
        try {
            empty.getPropertiesWithGivenReferencePrefix( "" );
            return fail( "an empty prefix must be rejected" );
        }
        catch ( final IllegalArgumentException expected ) {
            // as documented
        }
        final PropertiesList one = new PropertiesList();
        one.addProperty( p( "only:one", "v" ) );
        if ( ( one.size() != 1 ) || !"only:one".equals( refs( one ) ) ) {
            return fail( "a single property is kept as-is" );
        }
        return true;
    }
}

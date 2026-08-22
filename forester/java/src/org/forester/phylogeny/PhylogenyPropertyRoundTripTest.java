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

package org.forester.phylogeny;

import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Headless round-trip test for PHYLOGENY-level {@code <property applies_to="phylogeny">}: writing a tree-level property
 * with {@link org.forester.io.writers.PhylogenyWriter} and reading it back with the phyloXML parser
 * ({@link org.forester.io.parsers.phyloxml.PhyloXmlHandler}). Verifies the property returns at the phylogeny level and
 * on NO node, that {@link Phylogeny#copy()} preserves it (the undo-snapshot path), that the element is written AFTER
 * the {@code <clade>} (the phyloXML schema order), and that a tree with no tree-level property emits none.
 */
public final class PhylogenyPropertyRoundTripTest {

    private static final String REF   = "aptx:time_axis";
    private static final String VALUE = "v1;GEOLOGIC;0.0;0.0;1;0";

    public static void main( final String[] args ) {
        System.out.println( "PhylogenyPropertyRoundTrip: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            final Phylogeny phy = Phylogeny.createInstanceFromNhxString( "((A,B),C);" );
            phy.setProperties( new PropertiesList() );
            phy.getProperties().addProperty( new Property( REF, VALUE, "", "xsd:string", AppliesTo.PHYLOGENY ) );

            // --- write ---
            final String xml = phy.toPhyloXML( 0 ).toString();
            final int prop_pos = xml.indexOf( "applies_to=\"phylogeny\"" );
            if ( prop_pos < 0 ) {
                return fail( "the phylogeny-level property was not written:\n" + xml );
            }
            // SCHEMA ORDER: in phyloXML 1.10 the <phylogeny> sequence is ...clade?, ..., property*, so a
            // <property applies_to="phylogeny"> must be emitted AFTER the closing </clade>.
            final int last_clade_close = xml.lastIndexOf( "</clade>" );
            if ( ( last_clade_close < 0 ) || ( last_clade_close > prop_pos ) ) {
                return fail( "schema order: <property applies_to=\"phylogeny\"> must come AFTER </clade>:\n" + xml );
            }

            // --- re-parse ---
            final Phylogeny back = parse( xml );
            if ( !back.isHasProperties() ) {
                return fail( "the phylogeny-level property did not round-trip" );
            }
            final List<Property> props = back.getProperties().getProperties();
            if ( ( props.size() != 1 ) || !VALUE.equals( props.get( 0 ).getValue() )
                    || !REF.equals( props.get( 0 ).getRef() )
                    || ( AppliesTo.PHYLOGENY != props.get( 0 ).getAppliesTo() ) ) {
                return fail( "round-tripped property wrong: " + props );
            }
            // ...and it must be attached to NO node (it is tree-level, not node data)
            for( final PhylogenyNodeIterator it = back.iteratorPreorder(); it.hasNext(); ) {
                final PropertiesList node_props = it.next().getNodeData().getProperties();
                if ( ( node_props != null ) && !node_props.getProperties( REF ).isEmpty() ) {
                    return fail( "the tree-level property leaked onto a NODE" );
                }
            }

            // --- copy() preserves it (undo-snapshot path) and the property LIST is independent ---
            final Phylogeny copy = back.copy();
            if ( !copy.isHasProperties() || !VALUE.equals( copy.getProperties().getProperties().get( 0 ).getValue() ) ) {
                return fail( "phy.copy() dropped the tree-level property (undo would lose it)" );
            }
            final int before = back.getProperties().size();
            copy.getProperties()
                    .addProperty( new Property( "aptx:other", "x", "", "xsd:string", AppliesTo.PHYLOGENY ) );
            if ( back.getProperties().size() != before ) {
                return fail( "the copy's property list is not independent of the original (undo would corrupt)" );
            }

            // --- a tree with NO tree-level property emits none ---
            final Phylogeny bare = Phylogeny.createInstanceFromNhxString( "(A,B);" );
            if ( bare.toPhyloXML( 0 ).toString().contains( "applies_to=\"phylogeny\"" ) ) {
                return fail( "a tree without a tree-level property must not emit one" );
            }

            // --- clearing via setProperties(null) then writing emits none ---
            back.setProperties( null );
            if ( back.toPhyloXML( 0 ).toString().contains( "applies_to=\"phylogeny\"" ) ) {
                return fail( "after setProperties(null) no phylogeny-level property should be written" );
            }
            return true;
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            return fail( "unexpected exception: " + e.getMessage() );
        }
    }

    private static Phylogeny parse( final String xml ) throws Exception {
        final PhyloXmlParser p = PhyloXmlParser.createPhyloXmlParser();
        final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance().create( new StringBuffer( xml ), p );
        if ( p.getErrorCount() > 0 ) {
            throw new IllegalStateException( "phyloXML parse errors: " + p.getErrorMessages() );
        }
        return phys[ 0 ];
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [PhylogenyPropertyRoundTripTest] " + msg );
        return false;
    }

    private PhylogenyPropertyRoundTripTest() {
    }
}

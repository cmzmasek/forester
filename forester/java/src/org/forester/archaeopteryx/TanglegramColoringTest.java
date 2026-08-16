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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.forester.archaeopteryx.TanglegramColoring.Field;
import org.forester.archaeopteryx.TanglegramLinker.Link;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headless tests for {@link TanglegramColoring}: which categorical fields are offered (taxonomy + categorical
 * properties, but NOT numeric ones), and that a field maps its distinct values to distinct colours.
 */
public final class TanglegramColoringTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramColoring: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        return availableFieldsOk() && colorMapOk();
    }

    private static boolean availableFieldsOk() {
        // 4 tips, 3 distinct clades (< leaves) so 'clade' is colorable; 4 distinct numeric values -> a gradient
        final Phylogeny tree = tree( clade( tip( "Homo sapiens", "Primates", 1 ),
                                            tip( "Pan troglodytes", "Primates", 2 ),
                                            tip( "Mus musculus", "Rodentia", 3 ),
                                            tip( "Felis catus", "Carnivora", 4 ) ) );
        final List<String> labels = new ArrayList<>();
        for( final Field f : TanglegramColoring.availableFields( tree ) ) {
            labels.add( f.label() );
        }
        if ( !labels.contains( "Taxonomy: Scientific Name" ) ) {
            return fail( "expected a scientific-name field, got " + labels );
        }
        if ( !labels.contains( "data:clade" ) ) {
            return fail( "expected the categorical data:clade field, got " + labels );
        }
        if ( labels.contains( "data:num" ) ) {
            return fail( "the numeric data:num should NOT be offered as a colour field, got " + labels );
        }
        return true;
    }

    private static boolean colorMapOk() {
        final PhylogenyNode a = tip( "Homo sapiens", "Primates", 1 );
        final PhylogenyNode b = tip( "Mus musculus", "Rodentia", 2 );
        final PhylogenyNode c = tip( "Pan troglodytes", "Primates", 3 );
        final Phylogeny tree = tree( clade( a, b, c ) );
        final Field clade = fieldFor( tree, "data:clade" );
        if ( clade == null ) {
            return fail( "no data:clade field available" );
        }
        final List<Link> links = new ArrayList<>();
        links.add( new Link( a, a, "" ) );
        links.add( new Link( b, b, "" ) );
        links.add( new Link( c, c, "" ) );
        final Map<String, Color> map = TanglegramColoring.colorMap( links, clade );
        if ( map.size() != 2 ) {
            return fail( "expected 2 distinct clade colours (Primates, Rodentia), got " + map.size() );
        }
        if ( ( map.get( "Primates" ) == null ) || ( map.get( "Rodentia" ) == null ) ) {
            return fail( "a clade colour is missing: " + map.keySet() );
        }
        if ( map.get( "Primates" ).equals( map.get( "Rodentia" ) ) ) {
            return fail( "distinct clades must get distinct colours" );
        }
        return true;
    }

    private static Field fieldFor( final Phylogeny tree, final String label ) {
        for( final Field f : TanglegramColoring.availableFields( tree ) ) {
            if ( label.equals( f.label() ) ) {
                return f;
            }
        }
        return null;
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramColoring test failed: " + message );
        return false;
    }

    private static PhylogenyNode tip( final String scientific_name, final String clade, final int num ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( scientific_name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        n.getNodeData().setTaxonomy( t );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "data:clade", clade, "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:num", Integer.toString( num ), "", "xsd:decimal", AppliesTo.NODE ) );
        n.getNodeData().setProperties( pl );
        return n;
    }

    private static PhylogenyNode clade( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }
}

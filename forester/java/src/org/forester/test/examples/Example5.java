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

package org.forester.test.examples;

import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

public class Example5 {

    public static void main( final String[] args ) {
        // Creating a new rooted tree with two external nodes.
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode d1 = new PhylogenyNode();
        final PhylogenyNode d2 = new PhylogenyNode();
        // Setting of distances.
        d1.setDistanceToParent( 1.2 );
        d2.setDistanceToParent( 2.4 );
        // Adding species information.
        final Taxonomy t1 = new Taxonomy();
        t1.setScientificName( "Nematostella vectensis" );
        d1.getNodeData().addTaxonomy( t1 );
        final Taxonomy t2 = new Taxonomy();
        t2.setScientificName( "Monosiga brevicollis" );
        d2.getNodeData().addTaxonomy( t2 );
        // Adding gene names.
        final Sequence s1 = new Sequence();
        s1.setName( "Bcl-2" );
        d1.getNodeData().addSequence( s1 );
        final Sequence s2 = new Sequence();
        s2.setName( "Bcl-2" );
        d2.getNodeData().addSequence( s2 );
        // Root is a speciation.
        final Event ev = new Event();
        ev.setSpeciations( 1 );
        ev.setDuplications( 0 );
        root.getNodeData().setEvent( ev );
        // Putting the tree together.
        root.addAsChild( d1 );
        root.addAsChild( d2 );
        phy.setRoot( root );
        phy.setRooted( true );
        // Displaying the newly created tree with Archaeopteryx.
        Archaeopteryx.createApplication( phy );
    }
}

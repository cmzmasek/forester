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

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class Example3 {

    public static void main( final String[] args ) {
        // Creating a new rooted tree with four external nodes.
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        final PhylogenyNode d1 = new PhylogenyNode();
        final PhylogenyNode d2 = new PhylogenyNode();
        final PhylogenyNode d11 = new PhylogenyNode();
        final PhylogenyNode d12 = new PhylogenyNode();
        root.setName( "root" );
        d1.setName( "1" );
        d2.setName( "2" );
        d11.setName( "1-1" );
        d12.setName( "1-2" );
        root.addAsChild( d1 );
        root.addAsChild( d2 );
        d2.addAsChild( d11 );
        d2.addAsChild( d12 );
        phy.setRoot( root );
        phy.setRooted( true );
        // Using a variety of iterators to visit the nodes of the newly created tree.
        System.out.println( "post-order:" );
        for( final PhylogenyNodeIterator it = phy.iteratorPostorder(); it.hasNext(); ) {
            System.out.println( it.next().getName() );
        }
        System.out.println( "pre-order:" );
        for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
            System.out.println( it.next().getName() );
        }
        System.out.println( "level-order:" );
        for( final PhylogenyNodeIterator it = phy.iteratorLevelOrder(); it.hasNext(); ) {
            System.out.println( it.next().getName() );
        }
        System.out.println( "external nodes only:" );
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            System.out.println( it.next().getName() );
        }
    }
}

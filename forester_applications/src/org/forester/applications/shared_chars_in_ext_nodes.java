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

package org.forester.applications;

import java.io.File;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class shared_chars_in_ext_nodes {

    final static boolean SIMPLE = true;

    public static void main( final String args[] ) {
        if ( args.length != 2 ) {
            System.err.println();
            System.err.println( "shared_chars_in_ext_nodes: wrong number of arguments" );
            System.err.println( "Usage: \"shared_chars_in_ext_nodes <intree> <node name>" );
            System.err.println();
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        final String node_name = args[ 1 ];
        Phylogeny phy = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            phy = factory.create( infile, org.forester.io.parsers.util.ParserUtils
                                  .createParserDependingOnFileType( infile, true ) )[ 0 ];
        }
        catch ( final Exception e ) {
            System.err.println( e + "\nCould not read " + infile + "\n" );
            System.exit( -1 );
        }
        final SortedSet<String> a = phy.getNode( node_name ).getNodeData().getBinaryCharacters().getGainedCharacters();
        for( final PhylogenyNodeIterator it = phy.iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final SortedSet<String> b = n.getNodeData().getBinaryCharacters().getGainedCharacters();
            final SortedSet<String> a_copy = copy( a );
            a_copy.retainAll( b );
            final double ratio = ( double ) a_copy.size() / b.size();
            System.out.println( n.getName() + "\t\"" + a_copy.size() + "/" + b.size() + "\"\t" + ratio );
        }
    }

    private static SortedSet<String> copy( final SortedSet<String> set ) {
        final SortedSet<String> copy = new TreeSet<String>();
        for( final String i : set ) {
            copy.add( i );
        }
        return copy;
    }
}

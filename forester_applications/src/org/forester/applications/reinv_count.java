
package org.forester.applications;

// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2011 Christian M. Zmasek
// Copyright (C) 2008-2011 Burnham Institute for Medical Research
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: phylosoft @ gmail . com
// WWW: www.phylosoft.org/forester
// javac -cp ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/reinv_count.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.reinv_count
import java.io.File;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

public class reinv_count {

    public static void main( final String args[] ) {
        if ( args.length != 2 ) {
            System.err.println();
            System.err.println( "reinv_count: wrong number of arguments" );
            System.err.println( "Usage: \"reinv_count <intree> <name>" );
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
        for( final PhylogenyNodeIterator ite = phy.iteratorExternalForward(); ite.hasNext(); ) {
            final PhylogenyNode target_node = ite.next();
            final SortedSet<String> target_dcs = getAllExternalPresentAndGainedCharacters( target_node );
            //System.out.println( "Target DCs:" + target_dcs.size() );
            int counter = 0;
            final SortedSet<String> reinv = new TreeSet<String>();
            for( final String target_dc : target_dcs ) {
                int c = 0;
                for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                    final PhylogenyNode n = it.next();
                    final SortedSet<String> n_gained_dcs = n.getNodeData().getBinaryCharacters().getGainedCharacters();
                    if ( n_gained_dcs.contains( target_dc ) ) {
                        c++;
                    }
                }
                if ( c > 1 ) {
                    counter++;
                    reinv.add( target_dc );
                }
            }
            // System.out.println();
            //System.out.println( "reinv:" + reinv );
            //System.out.println();
            // System.out.println( "Target DCs:" + target_dcs.size() );
            // System.out.println( "reinv size:" + reinv.size() );
            // System.out.println( ">1:" + counter );
            final double ratio = ( double ) counter / target_dcs.size();
            System.out.println( target_node.getName() + "\t" + counter + "/" + target_dcs.size() + "\t" + ratio );
        }
    }

    private static SortedSet<String> getAllExternalPresentAndGainedCharacters( final PhylogenyNode node ) {
        final SortedSet<String> chars = new TreeSet<String>();
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        for( final PhylogenyNode desc : descs ) {
            chars.addAll( desc.getNodeData().getBinaryCharacters().getGainedCharacters() );
            chars.addAll( desc.getNodeData().getBinaryCharacters().getPresentCharacters() );
        }
        return chars;
    }
}

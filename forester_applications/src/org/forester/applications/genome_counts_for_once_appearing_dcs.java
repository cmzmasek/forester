
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
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/genome_counts_for_once_appearing_dcs.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.genome_counts_for_once_appearing_dcs
import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

public class genome_counts_for_once_appearing_dcs {

    public static void main( final String args[] ) {
        if ( args.length != 1 ) {
            System.err.println();
            System.err.println( "genome_counts_for_once_appearing_dcs: wrong number of arguments" );
            System.err.println( "Usage: \"genome_counts_for_once_appearing_dcs <intree>" );
            System.err.println();
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
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
        final SortedSet<String> all_dcs = getAllExternalPresentAndGainedCharacters( phy.getRoot() );
        final SortedSet<String> appearing_once_dcs = new TreeSet<String>();
        System.out.println( "All DCs: " + all_dcs.size() );
        for( final String dc : all_dcs ) {
            int reappearing_count = 0;
            for( final PhylogenyNodeIterator it = phy.iteratorPreorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                SortedSet<String> n_gained_dcs = null;
                if ( n.isRoot() ) {
                    n_gained_dcs = n.getNodeData().getBinaryCharacters().getPresentCharacters();
                }
                else {
                    n_gained_dcs = n.getNodeData().getBinaryCharacters().getGainedCharacters();
                }
                if ( n_gained_dcs.contains( dc ) ) {
                    reappearing_count++;
                }
            }
            if ( reappearing_count < 1 ) {
                System.out.println( "error: " + dc );
                System.exit( -1 );
            }
            if ( reappearing_count == 1 ) {
                appearing_once_dcs.add( dc );
            }
        }
        System.out.println( "Appearing once DCs: " + appearing_once_dcs.size() );
        final DescriptiveStatistics stats = new BasicDescriptiveStatistics();
        final Map<String, Set<String>> node_to_chars = new HashMap<String, Set<String>>();
        final SortedSet<String> appearing_in_all_dcs = new TreeSet<String>();
        for( final String appearing_once_dc : appearing_once_dcs ) {
            int count = 0;
            for( final PhylogenyNodeIterator ite = phy.iteratorExternalForward(); ite.hasNext(); ) {
                final PhylogenyNode ext_node = ite.next();
                if ( !node_to_chars.containsKey( ext_node.getName() ) ) {
                    node_to_chars.put( ext_node.getName(), getAllExternalPresentAndGainedCharacters( ext_node ) );
                }
                if ( node_to_chars.get( ext_node.getName() ).contains( appearing_once_dc ) ) {
                    count++;
                }
            }
            if ( count < 1 ) {
                System.out.println( "error, count is <1" );
                System.exit( -1 );
            }
            if ( count == phy.getNumberOfExternalNodes() ) {
                appearing_in_all_dcs.add( appearing_once_dc );
            }
            stats.addValue( count );
        }
        System.out.println();
        System.out.println( stats.toString() );
        System.out.println();
        final int[] bins = BasicDescriptiveStatistics.performBinning( stats.getDataAsDoubleArray(), 1, 172, 172 );
        for( int i = 0; i < bins.length; i++ ) {
            System.out.println( ( i + 1 ) + "\t" + bins[ i ] );
        }
        System.out.println();
        System.out.println( "appearing in all:" );
        for( final String i : appearing_in_all_dcs ) {
            System.out.println( i );
        }
        System.out.println();
        for( final String dc : appearing_once_dcs ) {
            System.out.println( "1\t" + dc );
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

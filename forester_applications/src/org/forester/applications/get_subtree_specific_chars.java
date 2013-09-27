
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
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/get_subtree_specific_chars.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.get_subtree_specific_chars
import java.io.File;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class get_subtree_specific_chars {

    final static boolean SIMPLE = true;

    public static void main( final String args[] ) {
        if ( args.length != 1 ) {
            System.err.println();
            System.err.println( "get_subtree_specific_chars: wrong number of arguments" );
            System.err.println( "Usage: \"get_subtree_specific_chars <intree>" );
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
        final SortedSet<Long> all_external_ids = getAllExternalDescendantsNodeIds( phy.getRoot() );
        final SortedSet<String> all_chars = getAllExternalPresentAndGainedCharacters( phy.getRoot() );
        System.out.println( "Sum of all external characters:\t" + all_chars.size() );
        System.out.println();
        for( final PhylogenyNodeIterator iter = phy.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( !SIMPLE && node.isExternal() ) {
                continue;
            }
            if ( !node.isRoot() ) {
                // System.out.println();
                if ( node.getNodeData().isHasTaxonomy()
                        && !ForesterUtil.isEmpty( node.getNodeData().getTaxonomy().getScientificName() ) ) {
                    System.out.print( node.getNodeData().getTaxonomy().getScientificName() );
                }
                else {
                    System.out.print( node.getName() );
                }
                // System.out.println( ":" );
                System.out.print( "\t" );
                final SortedSet<Long> external_ids = getAllExternalDescendantsNodeIds( node );
                final SortedSet<Long> not_external_ids = copy( all_external_ids );
                not_external_ids.removeAll( external_ids );
                final SortedSet<String> not_node_chars = new TreeSet<String>();
                for( final Long id : not_external_ids ) {
                    not_node_chars.addAll( getAllExternalPresentAndGainedCharacters( phy.getNode( id ) ) );
                }
                final SortedSet<String> node_chars = getAllExternalPresentAndGainedCharacters( node );
                final SortedSet<String> unique_chars = new TreeSet<String>();
                for( final String node_char : node_chars ) {
                    if ( !not_node_chars.contains( node_char ) ) {
                        if ( SIMPLE ) {
                            unique_chars.add( node_char );
                        }
                        else {
                            boolean found = true;
                            for( final Long external_id : external_ids ) {
                                if ( !phy.getNode( external_id ).getNodeData().getBinaryCharacters()
                                        .getGainedCharacters().contains( node_char )
                                        && !phy.getNode( external_id ).getNodeData().getBinaryCharacters()
                                                .getPresentCharacters().contains( node_char ) ) {
                                    found = false;
                                    break;
                                }
                            }
                            if ( found ) {
                                unique_chars.add( node_char );
                            }
                        }
                    }
                }
                // System.out.println( "\tSUM:\t" + unique_chars.size() );
                // System.out.println( unique_chars.size() );
                int counter = 0;
                System.out.print( "\t" + unique_chars.size() );
                for( final String unique_char : unique_chars ) {
                    // System.out.println( "\t" + counter + ":\t" + unique_char
                    // );
                    // System.out.println( "\t" + counter + ":\t" + unique_char
                    // );
                    System.out.print( "\t" + unique_char );
                    ++counter;
                }
                System.out.println();
            }
        }
    }

    private static SortedSet<Long> copy( final SortedSet<Long> set ) {
        final SortedSet<Long> copy = new TreeSet<Long>();
        for( final Long i : set ) {
            copy.add( i );
        }
        return copy;
    }

    private static SortedSet<Long> getAllExternalDescendantsNodeIds( final PhylogenyNode node ) {
        final SortedSet<Long> ids = new TreeSet<Long>();
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        for( final PhylogenyNode desc : descs ) {
            ids.add( desc.getId() );
        }
        return ids;
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

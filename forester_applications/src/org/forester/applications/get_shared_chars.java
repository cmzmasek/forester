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
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

public class get_shared_chars {

    public final static boolean DEBUG = true;

    public static void main( final String args[] ) {
        if ( args.length < 2 ) {
            System.err.println();
            System.err.println( "get_subtree_specific_chars: wrong number of arguments" );
            System.err.println( "Usage: \"get_shared_chars <intree> <subtree 1> <subtree 2> ... <subtree n>" );
            System.err.println();
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        Phylogeny phy = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            phy = factory.create( infile, ParserUtils.createParserDependingOnFileType( infile, true ) )[ 0 ];
        }
        catch ( final Exception e ) {
            System.err.println( e + "\nCould not read " + infile + "\n" );
            System.exit( -1 );
        }
        final SortedSet<Long> outside_external_ids = getAllExternalDescendantsNodeIds( phy.getRoot() );
        final SortedSet<String> all_chars = getAllExternalPresentAndGainedCharacters( phy.getRoot() );
        System.out.println( "Sum of all external characters:\t" + all_chars.size() );
        final SortedSet<String> all_shared_chars = new TreeSet<String>();
        for( int i = 1; i < args.length; ++i ) {
            System.out.print( args[ i ] + "\t" );
            final PhylogenyNode current_node = phy.getNode( args[ i ] );
            if ( i == 1 ) {
                all_shared_chars.addAll( getAllExternalPresentAndGainedCharacters( current_node ) );
            }
            else {
                all_shared_chars.retainAll( getAllExternalPresentAndGainedCharacters( current_node ) );
            }
            outside_external_ids.removeAll( getAllExternalDescendantsNodeIds( current_node ) );
        }
        System.out.println();
        if ( DEBUG ) {
            System.out.println( "Number of outside nodes: " + outside_external_ids.size() );
        }
        final SortedSet<String> outside_chars = new TreeSet<String>();
        System.out.println( "All shared characters\t" + all_shared_chars.size() );
        for( final Long id : outside_external_ids ) {
            outside_chars.addAll( getAllExternalPresentAndGainedCharacters( phy.getNode( id ) ) );
        }
        final SortedSet<String> unique_shared_chars = copy( all_shared_chars );
        unique_shared_chars.removeAll( outside_chars );
        System.out.println( "Unique shared characters\t" + unique_shared_chars.size() );
        System.out.println();
        System.out.println( "Unique shared characters:" );
        for( final String unique_shared_char : unique_shared_chars ) {
            System.out.println( unique_shared_char );
        }
    }

    private static SortedSet<String> copy( final SortedSet<String> set ) {
        final SortedSet<String> copy = new TreeSet<String>();
        for( final String s : set ) {
            copy.add( s );
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

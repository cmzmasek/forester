
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
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/core_chars.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.core_chars
import java.io.File;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

public class core_chars {

    final static boolean SIMPLE = true;

    public static void main( final String args[] ) {
        if ( args.length != 1 ) {
            System.err.println();
            System.err.println( "core_chars: wrong number of arguments" );
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
        final SortedSet<String> a = getAllExternalPresentAndGainedCharacters( phy.getNode( "Opisthokonta" ) );
        final SortedSet<String> b = getAllExternalPresentAndGainedCharacters( phy.getNode( "THETRA" ) );
        final SortedSet<String> c = getAllExternalPresentAndGainedCharacters( phy.getNode( "Amoebozoa" ) );
        final SortedSet<String> d = getAllExternalPresentAndGainedCharacters( phy.getNode( "Archaeplastida" ) );
        final SortedSet<String> e = getAllExternalPresentAndGainedCharacters( phy.getNode( "Chromalveolate" ) );
        final SortedSet<String> f = getAllExternalPresentAndGainedCharacters( phy.getNode( "Excavata" ) );
        a.retainAll( b );
        a.retainAll( c );
        a.retainAll( d );
        a.retainAll( e );
        a.retainAll( f );
        System.out.println( a.size() );
        for( final String s : a ) {
            System.out.println( s );
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

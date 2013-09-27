// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/get_loss_nodes.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.get_loss_nodes

package org.forester.applications;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class get_loss_nodes {

    public static void main( final String args[] ) {
        if ( args.length != 2 ) {
            System.out.println( "get_loss_nodes: Wrong number of arguments" );
            System.out.println( "Usage: \"get_loss_nodes <phylogeny file> <file with characters>\"" );
            System.exit( -1 );
        }
        final File phylogeny_infile = new File( args[ 0 ] );
        Phylogeny p = null;
        try {
            final PhylogenyParser pp = org.forester.io.parsers.util.ParserUtils
                    .createParserDependingOnFileType( phylogeny_infile, true );
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            p = factory.create( phylogeny_infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        List<String> chars = null;
        try {
            chars = ForesterUtil.file2list( new File( args[ 1 ] ) );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
        for( final String c : chars ) {
            boolean found = false;
            for( final PhylogenyNodeIterator it = p.iteratorPostorder(); it.hasNext(); ) {
                final PhylogenyNode n = it.next();
                if ( n.getNodeData().getBinaryCharacters().getLostCharacters().contains( c ) ) {
                    if ( n.getNodeData().isHasTaxonomy()
                            && !ForesterUtil.isEmpty( n.getNodeData().getTaxonomy().getScientificName() ) ) {
                        System.out.println( c + "\t" + n.getNodeData().getTaxonomy().getScientificName() );
                    }
                    else {
                        System.out.println( c + "\t" + n.getName() );
                    }
                    found = true;
                }
            }
            if ( !found ) {
                System.out.println( c + "\t" + "never lost" );
            }
        }
    }
}
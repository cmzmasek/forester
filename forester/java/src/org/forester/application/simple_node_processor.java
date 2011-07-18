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
// WWW: www.phylosoft.org

package org.forester.application;

import java.io.File;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;

public class simple_node_processor {

    private final static String BASE = "b_";

    public static void main( final String args[] ) {
        File in = null;
        File out = null;
        if ( ( args.length != 2 ) ) {
           // System.exit( -1 );
            if ( ( args.length == 0 ) ) {
                in = new File("C:\\Users\\zma\\dollo.xml");
                out = null;
            }
        }
        try {
            System.out.println( "...");
            CommandLineArguments cla = null;
            cla = new CommandLineArguments( args );
           // in = cla.getFile( 0 );
           // out = cla.getFile( 1 );
           // if ( out.exists() ) {
          //      System.out.println( out + " already exists" );
          //      System.exit( -1 );
          //  }
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = new PhyloXmlParser();
            final Phylogeny[] phylogenies_0 = factory.create( in, xml_parser );
            final Phylogeny phylogeny_0 = phylogenies_0[ 0 ];
            final PhylogenyNodeIterator it = phylogeny_0.iteratorPostorder();
            int i = 0;
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                processNode( node, i );
                i++;
            }
            //  final PhylogenyWriter writer = new PhylogenyWriter();
            //  writer.toPhyloXML( out, phylogeny_0, 0 );
        }
        catch ( final Exception e ) {
            System.out.println( e.getLocalizedMessage() );
            e.printStackTrace();
            System.exit( -1 );
        }
    }

    //    private static void processNode( final PhylogenyNode node, final int i ) {
    //        node.setDistanceToParent( PhylogenyNode.DISTANCE_DEFAULT );
    //        if ( !node.isExternal() ) {
    //            if ( ( node.getName() == null ) || node.getName().isEmpty() ) {
    //                node.setName( BASE + i );
    //            }
    //        }
    //    }
    private static void processNode( final PhylogenyNode node, final int i ) {
        if ( node.isExternal() ) {
            final String c = "" + node.getNodeData().getBinaryCharacters().getPresentCount();
            final String s = node.getNodeData().getTaxonomy().getScientificName();
            System.out.println( s + "\t" + c );
        }
    }
}

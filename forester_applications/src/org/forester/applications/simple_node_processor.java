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
// javac -cp ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/simple_node_processor.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.simple_node_processor

package org.forester.applications;

import java.io.File;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class simple_node_processor {

    public static void main( final String args[] ) {
        File in = null;
        final File out = null;
        try {
            CommandLineArguments cla = null;
            cla = new CommandLineArguments( args );
            in = cla.getFile( 0 );
            //   in = new File( "");
            //out = cla.getFile( 1 );
            // if ( out.exists() ) {
            //      System.out.println( out + " already exists" );
            //      System.exit( -1 );
            //  }
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            final Phylogeny[] phylogenies_0 = factory.create( in, xml_parser );
            final Phylogeny phylogeny_0 = phylogenies_0[ 0 ];
            final PhylogenyNodeIterator it = phylogeny_0.iteratorPostorder();
            int i = 0;
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                processNode( node, i, phylogeny_0 );
                i++;
            }
            final PhylogenyWriter writer = new PhylogenyWriter();
            //writer.toPhyloXML( out, phylogeny_0, 0 );
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
    private static void processNode( final PhylogenyNode node, final int i, final Phylogeny phy ) {
        //if ( node.isExternal() ) {
        //    final String c = "" + node.getNodeData().getBinaryCharacters().getPresentCount();
        //    final String s = node.getNodeData().getTaxonomy().getScientificName();
        //    System.out.println( s + "\t" + c );
        //}
        //        if ( !node.isExternal() ) {
        //            if ( !node.getNodeData().isHasTaxonomy() ) {
        //                if ( !ForesterUtil.isEmpty( node.getName() ) ) {
        //                    if ( ( node.getName().indexOf( "_" ) < 0 ) && ( node.getName().indexOf( "&" ) < 0 )
        //                            && ( node.getName().indexOf( " " ) < 0 ) ) {
        //                        Taxonomy t = new Taxonomy();
        //                        t.setScientificName( node.getName() );
        //                        node.getNodeData().addTaxonomy( t );
        //                        node.setName( "" );
        //                    }
        //                }
        //            }
        //        }
        if ( node.isExternal() ) {
            //final Taxonomy t = node.getNodeData().getTaxonomy();
            //System.out.println( t.getTaxonomyCode() + "\t" + t.getScientificName() + "\t" + t.getCommonName()
            //        + "\t" + label );
            //            if ( node.getNodeData().isHasTaxonomy() ) {
            //                final Taxonomy t = node.getNodeData().getTaxonomy();
            //                if ( !ForesterUtil.isEmpty( t.getTaxonomyCode() ) && ( t.getTaxonomyCode().length() == 5 ) ) {
            //                    if ( node.getName().equalsIgnoreCase( t.getTaxonomyCode() ) ) {
            //                        node.setName( "" );
            //                    }
            //                }
            //            }
            if ( node.getNodeData().isHasTaxonomy() ) {
                final Taxonomy t = node.getNodeData().getTaxonomy();
                if ( !ForesterUtil.isEmpty( t.getTaxonomyCode() ) ) {
                    final String c = t.getTaxonomyCode();
                    if ( c.indexOf( "XX" ) == 3 ) {
                        System.out.println( "FAKE_CODE_TO_ID_MAP.put( \"" + c + "\", " + t.getIdentifier().getValue()
                                + ");" );
                    }
                    //   SurfacingUtil.obtainHexColorStringDependingOnTaxonomyGroup( t.getTaxonomyCode(), phy );
                }
            }
        }
    }
}

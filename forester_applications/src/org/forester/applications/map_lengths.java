// $Id:
//
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2011 Christian M Zmasek
// Copyright (C) 2011 Sanford-Burnham Medical Research Institute
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

package org.forester.applications;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.forester.archaeopteryx.Archaeopteryx;
import org.forester.io.parsers.FastaParser;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesMap;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.sequence.Sequence;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class map_lengths {

    final static private String PRG_NAME = "map_lengths";

    public static void main( final String[] args ) {
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            ;
            final Phylogeny[] phylogenies_0 = factory.create( cla.getFile( 0 ), xml_parser );
            final Phylogeny phy = phylogenies_0[ 0 ];
            for( int i = 1; i < cla.getNumberOfNames(); i++ ) {
                final String fasta_name = cla.getName( i );
                final List<Sequence> seqs = FastaParser.parse( new File( fasta_name ) );
                for( int s = 0; s < seqs.size(); s++ ) {
                    final Sequence seq = seqs.get( s );
                    final int actual_length = seq.getLength() - seq.getNumberOfGapResidues();
                    String node_name = "" + seq.getIdentifier();
                    node_name = node_name.substring( 0, node_name.indexOf( "/" ) );
                    final PhylogenyNode n = phy.getNode( node_name );
                    if ( n.getNodeData().getProperties() == null ) {
                        n.getNodeData().setProperties( new PropertiesMap() );
                    }
                    final PropertiesMap properties = n.getNodeData().getProperties();
                    final Property p = new Property( "r:" + i, "" + actual_length, "", "xsd:integer", AppliesTo.NODE );
                    properties.addProperty( p );
                }
            }
            Archaeopteryx.createApplication( phy );
        }
        catch ( final IOException e ) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }
}

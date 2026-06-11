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

package org.forester.application;

import java.io.File;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Annotation;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;

public class decoratorX {

    private static final int SEQ_NAME_COLUMN = 1;
    private static final int SPECIES_COLUMN  = 2;
    private static final int SEQ_COLUMN      = 3;
    private static final int TARGET_COLUMN   = 4;

    public static void main( final String args[] ) {
        File intree = null;
        File outtree1 = null;
        File outtree2 = null;
        File intable = null;
        try {
            CommandLineArguments cla = null;
            cla = new CommandLineArguments( args );
            intree = cla.getFile( 0 );
            intable = cla.getFile( 1 );
            outtree1 = cla.getFile( 2 );
            outtree2 = cla.getFile( 3 );
            if ( outtree1.exists() ) {
                System.out.println( outtree1 + " already exists" );
                System.exit( -1 );
            }
            if ( outtree2.exists() ) {
                System.out.println( outtree2 + " already exists" );
                System.exit( -1 );
            }
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhyloXmlParser xml_parser = PhyloXmlParser.createPhyloXmlParserXsdValidating();
            final Phylogeny phy = factory.create( intree, xml_parser )[ 0 ];
            final BasicTable<String> t = BasicTableParser.parse( intable, '\t' );
            final PhylogenyNodeIterator it = phy.iteratorExternalForward();
            int i = 0;
            while ( it.hasNext() ) {
                final PhylogenyNode node = it.next();
                processNode( node, t );
                i++;
            }
            final PhylogenyWriter writer1 = new PhylogenyWriter();
            writer1.toPhyloXML( outtree1, phy, 0 );
            final PhylogenyNodeIterator it2 = phy.iteratorExternalForward();
            while ( it2.hasNext() ) {
                final PhylogenyNode node = it2.next();
                processNode2( node, phy );
            }
            final PhylogenyWriter writer2 = new PhylogenyWriter();
            writer2.toPhyloXML( outtree2, phy, 0 );
        }
        catch ( final Exception e ) {
            System.out.println( e.getLocalizedMessage() );
            System.exit( -1 );
        }
    }

    private static void processNode( final PhylogenyNode node, final BasicTable<String> t ) throws Exception {
        final String node_seq = node.getNodeData().getSequence().getMolecularSequence().toUpperCase();
        boolean found = false;
        String found_row = "";
        String found_protein_name = "";
        String found_species = "";
        for( int row = 0; row < t.getNumberOfRows(); ++row ) {
            final String table_seq = t.getValueAsString( SEQ_COLUMN, row ).toUpperCase();
            if ( table_seq.contains( node_seq ) ) {
                if ( found ) {
                    if ( !found_protein_name.equals( t.getValueAsString( SEQ_NAME_COLUMN, row ) )
                            || !found_species.equals( t.getValueAsString( SPECIES_COLUMN, row ) ) ) {
                        throw new Exception( "Sequence from node " + node + " is not unique: " + node_seq + "\n"
                                + "Already found in row " + found_row );
                    }
                }
                else {
                    found = true;
                    found_row = t.getRowAsString( row, ", " );
                    found_protein_name = t.getValueAsString( SEQ_NAME_COLUMN, row );
                    found_species = t.getValueAsString( SPECIES_COLUMN, row );
                }
                final Annotation annotation = new Annotation( "target", t.getValueAsString( TARGET_COLUMN, row ) );
                node.getNodeData().getSequence().addAnnotation( annotation );
                System.out.println( node + "->" + annotation );
            }
        }
    }

    private static void processNode2( final PhylogenyNode node, final Phylogeny t ) {
        if ( ( node.getNodeData().getSequence().getAnnotations() == null )
                || node.getNodeData().getSequence().getAnnotations().isEmpty() ) {
            t.deleteSubtree( node, true );
        }
    }
}

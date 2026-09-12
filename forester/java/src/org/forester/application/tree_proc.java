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
import java.io.IOException;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.ForesterUtil;

public class tree_proc {

    private static final String PRG_DATE     = "2021-10-08";
    private static final String PRG_VERSION  = "0.0.1";
    private static final String PRG_NAME     = "tree_proc";
    private static final String XSD_STRING   = "xsd:string";
    private static final String VIPR_HOST    = "vipr:Hosts";
    private static final String VIPR_YEAR    = "vipr:Year";
    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_VERSION, PRG_DATE );
        System.out.println();
        if ( ( args.length != 2 ) ) {
            System.out.println( PRG_NAME + ": Wrong number of arguments.\n" );
            System.out.println( "Usage: " + PRG_NAME + " <in-tree> <out-tree> \n" );
            System.exit( -1 );
        }
        final File intree = new File( args[ 0 ] );
        final File outtree = new File( args[ 1 ] );
        final String e0 = ForesterUtil.isWritableFile( outtree );
        if ( !ForesterUtil.isEmpty( e0 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e0 );
        }
        final String e1 = ForesterUtil.isReadableFile( intree );
        if ( !ForesterUtil.isEmpty( e1 ) ) {
            ForesterUtil.fatalError( PRG_NAME, e1 );
        }
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intree, true );
            p = factory.create( intree, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, "Could not read \"" + intree + "\" [" + e.getMessage() + "]" );
        }
        ForesterUtil
                .programMessage( PRG_NAME,
                                 "Successfully read in tree with " + p.getNumberOfExternalNodes() + " external nodes" );
        final int properties_added = 0;
        for( final PhylogenyNodeIterator iter = p.iteratorPostorder(); iter.hasNext(); ) {
            final PhylogenyNode node = iter.next();
            if ( node.isExternal() ) {
                node.getNodeData().setSequence(null);


            }


        }
        try {
            final PhylogenyWriter writer = new PhylogenyWriter();
            writer.toPhyloXML( p, 0, outtree );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage() );
        }
        System.out.println();
        System.out.println( "Sum of properties added: " + properties_added );
        System.out.println( "Wrote outtree to       : " + outtree );
        System.out.println();
    }
}

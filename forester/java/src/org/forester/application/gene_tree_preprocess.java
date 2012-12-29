// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2012 Christian M. Zmasek
// Copyright (C) 2008-2012 Burnham Institute for Medical Research
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.SequenceDbWsTools;

public class gene_tree_preprocess {

    final static private String HELP_OPTION_1           = "help";
    final static private String HELP_OPTION_2           = "h";
    final static private String PRG_NAME                = "gene_tree_preprocess";
    final static private String PRG_DESC                = "gene tree preprocessing for SDI analysis";
    final static private String PRG_VERSION             = "1.01";
    final static private String PRG_DATE                = "2012.06.07";
    final static private String E_MAIL                  = "phylosoft@gmail.com";
    final static private String WWW                     = "www.phylosoft.org/forester";
    private final static int    DEFAULT_LINES_TO_RETURN = 50;

    public static void main( final String[] args ) {
        try {
            final CommandLineArguments cla = new CommandLineArguments( args );
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) || ( args.length != 1 ) ) {
                printHelp();
                System.exit( 0 );
            }
            final File in = cla.getFile( 0 );
            Phylogeny phy = null;
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            try {
                phy = factory.create( in, ParserUtils.createParserDependingOnFileType( in, true ) )[ 0 ];
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "failed to read phylogeny from [" + in + "]: " + e.getLocalizedMessage() );
            }
            final File outtree = new File( ForesterUtil.removeSuffix( in.toString() )
                    + "_preprocessed_gene_tree.phylo.xml" );
            final File removed_nodes = new File( ForesterUtil.removeSuffix( in.toString() ) + "_removed_nodes.txt" );
            final File present_species = new File( ForesterUtil.removeSuffix( in.toString() ) + "_species_present.txt" );
            checkForOutputFileWriteability( outtree );
            checkForOutputFileWriteability( removed_nodes );
            checkForOutputFileWriteability( present_species );
            if ( phy.getNumberOfExternalNodes() < 2 ) {
                ForesterUtil.fatalError( PRG_NAME, "phylogeny has " + phy.getNumberOfExternalNodes()
                        + " external node(s), aborting" );
            }
            final SortedSet<String> not_found = SequenceDbWsTools.obtainSeqInformation( phy,
                                                                                        true,
                                                                                        false,
                                                                                        DEFAULT_LINES_TO_RETURN );
            for( final String remove_me : not_found ) {
                phy.deleteSubtree( phy.getNode( remove_me ), true );
            }
            phy.clearHashIdToNodeMap();
            phy.externalNodesHaveChanged();
            if ( phy.getNumberOfExternalNodes() < 2 ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "after removal of unresolvable external nodes, phylogeny has "
                                                 + phy.getNumberOfExternalNodes() + " external node(s), aborting" );
            }
            try {
                final PhylogenyWriter writer = new PhylogenyWriter();
                writer.toPhyloXML( phy, 0, outtree );
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME, "failed to write to [" + outtree + "]: " + e.getLocalizedMessage() );
            }
            ForesterUtil.programMessage( PRG_NAME, "wrote output phylogeny to: " + outtree );
            final SortedSet<String> species_set = new TreeSet<String>();
            for( final PhylogenyNodeIterator iter = phy.iteratorExternalForward(); iter.hasNext(); ) {
                final PhylogenyNode node = iter.next();
                if ( node.getNodeData().isHasTaxonomy() ) {
                    final String sn = node.getNodeData().getTaxonomy().getScientificName();
                    if ( !ForesterUtil.isEmpty( sn ) ) {
                        species_set.add( sn );
                    }
                }
            }
            try {
                final BufferedWriter out = new BufferedWriter( new FileWriter( present_species ) );
                for( final String species : species_set ) {
                    out.write( species );
                    out.newLine();
                }
                out.close();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "failed to write to [" + present_species + "]: " + e.getLocalizedMessage() );
            }
            ForesterUtil.programMessage( PRG_NAME, "wrote present species to: " + present_species );
            try {
                final BufferedWriter out = new BufferedWriter( new FileWriter( removed_nodes ) );
                for( final String remove_me : not_found ) {
                    out.write( remove_me );
                    out.newLine();
                }
                out.close();
            }
            catch ( final IOException e ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "failed to write to [" + removed_nodes + "]: " + e.getLocalizedMessage() );
            }
            ForesterUtil.programMessage( PRG_NAME, "wrote removed external nodes labels to: " + removed_nodes );
            ForesterUtil.programMessage( PRG_NAME, "OK" );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private static void checkForOutputFileWriteability( final File outfile ) {
        final String error = ForesterUtil.isWritableFile( outfile );
        if ( !ForesterUtil.isEmpty( error ) ) {
            ForesterUtil.fatalError( PRG_NAME, error );
        }
    }

    private static void printHelp() {
        ForesterUtil.printProgramInformation( PRG_NAME,
                                              PRG_DESC,
                                              PRG_VERSION,
                                              PRG_DATE,
                                              E_MAIL,
                                              WWW,
                                              ForesterUtil.getForesterLibraryInformation() );
        System.out.print( "Usage: " );
        System.out.println( PRG_NAME + " <input phylogeny file>" );
        System.out.println();
    }
}

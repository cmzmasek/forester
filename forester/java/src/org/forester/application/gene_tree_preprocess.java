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
// WWW: www.phylosoft.org/forester

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.SortedSet;

import org.forester.archaeopteryx.tools.SequenceDataRetriver;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class gene_tree_preprocess {

    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";
    final static private String FROM_OPTION   = "f";
    final static private String TO_OPTION     = "t";
    final static private String STEP_OPTION   = "s";
    final static private String WINDOW_OPTION = "w";
    final static private String PRG_NAME      = "gene_tree_preprocess";
    final static private String PRG_DESC      = "gene tree preprocessing for SDI analysis";
    final static private String PRG_VERSION   = "1.00";
    final static private String PRG_DATE      = "2012.06.07";
    final static private String E_MAIL        = "phylosoft@gmail.com";
    final static private String WWW           = "www.phylosoft.org/forester/";

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
                                         "failed to read target phylogenies from [" + in + "]: "
                                                 + e.getLocalizedMessage() );
            }
            final File outtree = new File( ForesterUtil.removeSuffix( in.toString() )
                    + "_preprocessed_gene_tree.phylo.xml" );
            final File removed_nodes = new File( ForesterUtil.removeSuffix( in.toString() ) + "_removed_nodes.txt" );
            checkForOutputFileWriteability( outtree );
            checkForOutputFileWriteability( removed_nodes );
            if ( phy.getNumberOfExternalNodes() < 2 ) {
                ForesterUtil.fatalError( PRG_NAME, "phylogeny has " + phy.getNumberOfExternalNodes()
                        + " external node(s), aborting" );
            }
            final SortedSet<String> not_found = SequenceDataRetriver.obtainSeqInformation( phy );
            for( final String remove_me : not_found ) {
                System.out.println( " not found: " + not_found );
                PhylogenyMethods.removeNode( phy.getNode( remove_me ), phy );
            }
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
            //  ForesterUtil.programMessage( PRG_NAME, "wrote output to: [" + outfile + "]" );
            ForesterUtil.programMessage( PRG_NAME, "OK" );
            System.out.println();
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    public static void checkForOutputFileWriteability( final File outfile ) {
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
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( PRG_NAME + " <options> <msa input file>" );
        System.out.println();
        System.out.println( " options: " );
        System.out.println();
        System.out.println( "   -" + FROM_OPTION + "=<integer>: from (msa column)" );
        System.out.println( "   -" + TO_OPTION + "=<integer>: to (msa column)" );
        System.out.println( "    or" );
        System.out.println( "   -" + WINDOW_OPTION + "=<integer>: window size (msa columns)" );
        System.out.println( "   -" + STEP_OPTION + "=<integer>: step size (msa columns)" );
        System.out.println();
        System.out.println();
        System.out.println();
    }
}

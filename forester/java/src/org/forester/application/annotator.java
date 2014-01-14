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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.forester.analysis.AncestralTaxonomyInference;
import org.forester.analysis.AncestralTaxonomyInferenceException;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;
import org.forester.ws.seqdb.SequenceDbWsTools;

public final class annotator {

    final static private String PRG_NAME    = "annotator";
    final static private String PRG_VERSION = "1.00";
    final static private String PRG_DATE    = "131122";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( annotator.PRG_NAME, annotator.PRG_VERSION, annotator.PRG_DATE );
        System.out.println();
        if ( ( args.length != 2 ) ) {
            annotator.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        final File indir = cla.getFile( 0 );
        final File outdir = cla.getFile( 1 );
        if ( !indir.isDirectory() ) {
            ForesterUtil.fatalError( PRG_NAME, indir + " is not a directory" );
        }
        if ( !outdir.isDirectory() ) {
            ForesterUtil.fatalError( PRG_NAME, outdir + " is not a directory" );
        }
        final File[] list_of_files = indir.listFiles();
        final List<File> infiles = new ArrayList<File>();
        for( final File file : list_of_files ) {
            if ( file.isFile() && file.canRead() && file.toString().toLowerCase().endsWith( ".xml" ) ) {
                infiles.add( file );
            }
        }
        Collections.sort( infiles );
        int c = 0;
        for( final File infile : infiles ) {
            System.out.println( ++c + "/" + infiles.size() + ": " + infile );
            final File outfile = new File( outdir.getAbsolutePath().toString() + "/" + infile.getName() );
            if ( outfile.exists() ) {
                System.out.println( outfile + " already exists" );
            }
            else {
                Phylogeny phy = null;
                try {
                    final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                    final Phylogeny[] phylogenies = factory.create( infile,
                                                                    PhyloXmlParser.createPhyloXmlParserXsdValidating() );
                    phy = phylogenies[ 0 ];
                }
                catch ( final Exception e ) {
                    ForesterUtil
                            .fatalError( PRG_NAME, "failed to read phylgenies from [" + infile + "] [" + e.getMessage()
                                    + "]" );
                }
                try {
                    obtainSeqInformation( phy );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
                }
                //                try {
                //                    inferTaxonomyFromDescendents( phy );
                //                }
                //                catch ( final IOException e ) {
                //                    ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
                //                }
                //                catch ( final AncestralTaxonomyInferenceException e ) {
                //                    ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
                //                }
                //phy.setRerootable( false );
                try {
                    final PhylogenyWriter w = new PhylogenyWriter();
                    w.toPhyloXML( phy, 0, outfile );
                }
                catch ( final IOException e ) {
                    ForesterUtil.fatalError( PRG_NAME, "failed to write output [" + e.getMessage() + "]" );
                }
            }
        }
    }

    private static void obtainSeqInformation( final Phylogeny phy ) throws IOException {
        SequenceDbWsTools.obtainSeqInformation( phy, true, true, SequenceDbWsTools.DEFAULT_LINES_TO_RETURN );
    }

    private static void inferTaxonomyFromDescendents( final Phylogeny phy ) throws IOException,
            AncestralTaxonomyInferenceException {
        AncestralTaxonomyInference.inferTaxonomyFromDescendents( phy );
    }

    private static void argumentsError() {
        System.out.println( annotator.PRG_NAME + " <indir> <outdir>" );
        System.out.println();
        System.exit( -1 );
    }
}

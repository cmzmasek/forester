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
import java.util.List;

import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

public class phylostrip {

    public static void main( final String args[] ) {
        if ( args.length < 4 ) {
            System.out.println( "\nstrip: Wrong number of arguments.\n" );
            System.out
            .println( "Usage: \"phylostrip <in-tree> <out-tree> <options> [name1] [name2] ... OR [ref-tree]\"\n" );
            System.out.println( " Options: -knn to keep listed nodes" );
            System.out.println( "          -rnn to remove listed nodes" );
            System.out.println( "          -knnp to keep nodes found in [ref-tree]" );
            System.out.println( "          -rnnp to remove nodes found in [ref-tree]" );
            System.out.println( "          -ktc to keep only nodes from listed taxonomy codes\n" );
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        final File outfile = new File( args[ 1 ] );
        final String options = args[ 2 ];
        Phylogeny p = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( infile, true );
            p = factory.create( infile, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + infile + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        boolean keep = false;
        boolean from_p0 = false;
        boolean ktc = false;
        if ( options.trim().toLowerCase().equals( "-knn" ) ) {
            keep = true;
        }
        else if ( options.trim().toLowerCase().equals( "-knnp" ) ) {
            keep = true;
            from_p0 = true;
        }
        else if ( options.trim().toLowerCase().equals( "-rnnp" ) ) {
            from_p0 = true;
        }
        else if ( options.trim().toLowerCase().equals( "-ktc" ) ) {
            ktc = true;
        }
        else if ( !options.trim().toLowerCase().equals( "-rnn" ) ) {
            System.out.println( "\nUnknown option \"" + options + "\"\n" );
            System.exit( -1 );
        }
        String[] names = null;
        if ( from_p0 ) {
            names = phylostrip.readInNamesFromPhylogeny( args[ 3 ] );
        }
        else {
            names = new String[ args.length - 3 ];
            for( int i = 0; i < ( args.length - 3 ); ++i ) {
                names[ i ] = args[ i + 3 ];
            }
        }
        if ( ktc ) {
            final List<Taxonomy> taxonomies_to_keep = new ArrayList<Taxonomy>();
            for( final String n : names ) {
                final Taxonomy t = new Taxonomy();
                try {
                    t.setTaxonomyCode( n );
                }
                catch ( final PhyloXmlDataFormatException e ) {
                    System.out.println( e.getMessage() );
                    System.exit( -1 );
                }
                taxonomies_to_keep.add( t );
            }
            PhylogenyMethods.deleteExternalNodesPositiveSelectionT( taxonomies_to_keep, p );
        }
        else if ( keep ) {
            PhylogenyMethods.deleteExternalNodesPositiveSelection( names, p );
        }
        else {
            PhylogenyMethods.deleteExternalNodesNegativeSelection( names, p );
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( outfile, p, 0 );
        }
        catch ( final IOException e ) {
            System.out.println( "\nFailure to write output [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
    }

    private static String[] readInNamesFromPhylogeny( final String file ) {
        Phylogeny p0 = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final File f = new File( file );
            final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( f, true );
            p0 = factory.create( f, pp )[ 0 ];
        }
        catch ( final Exception e ) {
            System.out.println( "\nCould not read \"" + file + "\" [" + e.getMessage() + "]\n" );
            System.exit( -1 );
        }
        return p0.getAllExternalNodeNames();
    }
}

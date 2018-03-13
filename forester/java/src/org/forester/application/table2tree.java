// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2018 Christian M. Zmasek
// Copyright (C) 2018 J. Craig Venter Institute
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
// WWW: https://sites.google.com/site/cmzmasek/home/software/forester

package org.forester.application;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Identifier;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public final class table2tree {

    final static private String PRG_NAME    = "table2tree";
    final static private String PRG_VERSION = "1.00";
    final static private String PRG_DATE    = "180312";

    public static void main( final String args[] ) {
        ForesterUtil.printProgramInformation( table2tree.PRG_NAME, table2tree.PRG_VERSION, table2tree.PRG_DATE );
        System.out.println();
        if ( ( args.length != 2 ) ) {
            table2tree.argumentsError();
        }
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        final File intable = cla.getFile( 0 );
        final File outfile = cla.getFile( 1 );
        if ( outfile.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, outfile + " already exists" );
        }
        if ( !intable.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, intable + " does not exist" );
        }
        BasicTable<String> t = null;
        try {
            t = BasicTableParser.parse( intable, '\t' );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        System.out.println( "Number of rows: " + t.getNumberOfRows() );
        final Phylogeny p = new Phylogeny();
        p.setRooted( true );
        p.setRerootable( false );
        p.setDescription( "based on " + intable );
        final PhylogenyNode root = new PhylogenyNode();
        root.setName( "r" );
        p.setRoot( root );
        for( int r = 0; r < t.getNumberOfRows(); ++r ) {
            final String sn = t.getValueAsString( 0, r );
            final String id = t.getValueAsString( 1, r );
            final String code = t.getValueAsString( 2, r );
            final String level1 = t.getValueAsString( 3, r );
            final String level2 = t.getValueAsString( 4, r );
            if ( level1.startsWith( "unclassified" ) || level1.startsWith( "unassigned" )
                    || level2.startsWith( "unclassified" ) || level2.startsWith( "unassigned" ) ) {
                continue;
            }
            final List<PhylogenyNode> nodes1 = p.getNodesViaScientificName( level1 );
            if ( nodes1.size() < 1 ) {
                final PhylogenyNode n = new PhylogenyNode();
                final Taxonomy tax = new Taxonomy();
                tax.setScientificName( level1 );
                n.getNodeData().addTaxonomy( tax );
                p.getRoot().addAsChild( n );
                nodes1.add( n );
            }
            else if ( nodes1.size() > 1 ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "This should not have happened: " + level1 + " is somehow unspecific" );
            }
            final List<PhylogenyNode> nodes2 = p.getNodesViaScientificName( level2 );
            if ( nodes2.size() < 1 ) {
                final PhylogenyNode n = new PhylogenyNode();
                final Taxonomy tax = new Taxonomy();
                tax.setScientificName( level2 );
                n.getNodeData().addTaxonomy( tax );
                nodes1.get( 0 ).addAsChild( n );
                nodes2.add( n );
            }
            else if ( nodes2.size() > 1 ) {
                ForesterUtil.fatalError( PRG_NAME,
                                         "This should not have happened: " + level2 + " is somehow unspecific" );
            }
            final PhylogenyNode n = new PhylogenyNode();
            final Taxonomy tax = new Taxonomy();
            tax.setScientificName( sn );
            try {
                tax.setTaxonomyCode( code );
            }
            catch ( final PhyloXmlDataFormatException e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
            tax.setIdentifier( new Identifier( id, "uniprot" ) );
            n.getNodeData().addTaxonomy( tax );
            nodes2.get( 0 ).addAsChild( n );
        }
        try {
            final PhylogenyWriter w = new PhylogenyWriter();
            w.toPhyloXML( p, 0, outfile );
        }
        catch ( final IOException e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        System.out.println( "Successfully wrote: " + outfile );
        System.out.println();
    }

    private static void argumentsError() {
        System.out.println( PRG_NAME + " <infile> <outfile>" );
        System.out.println();
        System.exit( -1 );
    }
}

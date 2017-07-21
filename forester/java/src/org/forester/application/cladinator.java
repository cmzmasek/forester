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
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

public class cladinator {

    final static private String PRG_NAME      = "cladinator";
    final static private String PRG_VERSION   = "0.100";
    final static private String PRG_DATE      = "170721";
    final static private String PRG_DESC      = "clades within clades";
    final static private String E_MAIL        = "phyloxml@gmail.com";
    final static private String WWW           = "https://sites.google.com/site/cmzmasek/home/software/forester";
    final static private String HELP_OPTION_1 = "help";
    final static private String HELP_OPTION_2 = "h";

    public static void main( final String args[] ) {
        try {
            ForesterUtil.printProgramInformation( PRG_NAME,
                                                  PRG_DESC,
                                                  PRG_VERSION,
                                                  PRG_DATE,
                                                  E_MAIL,
                                                  WWW,
                                                  ForesterUtil.getForesterLibraryInformation() );
            CommandLineArguments cla = null;
            try {
                cla = new CommandLineArguments( args );
            }
            catch ( final Exception e ) {
                ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
            }
            if ( cla.isOptionSet( HELP_OPTION_1 ) || cla.isOptionSet( HELP_OPTION_2 ) ) {
                System.out.println();
                print_help();
                System.exit( 0 );
            }
            else if ( ( args.length != 2 ) ) {
                System.out.println();
                System.out.println( "Wrong number of arguments." );
                System.out.println();
                print_help();
                System.exit( -1 );
            }
            final List<String> allowed_options = new ArrayList<String>();
            final File intreefile = cla.getFile( 0 );
            final String query = cla.getName( 1 );
            System.out.println( "Input tree: " + intreefile );
            System.out.println( "Query:      " + query );
            Phylogeny p = null;
            try {
                final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
                final PhylogenyParser pp = ParserUtils.createParserDependingOnFileType( intreefile, true );
                p = factory.create( intreefile, pp )[ 0 ];
            }
            catch ( final Exception e ) {
                System.out.println( "\nCould not read \"" + intreefile + "\" [" + e.getMessage() + "]\n" );
                System.exit( -1 );
            }
            execute( p, query );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
    }

    private static void execute( final Phylogeny p, final String query ) {
        final PhylogenyNode qnode = p.getNode( query );
        if ( qnode.isRoot() ) {
            throw new IllegalStateException( "Unexpected error: Query " + query
                    + " is root. This should have never happened" );
        }
        if ( qnode.getParent().isRoot() ) {
            throw new IllegalStateException( "Unexpected error: Parent of query " + query
                    + " is root. This should have never happened" );
        }
        final PhylogenyNode qnode_pp = qnode.getParent().getParent();
        final List<PhylogenyNode> qnode_ext_nodes = qnode_pp.getAllExternalDescendants();
        final int lec_ext_nodes = qnode_ext_nodes.size() - 1;
        final int p_ext_nodes = p.getNumberOfExternalNodes() - 1;
        final double lec_ratio = ( 100.0 * lec_ext_nodes ) / p_ext_nodes;
        final List<String> qnode_ext_nodes_names = new ArrayList<String>();
        for( final PhylogenyNode qnode_ext_node : qnode_ext_nodes ) {
            String name = qnode_ext_node.getName();
            if ( ForesterUtil.isEmptyTrimmed( name ) ) {
                throw new IllegalArgumentException( "external node(s) with empty names found" );
            }
            name = name.trim();
            if ( !name.equals( query ) ) {
                qnode_ext_nodes_names.add( name );
            }
        }
        final String greatest_common_prefix = ForesterUtil.greatestCommonPrefix( qnode_ext_nodes_names );
        System.out.println( );
        System.out.println( "Results:");
        if ( greatest_common_prefix.length() < 1 ) {
            System.out.println( "WARNING: No greatest common prefix" );
        }
        else {
            System.out.println( "Greatest common prefix: " + greatest_common_prefix );
        }
        if ( qnode_pp.isRoot() ) {
            System.out.println( "WARNING: Least Encompassing Clade is entire tree" );
        }
        System.out.println( "Least Encompassing Clade has " + lec_ext_nodes + " external nodes (" +lec_ratio + "% of a total of "+ p_ext_nodes +")" );
    }

    private final static void print_help() {
        System.out.println( "Usage: " + PRG_NAME
                + " <gene tree file> <query>" );
        System.out.println();
    }
}

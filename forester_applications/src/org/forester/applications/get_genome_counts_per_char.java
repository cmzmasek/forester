
package org.forester.applications;

// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2011 Christian M. Zmasek
// Copyright (C) 2008-2011 Burnham Institute for Medical Research
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
// javac -cp ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/get_genome_counts_per_char.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.get_genome_counts_per_char
import java.io.File;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;

public class get_genome_counts_per_char {

    final static boolean SIMPLE = true;

    public static void main( final String args[] ) {
        if ( args.length != 1 ) {
            System.err.println();
            System.err.println( "get_genome_counts_per_char: wrong number of arguments" );
            System.err.println( "Usage: \"get_subtree_specific_chars <intree>" );
            System.err.println();
            System.exit( -1 );
        }
        final File infile = new File( args[ 0 ] );
        Phylogeny phy = null;
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            phy = factory.create( infile, org.forester.io.parsers.util.ParserUtils
                    .createParserDependingOnFileType( infile, true ) )[ 0 ];
        }
        catch ( final Exception e ) {
            System.err.println( e + "\nCould not read " + infile + "\n" );
            System.exit( -1 );
        }
        final SortedSet<String> all_chars = getAllExternalPresentAndGainedCharacters( phy.getRoot() );
        final SortedSet<String> human = getAllExternalPresentAndGainedCharacters( phy.getNode( "HUMAN" ) );
        final SortedSet<String> primates = getAllExternalPresentAndGainedCharacters( find( "Primates", phy ) );
        final SortedSet<String> mammalia = getAllExternalPresentAndGainedCharacters( find( "Mammalia", phy ) );
        final SortedSet<String> metazoa = getAllExternalPresentAndGainedCharacters( find( "Metazoa", phy ) );
        final SortedSet<String> fungi = getAllExternalPresentAndGainedCharacters( find( "Fungi", phy ) );
        final SortedSet<String> plants = getAllExternalPresentAndGainedCharacters( find( "Viridiplantae", phy ) );
        System.out.println( "Sum of all external characters:\t" + all_chars.size() );
        System.out.println();
        final List<PhylogenyNode> ext = phy.getRoot().getAllExternalDescendants();
        System.out.println( "genomes" + "\t" + ext.size() );
        for( final String c : all_chars ) {
            int count = 0;
            for( final PhylogenyNode e : ext ) {
                if ( e.getNodeData().getBinaryCharacters().getGainedCharacters().contains( c )
                        || e.getNodeData().getBinaryCharacters().getPresentCharacters().contains( c ) ) {
                    count++;
                }
            }
            if ( count < 1 ) {
                System.err.println( "error" );
                System.exit( -1 );
            }
            System.out.print( c + "\t" + count + "\t" );
            if ( human.contains( c ) ) {
                System.out.print( "HUMAN" + "\t" );
            }
            else {
                System.out.print( "" + "\t" );
            }
            if ( primates.contains( c ) ) {
                System.out.print( "PRIMATES" + "\t" );
            }
            else {
                System.out.print( "" + "\t" );
            }
            if ( mammalia.contains( c ) ) {
                System.out.print( "MAMMALS" + "\t" );
            }
            else {
                System.out.print( "" + "\t" );
            }
            if ( metazoa.contains( c ) ) {
                System.out.print( "METAZOA" + "\t" );
            }
            else {
                System.out.print( "" + "\t" );
            }
            if ( fungi.contains( c ) ) {
                System.out.print( "FUNGI" + "\t" );
            }
            else {
                System.out.print( "" + "\t" );
            }
            if ( plants.contains( c ) ) {
                System.out.print( "PLANTS" + "\t" );
            }
            else {
                System.out.print( "" + "\t" );
            }
            System.out.println();
        }
    }

    private static PhylogenyNode find( final String s, final Phylogeny phy ) {
        final List<PhylogenyNode> l = PhylogenyMethods.searchData( s, phy, true, false, false );
        if ( l.size() != 1 ) {
            System.err.println( "error: " + s );
            System.exit( -1 );
        }
        return l.get( 0 );
    }

    private static SortedSet<String> getAllExternalPresentAndGainedCharacters( final PhylogenyNode node ) {
        final SortedSet<String> chars = new TreeSet<String>();
        final List<PhylogenyNode> descs = node.getAllExternalDescendants();
        for( final PhylogenyNode desc : descs ) {
            chars.addAll( desc.getNodeData().getBinaryCharacters().getGainedCharacters() );
            chars.addAll( desc.getNodeData().getBinaryCharacters().getPresentCharacters() );
        }
        return chars;
    }
}

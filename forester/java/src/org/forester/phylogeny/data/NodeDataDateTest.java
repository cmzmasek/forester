// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
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
// Contact: czmasek at jcvi dot org

package org.forester.phylogeny.data;

import java.io.File;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * {@link NodeData#isHasDate()}: a date of exactly 0 is a date. Every contemporaneous BEAST tip sits at height 0, and
 * the present is 0 on any age scale; asking "is the number zero" where "is the number there" was meant dropped such
 * a date from every copy of the tree (undo, subtree) and from every saved phyloXML file. Tested where the loss
 * happened -- the copy and the written file -- and not only at the predicate.
 */
public final class NodeDataDateTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeData date of zero: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return predicateOk() && copyOk() && phyloXmlOk();
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    private static NodeData with( final Date d ) {
        final NodeData nd = new NodeData();
        nd.setDate( d );
        return nd;
    }

    private static boolean predicateOk() {
        // zero is a value, in each of the three numbers a date can state
        if ( !with( new Date( "", new BigDecimal( "0" ), null, null, "" ) ).isHasDate() ) {
            return fail( "a date VALUE of 0 is a date" );
        }
        if ( !with( new Date( "", new BigDecimal( "0.0" ), null, null, "" ) ).isHasDate() ) {
            return fail( "a date VALUE of 0.0 is a date" );
        }
        if ( !with( new Date( "", null, new BigDecimal( "0" ), null, "" ) ).isHasDate() ) {
            return fail( "a date MINIMUM of 0 is a date" );
        }
        if ( !with( new Date( "", null, null, new BigDecimal( "0" ), "" ) ).isHasDate() ) {
            return fail( "a date MAXIMUM of 0 is a date" );
        }
        if ( with( new Date( "", new BigDecimal( "0" ), null, null, "" ) ).isEmpty() ) {
            return fail( "node data stating a date of 0 is not empty" );
        }
        // what was a date before still is
        if ( !with( new Date( "", new BigDecimal( "2.1" ), null, null, "" ) ).isHasDate()
                || !with( new Date( "K-Pg", null, null, null, "" ) ).isHasDate()
                || !with( new Date( "", null, null, null, "mya" ) ).isHasDate() ) {
            return fail( "a value, a description or a unit alone is a date, as before" );
        }
        // deliberate non-behaviour: a date object that states NOTHING is still no date
        if ( new NodeData().isHasDate() ) {
            return fail( "no date object: no date" );
        }
        if ( with( new Date() ).isHasDate() || with( new Date( "" ) ).isHasDate() ) {
            return fail( "a date object stating nothing is no date" );
        }
        if ( !with( new Date() ).isEmpty() ) {
            return fail( "node data whose date states nothing is still empty" );
        }
        return true;
    }

    /** root (2.1) -> a (0.0), b (0); no unit, no description, no interval: the shape of a BEAST tree in heights. */
    private static Phylogeny heights() {
        final PhylogenyNode root = new PhylogenyNode();
        root.getNodeData().setDate( new Date( "", new BigDecimal( "2.1" ), null, null, "" ) );
        final PhylogenyNode a = new PhylogenyNode();
        a.setName( "a" );
        a.setDistanceToParent( 2.1 );
        a.getNodeData().setDate( new Date( "", new BigDecimal( "0.0" ), null, null, "" ) );
        final PhylogenyNode b = new PhylogenyNode();
        b.setName( "b" );
        b.setDistanceToParent( 2.1 );
        b.getNodeData().setDate( new Date( "", new BigDecimal( "0" ), null, null, "" ) );
        root.addAsChild( a );
        root.addAsChild( b );
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.setRooted( true );
        p.externalNodesHaveChanged();
        return p;
    }

    /** The date value of the node named {@code name} AS TEXT, or "none". */
    private static String dateOf( final Phylogeny p, final String name ) {
        final Date d = p.getNode( name ).getNodeData().getDate();
        return ( ( d == null ) || ( d.getValue() == null ) ) ? "none" : d.getValue().toPlainString();
    }

    private static boolean copyOk() {
        final NodeData copy = ( NodeData ) with( new Date( "", new BigDecimal( "0.0" ), null, null, "" ) ).copy();
        if ( ( copy.getDate() == null ) || ( copy.getDate().getValue() == null )
                || !"0.0".equals( copy.getDate().getValue().toPlainString() ) ) {
            return fail( "a copy of node data must keep a date of 0.0" );
        }
        final Phylogeny tree_copy = heights().copy();
        if ( !"0.0".equals( dateOf( tree_copy, "a" ) ) || !"0".equals( dateOf( tree_copy, "b" ) ) ) {
            return fail( "a copy of a tree must keep its tips' dates of 0; got a=" + dateOf( tree_copy, "a" ) + " b="
                    + dateOf( tree_copy, "b" ) );
        }
        return true;
    }

    private static boolean phyloXmlOk() throws Exception {
        final String xml = new PhylogenyWriter().toPhyloXML( heights(), 0 ).toString();
        // the written TEXT, not what a reader makes of it: three dates, two of them zero
        if ( count( xml, "<date>" ) != 3 ) {
            return fail( "all three dates must be written, got " + count( xml, "<date>" ) + " in\n" + xml );
        }
        if ( ( count( xml, "<value>0.0</value>" ) != 1 ) || ( count( xml, "<value>0</value>" ) != 1 ) ) {
            return fail( "a date of 0 must be written as it was stated, in\n" + xml );
        }
        final File f = File.createTempFile( "node_data_date", ".xml" );
        f.deleteOnExit();
        Files.write( f.toPath(), xml.getBytes( StandardCharsets.UTF_8 ) );
        final Phylogeny[] back = ParserBasedPhylogenyFactory.getInstance()
                .create( f, PhyloXmlParser.createPhyloXmlParserXsdValidating() );
        if ( ( back == null ) || ( back.length != 1 ) ) {
            return fail( "the written file must read back as one tree" );
        }
        if ( !"0.0".equals( dateOf( back[ 0 ], "a" ) ) || !"0".equals( dateOf( back[ 0 ], "b" ) ) ) {
            return fail( "a saved and reopened tree must keep its tips' dates of 0; got a=" + dateOf( back[ 0 ], "a" )
                    + " b=" + dateOf( back[ 0 ], "b" ) );
        }
        if ( !back[ 0 ].getNode( "a" ).getNodeData().isHasDate() ) {
            return fail( "the reopened tip must say it has a date" );
        }
        return true;
    }

    private static int count( final String s, final String what ) {
        int n = 0;
        for( int i = s.indexOf( what ); i >= 0; i = s.indexOf( what, i + what.length() ) ) {
            ++n;
        }
        return n;
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [NodeDataDateTest] " + m );
        return false;
    }

    private NodeDataDateTest() {
    }
}

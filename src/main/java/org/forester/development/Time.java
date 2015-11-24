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

package org.forester.development;

import java.io.File;
import java.io.IOException;
import java.util.Date;

import org.forester.io.parsers.nhx.NHXParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;
import org.forester.util.BasicDescriptiveStatistics;
import org.forester.util.DescriptiveStatistics;

@SuppressWarnings( "unused")
public final class Time {

    public static void main( final String[] args ) {
        try {
            final DescriptiveStatistics parse_stats = new BasicDescriptiveStatistics();
            final DescriptiveStatistics post_stats = new BasicDescriptiveStatistics();
            final DescriptiveStatistics pre_stats = new BasicDescriptiveStatistics();
            final File f = new File( args[ 0 ] );
            Phylogeny phy = null;
            for( int i = 0; i < 10; i++ ) {
                final long start_time = new Date().getTime();
                phy = ParserBasedPhylogenyFactory.getInstance().create( f, new NHXParser() )[ 0 ];
                System.out.println( phy.getNumberOfExternalNodes() );
                parse_stats.addValue( new Date().getTime() - start_time );
                //
                PhylogenyNode n = null;
                final long start_time_post = new Date().getTime();
                final PhylogenyNodeIterator post = phy.iteratorPostorder();
                while ( post.hasNext() ) {
                    n = post.next();
                }
                post_stats.addValue( new Date().getTime() - start_time_post );
                //
                final long start_time_pre = new Date().getTime();
                final PhylogenyNodeIterator pre = phy.iteratorPreorder();
                while ( pre.hasNext() ) {
                    n = pre.next();
                }
                pre_stats.addValue( new Date().getTime() - start_time_pre );
            }
            System.out.println( "Parsing [ms]:" );
            System.out.println( parse_stats.toString() );
            System.out.println( "Post-order [ms]:" );
            System.out.println( post_stats.toString() );
            System.out.println( "Pre-order [ms]:" );
            System.out.println( pre_stats.toString() );
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}

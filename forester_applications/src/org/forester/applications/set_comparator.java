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
// WWW: www.phylosoft.org/forester
// javac -cp ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// ~/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/org/forester/applications/set_comparator.java
// java -Xmx2048m -cp
// /home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester_applications/src/:/home/czmasek/SOFTWARE_DEV/ECLIPSE_WORKSPACE/forester/java/forester.jar
// org.forester.applications.set_comparator

package org.forester.applications;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.forester.util.ForesterUtil;

public class set_comparator {

    public static void main( final String args[] ) {
        try {
            if ( args.length != 2 ) {
                System.out.println( "Usage: set_comparator <set A> <set B>" );
                System.exit( -1 );
            }
            Set<String> set_a = ForesterUtil.file2set( new File( args[ 0 ] ) );
            final Set<String> set_b = ForesterUtil.file2set( new File( args[ 1 ] ) );
            System.out.println( "# A SIZE: " + set_a.size() );
            System.out.println( "# B SIZE: " + set_b.size() );
            set_a.retainAll( set_b );
            System.out.println( "# INTERSECTION (" + set_a.size() + "):" );
            for( final String s : set_a ) {
                System.out.println( s );
            }
            set_a = ForesterUtil.file2set( new File( args[ 0 ] ) );
            System.out.println();
            set_a.removeAll( set_b );
            System.out.println( "# A ONLY (" + set_a.size() + "):" );
            for( final String s : set_a ) {
                System.out.println( s );
            }
            set_a = ForesterUtil.file2set( new File( args[ 0 ] ) );
            System.out.println();
            set_b.removeAll( set_a );
            System.out.println( "# B ONLY (" + set_b.size() + "):" );
            for( final String s : set_b ) {
                System.out.println( s );
            }
        }
        catch ( final IOException e ) {
            e.printStackTrace();
            System.exit( -1 );
        }
    }
}
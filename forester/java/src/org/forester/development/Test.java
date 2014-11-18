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
import java.util.Date;
import java.util.Locale;

import org.forester.util.ForesterUtil;

/*
 * *
 */
public class Test {

    private final static String PATH_TO_TEST_DATA = System.getProperty( "user.dir" ) + ForesterUtil.getFileSeparator()
            + "test_data" + ForesterUtil.getFileSeparator();

    public static void main( final String[] args ) {
        System.out.println( "[Java version: " + ForesterUtil.JAVA_VERSION + " " + ForesterUtil.JAVA_VENDOR + "]" );
        System.out.println( "[OS: " + ForesterUtil.OS_NAME + " " + ForesterUtil.OS_ARCH + " " + ForesterUtil.OS_VERSION
                            + "]" );
        Locale.setDefault( Locale.US );
        System.out.println( "[Locale: " + Locale.getDefault() + "]" );
        final int failed = 0;
        final int succeeded = 0;
        System.out.print( "[Test if directory with files for testing exists/is readable: " );
        if ( Test.testDir( PATH_TO_TEST_DATA ) ) {
            System.out.println( "OK.]" );
        }
        else {
            System.out.println( "could not find/read from directory \"" + PATH_TO_TEST_DATA + "\".]" );
            System.out.println( "Testing aborted." );
            System.exit( -1 );
        }
        final long start_time = new Date().getTime();
        System.out.println( "\nTime requirement:  " + ( new Date().getTime() - start_time ) + "ms." );
        System.out.println();
        System.out.println( "Successful tests: " + succeeded );
        System.out.println( "Failed     tests: " + failed );
        System.out.println();
        if ( failed < 1 ) {
            System.out.println( "OK." );
        }
        else {
            System.out.println( "Not OK." );
        }
    }

    private static boolean testDir( final String file ) {
        try {
            final File f = new File( file );
            if ( !f.exists() ) {
                return false;
            }
            if ( !f.isDirectory() ) {
                return false;
            }
            if ( !f.canRead() ) {
                return false;
            }
        }
        catch ( final Exception e ) {
            return false;
        }
        return true;
    }
}

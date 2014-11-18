// $Id:
// forester -- software libraries and applications
// for genomics and evolutionary biology research.
//
// Copyright (C) 2010 Christian M Zmasek
// Copyright (C) 2010 Sanford-Burnham Medical Research Institute
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

package org.forester.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class CommandProcessBuilder {

    public static Process execute( final List<String> command, final File working_dir ) throws InterruptedException,
    IOException {
        final ProcessBuilder builder = new ProcessBuilder( command );
        if ( working_dir != null ) {
            if ( !working_dir.exists() ) {
                throw new IllegalArgumentException( "directory [" + working_dir.getAbsolutePath() + "] does not exist" );
            }
            if ( !working_dir.isDirectory() ) {
                throw new IllegalArgumentException( "[" + working_dir.getAbsolutePath() + "] is not a directory" );
            }
            if ( !working_dir.canWrite() ) {
                throw new IllegalArgumentException( "cannot write to [" + working_dir.getAbsolutePath() + "]" );
            }
            builder.directory( working_dir );
        }
        final Process process = builder.start();
        return process;
    }

    public static void main( final String args[] ) {
        final List<String> command = new ArrayList<String>();
        command.add( System.getenv( "windir" ) + "\\system32\\" + "tree.com" );
        command.add( "/A" );
        Process p;
        System.out.println( "Directory : " + System.getenv( "temp" ) );
        try {
            p = CommandProcessBuilder.execute( command, new File( System.getenv( "temp" ) ) );
            final InputStream is = p.getInputStream();
            final InputStreamReader isr = new InputStreamReader( is );
            final BufferedReader br = new BufferedReader( isr );
            String line;
            while ( ( line = br.readLine() ) != null ) {
                System.out.println( line );
            }
            System.out.println( "OK." );
        }
        catch ( final InterruptedException e ) {
            e.printStackTrace();
        }
        catch ( final IOException e ) {
            e.printStackTrace();
        }
    }
}

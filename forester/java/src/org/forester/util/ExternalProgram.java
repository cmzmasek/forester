// forester -- software libraries and applications
// for evolutionary biology and genomics.
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: czmasek at jcvi dot org

package org.forester.util;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class ExternalProgram {

    public static boolean isExecuteableFile( final File path_to_cmd_f ) {
        if ( !path_to_cmd_f.exists() ) {
            return false;
        }
        else if ( path_to_cmd_f.isDirectory() ) {
            return false;
        }
        else if ( !path_to_cmd_f.canExecute() ) {
            return false;
        }
        return true;
    }
    private Process      _process;
    private final String _path_to_cmd;

    public ExternalProgram( final String path_to_cmd ) {
        final File path_to_cmd_f = new File( path_to_cmd );
        checkCmdFile( path_to_cmd_f );
        _path_to_cmd = path_to_cmd_f.getAbsolutePath();
    }

    private void checkCmdFile( final File path_to_cmd_f ) {
        if ( !path_to_cmd_f.exists() ) {
            throw new IllegalArgumentException( "[" + path_to_cmd_f.getAbsolutePath() + "] does not exist" );
        }
        else if ( path_to_cmd_f.isDirectory() ) {
            throw new IllegalArgumentException( "[" + path_to_cmd_f.getAbsolutePath() + "] is a directory" );
        }
        else if ( !path_to_cmd_f.canExecute() ) {
            throw new IllegalArgumentException( "[" + path_to_cmd_f.getAbsolutePath() + "] is not executeable" );
        }
    }

    public InputStream getErrorStream() {
        return getProcess().getErrorStream();
    }

    public InputStream getInputStream() {
        return getProcess().getInputStream();
    }

    public OutputStream getOutputStream() {
        return getProcess().getOutputStream();
    }

    private String getPathToCmd() {
        return _path_to_cmd;
    }

    private Process getProcess() {
        return _process;
    }

    public Process launch( final String[] opts ) throws IOException, InterruptedException {
        String[] cmd;
        if ( ( opts == null ) || ( opts.length < 1 ) ) {
            cmd = new String[ 1 ];
        }
        else {
            cmd = new String[ opts.length + 1 ];
            for( int i = 0; i < opts.length; i++ ) {
                cmd[ i + 1 ] = opts[ i ];
            }
        }
        cmd[ 0 ] = getPathToCmd();
        System.out.println();
        for( final String element : cmd ) {
            System.out.print( element + " " );
        }
        System.out.println();
        setProcess( Runtime.getRuntime().exec( cmd ) );
        return getProcess();
    }

    private void setProcess( final Process process ) {
        _process = process;
    }

    public int waitFor() {
        try {
            return getProcess().waitFor();
        }
        catch ( final InterruptedException e ) {
            // TODO Auto-generated catch block
            getProcess().destroy();
            e.printStackTrace();
            return -1;
        }
    }
}

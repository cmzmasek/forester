// $Id:
/**
 * This class is intended to be used with the SystemCommandExecutor class to let
 * users execute system commands from Java applications.
 *
 * This class is based on work that was shared in a JavaWorld article named
 * "When System.exec() won't". That article is available at this url:
 *
 * http://www.javaworld.com/javaworld/jw-12-2000/jw-1229-traps.html
 *
 * Documentation for this class is available at this URL:
 *
 * http://devdaily.com/java/java-processbuilder-process-system-exec
 *
 *
 * Copyright 2010 alvin j. alexander, devdaily.com.
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please ee the following page for the LGPL license:
 * http://www.gnu.org/licenses/lgpl.txt
 *
 */

package org.forester.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;

class ThreadedStreamHandler extends Thread {

    InputStream     inputStream;
    String          adminPassword;
    OutputStream    outputStream;
    PrintWriter     printWriter;
    StringBuilder   outputBuffer    = new StringBuilder( 65536 );
    private boolean sudoIsRequested = false;

    /**
     * A simple constructor for when the sudo command is not necessary.
     * This constructor will just run the command you provide, without
     * running sudo before the command, and without expecting a password.
     *
     * @param inputStream
     * @param streamType
     */
    ThreadedStreamHandler( final InputStream inputStream ) {
        this.inputStream = inputStream;
    }

    /**
     * Use this constructor when you want to invoke the 'sudo' command.
     * The outputStream must not be null. If it is, you'll regret it. :)
     *
     * TODO this currently hangs if the admin password given for the sudo command is wrong.
     *
     * @param inputStream
     * @param streamType
     * @param outputStream
     * @param adminPassword
     */
    ThreadedStreamHandler( final InputStream inputStream, final OutputStream outputStream, final String adminPassword ) {
        this.inputStream = inputStream;
        this.outputStream = outputStream;
        printWriter = new PrintWriter( outputStream );
        this.adminPassword = adminPassword;
        sudoIsRequested = true;
    }

    ThreadedStreamHandler( final InputStream inputStream, final OutputStream outputStream ) {
        this.inputStream = inputStream;
        this.outputStream = outputStream;
        printWriter = new PrintWriter( outputStream );
        sudoIsRequested = false;
    }

    private void doSleep( final long millis ) {
        try {
            Thread.sleep( millis );
        }
        catch ( final InterruptedException e ) {
            // ignore
        }
    }

    public StringBuilder getOutputBuffer() {
        return outputBuffer;
    }

    @Override
    public void run() {
        // on mac os x 10.5.x, when i run a 'sudo' command, i need to write
        // the admin password out immediately; that's why this code is
        // here.
        if ( sudoIsRequested ) {
            //doSleep(500);
            printWriter.println( adminPassword );
            printWriter.flush();
        }
        BufferedReader bufferedReader = null;
        final String newline = ForesterUtil.LINE_SEPARATOR;
        try {
            bufferedReader = new BufferedReader( new InputStreamReader( inputStream ) );
            String line = null;
            while ( ( line = bufferedReader.readLine() ) != null ) {
                // outputBuffer.append( line + "\n" ); // CMZ change
                outputBuffer.append( line );
                outputBuffer.append( newline );
            }
        }
        catch ( final IOException ioe ) {
            // TODO handle this better
            ioe.printStackTrace();
        }
        catch ( final Throwable t ) {
            // TODO handle this better
            t.printStackTrace();
        }
        finally {
            try {
                bufferedReader.close();
            }
            catch ( final IOException e ) {
                // ignore this one
            }
        }
    }
}


package org.forester.util;

import java.io.BufferedWriter;
import java.io.IOException;

public final class EasyWriter extends BufferedWriter {

    private final static String LINE_SEPARATOR = ForesterUtil.LINE_SEPARATOR;

    public EasyWriter( final BufferedWriter out ) {
        super( out );
    }

    public void println( final String s ) throws IOException {
        write( s );
        write( LINE_SEPARATOR );
    }

    public void println() throws IOException {
        write( LINE_SEPARATOR );
    }

    public void print( final String s ) throws IOException {
        write( s );
    }
}

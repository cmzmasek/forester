
package org.forester.ws.seqdb;

public class DatabaseTools {

    static String extract( final String target, final String a, final String b ) {
        final int i_a = target.indexOf( a );
        final int i_b = target.indexOf( b );
        if ( ( i_a < 0 ) || ( i_b < i_a ) ) {
            throw new IllegalArgumentException( "attempt to extract from [" + target + "] between [" + a + "] and ["
                    + b + "]" );
        }
        return target.substring( i_a + a.length(), i_b ).trim();
    }

    static String extract( final String target, final String a ) {
        final int i_a = target.indexOf( a );
        return target.substring( i_a + a.length() ).trim();
    }
}

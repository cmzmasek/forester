
package org.forester.development;

import java.util.Random;

public class RandomThing {

    public static void main( final String[] args ) {
        final Random rg = new Random();
        for( int i = 0; i < 200; i++ ) {
            for( int j = 0; j < 30; j++ ) {
                System.out.print( "\t" + rg.nextFloat() );
            }
            System.out.println();
        }
    }
}

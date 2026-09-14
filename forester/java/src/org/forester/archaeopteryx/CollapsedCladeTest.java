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

package org.forester.archaeopteryx;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * The collapsed-clade rules shared with Archaeopteryx.js (archaeopteryx.js 407f204): rows and layout weight, wedge
 * height, the name (node name, else a Color-by value 95% of the tips share, else the common name prefix), the label
 * (" · N tips", " [found/total]"), the dominant colour and the fully-marked state. Headless.
 */
public final class CollapsedCladeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CollapsedClade: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            rowsAndWeight();
            height();
            dominantValue();
            prefixName();
            nameOrder();
            label();
            mostFrequent();
            fullyMarked();
            return true;
        }
        catch ( final AssertionError e ) {
            System.out.println( "  [CollapsedCladeTest] " + e.getMessage() );
            return false;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void rowsAndWeight() {
        // 1 + log2(tips) / 4: a pair takes 1.25 rows, 16 tips 2, 64 tips reach the 2.5 cap, and stay there
        ck( near( CollapsedClade.rows( 2 ), 1.25 ), "a pair takes 1.25 rows, got " + CollapsedClade.rows( 2 ) );
        ck( near( CollapsedClade.rows( 1 ), 1.25 ), "a single tip counts as a pair" );
        ck( near( CollapsedClade.rows( 0 ), 1.25 ), "no tips counts as a pair" );
        ck( near( CollapsedClade.rows( 16 ), 2.0 ), "16 tips take 2 rows, got " + CollapsedClade.rows( 16 ) );
        ck( near( CollapsedClade.rows( 64 ), 2.5 ), "64 tips reach the cap of 2.5 rows" );
        ck( near( CollapsedClade.rows( 100000 ), 2.5 ), "the rows never pass 2.5" );
        ck( CollapsedClade.rows( 5 ) > CollapsedClade.rows( 4 ), "more tips, more rows (below the cap)" );
        // the layout weight against a tip's 1 is 2 rows - 1, so the gap to a neighbour tip is (1 + weight) / 2 rows
        ck( near( CollapsedClade.rowWeight( 2 ), 1.5 ), "a pair weighs 1.5, got " + CollapsedClade.rowWeight( 2 ) );
        ck( near( CollapsedClade.rowWeight( 16 ), 3.0 ), "16 tips weigh 3" );
        ck( near( CollapsedClade.rowWeight( 1000 ), 4.0 ), "the weight never passes 4" );
    }

    private static void height() {
        // rows x row unit x 0.82, never under 6 px
        ck( near( CollapsedClade.height( 16, 20 ), 2.0 * 20 * 0.82 ), "16 tips at a 20 px row: 32.8 px, got "
                + CollapsedClade.height( 16, 20 ) );
        ck( near( CollapsedClade.height( 2, 1 ), 6 ), "a thin row still draws 6 px" );
        ck( near( CollapsedClade.height( 64, 10 ), 2.5 * 10 * 0.82 ), "the cap applies to the height too" );
    }

    private static void dominantValue() {
        final List<String> v = new ArrayList<>();
        for( int i = 0; i < 19; ++i ) {
            v.add( "Bovine" );
        }
        v.add( "Swine" );
        ck( "Bovine".equals( CollapsedClade.dominantValue( v, 20 ) ), "19 of 20 is 95%: the value names the clade" );
        v.set( 0, "Swine" );
        ck( CollapsedClade.dominantValue( v, 20 ) == null, "18 of 20 is 90%: no name" );
        // a tip WITHOUT a value still counts against the share
        final List<String> w = new ArrayList<>( Collections.nCopies( 19, "Bovine" ) );
        w.add( null );
        ck( "Bovine".equals( CollapsedClade.dominantValue( w, 20 ) ), "19 valued of 20 tips is 95%" );
        final List<String> x = new ArrayList<>( Collections.nCopies( 18, "Bovine" ) );
        x.add( null );
        x.add( null );
        ck( CollapsedClade.dominantValue( x, 20 ) == null, "the unvalued tips are in the denominator: 18 of 20 is 90%" );
        ck( CollapsedClade.dominantValue( Arrays.asList( null, null ), 2 ) == null, "no values, no name" );
        ck( CollapsedClade.dominantValue( null, 3 ) == null, "no Color-by, no name" );
    }

    private static void prefixName() {
        ck( "SARS_CoV_2/human".equals( CollapsedClade.prefixName( "SARS_CoV_2/human/" ) ), "a trailing separator goes" );
        ck( "Influenza A virus".equals( CollapsedClade.prefixName( "Influenza A virus " ) ), "trailing space goes" );
        ck( "abc".equals( CollapsedClade.prefixName( "abc_-.:|/ " ) ), "every trailing separator kind goes" );
        ck( "a_b".equals( CollapsedClade.prefixName( "a_b" ) ), "an inner separator stays" );
        ck( "".equals( CollapsedClade.prefixName( "A_" ) ), "under 2 characters left: no name" );
        ck( "AB".equals( CollapsedClade.prefixName( "AB|" ) ), "2 characters left is enough" );
        ck( "".equals( CollapsedClade.prefixName( "" ) ) && "".equals( CollapsedClade.prefixName( null ) ), "no prefix" );
    }

    private static void nameOrder() {
        ck( "Clade A".equals( CollapsedClade.name( " Clade A ", "Bovine", "SARS_" ) ), "the node's own name first, trimmed" );
        ck( "Bovine".equals( CollapsedClade.name( "  ", "Bovine", "SARS_CoV_2" ) ),
            "a blank node name falls to the Color-by value" );
        ck( "SARS_CoV_2".equals( CollapsedClade.name( null, null, "SARS_CoV_2/" ) ), "then the prefix, trimmed" );
        ck( "".equals( CollapsedClade.name( null, null, "" ) ), "nothing names it" );
    }

    private static void label() {
        ck( "Bovine · 12 tips".equals( CollapsedClade.label( "Bovine", 12, 0 ) ), "name, middle dot, tip count: "
                + CollapsedClade.label( "Bovine", 12, 0 ) );
        ck( "12 tips".equals( CollapsedClade.label( "", 12, 0 ) ), "no name: just the count" );
        ck( "1 tip".equals( CollapsedClade.label( null, 1, 0 ) ), "one tip is singular" );
        ck( "12 tips [3/12]".equals( CollapsedClade.label( "", 12, 3 ) ), "hits inside add [found/total]" );
        ck( "X · 2 tips [2/2]".equals( CollapsedClade.label( "X", 2, 2 ) ), "all found" );
    }

    private static void mostFrequent() {
        final Color r = Color.RED, b = Color.BLUE;
        ck( b.equals( CollapsedClade.mostFrequent( Arrays.asList( r, b, b, null, null, null ) ) ),
            "the most frequent vote wins, nulls cast none" );
        ck( r.equals( CollapsedClade.mostFrequent( Arrays.asList( r, r, b, b ) ) ),
            "a tie goes to the vote that reached the top count first (red reached 2 first)" );
        ck( b.equals( CollapsedClade.mostFrequent( Arrays.asList( r, b, b, r ) ) ),
            "...and a later tie does not take it over (blue reached 2 first; the JS running count keeps it)" );
        ck( b.equals( CollapsedClade.mostFrequent( Arrays.asList( r, b, b, r, b ) ) ), "a later majority takes over" );
        ck( CollapsedClade.mostFrequent( Arrays.asList( null, null ) ) == null, "no votes, no colour" );
    }

    private static void fullyMarked() {
        ck( CollapsedClade.fullyMarked( 4, 4 ), "every tip found is fully marked" );
        ck( !CollapsedClade.fullyMarked( 3, 4 ), "some tips found is not" );
        ck( !CollapsedClade.fullyMarked( 0, 0 ), "an empty clade is not" );
    }

    private static boolean near( final double a, final double b ) {
        return Math.abs( a - b ) < 1e-9;
    }

    private static void ck( final boolean cond, final String msg ) {
        if ( !cond ) {
            throw new AssertionError( msg );
        }
    }

    private CollapsedCladeTest() {
    }
}

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

package org.forester.io.parsers.json;

import java.io.IOException;
import java.util.List;
import java.util.Map;

/** Headless tests for {@link JsonParser}: the value model, escapes, numbers, whitespace, order, and malformation. */
public final class JsonParserTest {

    public static void main( final String[] args ) {
        System.out.println( "JsonParser: " + ( test() ? "OK." : "FAILED." ) );
    }

    @SuppressWarnings( "unchecked" )
    public static boolean test() {
        try {
            // primitives
            if ( !Double.valueOf( 42 ).equals( JsonParser.parse( "42" ) ) ) {
                return fail( "integer number" );
            }
            if ( !Double.valueOf( -1.5e3 ).equals( JsonParser.parse( "  -1.5e3 " ) ) ) {
                return fail( "signed exponent number + surrounding whitespace" );
            }
            if ( !"hi".equals( JsonParser.parse( "\"hi\"" ) ) ) {
                return fail( "string" );
            }
            if ( !Boolean.TRUE.equals( JsonParser.parse( "true" ) ) || !Boolean.FALSE.equals( JsonParser.parse( "false" ) ) ) {
                return fail( "booleans" );
            }
            if ( JsonParser.parse( "null" ) != null ) {
                return fail( "null" );
            }
            // string escapes, incl. \\uXXXX
            if ( !"a\"b\\c/d\n\te\u00e9".equals( JsonParser.parse( "\"a\\\"b\\\\c\\/d\\n\\te\\u00e9\"" ) ) ) {
                return fail( "string escapes" );
            }
            // object: insertion order preserved, mixed value types, nesting
            final Map<String, Object> o = (Map<String, Object>) JsonParser
                    .parse( "{ \"z\": 1, \"a\": [true, null, \"x\"], \"n\": {\"k\": 2.0} }" );
            if ( o.size() != 3 ) {
                return fail( "object size" );
            }
            if ( !new java.util.ArrayList<>( o.keySet() ).equals( java.util.Arrays.asList( "z", "a", "n" ) ) ) {
                return fail( "object must preserve key insertion order" );
            }
            final List<Object> a = (List<Object>) o.get( "a" );
            if ( ( a.size() != 3 ) || !Boolean.TRUE.equals( a.get( 0 ) ) || ( a.get( 1 ) != null )
                    || !"x".equals( a.get( 2 ) ) ) {
                return fail( "nested array" );
            }
            if ( !Double.valueOf( 2.0 ).equals( ( (Map<String, Object>) o.get( "n" ) ).get( "k" ) ) ) {
                return fail( "nested object" );
            }
            // empty object / array
            if ( !( (Map<?, ?>) JsonParser.parse( "{}" ) ).isEmpty() || !( (List<?>) JsonParser.parse( "[ ]" ) ).isEmpty() ) {
                return fail( "empty object/array" );
            }
            // malformation must throw IOException, never a runtime exception
            final String[] bad = { "", "  ", "{", "[1,]", "{\"a\":}", "{\"a\" 1}", "\"unterminated", "truX",
                    "1 2", "{\"a\":1}x", "01a", "\"\\x\"", "NaN", "Infinity", "[1,2", "{,}" };
            for ( final String b : bad ) {
                boolean threw = false;
                try {
                    JsonParser.parse( b );
                }
                catch ( final IOException e ) {
                    threw = true;
                }
                if ( !threw ) {
                    return fail( "malformed JSON must throw: <" + b + ">" );
                }
            }
            // a null argument throws (not NPE)
            try {
                JsonParser.parse( null );
                return fail( "null input must throw IOException" );
            }
            catch ( final IOException expected ) {
                // ok
            }
            return true;
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [JsonParserTest] " + msg );
        return false;
    }

    private JsonParserTest() {
    }
}

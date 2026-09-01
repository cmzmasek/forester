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

/**
 * Escaping for text stored inside a phyloXML {@code <property>} value.
 * <p>
 * Two hazards make this necessary. {@code XmlElement.getValueAsString()} COLLAPSES whitespace runs to a single
 * space and trims on reload, so a value containing a double, leading or trailing space comes back corrupted; and
 * anything used as a field separator has to be unambiguous inside the fields themselves.
 * <p>
 * The escape set is the annotation-import profile's, which proved the approach, PLUS {@code =} -- the figure spec
 * stores {@code key=value} pairs, where the profile uses fixed tab-separated fields. That extra escape is why the
 * profile is deliberately NOT migrated onto this class: doing so would change an already-shipped encoding for no
 * benefit. Anything NEW that stores text in a property should use this one.
 */
final class PropertyTextCodec {

    /** Escapes the separators AND every whitespace character (including the plain space). */
    static String esc( final String s ) {
        if ( s == null ) {
            return "";
        }
        final StringBuilder b = new StringBuilder();
        for( int i = 0; i < s.length(); i++ ) {
            final char c = s.charAt( i );
            switch ( c ) {
                case '\\': b.append( "\\\\" ); break;
                case '\t': b.append( "\\t" ); break;
                case '\n': b.append( "\\n" ); break;
                case '\r': b.append( "\\r" ); break;
                case ' ':  b.append( "\\_" ); break;
                case '|':  b.append( "\\p" ); break;
                case '~':  b.append( "\\s" ); break;
                case ';':  b.append( "\\c" ); break;
                case '=':  b.append( "\\e" ); break;
                default:   b.append( c );
            }
        }
        return b.toString();
    }

    /** The inverse of {@link #esc}. An unknown escape is kept verbatim rather than dropped. */
    static String unesc( final String s ) {
        if ( s == null ) {
            return "";
        }
        final StringBuilder b = new StringBuilder();
        for( int i = 0; i < s.length(); i++ ) {
            final char c = s.charAt( i );
            if ( ( c == '\\' ) && ( ( i + 1 ) < s.length() ) ) {
                final char n = s.charAt( ++i );
                switch ( n ) {
                    case 't':  b.append( '\t' ); break;
                    case 'n':  b.append( '\n' ); break;
                    case 'r':  b.append( '\r' ); break;
                    case '_':  b.append( ' ' ); break;
                    case 'p':  b.append( '|' ); break;
                    case 's':  b.append( '~' ); break;
                    case 'c':  b.append( ';' ); break;
                    case 'e':  b.append( '=' ); break;
                    case '\\': b.append( '\\' ); break;
                    default:   b.append( '\\' ).append( n ); break;
                }
            }
            else {
                b.append( c );
            }
        }
        return b.toString();
    }

    private PropertyTextCodec() {
    }
}

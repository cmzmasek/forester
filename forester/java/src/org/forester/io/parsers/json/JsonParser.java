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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * A minimal, dependency-free JSON parser (RFC 8259) producing a tree of plain Java values: a JSON object becomes a
 * {@link LinkedHashMap}&lt;String,Object&gt; (insertion order preserved), an array an {@link ArrayList}&lt;Object&gt;,
 * a string a {@link String}, a number a {@link Double}, {@code true}/{@code false} a {@link Boolean}, and {@code null}
 * a Java {@code null}. forester hand-writes its parsers and carries no JSON dependency; this is the first JSON forester
 * reads (for Auspice / Nextstrain {@code dataset.json}). Read-only -- not a serializer.
 */
public final class JsonParser {

    private static final int MAX_DEPTH = 10000; // guard against a StackOverflow on pathologically deep nesting

    private final String     _s;
    private int              _pos;
    private int              _depth;

    private JsonParser( final String s ) {
        _s = s;
        _pos = 0;
        _depth = 0;
    }

    /** Parse a complete JSON document to a tree of plain Java values (see the class doc). Throws {@link IOException}
     *  on any malformation (with the position), never a runtime exception. */
    public static Object parse( final String json ) throws IOException {
        if ( json == null ) {
            throw new IOException( "JSON input is null" );
        }
        // a UTF-8 BOM (U+FEFF) prepended by some editors/exporters is not whitespace -- strip a leading one
        final boolean bom = !json.isEmpty() && ( json.charAt( 0 ) == 0xFEFF );
        final JsonParser p = new JsonParser( bom ? json.substring( 1 ) : json );
        p.skipWhitespace();
        final Object value = p.parseValue();
        p.skipWhitespace();
        if ( p._pos != p._s.length() ) {
            throw p.error( "unexpected trailing content" );
        }
        return value;
    }

    private Object parseValue() throws IOException {
        if ( _pos >= _s.length() ) {
            throw error( "unexpected end of input" );
        }
        final char c = _s.charAt( _pos );
        switch ( c ) {
            case '{':
                return parseObject();
            case '[':
                return parseArray();
            case '"':
                return parseString();
            case 't':
                return parseLiteral( "true", Boolean.TRUE );
            case 'f':
                return parseLiteral( "false", Boolean.FALSE );
            case 'n':
                return parseLiteral( "null", null );
            default:
                if ( ( c == '-' ) || ( ( c >= '0' ) && ( c <= '9' ) ) ) {
                    return parseNumber();
                }
                throw error( "unexpected character '" + c + "'" );
        }
    }

    private Map<String, Object> parseObject() throws IOException {
        enter();
        final Map<String, Object> obj = new LinkedHashMap<>();
        _pos++; // consume '{'
        skipWhitespace();
        if ( peek() == '}' ) {
            _pos++;
            leave();
            return obj;
        }
        while ( true ) {
            skipWhitespace();
            if ( peek() != '"' ) {
                throw error( "expected a string key" );
            }
            final String key = parseString();
            skipWhitespace();
            if ( peek() != ':' ) {
                throw error( "expected ':' after a key" );
            }
            _pos++; // consume ':'
            skipWhitespace();
            obj.put( key, parseValue() );
            skipWhitespace();
            final char c = peek();
            if ( c == ',' ) {
                _pos++;
                continue;
            }
            if ( c == '}' ) {
                _pos++;
                leave();
                return obj;
            }
            throw error( "expected ',' or '}'" );
        }
    }

    private List<Object> parseArray() throws IOException {
        enter();
        final List<Object> arr = new ArrayList<>();
        _pos++; // consume '['
        skipWhitespace();
        if ( peek() == ']' ) {
            _pos++;
            leave();
            return arr;
        }
        while ( true ) {
            skipWhitespace();
            arr.add( parseValue() );
            skipWhitespace();
            final char c = peek();
            if ( c == ',' ) {
                _pos++;
                continue;
            }
            if ( c == ']' ) {
                _pos++;
                leave();
                return arr;
            }
            throw error( "expected ',' or ']'" );
        }
    }

    private String parseString() throws IOException {
        final StringBuilder sb = new StringBuilder();
        _pos++; // consume opening '"'
        while ( true ) {
            if ( _pos >= _s.length() ) {
                throw error( "unterminated string" );
            }
            final char c = _s.charAt( _pos++ );
            if ( c == '"' ) {
                return sb.toString();
            }
            if ( c == '\\' ) {
                if ( _pos >= _s.length() ) {
                    throw error( "unterminated escape" );
                }
                final char e = _s.charAt( _pos++ );
                switch ( e ) {
                    case '"':
                        sb.append( '"' );
                        break;
                    case '\\':
                        sb.append( '\\' );
                        break;
                    case '/':
                        sb.append( '/' );
                        break;
                    case 'b':
                        sb.append( '\b' );
                        break;
                    case 'f':
                        sb.append( '\f' );
                        break;
                    case 'n':
                        sb.append( '\n' );
                        break;
                    case 'r':
                        sb.append( '\r' );
                        break;
                    case 't':
                        sb.append( '\t' );
                        break;
                    case 'u':
                        if ( ( _pos + 4 ) > _s.length() ) {
                            throw error( "truncated \\u escape" );
                        }
                        try {
                            sb.append( (char) Integer.parseInt( _s.substring( _pos, _pos + 4 ), 16 ) );
                        }
                        catch ( final NumberFormatException nfe ) {
                            throw error( "invalid \\u escape" );
                        }
                        _pos += 4;
                        break;
                    default:
                        throw error( "invalid escape '\\" + e + "'" );
                }
            }
            else if ( c < 0x20 ) {
                throw error( "unescaped control character in a string" );
            }
            else {
                sb.append( c );
            }
        }
    }

    private Double parseNumber() throws IOException {
        final int start = _pos;
        if ( peek() == '-' ) {
            _pos++;
        }
        while ( _pos < _s.length() ) {
            final char c = _s.charAt( _pos );
            if ( ( ( c >= '0' ) && ( c <= '9' ) ) || ( c == '.' ) || ( c == 'e' ) || ( c == 'E' ) || ( c == '+' )
                    || ( c == '-' ) ) {
                _pos++;
            }
            else {
                break;
            }
        }
        final String num = _s.substring( start, _pos );
        try {
            final double d = Double.parseDouble( num );
            if ( Double.isNaN( d ) || Double.isInfinite( d ) ) {
                throw error( "non-finite number '" + num + "'" ); // JSON has no NaN/Infinity
            }
            return Double.valueOf( d );
        }
        catch ( final NumberFormatException nfe ) {
            throw error( "invalid number '" + num + "'" );
        }
    }

    private Object parseLiteral( final String literal, final Object value ) throws IOException {
        if ( !_s.startsWith( literal, _pos ) ) {
            throw error( "invalid literal (expected '" + literal + "')" );
        }
        _pos += literal.length();
        return value;
    }

    private char peek() {
        return ( _pos < _s.length() ) ? _s.charAt( _pos ) : '\0';
    }

    private void skipWhitespace() {
        while ( _pos < _s.length() ) {
            final char c = _s.charAt( _pos );
            if ( ( c == ' ' ) || ( c == '\t' ) || ( c == '\n' ) || ( c == '\r' ) ) {
                _pos++;
            }
            else {
                break;
            }
        }
    }

    private void enter() throws IOException {
        if ( ++_depth > MAX_DEPTH ) {
            throw error( "JSON nesting too deep" );
        }
    }

    private void leave() {
        _depth--;
    }

    private IOException error( final String msg ) {
        return new IOException( "JSON parse error at position " + _pos + ": " + msg );
    }
}

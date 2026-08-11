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

package org.forester.io.parsers.nhx;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.util.ForesterUtil;

/**
 * Parses a BEAST / BEAST X / BEAST 2 / TreeAnnotator FigTree-style {@code [&key=value,...]} annotation blob
 * (the content between {@code [} and {@code ]}) onto a {@link PhylogenyNode}, mapping each field to the phyloXML
 * structure the app's existing display features already consume:
 * <ul>
 * <li>{@code posterior} &rarr; a {@link org.forester.phylogeny.data.Confidence} of type "posterior" (Confidence
 *     Values display + support colour/symbol);</li>
 * <li>node age {@code height}/{@code height_mean}/{@code height_median} + {@code height_95%_HPD={lo,hi}} (or
 *     {@code height_range}) + {@code date} &rarr; a {@link Date} with value/min/max/desc (the Node Age Bars (HPD)
 *     feature draws the interval);</li>
 * <li>every other field ({@code rate}, {@code length_*}, discrete traits, {@code location}, ...) &rarr; a node
 *     {@link Property} {@code beast:<key>} (numeric &rarr; {@code xsd:decimal} so Color-by / Size-by / Annotation
 *     Columns pick it up; otherwise {@code xsd:string}).</li>
 * </ul>
 * The branch length lives on the Newick {@code :length} and is left untouched. Every {@code key=value} field is
 * preserved -- recognised ones as native structures, the rest as {@code beast:*} node properties. Robust: a
 * malformed field is skipped, never aborting the parse.
 */
public final class BeastAnnotationParser {

    private static final String HPD_SUFFIX = "_95%_hpd";

    private BeastAnnotationParser() {
        // pure utility
    }

    /** Apply the annotation content {@code inner_blob} (e.g. {@code &posterior=0.99,height_95%_HPD={1.4,1.5}}) to
     *  {@code node}. A leading {@code &} is optional. No-op on null/empty input. */
    public static void apply( final String inner_blob, final PhylogenyNode node ) {
        if ( ( inner_blob == null ) || ( node == null ) ) {
            return;
        }
        String blob = inner_blob.trim();
        if ( blob.startsWith( "&" ) ) {
            blob = blob.substring( 1 );
        }
        if ( blob.isEmpty() ) {
            return;
        }
        String height_median = null;
        String height_mean = null;
        String height = null;
        String[] hpd = null;
        String[] range = null;
        String date_desc = null;
        for( final String token : splitTopLevel( blob ) ) {
            final int eq = token.indexOf( '=' );
            if ( eq <= 0 ) {
                continue; // not a key=value field
            }
            final String key = token.substring( 0, eq ).trim();
            // NHXParser's streaming pre-parse maps a ':' inside a quoted value to the BELL char (7); restore it
            final String value = stripQuotes( token.substring( eq + 1 ).trim() ).replace( '\u0007', ':' );
            if ( ForesterUtil.isEmpty( key ) || ForesterUtil.isEmpty( value ) ) {
                continue;
            }
            final String kl = key.toLowerCase( Locale.ROOT );
            try {
                if ( kl.equals( "posterior" ) ) {
                    final Double d = parseNumber( value );
                    if ( d != null ) {
                        PhylogenyMethods.setConfidence( node, d.doubleValue(), "posterior" );
                    }
                }
                else if ( kl.equals( "height_median" ) ) {
                    height_median = value;
                }
                else if ( kl.equals( "height_mean" ) ) {
                    height_mean = value;
                }
                else if ( kl.equals( "height" ) ) {
                    height = value;
                }
                else if ( kl.equals( "height" + HPD_SUFFIX ) ) {
                    hpd = parseInterval( value );
                }
                else if ( kl.equals( "height_range" ) ) {
                    range = parseInterval( value );
                }
                else if ( kl.equals( "date" ) ) {
                    date_desc = value;
                }
                else {
                    addProperty( node, key, value );
                }
            }
            catch ( final Exception e ) {
                // a single malformed field must not abort the whole parse
            }
        }
        applyDate( node, firstNonEmpty( height_median, height_mean, height ),
                   ( hpd != null ) ? hpd : range, date_desc );
    }

    /** Attach a {@link Date} (age point value + younger/older HPD bounds + optional calendar desc) when any age
     *  information was present. min = lower/younger bound, max = upper/older bound (what the HPD bars expect). */
    private static void applyDate( final PhylogenyNode node, final String value, final String[] interval,
                                   final String date_desc ) {
        final boolean has_interval = ( interval != null ) && ( interval[ 0 ] != null ) && ( interval[ 1 ] != null );
        if ( ForesterUtil.isEmpty( value ) && !has_interval && ForesterUtil.isEmpty( date_desc ) ) {
            return;
        }
        // parse each piece independently, so an unparseable point value never discards a valid {lo,hi} interval
        final BigDecimal v = toBigDecimal( value );
        final BigDecimal min = has_interval ? toBigDecimal( interval[ 0 ] ) : null;
        final BigDecimal max = has_interval ? toBigDecimal( interval[ 1 ] ) : null;
        if ( ( v == null ) && ( min == null ) && ( max == null ) && ForesterUtil.isEmpty( date_desc ) ) {
            return; // nothing usable survived
        }
        node.getNodeData().setDate( new Date( ForesterUtil.isEmpty( date_desc ) ? "" : date_desc, v, min, max, "" ) );
    }

    private static BigDecimal toBigDecimal( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            return null;
        }
        try {
            return new BigDecimal( s );
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
    }

    /** Add a {@code beast:<key>} node property; numeric values are typed {@code xsd:decimal} (so Color-by picks
     *  them up), everything else {@code xsd:string}. */
    private static void addProperty( final PhylogenyNode node, final String key, final String value ) {
        PropertiesList pl = node.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            node.getNodeData().setProperties( pl );
        }
        final String datatype = ( parseNumber( value ) != null ) ? "xsd:decimal" : "xsd:string";
        pl.addProperty( new Property( "beast:" + refKey( key ), value, "", datatype, AppliesTo.NODE ) );
    }

    /** Split on TOP-LEVEL commas only -- commas inside {@code {...}} sets or {@code "..."}/{@code '...'} quotes
     *  are NOT separators (so {@code height_95%_HPD={1.4,1.5}} stays one token). */
    static List<String> splitTopLevel( final String s ) {
        final List<String> out = new ArrayList<String>();
        int depth = 0;
        boolean in_dq = false;
        boolean in_sq = false;
        final StringBuilder cur = new StringBuilder();
        for( int i = 0; i < s.length(); i++ ) {
            final char c = s.charAt( i );
            if ( in_dq ) {
                if ( c == '"' ) {
                    in_dq = false;
                }
                cur.append( c );
            }
            else if ( in_sq ) {
                if ( c == '\'' ) {
                    in_sq = false;
                }
                cur.append( c );
            }
            else if ( c == '"' ) {
                in_dq = true;
                cur.append( c );
            }
            else if ( c == '\'' ) {
                in_sq = true;
                cur.append( c );
            }
            else if ( ( c == '{' ) || c == '[' ) {
                depth++;
                cur.append( c );
            }
            else if ( ( c == '}' ) || c == ']' ) {
                if ( depth > 0 ) {
                    depth--;
                }
                cur.append( c );
            }
            else if ( ( c == ',' ) && ( depth == 0 ) ) {
                out.add( cur.toString() );
                cur.setLength( 0 );
            }
            else {
                cur.append( c );
            }
        }
        if ( cur.length() > 0 ) {
            out.add( cur.toString() );
        }
        return out;
    }

    /** Parse a two-value BEAST set {@code {lo,hi}} (or {@code [lo,hi]}) into the two raw numeric strings, or null
     *  if it is not a two-number set. Raw strings preserve full precision (no double round-trip). */
    static String[] parseInterval( final String v ) {
        String s = v.trim();
        if ( ( s.length() < 3 ) || !( ( s.charAt( 0 ) == '{' ) || ( s.charAt( 0 ) == '[' ) ) ) {
            return null;
        }
        s = s.substring( 1, s.length() - 1 ); // drop the braces/brackets
        final List<String> parts = splitTopLevel( s );
        if ( parts.size() != 2 ) {
            return null;
        }
        final String lo = parts.get( 0 ).trim();
        final String hi = parts.get( 1 ).trim();
        if ( ( parseNumber( lo ) == null ) || ( parseNumber( hi ) == null ) ) {
            return null;
        }
        return new String[] { lo, hi };
    }

    /** {@link Double} value iff {@code v} parses as a finite number, else null. */
    static Double parseNumber( final String v ) {
        try {
            final double d = Double.parseDouble( v );
            return ( Double.isNaN( d ) || Double.isInfinite( d ) ) ? null : Double.valueOf( d );
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
    }

    private static String stripQuotes( final String v ) {
        if ( ( v.length() >= 2 )
                && ( ( ( v.charAt( 0 ) == '"' ) && ( v.charAt( v.length() - 1 ) == '"' ) )
                        || ( ( v.charAt( 0 ) == '\'' ) && ( v.charAt( v.length() - 1 ) == '\'' ) ) ) ) {
            return v.substring( 1, v.length() - 1 );
        }
        return v;
    }

    /** A property-ref-safe rendering of a BEAST key: keep letters/digits, collapse every other run to a single
     *  underscore, and drop a trailing one (so {@code rate_95%_HPD} yields a clean {@code beast:rate_95_HPD} ref). */
    private static String refKey( final String key ) {
        final StringBuilder sb = new StringBuilder( key.length() );
        boolean prev_underscore = false;
        for( int i = 0; i < key.length(); i++ ) {
            final char c = key.charAt( i );
            if ( Character.isLetterOrDigit( c ) ) {
                sb.append( c );
                prev_underscore = false;
            }
            else if ( !prev_underscore ) {
                sb.append( '_' );
                prev_underscore = true;
            }
        }
        if ( ( sb.length() > 0 ) && ( sb.charAt( sb.length() - 1 ) == '_' ) ) {
            sb.setLength( sb.length() - 1 );
        }
        return sb.toString();
    }

    private static String firstNonEmpty( final String... vals ) {
        for( final String v : vals ) {
            if ( !ForesterUtil.isEmpty( v ) ) {
                return v;
            }
        }
        return null;
    }
}

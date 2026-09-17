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

import java.awt.Color;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.forester.io.parsers.json.AuspiceJsonParser;
import org.forester.phylogeny.PhylogenyMethods;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.BranchColor;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.NodeVisualData;
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
 * <li>the Auspice / Nextstrain "download Nexus" vocabulary lands exactly where {@link AuspiceJsonParser} puts the
 *     same dataset's JSON: {@code num_date} &rarr; the {@link Date} value with unit {@code "year"} (the unit, not the
 *     value, is what derives the calendar axis) plus a numeric {@code nextstrain:num_date} property;
 *     {@code num_date_CI={lo,hi}} &rarr; the date's min/max; {@code div} &rarr; {@code nextstrain:div} (the Time | Div
 *     toggle). A {@code num_date} outranks every {@code height*} as the date value, and only it carries the unit.</li>
 * <li>TreeTime's {@code mutations} / {@code mcc} are always text ({@code xsd:string}), never a number.</li>
 * <li>MrBayes' {@code prob} (+ {@code prob_stddev}) &rarr; a confidence of type "posterior probability" with its
 *     standard deviation; {@code bootstrap} &rarr; a confidence of type "bootstrap"; FigTree's {@code !color} &rarr;
 *     the branch colour, written either as {@code #RRGGBB} or -- what FigTree really writes -- as {@code #} + Java's
 *     SIGNED {@code Color.getRGB()} int ({@code #-8381639}). These are read wherever they stand in the blob: a
 *     FigTree-coloured BEAST tree LEADS with {@code !color} and still carries its posterior, height and rates.</li>
 * </ul>
 * Keys are compared lower-cased and without {@code %} and {@code _}, so {@code height_95%_HPD},
 * {@code height_95%HPD} and {@code height95%HPD} are one key -- while {@code length_95%HPD}, a BRANCH-LENGTH interval,
 * stays a different one and never becomes a date's interval.
 * <p>
 * {@code beast:} here means "came from a bracket annotation", not "produced by BEAST".
 * <p>
 * The branch length lives on the Newick {@code :length} and is left untouched. Every {@code key=value} field is
 * preserved -- recognised ones as native structures, the rest as {@code beast:*} node properties. Robust: a
 * malformed field is skipped, never aborting the parse.
 */
public final class BeastAnnotationParser {

    /** The type MrBayes' {@code prob} has always been given here. */
    static final String         MRBAYES_CONFIDENCE_TYPE = "posterior probability";
    private static final String BEAST_PREFIX = "beast:";
    /** The unit of a {@code num_date}: decimal calendar years, as in {@link AuspiceJsonParser}. */
    static final String         YEAR_UNIT  = "year";

    private BeastAnnotationParser() {
        // pure utility
    }

    /** Apply the annotation content {@code inner_blob} (e.g. {@code &posterior=0.99,height_95%_HPD={1.4,1.5}}) to
     *  {@code node}. A leading {@code &} is optional. No-op on null/empty input. */
    public static void apply( final String inner_blob, final PhylogenyNode node ) {
        apply( inner_blob, node, false );
    }

    /** Apply the annotation a Nexus TAXLABELS entry carries ({@code 'NewYork_454'[&!color=#-8381639]}) to the tip of
     *  that name. It is read like any other blob, with ONE difference: FigTree writes a taxon's colour there, and a
     *  coloured TAXON is a coloured LABEL -- so {@code !color} becomes the tip's label (font) colour, not its branch
     *  colour (which FigTree writes into the tree string). */
    public static void applyTaxlabelAnnotation( final String inner_blob, final PhylogenyNode tip ) {
        apply( inner_blob, tip, true );
    }

    private static void apply( final String inner_blob, final PhylogenyNode node, final boolean from_taxlabel ) {
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
        String hpd_raw = null;
        String range_raw = null;
        String date_desc = null;
        String num_date = null;
        String num_date_ci_raw = null;
        String[] num_date_ci = null;
        Double prob = null;
        Double prob_stddev = null;
        for( final String token : splitTopLevel( blob ) ) {
            final int eq = token.indexOf( '=' );
            if ( eq <= 0 ) {
                continue; // not a key=value field
            }
            String key = token.substring( 0, eq ).trim();
            if ( key.startsWith( "&" ) ) {
                key = key.substring( 1 ).trim(); // "[&bootstrap=69,&!color=#FFFFFF]": some writers lead EVERY field with '&'
            }
            // NHXParser's streaming pre-parse maps a ':' inside a quoted value to the BELL char (7); restore it
            final String value = stripQuotes( token.substring( eq + 1 ).trim() ).replace( '\u0007', ':' )
                    .replace( NHXParser.BLOB_OPEN_BRACKET, '[' ).replace( NHXParser.BLOB_CLOSE_BRACKET, ']' );
            if ( ForesterUtil.isEmpty( key ) || ForesterUtil.isEmpty( value ) ) {
                continue;
            }
            final String kl = normalizedKey( key );
            try {
                if ( kl.equals( "posterior" ) ) {
                    final Double d = parseNumber( value );
                    if ( d != null ) {
                        PhylogenyMethods.setConfidence( node, d.doubleValue(), "posterior" );
                    }
                }
                else if ( kl.equals( "heightmedian" ) ) {
                    height_median = value;
                }
                else if ( kl.equals( "heightmean" ) ) {
                    height_mean = value;
                }
                else if ( kl.equals( "height" ) ) {
                    height = value;
                }
                else if ( kl.equals( "height95hpd" ) ) {
                    hpd = parseInterval( value );
                    hpd_raw = value;
                }
                else if ( kl.equals( "heightrange" ) ) {
                    range = parseInterval( value );
                    range_raw = value;
                }
                else if ( kl.equals( "date" ) ) {
                    date_desc = value;
                }
                else if ( kl.equals( "numdate" ) && ( parseNumber( value ) != null ) ) {
                    num_date = value;
                    addProperty( node, AuspiceJsonParser.PREFIX + "num_date", value, false );
                }
                else if ( kl.equals( "numdateci" ) && ( parseInterval( value ) != null ) ) {
                    num_date_ci = parseInterval( value );
                    num_date_ci_raw = value;
                }
                else if ( kl.equals( "div" ) && ( parseNumber( value ) != null ) ) {
                    addProperty( node, AuspiceJsonParser.PREFIX + "div", value, false );
                }
                else if ( kl.equals( "prob" ) && ( parseNumber( value ) != null ) ) {
                    prob = parseNumber( value );
                }
                else if ( kl.equals( "probstddev" ) && ( parseNumber( value ) != null ) ) {
                    prob_stddev = parseNumber( value );
                }
                else if ( ( kl.equals( "bootstrap" ) || kl.equals( "boot" ) ) && ( parseNumber( value ) != null )
                        && ( parseNumber( value ).doubleValue() >= 0 ) ) {
                    node.getBranchData().addConfidence( new Confidence( parseNumber( value ).doubleValue(),
                                                                        "bootstrap" ) );
                }
                else if ( ( kl.equals( "!color" ) || kl.equals( "!colour" ) || kl.equals( "color" )
                        || kl.equals( "colour" ) ) && ( parseColor( value ) != null ) ) {
                    if ( from_taxlabel ) {
                        if ( node.getNodeData().getNodeVisualData() == null ) {
                            node.getNodeData().setNodeVisualData( new NodeVisualData() );
                        }
                        node.getNodeData().getNodeVisualData().setFontColor( parseColor( value ) );
                    }
                    else {
                        node.getBranchData().setBranchColor( new BranchColor( parseColor( value ) ) );
                    }
                }
                else if ( kl.equals( "mutations" ) || kl.equals( "mcc" ) ) {
                    addProperty( node, BEAST_PREFIX + refKey( key ), value, true );
                }
                else {
                    // A "!" key is one of FigTree's display DIRECTIVES (!color, !rotate, !collapse, ...), never a
                    // measurement: kept, but as text -- a refused "!color=-8381639" typed as a decimal would be
                    // offered by Color-by as a numeric trait with a gradient of its own.
                    addProperty( node, BEAST_PREFIX + refKey( key ), value, key.startsWith( "!" ) );
                }
            }
            catch ( final Exception e ) {
                // a single malformed field must not abort the whole parse
            }
        }
        if ( ( prob != null ) && ( prob.doubleValue() >= 0 ) ) {
            node.getBranchData().addConfidence( ( ( prob_stddev != null ) && ( prob_stddev.doubleValue() >= 0 ) )
                    ? new Confidence( prob.doubleValue(), MRBAYES_CONFIDENCE_TYPE, prob_stddev.doubleValue() )
                    : new Confidence( prob.doubleValue(), MRBAYES_CONFIDENCE_TYPE ) );
        }
        else if ( prob_stddev != null ) {
            // a deviation with no probability to qualify: keep it as data rather than lose it
            addProperty( node, BEAST_PREFIX + "prob_stddev", String.valueOf( prob_stddev ), false );
        }
        if ( num_date != null ) {
            // a calendar date outranks every age, and brings its OWN interval: a height HPD is on the age scale
            applyDate( node, num_date, num_date_ci, date_desc, YEAR_UNIT );
            // the heights it outranked are kept as what they were written as, rather than dropped (no real file
            // carries both, so nothing here is lost to a guess)
            final String[][] outranked = { { "height_median", height_median }, { "height_mean", height_mean },
                    { "height", height }, { "height_95_HPD", hpd_raw }, { "height_range", range_raw } };
            for( final String[] o : outranked ) {
                if ( o[ 1 ] != null ) {
                    addProperty( node, BEAST_PREFIX + o[ 0 ], o[ 1 ], false );
                }
            }
        }
        else {
            if ( num_date_ci_raw != null ) {
                // an interval with no point date to hang on: keep it as data rather than lose it
                addProperty( node, AuspiceJsonParser.PREFIX + "num_date_CI", num_date_ci_raw, true );
            }
            applyDate( node, firstNonEmpty( height_median, height_mean, height ), ( hpd != null ) ? hpd : range,
                       date_desc, "" );
        }
    }

    /** Attach a {@link Date} (age point value + younger/older HPD bounds + optional calendar desc) when any age
     *  information was present. min = lower/younger bound, max = upper/older bound (what the HPD bars expect). */
    private static void applyDate( final PhylogenyNode node, final String value, final String[] interval,
                                   final String date_desc, final String unit ) {
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
        // the unit belongs to the point VALUE: without one that parsed there is nothing for it to describe
        node.getNodeData().setDate( new Date( ForesterUtil.isEmpty( date_desc ) ? "" : date_desc, v, min, max,
                                              ( v != null ) ? unit : "" ) );
    }

    /** A key as it is COMPARED: lower-cased, without {@code %} and {@code _} (the HPD keys are spelled three ways
     *  across BEAST 1 / BEAST 2 / MrBayes). The stored property ref keeps the key as written. */
    static String normalizedKey( final String key ) {
        return key.toLowerCase( Locale.ROOT ).replace( "%", "" ).replace( "_", "" );
    }

    /** A FigTree colour: {@code #RRGGBB}, or {@code #} + the signed int of Java's {@code Color.getRGB()} -- which is
     *  what FigTree writes ({@code #-8381639} is 0xFF801B39, i.e. opaque 0x801B39). A branch colour has no alpha, so
     *  the alpha byte is dropped. Null when it is neither. */
    static Color parseColor( final String v ) {
        if ( ( v == null ) || ( v.length() < 2 ) || ( v.charAt( 0 ) != '#' ) ) {
            return null;
        }
        final String body = v.substring( 1 );
        try {
            if ( body.matches( "[0-9a-fA-F]{6}" ) ) {
                return new Color( Integer.parseInt( body, 16 ) );
            }
            if ( body.matches( "-?[0-9]+" ) ) {
                return new Color( Integer.parseInt( body ) ); // new Color(int) is opaque: the alpha byte is dropped
            }
        }
        catch ( final NumberFormatException e ) {
            // beyond an int: not a colour
        }
        return null;
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

    /** Add a node property under the complete ref {@code ref}; numeric values are typed {@code xsd:decimal} (so
     *  Color-by picks them up) unless {@code always_text}, everything else {@code xsd:string}. */
    private static void addProperty( final PhylogenyNode node, final String ref, final String value,
                                     final boolean always_text ) {
        PropertiesList pl = node.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            node.getNodeData().setProperties( pl );
        }
        final String datatype = ( !always_text && ( parseNumber( value ) != null ) ) ? "xsd:decimal" : "xsd:string";
        pl.addProperty( new Property( ref, value, "", datatype, AppliesTo.NODE ) );
    }

    /** Split on TOP-LEVEL commas only -- commas inside {@code {...}} sets or inside a quoted VALUE are NOT separators
     *  (so {@code height_95%_HPD={1.4,1.5}} stays one token). A quote that does not open a value is itself data:
     *  {@code country=Côte d'Ivoire,region=Africa} is two fields. The same rule as NHXParser's streaming scanner, and
     *  JOINT with Archaeopteryx.js (splitTopLevelCommas / opensBlobQuote / blobQuoteClose). */
    static List<String> splitTopLevel( final String s ) {
        final List<String> out = new ArrayList<String>();
        int depth = 0;
        final StringBuilder cur = new StringBuilder();
        for( int i = 0; i < s.length(); i++ ) {
            final char c = s.charAt( i );
            if ( ( ( c == '"' ) || ( c == '\'' ) ) && opensBlobQuote( s, i ) ) {
                final int close = blobQuoteClose( s, i );
                if ( close > -1 ) {
                    cur.append( s, i, close + 1 ); // a quoted value, its commas data
                    i = close;
                    continue;
                }
            }
            if ( ( c == '{' ) || ( c == '[' ) ) {
                depth++;
                cur.append( c );
            }
            else if ( ( c == '}' ) || ( c == ']' ) ) {
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

    /** A quote opens a run only where a value can START: straight after '=', or after '{', '[' or ',' . */
    static boolean opensBlobQuote( final String s, final int p ) {
        int m = p - 1;
        while ( ( m >= 0 ) && Character.isWhitespace( s.charAt( m ) ) ) {
            --m;
        }
        return ( m >= 0 ) && ( "={[,".indexOf( s.charAt( m ) ) >= 0 );
    }

    /** The index of the quote closing the run opened at {@code p} -- the same character standing where a value can
     *  END (before ',', '}', ']' or the end) -- or -1. */
    static int blobQuoteClose( final String s, final int p ) {
        final char q = s.charAt( p );
        for( int k = p + 1; k < s.length(); ++k ) {
            if ( s.charAt( k ) != q ) {
                continue;
            }
            int m = k + 1;
            while ( ( m < s.length() ) && Character.isWhitespace( s.charAt( m ) ) ) {
                ++m;
            }
            if ( ( m >= s.length() ) || ( ",}]".indexOf( s.charAt( m ) ) >= 0 ) ) {
                return k;
            }
        }
        return -1;
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

    // A plain decimal with an optional exponent -- what Archaeopteryx.js reads as a number too. Java alone would
    // also take "3f", "1d" and hex floats, and type a clade called "3f" as xsd:decimal (invalid phyloXML).
    private static final java.util.regex.Pattern NUMBER_PATTERN = java.util.regex.Pattern
            .compile( "[-+]?(\\d+\\.?\\d*|\\.\\d+)([eE][-+]?\\d+)?" );

    /** {@link Double} value iff {@code v} parses as a finite number, else null. */
    static Double parseNumber( final String v ) {
        if ( ( v == null ) || !NUMBER_PATTERN.matcher( v ).matches() ) {
            return null;
        }
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

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

package org.forester.archaeopteryx.tools;

import java.math.BigDecimal;
import java.time.LocalDate;
import java.time.Year;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.util.ForesterUtil;

/**
 * Pulls a sampling DATE out of each external-node LABEL and (via the caller) writes it into the tip's phyloXML
 * {@code <date>} (value = decimal year, unit "year") plus a numeric {@code data:date} property -- so a tip-dated
 * Newick / phyloXML / Nexus tree (the common way molecular-epidemiology trees are shared: BEAST / TreeTime / augur /
 * GISAID Newick, where the date is baked into the strain name) lights up the same Calendar axis, Color-by-date and
 * time-tree machinery an Auspice dataset gets for free. This is the date sibling of {@link LabelDataExtractor}.
 *
 * <p>Recognised formats, tried most-specific first (the rightmost match within a format wins, since the date usually
 * trails the label):
 * <ol>
 * <li>ISO date {@code 2021-03-15};</li>
 * <li>month-name date {@code 15-Mar-2021} (the GenBank {@code collection_date} style);</li>
 * <li>numeric slash/dot date {@code 15/03/2021}, {@code 2021/03/15}, {@code 2021.03.15};</li>
 * <li>ISO year-month {@code 2021-03};</li>
 * <li>month-name year {@code Mar-2021};</li>
 * <li>decimal year {@code 2021.37};</li>
 * <li>a bare 4-digit year token {@code A/Texas/50/2012}.</li>
 * </ol>
 *
 * <p><b>Conventions</b> (chosen with the user): an incomplete date maps to the MIDPOINT of its interval ({@code 2021}
 * &rarr; 2021.5, {@code 2021-03} &rarr; mid-March), matching BEAST/TreeTime. An ambiguous numeric date (both non-year
 * fields &le; 12, e.g. {@code 05/03/2021}) is read <b>day-first</b> by default (rest-of-world), overridable per run.
 *
 * <p>The parsers are pure and the node-writer takes a single node, so the whole thing is unit-testable with no GUI.
 */
public final class TipDateExtractor {

    /** The numeric decimal-year property written alongside the {@code <date>}, so the sampling date can drive
     *  "Color by" (a date gradient), search, and annotation columns -- the generic twin of {@code nextstrain:num_date}. */
    public static final String DATE_PROPERTY_REF = "data:date";

    // plausible calendar-year window for a BARE year token / decimal year (guards against matching a random 4-digit
    // strain number); an explicit ISO / month-name / slash date is strong evidence and gets a looser sanity window
    private static final int   MIN_BARE_YEAR     = 1900;
    private static final int   MAX_BARE_YEAR     = 2100;
    private static final int   MIN_YEAR          = 1000;
    private static final int   MAX_YEAR          = 2999;

    private static final Pattern ISO_FULL      = Pattern.compile( "(?<![0-9])(\\d{4})-(\\d{1,2})-(\\d{1,2})(?![0-9])" );
    private static final Pattern MONTH_FULL    = Pattern
            .compile( "(?<![A-Za-z0-9])(\\d{1,2})[-\\s]([A-Za-z]{3,9})[-\\s](\\d{4})(?![0-9])" );
    private static final Pattern SLASH_FULL    = Pattern
            .compile( "(?<![0-9])(\\d{1,4})[/.](\\d{1,2})[/.](\\d{1,4})(?![0-9])" );
    private static final Pattern ISO_PARTIAL   = Pattern.compile( "(?<![0-9])(\\d{4})-(\\d{1,2})(?![0-9-])" );
    private static final Pattern MONTH_PARTIAL = Pattern.compile( "(?<![A-Za-z0-9])([A-Za-z]{3,9})[-\\s](\\d{4})(?![0-9])" );
    private static final Pattern DECIMAL_YEAR  = Pattern.compile( "(?<![0-9.])(\\d{4}\\.\\d+)(?![0-9])" );
    private static final Pattern BARE_YEAR     = Pattern.compile( "(?<![0-9.])(\\d{4})(?![0-9.])" );

    private static final String[] MONTHS       = { "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct",
            "nov", "dec" };

    /** The order of the two non-year fields in an ambiguous numeric date (e.g. {@code 05/03/2021}). */
    public enum DayMonthOrder {
        DAY_FIRST, MONTH_FIRST
    }

    /** How much of a date the label carried. */
    public enum Precision {
        YEAR, MONTH, DAY
    }

    /** A date recognised in a label: the decimal-year value, the exact substring matched, its precision, a human label
     *  of the format, and whether the numeric day/month order had to be guessed (ambiguous). */
    public record DateMatch(double decimalYear, String matchedText, Precision precision, String formatLabel,
            boolean ambiguous) {
    }

    /** One tip's preview row (drives the dialog table + counts): the node, its label, the date found (null = none), and
     *  whether the tip already had a {@code <date>} (so it can be skipped). */
    public record TipDate(PhylogenyNode node, String label, DateMatch match, boolean alreadyDated) {
    }

    /** Parse the first (most-specific, then rightmost) date out of a label; null if none is found. */
    public static DateMatch parse( final String label, final DayMonthOrder order ) {
        if ( ForesterUtil.isEmpty( label ) ) {
            return null;
        }
        DateMatch m;
        if ( ( m = matchIsoFull( label ) ) != null ) {
            return m;
        }
        if ( ( m = matchMonthFull( label ) ) != null ) {
            return m;
        }
        if ( ( m = matchSlashFull( label, order ) ) != null ) {
            return m;
        }
        if ( ( m = matchIsoPartial( label ) ) != null ) {
            return m;
        }
        if ( ( m = matchMonthPartial( label ) ) != null ) {
            return m;
        }
        if ( ( m = matchDecimalYear( label ) ) != null ) {
            return m;
        }
        return matchBareYear( label );
    }

    // ---- per-format matchers (each returns the RIGHTMOST valid match, or null) --------------------------------

    private static DateMatch matchIsoFull( final String s ) {
        final Matcher mat = ISO_FULL.matcher( s );
        DateMatch found = null;
        while ( mat.find() ) {
            final DateMatch m = fromYmd( Integer.parseInt( mat.group( 1 ) ), Integer.parseInt( mat.group( 2 ) ),
                                         Integer.parseInt( mat.group( 3 ) ), mat.group(), "ISO date (YYYY-MM-DD)",
                                         false );
            if ( m != null ) {
                found = m; // keep the last (rightmost)
            }
        }
        return found;
    }

    private static DateMatch matchMonthFull( final String s ) {
        final Matcher mat = MONTH_FULL.matcher( s );
        DateMatch found = null;
        while ( mat.find() ) {
            final int mon = monthNumber( mat.group( 2 ) );
            if ( mon > 0 ) {
                final DateMatch m = fromYmd( Integer.parseInt( mat.group( 3 ) ), mon, Integer.parseInt( mat.group( 1 ) ),
                                             mat.group(), "month-name date", false );
                if ( m != null ) {
                    found = m; // keep the last (rightmost)
                }
            }
        }
        return found;
    }

    private static DateMatch matchSlashFull( final String s, final DayMonthOrder order ) {
        final Matcher mat = SLASH_FULL.matcher( s );
        DateMatch found = null;
        while ( mat.find() ) {
            final int a = Integer.parseInt( mat.group( 1 ) );
            final int b = Integer.parseInt( mat.group( 2 ) );
            final int c = Integer.parseInt( mat.group( 3 ) );
            int year = -1, month = -1, day = -1;
            boolean ambiguous = false;
            if ( mat.group( 1 ).length() == 4 ) { // year-first (YYYY/MM/DD) -- unambiguous
                year = a;
                month = b;
                day = c;
            }
            else if ( mat.group( 3 ).length() == 4 ) { // year-last (D/M/YYYY or M/D/YYYY)
                year = c;
                if ( a > 12 ) {
                    day = a;
                    month = b;
                }
                else if ( b > 12 ) {
                    month = a;
                    day = b;
                }
                else { // both <= 12 -> genuinely ambiguous, use the chosen order
                    ambiguous = true;
                    if ( order == DayMonthOrder.MONTH_FIRST ) {
                        month = a;
                        day = b;
                    }
                    else {
                        day = a;
                        month = b;
                    }
                }
            }
            if ( year > 0 ) {
                final DateMatch m = fromYmd( year, month, day, mat.group(), "numeric date", ambiguous );
                if ( m != null ) {
                    found = m;
                }
            }
        }
        return found;
    }

    private static DateMatch matchIsoPartial( final String s ) {
        final Matcher mat = ISO_PARTIAL.matcher( s );
        DateMatch found = null;
        while ( mat.find() ) {
            final DateMatch m = fromYm( Integer.parseInt( mat.group( 1 ) ), Integer.parseInt( mat.group( 2 ) ),
                                        mat.group(), "ISO year-month (YYYY-MM)" );
            if ( m != null ) {
                found = m;
            }
        }
        return found;
    }

    private static DateMatch matchMonthPartial( final String s ) {
        final Matcher mat = MONTH_PARTIAL.matcher( s );
        DateMatch found = null;
        while ( mat.find() ) {
            final int mon = monthNumber( mat.group( 1 ) );
            if ( mon > 0 ) {
                final DateMatch m = fromYm( Integer.parseInt( mat.group( 2 ) ), mon, mat.group(), "month-name year" );
                if ( m != null ) {
                    found = m;
                }
            }
        }
        return found;
    }

    private static DateMatch matchDecimalYear( final String s ) {
        final Matcher mat = DECIMAL_YEAR.matcher( s );
        DateMatch found = null;
        while ( mat.find() ) {
            final double v = Double.parseDouble( mat.group( 1 ) );
            final int y = (int) Math.floor( v );
            if ( ( y >= MIN_BARE_YEAR ) && ( y <= MAX_BARE_YEAR ) ) {
                found = new DateMatch( v, mat.group( 1 ), Precision.DAY, "decimal year", false );
            }
        }
        return found;
    }

    private static DateMatch matchBareYear( final String s ) {
        final Matcher mat = BARE_YEAR.matcher( s );
        DateMatch found = null;
        while ( mat.find() ) {
            final int y = Integer.parseInt( mat.group( 1 ) );
            if ( ( y >= MIN_BARE_YEAR ) && ( y <= MAX_BARE_YEAR ) ) {
                found = new DateMatch( y + 0.5, mat.group( 1 ), Precision.YEAR, "year", false );
            }
        }
        return found;
    }

    // ---- decimal-year conversion (leap-year aware; incomplete -> interval midpoint) ---------------------------

    /** Decimal year for a (year, month, day) date; month or day &le; 0 means "not given" and maps to the interval
     *  midpoint (year -> .5; month -> mid-month). Leap-year aware. */
    public static double toDecimalYear( final int year, final int month, final int day ) {
        if ( month <= 0 ) {
            return year + 0.5; // year only -> midpoint of the year
        }
        if ( day <= 0 ) {
            final LocalDate first = LocalDate.of( year, month, 1 );
            final double mid_doy = ( first.getDayOfYear() - 1 ) + ( first.lengthOfMonth() / 2.0 ); // mid-month
            return year + ( mid_doy / Year.of( year ).length() );
        }
        final LocalDate d = LocalDate.of( year, month, day );
        return year + ( ( d.getDayOfYear() - 0.5 ) / Year.of( year ).length() ); // mid-day
    }

    // ---- whole-tree helpers -----------------------------------------------------------------------------------

    /** One preview row per external node: {@code {node, label, date-or-null, already-dated}}. Pure, GUI-free. */
    public static List<TipDate> preview( final Phylogeny phy, final DayMonthOrder order ) {
        final List<TipDate> rows = new ArrayList<>();
        if ( ( phy != null ) && !phy.isEmpty() ) {
            for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
                rows.add( new TipDate( ext, ext.getName(), parse( ext.getName(), order ), ext.getNodeData().isHasDate() ) );
            }
        }
        return rows;
    }

    /** A roll-up of a preview for the dialog's summary line: how many tips matched / did not, how many needed a
     *  day/month guess, the dominant format, and the matched date range (decimal years; 0..0 when none matched). */
    public record Summary(int matched, int unmatched, int ambiguous, String dominantFormat, double minYear,
            double maxYear) {
    }

    /** Roll a preview up into a {@link Summary} (pure; drives the dialog header). */
    public static Summary summarize( final List<TipDate> rows ) {
        int matched = 0, unmatched = 0, ambiguous = 0;
        final java.util.Map<String, Integer> formats = new java.util.LinkedHashMap<>();
        double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;
        for ( final TipDate t : rows ) {
            if ( t.match() == null ) {
                unmatched++;
                continue;
            }
            matched++;
            if ( t.match().ambiguous() ) {
                ambiguous++;
            }
            formats.merge( t.match().formatLabel(), 1, Integer::sum );
            min = Math.min( min, t.match().decimalYear() );
            max = Math.max( max, t.match().decimalYear() );
        }
        String dom = "";
        int best = -1;
        for ( final java.util.Map.Entry<String, Integer> e : formats.entrySet() ) {
            if ( e.getValue() > best ) {
                best = e.getValue();
                dom = e.getKey();
            }
        }
        return new Summary( matched, unmatched, ambiguous, dom, matched > 0 ? min : 0, matched > 0 ? max : 0 );
    }

    /** Whether a STRICT MAJORITY of the tips carry a parseable date -- the gate for the load-time auto-offer. */
    public static boolean mostLabelsHaveDates( final Phylogeny phy ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return false;
        }
        int tips = 0, dated = 0;
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            tips++;
            if ( parse( ext.getName(), DayMonthOrder.DAY_FIRST ) != null ) {
                dated++;
            }
        }
        return ( tips > 0 ) && ( ( dated * 2 ) > tips );
    }

    /** Whether any tip's date is a numeric date whose day/month order had to be guessed -- so the day/month toggle is
     *  worth showing (adaptive UI). */
    public static boolean hasAmbiguousDates( final Phylogeny phy, final DayMonthOrder order ) {
        if ( ( phy == null ) || phy.isEmpty() ) {
            return false;
        }
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            final DateMatch m = parse( ext.getName(), order );
            if ( ( m != null ) && m.ambiguous() ) {
                return true;
            }
        }
        return false;
    }

    /** Write a matched date onto a node: the {@code <date>} (value = rounded decimal year, unit "year") + the numeric
     *  {@code data:date} property (for Color-by-date). Returns false for a null match. */
    public static boolean applyToNode( final PhylogenyNode node, final DateMatch m ) {
        if ( ( node == null ) || ( m == null ) ) {
            return false;
        }
        final String year_str = decimalYearString( m.decimalYear() );
        node.getNodeData().setDate( new Date( "", new BigDecimal( year_str ), null, null, "year" ) );
        PropertiesList pl = node.getNodeData().getProperties();
        if ( pl == null ) {
            pl = new PropertiesList();
            node.getNodeData().setProperties( pl );
        }
        // do not duplicate the property on a re-run
        for ( final Property p : pl.getProperties() ) {
            if ( DATE_PROPERTY_REF.equals( p.getRef() ) ) {
                pl.getProperties().remove( p );
                break;
            }
        }
        pl.addProperty( new Property( DATE_PROPERTY_REF, year_str, "", "xsd:decimal", AppliesTo.NODE ) );
        return true;
    }

    /** A provenance sentence for the tree description (per the tree-mutation-provenance rule). */
    public static String provenanceSentence( final int matched, final int total, final String format ) {
        return "Extracted sampling dates from the tip labels (" + format + ") for " + matched + " of " + total
                + " tips.";
    }

    /** The rounded (5-decimal) decimal year as a plain decimal string (shared by the {@code <date>} value and the
     *  {@code data:date} property, so they always agree). */
    public static String decimalYearString( final double decimal_year ) {
        final double rounded = Math.round( decimal_year * 100000.0 ) / 100000.0;
        return BigDecimal.valueOf( rounded ).stripTrailingZeros().toPlainString();
    }

    // ---- internals --------------------------------------------------------------------------------------------

    /** Build a full-date match, validating the calendar date (returns null on an impossible date like 2021-13-40). */
    private static DateMatch fromYmd( final int year, final int month, final int day, final String matched,
                                      final String format, final boolean ambiguous ) {
        if ( ( year < MIN_YEAR ) || ( year > MAX_YEAR ) || ( month < 1 ) || ( month > 12 ) || ( day < 1 )
                || ( day > 31 ) ) {
            return null;
        }
        try {
            LocalDate.of( year, month, day ); // rejects e.g. Feb 30
        }
        catch ( final RuntimeException e ) {
            return null;
        }
        return new DateMatch( toDecimalYear( year, month, day ), matched, Precision.DAY, format, ambiguous );
    }

    /** Build a year-month match (day unknown -> mid-month). */
    private static DateMatch fromYm( final int year, final int month, final String matched, final String format ) {
        if ( ( year < MIN_YEAR ) || ( year > MAX_YEAR ) || ( month < 1 ) || ( month > 12 ) ) {
            return null;
        }
        return new DateMatch( toDecimalYear( year, month, 0 ), matched, Precision.MONTH, format, false );
    }

    /** 1-12 for a (possibly long) English month name via its first three letters, else 0. */
    private static int monthNumber( final String name ) {
        if ( ( name == null ) || ( name.length() < 3 ) ) {
            return 0;
        }
        final String key = name.substring( 0, 3 ).toLowerCase( Locale.ROOT );
        for ( int i = 0; i < MONTHS.length; ++i ) {
            if ( MONTHS[ i ].equals( key ) ) {
                return i + 1;
            }
        }
        return 0;
    }

    private TipDateExtractor() {
    }
}

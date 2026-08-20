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

import org.forester.archaeopteryx.tools.TipDateExtractor.DateMatch;
import org.forester.archaeopteryx.tools.TipDateExtractor.DayMonthOrder;
import org.forester.archaeopteryx.tools.TipDateExtractor.Precision;
import org.forester.archaeopteryx.tools.TipDateExtractor.Summary;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.io.parsers.nhx.NHXParser;

/** Headless tests for {@link TipDateExtractor}: every date format, the decimal-year conversion (leap-year + midpoint of
 *  an incomplete date), the day-first/month-first ambiguity, specificity + rightmost preference, non-matches, the
 *  node-writer (date + data:date property, no re-run duplicate), and the whole-tree helpers. */
public final class TipDateExtractorTest {

    private static final DayMonthOrder DF = DayMonthOrder.DAY_FIRST;

    public static void main( final String[] args ) {
        System.out.println( "TipDateExtractor: " + ( test() ? "OK." : "FAILED." ) );
    }

    public static boolean test() {
        try {
            // --- formats ---
            if ( !yearEq( "strain|2021-03-15", 2021.201, Precision.DAY ) ) {
                return fail( "ISO full date" );
            }
            if ( !yearEq( "GISAID/01-Dec-2015", 2015.916, Precision.DAY ) ) {
                return fail( "month-name date (01-Dec-2015)" );
            }
            if ( !yearEq( "x/15-March-2021", 2021.201, Precision.DAY ) ) {
                return fail( "long month name (March)" );
            }
            if ( !yearEq( "seq/2021/03/15", 2021.201, Precision.DAY ) ) {
                return fail( "year-first slash date" );
            }
            if ( !yearEq( "seq/15/03/2021", 2021.201, Precision.DAY ) ) {
                return fail( "day-first slash date, disambiguated (15>12)" );
            }
            if ( !yearEq( "seq_2021.03.15", 2021.201, Precision.DAY ) ) {
                return fail( "dot date must beat decimal-year (2021.03.15)" );
            }
            if ( !yearEq( "iso/2021-03", 2021.204, Precision.MONTH ) ) {
                return fail( "ISO year-month -> mid-month" );
            }
            if ( !yearEq( "May-2021", 2021.371, Precision.MONTH ) ) {
                return fail( "month-name year -> mid-month" );
            }
            if ( !yearEq( "sample_2021.37", 2021.37, Precision.DAY ) ) {
                return fail( "decimal year (2021.37)" );
            }
            if ( !yearEq( "A/Texas/50/2012", 2012.5, Precision.YEAR ) ) {
                return fail( "bare 4-digit year token -> midpoint" );
            }
            // --- leap year ---
            if ( !yearEq( "s/2020-03-15", 2020.204, Precision.DAY ) ) {
                return fail( "leap-year full date (2020)" );
            }
            // --- day/month ambiguity (both <= 12) ---
            final DateMatch amb_df = TipDateExtractor.parse( "s/05/03/2021", DayMonthOrder.DAY_FIRST );
            final DateMatch amb_mf = TipDateExtractor.parse( "s/05/03/2021", DayMonthOrder.MONTH_FIRST );
            if ( ( amb_df == null ) || !amb_df.ambiguous() || ( amb_mf == null ) || !amb_mf.ambiguous() ) {
                return fail( "an all-<=12 numeric date must be flagged ambiguous" );
            }
            // day-first reads 05 as day (March) ~2021.174; month-first reads 05 as month (May) ~2021.336
            if ( ( Math.abs( amb_df.decimalYear() - 2021.174 ) > 2e-3 )
                    || ( Math.abs( amb_mf.decimalYear() - 2021.336 ) > 2e-3 )
                    || ( amb_df.decimalYear() >= amb_mf.decimalYear() ) ) {
                return fail( "day-first vs month-first must differ: df=" + amb_df.decimalYear() + " mf="
                        + amb_mf.decimalYear() );
            }
            // --- specificity + rightmost ---
            if ( !yearEq( "hCoV-19/USA/CA-1234/2020|EPI_ISL_9/2021-03-15", 2021.201, Precision.DAY ) ) {
                return fail( "an explicit ISO date must win over earlier/trailing bare years" );
            }
            // --- non-matches ---
            if ( TipDateExtractor.parse( "just_a_name", DF ) != null ) {
                return fail( "a label with no date must return null" );
            }
            if ( TipDateExtractor.parse( "sample_1500", DF ) != null ) {
                return fail( "a 4-digit number below 1900 must not be read as a year" );
            }
            if ( TipDateExtractor.parse( "s/05/03/21", DF ) != null ) {
                return fail( "a 2-digit-year numeric date is deferred (no 4-digit year -> no match)" );
            }
            if ( TipDateExtractor.parse( "", DF ) != null ) {
                return fail( "empty label -> null" );
            }
            // --- toDecimalYear directly ---
            if ( ( Math.abs( TipDateExtractor.toDecimalYear( 2020, 0, 0 ) - 2020.5 ) > 1e-9 )
                    || ( Math.abs( TipDateExtractor.toDecimalYear( 2021, 3, 15 ) - 2021.201 ) > 1e-3 ) ) {
                return fail( "toDecimalYear midpoint / full date" );
            }
            // --- decimalYearString rounding ---
            if ( !"2012.5".equals( TipDateExtractor.decimalYearString( 2012.5 ) )
                    || !"2020".equals( TipDateExtractor.decimalYearString( 2020.0 ) ) ) {
                return fail( "decimalYearString: got " + TipDateExtractor.decimalYearString( 2020.0 ) );
            }
            if ( !applyToNodeOk() ) {
                return false;
            }
            if ( !wholeTreeHelpersOk() ) {
                return false;
            }
            return true;
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    /** applyToNode writes a <date> (year unit) + a numeric data:date property, and a re-run does not duplicate it. */
    private static boolean applyToNodeOk() {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( "seq/2021-03-15" );
        final DateMatch m = TipDateExtractor.parse( n.getName(), DF );
        if ( !TipDateExtractor.applyToNode( n, m ) ) {
            return fail( "applyToNode must succeed for a match" );
        }
        if ( !n.getNodeData().isHasDate() || ( n.getNodeData().getDate().getValue() == null )
                || ( Math.abs( n.getNodeData().getDate().getValue().doubleValue() - 2021.201 ) > 1e-3 )
                || !"year".equals( n.getNodeData().getDate().getUnit() ) ) {
            return fail( "applyToNode must set a <date> (year unit) at the decimal year" );
        }
        if ( !"2021.20137".equals( propVal( n, TipDateExtractor.DATE_PROPERTY_REF ) ) ) {
            return fail( "applyToNode must set data:date, got " + propVal( n, TipDateExtractor.DATE_PROPERTY_REF ) );
        }
        // re-run: still exactly one data:date property
        TipDateExtractor.applyToNode( n, m );
        int count = 0;
        for ( final Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( TipDateExtractor.DATE_PROPERTY_REF.equals( p.getRef() ) ) {
                count++;
            }
        }
        if ( count != 1 ) {
            return fail( "a re-run must not duplicate data:date, count=" + count );
        }
        // a null match is a no-op
        if ( TipDateExtractor.applyToNode( n, null ) ) {
            return fail( "applyToNode(null) must be a no-op" );
        }
        return true;
    }

    /** mostLabelsHaveDates (majority gate for the load offer) + hasAmbiguousDates (adaptive toggle) + preview. */
    private static boolean wholeTreeHelpersOk() throws Exception {
        // 3 dated tips + 1 undated -> majority
        final Phylogeny mostly = parse( "((a_2019-01-15:1,b_2020-06-02:1):1,(c_2021-03-15:1,d_no_date:1):1);" );
        if ( !TipDateExtractor.mostLabelsHaveDates( mostly ) ) {
            return fail( "mostLabelsHaveDates must be true for a majority-dated tree" );
        }
        final Phylogeny few = parse( "((a_x:1,b_y:1):1,(c_z:1,d_2021:1):1);" ); // only 1/4 dated
        if ( TipDateExtractor.mostLabelsHaveDates( few ) ) {
            return fail( "mostLabelsHaveDates must be false when few tips are dated" );
        }
        if ( TipDateExtractor.preview( mostly, DF ).size() != 4 ) {
            return fail( "preview must have one row per tip" );
        }
        final Summary sum = TipDateExtractor.summarize( TipDateExtractor.preview( mostly, DF ) );
        if ( ( sum.matched() != 3 ) || ( sum.unmatched() != 1 ) || ( sum.ambiguous() != 0 )
                || !sum.dominantFormat().contains( "ISO" ) || ( sum.minYear() >= sum.maxYear() ) ) {
            return fail( "summarize: " + sum );
        }
        // an all-<=12 numeric date (dots, safe in a Newick name) makes hasAmbiguousDates true
        final Phylogeny amb = parse( "(a_05.03.2021:1,b_2020-01-01:1);" );
        if ( !TipDateExtractor.hasAmbiguousDates( amb, DF ) ) {
            return fail( "hasAmbiguousDates must be true when an ambiguous numeric date is present" );
        }
        final Phylogeny unamb = parse( "(a_2021-03-15:1,b_2020-01-01:1);" );
        if ( TipDateExtractor.hasAmbiguousDates( unamb, DF ) ) {
            return fail( "hasAmbiguousDates must be false for unambiguous dates" );
        }
        return true;
    }

    private static boolean yearEq( final String label, final double expected, final Precision prec ) {
        final DateMatch m = TipDateExtractor.parse( label, DF );
        if ( m == null ) {
            return note( "no date found in '" + label + "'" );
        }
        if ( Math.abs( m.decimalYear() - expected ) > 2e-3 ) {
            return note( "'" + label + "' -> " + m.decimalYear() + ", expected ~" + expected );
        }
        if ( m.precision() != prec ) {
            return note( "'" + label + "' precision " + m.precision() + ", expected " + prec );
        }
        return true;
    }

    private static String propVal( final PhylogenyNode n, final String ref ) {
        if ( n.getNodeData().getProperties() == null ) {
            return null;
        }
        for ( final Property p : n.getNodeData().getProperties().getProperties() ) {
            if ( ref.equals( p.getRef() ) ) {
                return p.getValue();
            }
        }
        return null;
    }

    private static Phylogeny parse( final String nh ) throws Exception {
        return ParserBasedPhylogenyFactory.getInstance().create( nh, new NHXParser() )[ 0 ];
    }

    private static boolean note( final String msg ) {
        System.out.println( "  [TipDateExtractorTest] " + msg );
        return false;
    }

    private static boolean fail( final String msg ) {
        return note( msg );
    }

    private TipDateExtractorTest() {
    }
}

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

import java.util.ArrayList;
import java.util.List;

import org.forester.phylogeny.data.Confidence;

/**
 * Unit tests for {@link TreePanel#confidenceLabel}: the MAD support value is shown first and only
 * when its display option is on, is never hidden by the "min. confidence" threshold (low MAD is
 * good), while regular confidence values (e.g. bootstrap; high is good) are shown only at or above
 * the threshold -- so a MAD-rooted tree reads as e.g. "0/90". In {@code org.forester.archaeopteryx}
 * so it can reach the package-private method; run via the suite ({@link #test()}).
 */
public final class ConfidenceLabelTest {

    public static void main( final String[] args ) {
        System.out.println( "ConfidenceLabel: " + ( test() ? "OK." : "FAILED." ) );
        System.exit( test() ? 0 : 1 );
    }

    public static boolean test() {
        final List<Confidence> mad_boot = list( new Confidence( 0.0, "MAD" ), new Confidence( 90.0, "" ) );
        // MAD off (the default): only the regular value shows (>= threshold)
        if ( !"90".equals( TreePanel.confidenceLabel( mad_boot, false, 50.0, false, 0 ) ) ) {
            return fail( "MAD off should show only the bootstrap (90)" );
        }
        // MAD on: "MAD/regular"; the MAD value (0) is NOT hidden by the 50 threshold
        if ( !"0/90".equals( TreePanel.confidenceLabel( mad_boot, true, 50.0, false, 0 ) ) ) {
            return fail( "MAD on should show '0/90'" );
        }
        // a non-zero MAD keeps its own (up to 2-decimal) precision even though regular uses 0 digits
        if ( !"0.05/90".equals( TreePanel.confidenceLabel( list( new Confidence( 0.05, "MAD" ),
                                                                 new Confidence( 90.0, "" ) ), true, 50.0, false, 0 ) ) ) {
            return fail( "MAD should keep 2 decimals -> '0.05/90'" );
        }
        // a regular value below the threshold is hidden, but MAD still shows
        final List<Confidence> mad_low = list( new Confidence( 0.0, "MAD" ), new Confidence( 30.0, "" ) );
        if ( !"0".equals( TreePanel.confidenceLabel( mad_low, true, 50.0, false, 0 ) ) ) {
            return fail( "low bootstrap hidden, MAD shown -> '0'" );
        }
        if ( !"".equals( TreePanel.confidenceLabel( mad_low, false, 50.0, false, 0 ) ) ) {
            return fail( "MAD off and bootstrap below threshold -> empty" );
        }
        // no MAD confidence present at all
        if ( !"90".equals( TreePanel.confidenceLabel( list( new Confidence( 90.0, "" ) ), true, 50.0, false, 0 ) ) ) {
            return fail( "no MAD present -> '90'" );
        }
        // threshold 0 lets a low bootstrap through; MAD still off here
        if ( !"30".equals( TreePanel.confidenceLabel( mad_low, false, 0.0, false, 0 ) ) ) {
            return fail( "threshold 0 shows the low bootstrap (30)" );
        }
        // The "min. confidence shown" is a FRACTION of the tree's detected scale (paintConfidenceValues passes
        // fraction * confidenceScaleMaxFor(tree) as the absolute threshold below). Default 0.5 = hide below
        // half-scale, which is scale-correct where the old fixed 50 silently hid every 0-1 posterior value.
        if ( Options.MIN_CONFIDENCE_FRACTION_DEFAULT != 0.5 ) {
            return fail( "default min-confidence fraction should be 0.5, was " + Options.MIN_CONFIDENCE_FRACTION_DEFAULT );
        }
        final double frac = Options.MIN_CONFIDENCE_FRACTION_DEFAULT;
        // bootstrap scale (0-100): ceiling 100 -> cutoff 50 -> hide 30, show 70 (unchanged from the old behavior)
        final double boot_cut = frac * TreePanelUtil.confidenceScaleMaxFor( 90.0 );
        if ( !"".equals( TreePanel.confidenceLabel( list( new Confidence( 30.0, "" ) ), false, boot_cut, false, 0 ) )
                || !"70".equals( TreePanel.confidenceLabel( list( new Confidence( 70.0, "" ) ), false, boot_cut, false,
                                                            0 ) ) ) {
            return fail( "bootstrap default (cutoff 50) should hide 30 and show 70" );
        }
        // posterior scale (0-1): ceiling 1 -> cutoff 0.5 -> hide 0.3, show 0.9 (invisible under the old fixed 50)
        final double post_cut = frac * TreePanelUtil.confidenceScaleMaxFor( 0.98 );
        if ( !"".equals( TreePanel.confidenceLabel( list( new Confidence( 0.3, "" ) ), false, post_cut, false, 2 ) )
                || !"0.9".equals( TreePanel.confidenceLabel( list( new Confidence( 0.9, "" ) ), false, post_cut, false,
                                                             2 ) ) ) {
            return fail( "posterior default (cutoff 0.5) should hide 0.3 and show 0.9" );
        }
        return true;
    }

    private static List<Confidence> list( final Confidence... cs ) {
        final List<Confidence> l = new ArrayList<Confidence>();
        for ( final Confidence c : cs ) {
            l.add( c );
        }
        return l;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [ConfidenceLabelTest] " + message );
        return false;
    }

    private ConfidenceLabelTest() {
    }
}

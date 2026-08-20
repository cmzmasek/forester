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

import java.awt.GraphicsEnvironment;
import java.io.File;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.TIME_AXIS_TYPE;
import org.forester.archaeopteryx.tools.TipDateExtractor;
import org.forester.archaeopteryx.tools.TipDateExtractor.DayMonthOrder;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * The "Extract Dates from Labels" tool (Increment 1). Dogfoods the {@code date-in-labels.xml} demo: the preview dialog
 * reports 10 matched ISO-dominant tips; applying writes a {@code <date>} (year) + a numeric {@code data:date} property
 * to each tip, so the Calendar axis auto-derives and {@code data:date} is offered as a numeric Color-by ref; a
 * provenance sentence is appended; and Undo restores the undated tree. Headful; a green no-op when headless.
 */
public final class TipDateToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Tip date tool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = parse( "date-in-labels.xml" );
            if ( phy == null ) {
                return fail( "date-in-labels.xml demo missing" );
            }
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "td" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Phylogeny p = tp.getPhylogeny();
                    if ( p.getNumberOfExternalNodes() != 10 ) {
                        fail( ok, "demo must have 10 tips, has " + p.getNumberOfExternalNodes() );
                    }
                    if ( anyTipDated( p ) ) {
                        fail( ok, "the demo tips must start WITHOUT a <date> (the whole point is extraction)" );
                    }
                    // the load-time auto-offer (Inc 2) would fire for this tree: date-bearing labels, no dates yet
                    if ( !TipDateExtractor.shouldOffer( p ) ) {
                        fail( ok, "the load-time offer must fire for the demo (date-bearing labels, undated)" );
                    }
                    // the preview dialog: 10 matched, ISO-dominant, no ambiguity (so day/month order irrelevant)
                    final TipDateExtractionDialog dlg = new TipDateExtractionDialog( frame, p );
                    final String summary = dlg.summaryTextForTest();
                    dlg.dispose();
                    if ( ( summary == null ) || !summary.contains( "10" ) || !summary.contains( "ISO" ) ) {
                        fail( ok, "the preview summary must report 10 ISO matches, got: " + summary );
                    }

                    // apply: writes <date> (year) + data:date to all 10 tips
                    final int n = frame.applyTipDatesAndRefit( tp, DayMonthOrder.DAY_FIRST, true );
                    if ( n != 10 ) {
                        fail( ok, "must date all 10 tips, dated " + n );
                    }
                    for ( final PhylogenyNode ext : p.getExternalNodes() ) {
                        if ( !ext.getNodeData().isHasDate() || ( ext.getNodeData().getDate().getValue() == null )
                                || !"year".equals( ext.getNodeData().getDate().getUnit() ) ) {
                            fail( ok, "every tip must get a <date> (year), missing on " + ext.getName() );
                        }
                    }
                    if ( !hasTipProperty( p, TipDateExtractor.DATE_PROPERTY_REF )
                            || !PropertyColorScheme.numericRefs( p ).contains( TipDateExtractor.DATE_PROPERTY_REF ) ) {
                        fail( ok, "data:date must be a numeric Color-by ref (date gradient)" );
                    }
                    // the Calendar axis auto-derives from the new "year" dates
                    if ( tp.effectiveTimeAxisType() != TIME_AXIS_TYPE.CALENDAR ) {
                        fail( ok, "extracted year dates must auto-derive the CALENDAR axis, got "
                                + tp.effectiveTimeAxisType() );
                    }
                    // the tip-dated tree is now detected DATED -> the "Time tree" badge shows (Inc 2)
                    if ( AptxUtil.detectTimeTree( p ) != AptxUtil.TIME_TREE_KIND.DATED ) {
                        fail( ok, "a tip-dated tree must be detected DATED (time-tree badge), got "
                                + AptxUtil.detectTimeTree( p ) );
                    }
                    // provenance appended
                    if ( ( p.getDescription() == null ) || !p.getDescription().contains( "Extracted sampling dates" ) ) {
                        fail( ok, "a provenance sentence must be appended, got: " + p.getDescription() );
                    }

                    // re-running with skip-existing must write 0 (all tips are now dated) -- so the tool reports a
                    // no-op instead of a false "wrote dates" success, and pushes no undo checkpoint
                    if ( frame.applyTipDatesAndRefit( tp, DayMonthOrder.DAY_FIRST, true ) != 0 ) {
                        fail( ok, "a second run with skip-existing must write 0 (all tips already dated)" );
                    }

                    // Undo restores the undated tree
                    frame.undo();
                    if ( anyTipDated( frame.getMainPanel().getCurrentTreePanel().getPhylogeny() ) ) {
                        fail( ok, "Undo must restore the undated tree" );
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace();
                    fail( ok, "exception: " + t );
                }
                finally {
                    frame.dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return fail( "exception: " + t );
        }
    }

    private static boolean anyTipDated( final Phylogeny phy ) {
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( ext.getNodeData().isHasDate() ) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasTipProperty( final Phylogeny phy, final String ref ) {
        for ( final PhylogenyNode ext : phy.getExternalNodes() ) {
            if ( ext.getNodeData().getProperties() != null ) {
                for ( final org.forester.phylogeny.data.Property pr : ext.getNodeData().getProperties()
                        .getProperties() ) {
                    if ( ref.equals( pr.getRef() ) ) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    private static Phylogeny parse( final String demo ) {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
        if ( !f.exists() ) {
            return null;
        }
        try {
            return ParserBasedPhylogenyFactory.getInstance().create( f, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
        }
        catch ( final Exception e ) {
            return null;
        }
    }

    private static boolean fail( final String m ) {
        System.out.println( "  [TipDateToolTest] " + m );
        return false;
    }

    private static void fail( final boolean[] ok, final String m ) {
        ok[ 0 ] = false;
        fail( m );
    }

    private TipDateToolTest() {
    }
}

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
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.nio.file.Files;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.tools.NodeDataImporter;
import org.forester.phylogeny.Phylogeny;

/**
 * <b>Settings &gt; Layout &gt; Legend in Its Own Column</b>, on the case it exists for: a root-left clustergram whose
 * matrix runs to the right edge, where the legend's default corner sits on the annotation-column HEADERS and hides
 * the last of the gene names.
 * <p>
 * The two settings are asked the SAME question and must answer it the other way round. Both sides of the question are
 * production geometry, never this test's own idea of where either thing should be: the legend's real drawn bounds
 * ({@link TreePanel#getPropertyLegendBounds}, recorded by the paint itself), and the last column's header band -- its
 * x from {@link TreePanel#annotationColumnGrabPointForTest} (the point a click and a drag hit-test against), its
 * height from the layout's own {@code annotationHeaderTopReserve}.
 * <ul>
 * <li>OFF: the legend box overlaps that band. That is the bug, and it is asserted rather than assumed, so a layout
 * change that stopped the two colliding would fail here instead of quietly turning the check below into a no-op.</li>
 * <li>ON: it no longer does, the whole matrix ends left of the legend, and every column has moved left.</li>
 * </ul>
 * The fixture is the sparse accessory genome demo, joined with its table exactly as File &gt; Demo Trees does.
 * Headful; a green no-op when headless.
 */
public final class LegendColumnTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LegendColumn: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return legendColumnOk();
    }

    private static boolean fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [LegendColumnTest] " + msg );
        ok[ 0 ] = false;
        return false;
    }

    /** The demo tree with its table joined on, as File &gt; Demo Trees leaves it. */
    private static Phylogeny demo() throws Exception {
        final File dir = new File( System.getProperty( "user.dir" ), "forester/demo" );
        final Phylogeny phy = FigureRenderer.readTrees( new File( dir, "sparse-accessory-genome.xml" ) )[ 0 ];
        final NodeDataImporter.Table table = NodeDataImporter
                .parseTable( Files.readString( new File( dir, "sparse-accessory-genome.tsv" ).toPath() ) );
        NodeDataImporter.apply( phy, table, table.defaultKeyColumn(), NodeDataImporter.MatchBy.TIP_NAME );
        return phy;
    }

    private static boolean legendColumnOk() {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        try {
            final Phylogeny phy = demo();
            SwingUtilities.invokeAndWait( () -> {
                mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                               "legend column" );
                ( ( JFrame ) mf[ 0 ] ).setLocation( -32000, -32000 );
                ( ( JFrame ) mf[ 0 ] ).setSize( 1200, 800 );
                ( ( JFrame ) mf[ 0 ] ).setVisible( true );
            } );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    mf[ 0 ].applyClustergramPreset(); // rectangular, root on the left, the matrix at the right
                    final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                    mf[ 0 ].getMainPanel().getControlPanel().showWhole(); // one fit to start from; the clicks re-fit
                    final int columns = MatrixColumnOrder.matrixRefs( tp.getAnnotationColumnSpecs() ).size();
                    if ( columns != 18 ) {
                        fail( ok, "the fixture should give 18 matrix columns, got " + columns );
                        return;
                    }
                    // ON by default: the window agrees with the figure unless the user says otherwise
                    if ( !mf[ 0 ].getOptions().isReserveLegendColumn() ) {
                        fail( ok, "the setting must be ON by default" );
                    }
                    if ( !Options.createDefaultInstance().isReserveLegendColumn() ) {
                        fail( ok, "a fresh Options must default to ON (what Reset to Defaults restores)" );
                    }
                    // the backing menu item drives Options, and applyOptionsToMenuStates is its inverse
                    mf[ 0 ]._legend_column_cbmi.doClick();
                    if ( mf[ 0 ].getOptions().isReserveLegendColumn() ) {
                        fail( ok, "clicking the backing menu item must write the setting through to Options" );
                    }
                    // ...and the inverse must CORRECT an item that is out of step. Putting it back from where the
                    // click already left it would pass even if applyOptionsToMenuStates ignored this item entirely,
                    // so drive the item (not Options: setSelected fires no action) the WRONG way first.
                    mf[ 0 ]._legend_column_cbmi.setSelected( true );
                    mf[ 0 ].applyOptionsToMenuStates( mf[ 0 ].getOptions() );
                    if ( mf[ 0 ]._legend_column_cbmi.isSelected() ) {
                        fail( ok, "applyOptionsToMenuStates must put the menu item back where Options says" );
                    }
                    final Rectangle off_legend = legendAfterPaint( mf[ 0 ], tp, false );
                    final Point off_header = tp.annotationColumnGrabPointForTest( columns - 1 );
                    final int off_reserve = tp.legendColumnReserve();
                    final Rectangle on_legend = legendAfterPaint( mf[ 0 ], tp, true );
                    final Point on_header = tp.annotationColumnGrabPointForTest( columns - 1 );
                    final int on_reserve = tp.legendColumnReserve();
                    if ( ( off_legend == null ) || ( on_legend == null ) ) {
                        fail( ok, "a clustergram always draws the shared matrix legend; got " + off_legend + " / "
                                + on_legend );
                        return;
                    }
                    if ( ( off_header == null ) || ( on_header == null ) ) {
                        fail( ok, "the last column must have a header to hit-test; got " + off_header + " / "
                                + on_header );
                        return;
                    }
                    // the header band: the strip above the first tip that the rotated column names are drawn in,
                    // its height straight from the layout (annotationHeaderTopReserve), at the last column's x
                    final int band = tp.annotationHeaderTopReserveForTest();
                    if ( band <= 0 ) {
                        fail( ok, "the clustergram must reserve a band for its column headers, got " + band );
                        return;
                    }
                    final Rectangle off_band = new Rectangle( off_header.x, 0, 1, band );
                    final Rectangle on_band = new Rectangle( on_header.x, 0, 1, band );
                    // OFF: the bug. Assert it, or the ON case below proves nothing.
                    if ( off_reserve != 0 ) {
                        fail( ok, "with the setting OFF nothing may be reserved, got " + off_reserve );
                    }
                    if ( !off_legend.intersects( off_band ) ) {
                        fail( ok, "fixture: with the setting OFF the legend " + off_legend + " must cover the last"
                                + " column's header band " + off_band + " -- if it no longer does, this test stopped"
                                + " measuring the thing the setting fixes" );
                    }
                    // ON: the same header, clear of the legend, in a column of its own
                    if ( on_reserve <= 0 ) {
                        fail( ok, "with the setting ON the legend must get a column, reserve was " + on_reserve );
                    }
                    if ( on_legend.intersects( on_band ) ) {
                        fail( ok, "with the setting ON the legend " + on_legend + " still covers the last column's"
                                + " header band " + on_band );
                    }
                    if ( on_header.x >= on_legend.x ) {
                        fail( ok, "the last column must end left of the legend box: header x=" + on_header.x
                                + ", legend starts at x=" + on_legend.x );
                    }
                    // the column is what moved it: the header sits further left than it did with the setting off
                    if ( on_header.x >= off_header.x ) {
                        fail( ok, "reserving the column must pull the matrix left: header x went from " + off_header.x
                                + " to " + on_header.x );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                    t.printStackTrace();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
        finally {
            if ( mf[ 0 ] != null ) {
                try {
                    SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
                }
                catch ( final Exception ignored ) {
                    // disposing a frame that never came up is not a test failure
                }
            }
        }
        return ok[ 0 ];
    }

    /**
     * Puts the setting where {@code reserve} wants it THROUGH THE UI -- a click on the backing menu item, the same
     * event the Settings checkbox sends -- and paints. Nothing here re-fits the tree: the click's own handler has to,
     * or the new reserve would not reach the layout. The legend bounds this test reads are recorded by the paint, so
     * a measurement taken without one would be the previous layout's.
     */
    private static Rectangle legendAfterPaint( final MainFrame mf, final TreePanel tp, final boolean reserve ) {
        if ( mf.getOptions().isReserveLegendColumn() != reserve ) {
            mf._legend_column_cbmi.doClick();
        }
        tp.validate();
        final BufferedImage img = new BufferedImage( Math.max( tp.getWidth(), 100 ), Math.max( tp.getHeight(), 100 ),
                                                     BufferedImage.TYPE_INT_RGB );
        tp.printAll( img.getGraphics() );
        img.flush();
        final Rectangle r = tp.getPropertyLegendBounds();
        return ( r == null ) ? null : new Rectangle( r ); // a copy: the next paint overwrites the panel's own
    }
}

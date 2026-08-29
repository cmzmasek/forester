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
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.TreePanel.CLADE_LABEL_ANGLE;
import org.forester.archaeopteryx.TreePanel.CLADE_VIS;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * With several ranks annotated at once, the clade legend must say WHICH rank each row belongs to -- one
 * alphabetical heap of families, subfamilies and genera is unreadable, and was the state this increment fixed.
 * Asserts that a multi-level legend is split into one titled block per rank, that a single-level legend is left
 * exactly as it was (flat, no headings), and that the default legend corner moves off the right edge, where the
 * clade columns are drawn -- an annotated figure used to open with its legend sitting on top of the very marks
 * it describes, hiding whole ranks of them. Guarded to a no-op on a headless box.
 */
public final class CladeLegendSectionsTest {

    private static final String DEMO = "forester/demo/bat-phylogeny.xml";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "CladeLegendSections: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final File file = new File( System.getProperty( "user.dir" ), DEMO );
            if ( !file.exists() ) {
                System.out.println( "  [CladeLegendSectionsTest] missing demo tree " + file.getAbsolutePath() );
                return false;
            }
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( phys, new Configuration(), "clade-legend" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    exercise( mf[ 0 ], ok );
                }
                finally {
                    ( ( JFrame ) mf[ 0 ] ).dispose();
                }
            } );
            reachabilityOnAnimalTree( ok );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void exercise( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 1200, 900 );
        ( ( JFrame ) frame ).validate();
        final TreePanel tp = frame.getCurrentTreePanel();
        frame.getCurrentTreePanel().setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );

        // ---- one level: the flat legend, unchanged -----------------------------------------------------
        tp.setCladeLevels( specs( "family" ), CLADE_VIS.BARS, true );
        List<String> rows = drawAndReadRows( tp );
        if ( rows.isEmpty() ) {
            fail( ok, "a single-level clade legend should still list its taxa" );
            return;
        }
        if ( rows.contains( null ) ) {
            fail( ok, "a single-level legend must stay FLAT -- no rank headings, got " + rows );
        }

        // ---- three levels: one titled block per rank ---------------------------------------------------
        tp.setCladeLevels( specs( "genus", "subfamily", "family" ), CLADE_VIS.BARS, true );
        rows = drawAndReadRows( tp );
        int headings = 0;
        for ( final String row : rows ) {
            if ( row == null ) {
                ++headings;
            }
        }
        if ( headings != 3 ) {
            fail( ok, "three annotated ranks must give three legend headings, got " + headings + " in " + rows );
            return;
        }
        // the blocks run broadest-first (the way taxonomy is written), so the first block holds the families
        final List<String> first_block = blockAfter( rows, 0 );
        final List<String> last_block = blockAfter( rows, 2 );
        if ( !first_block.contains( "Vespertilionidae" ) ) {
            fail( ok, "the first legend block should be the FAMILY block, got " + first_block );
        }
        if ( !last_block.contains( "Myotis" ) ) {
            fail( ok, "the last legend block should be the GENUS block, got " + last_block );
        }
        if ( first_block.contains( "Myotis" ) || last_block.contains( "Vespertilionidae" ) ) {
            fail( ok, "a taxon must appear only under its own rank's heading" );
        }
        // every rank keeps rows: one shared cap would let a big rank crowd the others out entirely
        if ( blockAfter( rows, 1 ).isEmpty() ) {
            fail( ok, "the middle (subfamily) block must not be squeezed to nothing" );
        }

        // ---- the sort control must survive going multi-level -------------------------------------------
        // the sort order still drives the rows inside every block, so the control that changes it has to stay
        tp.setCladeLevels( specs( "genus", "subfamily", "family" ), CLADE_VIS.BARS, true );
        drawAndReadRows( tp );
        if ( tp.legendSortToggleBoundsForTest() == null ) {
            fail( ok, "a sectioned legend must still offer the [by count]/[A-Z] sort toggle" );
        }

        // ---- the default legend corner clears the clade columns ----------------------------------------
        if ( !tp.legendDefaultGoesLeftForTest() ) {
            fail( ok, "with bar columns on the right, the default legend corner must move off the right edge" );
        }
        tp.setCladeLevels( specs( "family" ), CLADE_VIS.BOXES, true );
        if ( tp.legendDefaultGoesLeftForTest() ) {
            fail( ok, "boxes draw no right-edge columns, so the legend corner must stay where it was" );
        }
        tp.setCladeLevels( specs( "genus", "family" ), CLADE_VIS.BARS, true );
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        if ( tp.legendDefaultGoesLeftForTest() ) {
            fail( ok, "the circular layout rings the tree rather than filling the right edge -- corner unchanged" );
        }
    }

    /**
     * "Show all" must not be a one-way trip. The expand/collapse control used to live at the BOTTOM of the legend
     * box, so expanding a long legend grew the box past the window and scrolled the way back off the screen --
     * and the sectioned renderer drew no control at all once nothing was hidden. The control now sits in the
     * TITLE row, and the box's top is clamped into view, so it stays reachable at any length. Checked for both
     * the flat (single-rank) and the sectioned (multi-rank) legend.
     */
    private static void expandControlStaysReachable( final TreePanel tp, final boolean[] ok ) {
        // a deliberately small viewport, so an expanded legend is certainly taller than it
        final Rectangle view = new Rectangle( 0, 0, 420, 260 );
        for ( final String what : new String[] { "flat", "sectioned" } ) {
            if ( "flat".equals( what ) ) {
                tp.setCladeLevels( specs( "order" ), CLADE_VIS.BARS, false );
            }
            else {
                tp.setCladeLevels( specs( "order", "class", "phylum" ), CLADE_VIS.BARS, false );
            }
            // the check only means something when entries are actually hidden; say so plainly if the fixture
            // ever stops providing enough, rather than failing further down with a confusing message
            if ( tp.cladeBandCount() <= 20 ) {
                fail( ok, what + ": this check needs a rank with more taxa than the default legend cap, got "
                        + tp.cladeBandCount() );
                continue;
            }
            drawLegend( tp, view );
            Rectangle chip = tp.legendMoreBoundsForTest();
            if ( chip == null ) {
                fail( ok, what + ": a legend with hidden entries must offer a way to show them" );
                continue;
            }
            if ( !view.contains( chip ) ) {
                fail( ok, what + ": the 'show all' control must be on screen, got " + chip + " in " + view );
                continue;
            }
            click( tp, chip ); // expand
            drawLegend( tp, view );
            chip = tp.legendMoreBoundsForTest();
            if ( chip == null ) {
                fail( ok, what + ": an EXPANDED legend must still offer a way back to the short list" );
                continue;
            }
            if ( !view.contains( chip ) ) {
                fail( ok, what + ": the collapse control must stay on screen when the legend is expanded, got "
                        + chip + " in " + view );
            }
            click( tp, chip ); // collapse again, leaving the panel as we found it
            drawLegend( tp, view );
            if ( tp.legendMoreBoundsForTest() == null ) {
                fail( ok, what + ": collapsing should restore the 'show all' control" );
            }
        }
    }

    /** Fires the real click handler at the centre of {@code r}. */
    private static void click( final TreePanel tp, final Rectangle r ) {
        tp.handleLegendClick( new java.awt.event.MouseEvent( tp, java.awt.event.MouseEvent.MOUSE_CLICKED,
                                                             System.currentTimeMillis(), 0, r.x + ( r.width / 2 ),
                                                             r.y + ( r.height / 2 ), 1, false ) );
    }

    private static void drawLegend( final TreePanel tp, final Rectangle bounds ) {
        final BufferedImage img = new BufferedImage( Math.max( 1, bounds.width ), Math.max( 1, bounds.height ),
                                                     BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        try {
            tp.drawLegendForTest( g, bounds, true );
        }
        finally {
            g.dispose();
        }
    }

    /**
     * Runs {@link #expandControlStaysReachable} on the ANIMAL demo rather than the bat one. Its phylum/class/order
     * names are explicit internal-node annotations, so the number of taxa at a rank is fixed by the file; the bat
     * tree's genera are resolved through the process-wide lineage-service cache, which other tests in the suite
     * populate, so its count is not something a test can rely on.
     */
    private static void reachabilityOnAnimalTree( final boolean[] ok ) {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/animal-tree-of-life.xml" );
        if ( !file.exists() ) {
            fail( ok, "missing demo tree " + file.getAbsolutePath() );
            return;
        }
        try {
            final Phylogeny[] phys = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( phys, new Configuration(), "clade-legend-reach" ) );
            SwingUtilities.invokeAndWait( () -> {
                try {
                    ( ( JFrame ) mf[ 0 ] ).setSize( 1000, 800 );
                    ( ( JFrame ) mf[ 0 ] ).validate();
                    expandControlStaysReachable( mf[ 0 ].getCurrentTreePanel(), ok );
                }
                finally {
                    ( ( JFrame ) mf[ 0 ] ).dispose();
                }
            } );
        }
        catch ( final Exception e ) {
            fail( ok, "reachability check threw " + e );
        }
    }

    /** Paints the legend into an off-screen image (which is what records the row list) and returns the rows. */
    private static List<String> drawAndReadRows( final TreePanel tp ) {
        final BufferedImage img = new BufferedImage( 600, 600, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        try {
            tp.drawLegendForTest( g, new Rectangle( 0, 0, 600, 600 ), true );
        }
        finally {
            g.dispose();
        }
        return new ArrayList<>( tp.legendRowLabelsForTest() );
    }

    /** The taxa listed under heading number {@code n} (0-based), i.e. up to the next heading. */
    private static List<String> blockAfter( final List<String> rows, final int n ) {
        final List<String> out = new ArrayList<>();
        int seen = -1;
        for ( final String row : rows ) {
            if ( row == null ) {
                ++seen;
                continue;
            }
            if ( seen == n ) {
                out.add( row );
            }
        }
        return out;
    }

    private static List<CladeLevel.Spec> specs( final String... ranks ) {
        final List<CladeLevel.Spec> out = new ArrayList<>();
        for ( final String r : ranks ) {
            out.add( new CladeLevel.Spec( r, CLADE_LABEL_ANGLE.VERTICAL ) );
        }
        return out;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [CladeLegendSectionsTest] " + msg );
        ok[ 0 ] = false;
    }

    private CladeLegendSectionsTest() {
    }
}

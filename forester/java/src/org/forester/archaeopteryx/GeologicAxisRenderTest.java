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

import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;

/**
 * Renders the dinosaur time-tree demo (forester/demo/dinosaur-time-tree.xml) with the geologic Time Axis ON and
 * asserts the coloured ICS period/epoch bands appear along the bottom -- and NONE when the axis is off. Also checks
 * the calibration (root age = the oldest {@code <date>} = 250 Ma) and the approved biological exception (the geologic
 * axis is a rectangular root-left overlay, so it does NOT apply in a circular layout). Headful; a green no-op when
 * headless. Dogfoods the demo tree (the "feature test loads its demo tree" convention).
 */
public final class GeologicAxisRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "GeologicAxisRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return axisRendersOk();
    }

    private static boolean axisRendersOk() {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/dinosaur-time-tree.xml" );
            if ( !file.exists() ) {
                return fail( "demo tree missing: " + file.getAbsolutePath() );
            }
            final PhyloXmlParser parser = PhyloXmlParser.createPhyloXmlParser();
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance().create( file, parser )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "geo" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                    final Options o = frame.getOptions();
                    o.setGraphicsExportWhiteBackground( false ); // keep our forced colors (no light-theme switch)
                    // force the tree to pure black-on-white so the ONLY saturated colors in the image are the
                    // official ICS band fills (green/cyan/purple/orange)
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BACKGROUND, new Color( 255, 255, 255 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH, new Color( 0, 0, 0 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.BRANCH_LENGTH, new Color( 0, 0, 0 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.SEQUENCE, new Color( 0, 0, 0 ) );
                    tp.getTreeColorSet().setColorforDefault( TreeColorSet.TAXONOMY, new Color( 0, 0, 0 ) );
                    tp.getTreeColorSet().setColorSchema( 0 );
                    // pin root-left rectangular: the developer's persisted orientation/graphics-type may be vertical or
                    // circular, and the first checks below verify the horizontal (root-left) axis specifically
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    o.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    final int w = 900, h = 560;
                    frame.showWhole();
                    tp.setSize( w, h );
                    tp.calcParametersForPainting( w, h );

                    // calibration: the root age is taken from the tree's oldest <date> (Archosauria at 250 Ma)
                    final double root_age = tp.timeAxisRootAgeMa();
                    if ( Math.abs( root_age - 250.0 ) > 0.5 ) {
                        fail( ok, "root age should be the oldest <date> (250 Ma), got " + root_age );
                    }
                    // an explicit "Set root age…" override applies to THIS tree, but must NOT leak onto a replacement
                    // tree (navigating into a subtree / undo / paste replaces the panel's tree)
                    tp.setTimeAxisRootAge( 500.0 );
                    if ( Math.abs( tp.timeAxisRootAgeMa() - 500.0 ) > 0.5 ) {
                        fail( ok, "an explicit root-age override should apply to the current tree, got "
                                + tp.timeAxisRootAgeMa() );
                    }
                    tp.setTree( phy.copy() ); // replace the panel's tree with a distinct object
                    if ( Math.abs( tp.timeAxisRootAgeMa() - 250.0 ) > 0.5 ) {
                        fail( ok, "the root-age override must not leak onto a replacement tree (should revert to the "
                                + "derived 250 Ma), got " + tp.timeAxisRootAgeMa() );
                    }
                    tp.setTimeAxisRootAge( 0 ); // clear for the rest of the checks (derived 250 Ma stands)
                    tp.calcParametersForPainting( w, h ); // re-lay-out the replacement tree before the render checks

                    // OFF: no geologic axis -> no saturated band pixels. Re-fit after each type change, exactly as the
                    // app does (SettingsDialog.maybeCalibrateAndFit -> reFitCurrentTree), so the layout reserve/R
                    // transform matches the option (an on-screen change without a re-fit is not a real usage path).
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    tp.calcParametersForPainting( w, h );
                    if ( tp.geologicAxisApplies() ) {
                        fail( ok, "the geologic axis must not apply when the Time Axis type is NONE" );
                    }
                    final int off = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );

                    // ON: the geologic axis draws the coloured ICS bands
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.GEOLOGIC );
                    tp.calcParametersForPainting( w, h );
                    if ( !tp.geologicAxisApplies() ) {
                        fail( ok, "the geologic axis must apply to a dated phylogram with a positive root age "
                                + "(isDrawPhylogram=" + tp.getControlPanel().isDrawPhylogram() + ", rootAge=" + root_age
                                + ")" );
                    }
                    final BufferedImage geo_img = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    final int on = countSaturated( geo_img );
                    if ( on <= ( off + 2000 ) ) {
                        fail( ok, "the geologic axis should add a lot of coloured band pixels (on=" + on + " off=" + off
                                + ")" );
                    }

                    // OPTIONAL grid lines (off by default): turning them on adds faint reference lines through the tree
                    o.setShowGeologicGridLines( false );
                    final BufferedImage grid_off = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    o.setShowGeologicGridLines( true );
                    final BufferedImage grid_on = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( diffPixels( grid_off, grid_on ) < 200 ) {
                        fail( ok, "the geologic grid lines should draw faint reference lines through the tree when on" );
                    }
                    o.setShowGeologicGridLines( false );

                    // OPTIONAL boundary ages (off by default): turning them on reserves an extra row (the tips compress
                    // -> smaller y-distance) AND draws the age labels, so the image changes
                    tp.calcParametersForPainting( w, h );
                    final float yd_no_ages = tp.getYdistance();
                    final BufferedImage ages_off = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    o.setShowGeologicBoundaryAges( true );
                    tp.calcParametersForPainting( w, h );
                    final float yd_ages = tp.getYdistance();
                    final BufferedImage ages_on = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                    if ( !( yd_ages < yd_no_ages ) ) {
                        fail( ok, "the boundary-age row must reserve extra space (compress the tips): y-distance no-ages="
                                + yd_no_ages + " ages=" + yd_ages );
                    }
                    if ( diffPixels( ages_off, ages_on ) < 200 ) {
                        fail( ok, "the geologic boundary ages should draw age labels when on" );
                    }
                    o.setShowGeologicBoundaryAges( false );
                    tp.calcParametersForPainting( w, h );

                    // the boundary-ages MENU TOGGLE must re-fit so the extra reserved row is applied (the checks above
                    // set the option + recalced directly; this exercises the actionPerformed re-fit a user's click
                    // triggers -- without it the age row would overlap the bottom tip). Geologic axis applies here.
                    if ( frame._show_geologic_ages_cbmi != null ) {
                        frame._show_geologic_ages_cbmi.setSelected( false );
                        o.setShowGeologicBoundaryAges( false );
                        frame.showWhole();
                        final float yd_off_click = tp.getYdistance();
                        frame._show_geologic_ages_cbmi.doClick(); // turns ON + re-fits via the real handler
                        final float yd_on_click = tp.getYdistance();
                        if ( !( yd_on_click < yd_off_click ) ) {
                            fail( ok, "toggling Geologic Boundary Ages via the menu must re-fit so the reserved row is "
                                    + "applied (tips compress): y-distance off=" + yd_off_click + " on=" + yd_on_click );
                        }
                        frame._show_geologic_ages_cbmi.setSelected( false );
                        o.setShowGeologicBoundaryAges( false );
                        frame.showWhole();
                        tp.setSize( w, h );
                        tp.calcParametersForPainting( w, h );
                    }

                    // the bottom band must reserve space -> the tip-spread is compressed (smaller y-distance) vs off
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    tp.calcParametersForPainting( w, h );
                    final float yd_off = tp.getYdistance();
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.GEOLOGIC );
                    tp.calcParametersForPainting( w, h );
                    final float yd_on = tp.getYdistance();
                    if ( !( yd_on < yd_off ) ) {
                        fail( ok, "the geologic axis must reserve a bottom band (compress the tip-spread): y-distance "
                                + "off=" + yd_off + " on=" + yd_on );
                    }
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );

                    // PARITY: the two-band geologic axis also renders in the root-TOP and root-BOTTOM orientations (it
                    // rides the canvas rotation R into a side band down the breadth edge)
                    for ( final Options.TREE_ORIENTATION ori : new Options.TREE_ORIENTATION[] {
                            Options.TREE_ORIENTATION.ROOT_TOP, Options.TREE_ORIENTATION.ROOT_BOTTOM } ) {
                        o.setTreeOrientation( ori );
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                        tp.calcParametersForPainting( w, h );
                        final int voff = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1,
                                false ) );
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.GEOLOGIC );
                        if ( !tp.geologicAxisApplies() ) {
                            fail( ok, "the geologic axis must apply in the " + ori + " orientation" );
                        }
                        tp.calcParametersForPainting( w, h );
                        final int von = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1,
                                false ) );
                        if ( von <= ( voff + 2000 ) ) {
                            fail( ok, "the geologic axis should add coloured band pixels in " + ori + " (on=" + von
                                    + " off=" + voff + ")" );
                        }
                    }
                    o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    tp.calcParametersForPainting( w, h );

                    // PARITY: the CIRCULAR layout shows the geologic scale as concentric coloured ICS rings
                    o.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    frame.showWhole();
                    final int coff = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.GEOLOGIC );
                    if ( !tp.geologicRingsApplyCircular() ) {
                        fail( ok, "the geologic rings must apply to a circular phylogram with a positive root age" );
                    }
                    if ( tp.geologicAxisApplies() ) {
                        fail( ok, "the rectangular two-band axis must NOT apply in a circular layout (rings apply there)" );
                    }
                    final int con = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    // a large delta: the translucent annuli fill most of the disc (tens of thousands of px), well beyond
                    // the few opaque period-name chips -- so this fails if the ring FILLS (not just the labels) go missing
                    if ( con <= ( coff + 15000 ) ) {
                        fail( ok, "the circular geologic rings should fill coloured annuli across the disc (on=" + con
                                + " off=" + coff + ")" );
                    }

                    // B&W circular export: the annuli render as light grey (geologicRingFill), and the period-name
                    // labels must use BLACK ink on that grey (not labelInkOn(vivid), which vanishes) -- exercise the
                    // path renders grey bands + does not crash
                    o.setPrintBlackAndWhite( true );
                    final int grey = countLightGrey( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                    if ( grey < 3000 ) {
                        fail( ok, "B&W circular geologic export should render light-grey annuli (got " + grey + " px)" );
                    }
                    o.setPrintBlackAndWhite( false );

                    // UNROOTED is N/A (an approved biological exception): neither predicate applies there
                    o.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                    if ( tp.geologicAxisApplies() || tp.geologicRingsApplyCircular() ) {
                        fail( ok, "the geologic axis is N/A in the UNROOTED layout (approved biological exception)" );
                    }

                    o.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    tp.calcParametersForPainting( w, h );

                    // DEEP TIME: a tree reaching into the Archean (oldest date > 2500 Ma) makes the axis adapt to the
                    // coarse Eon/Era band pair (epochs/periods run out) -- the Precambrian must be BANDED, not blank. Load
                    // the deep-time demo and assert the geologic bands still render (the coarse ranks paint through the
                    // same path). Guards the bandRanks() wiring in the paint methods.
                    final File deep_file = new File( System.getProperty( "user.dir" ),
                                                     "forester/demo/tree-of-life-deep-time.xml" );
                    if ( deep_file.exists() ) {
                        final Phylogeny deep = ParserBasedPhylogenyFactory.getInstance()
                                .create( deep_file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
                        tp.setTree( deep );
                        o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                        frame.showWhole(); // re-fit the replacement tree (setTree alone leaves a stale preferred size)
                        tp.setSize( w, h );
                        if ( tp.timeAxisRootAgeMa() < 2500.0 ) {
                            fail( ok, "the deep-time demo should reach into the Archean (> 2500 Ma), got "
                                    + tp.timeAxisRootAgeMa() );
                        }
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                        tp.calcParametersForPainting( w, h );
                        final int doff = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1,
                                false ) );
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.GEOLOGIC );
                        tp.calcParametersForPainting( w, h );
                        final int don = countSaturated( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1,
                                false ) );
                        if ( don <= ( doff + 2000 ) ) {
                            fail( ok, "the deep-time (Eon/Era) geologic axis should band the Precambrian (on=" + don
                                    + " off=" + doff + ")" );
                        }
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    }

                    // CALENDAR axis: a tip-dated tree (calendar-year dates) draws a labelled year ruler (not coloured
                    // bands). Load the SARS-CoV-2 demo; present = the most-recent tip (2022.5) is derived automatically.
                    final File cal_file = new File( System.getProperty( "user.dir" ),
                                                    "forester/demo/sars-cov-2-time-tree.xml" );
                    if ( cal_file.exists() ) {
                        final Phylogeny cal = ParserBasedPhylogenyFactory.getInstance()
                                .create( cal_file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
                        tp.setTree( cal );
                        o.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                        o.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                        frame.showWhole();
                        tp.setSize( w, h );
                        if ( Math.abs( tp.timeAxisPresentDate() - 2022.5 ) > 0.05 ) {
                            fail( ok, "the SARS demo present date should be the most-recent tip (2022.5), got "
                                    + tp.timeAxisPresentDate() );
                        }
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                        tp.calcParametersForPainting( w, h );
                        final BufferedImage cal_off = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.CALENDAR );
                        if ( !tp.calendarAxisApplies() ) {
                            fail( ok, "the calendar axis must apply to a tip-dated phylogram (present="
                                    + tp.timeAxisPresentDate() + ")" );
                        }
                        tp.calcParametersForPainting( w, h );
                        final BufferedImage cal_on = AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false );
                        // count DARK ink in the reserved bottom strip (the axis band, below the compressed tips): the
                        // year ticks + labels live there. Counting the strip -- not a whole-image diff -- guards the AXIS
                        // itself, not just the tip-compression the reserve would cause even if the axis weren't drawn.
                        final int strip = 26;
                        if ( countDarkInBottomStrip( cal_on, strip ) <= ( countDarkInBottomStrip( cal_off, strip )
                                + 40 ) ) {
                            fail( ok, "the calendar axis should draw a year ruler (ticks + labels) in the bottom band" );
                        }
                        // circular parity: the calendar rings apply in a circular phylogram; the rectangular axis does not
                        o.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                        if ( !tp.calendarRingsApplyCircular() || tp.calendarAxisApplies() ) {
                            fail( ok, "the calendar rings apply in circular; the rectangular calendar axis must not" );
                        }
                        o.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                        o.setTimeAxisType( Options.TIME_AXIS_TYPE.NONE );
                    }
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Count strongly-saturated colored pixels (max-min channel large, not near-black/near-white) -- with the tree
     *  forced black-on-white, these are exactly the official ICS geologic band fills. */
    private static int countSaturated( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) );
                final int min = Math.min( r, Math.min( g, b ) );
                // >22 catches both the opaque rectangular bands (diff ~120) AND the translucent circular ring fills
                // (diff ~28-46 over white); the black-on-white tree ink is neutral (diff 0), so it is never counted
                if ( ( ( max - min ) > 22 ) && ( max > 60 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Count near-light-grey pixels (the B&W geologic annulus fill ~235, and its epoch rings) -- neutral greys well
     *  above the black tree ink, below pure white. */
    private static int countLightGrey( final BufferedImage img ) {
        int n = 0;
        for ( int y = 0; y < img.getHeight(); ++y ) {
            for ( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                final int max = Math.max( r, Math.max( g, b ) );
                final int min = Math.min( r, Math.min( g, b ) );
                if ( ( ( max - min ) <= 8 ) && ( min >= 210 ) && ( max <= 245 ) ) { // neutral grey ~235, not white/black
                    ++n;
                }
            }
        }
        return n;
    }

    /** Count near-black (dark ink) pixels in the bottom {@code strip} rows -- where the calendar axis line/ticks/labels
     *  sit (below the reserve-compressed tips). Neutral dark, not colored. */
    private static int countDarkInBottomStrip( final BufferedImage img, final int strip ) {
        int n = 0;
        for ( int y = Math.max( 0, img.getHeight() - strip ); y < img.getHeight(); ++y ) {
            for ( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF, g = ( rgb >> 8 ) & 0xFF, b = rgb & 0xFF;
                if ( ( r < 110 ) && ( g < 110 ) && ( b < 110 ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Number of pixels that differ between two equal-size renders (a robust "this toggle changed the figure" check). */
    private static int diffPixels( final BufferedImage a, final BufferedImage b ) {
        if ( ( a.getWidth() != b.getWidth() ) || ( a.getHeight() != b.getHeight() ) ) {
            return Integer.MAX_VALUE;
        }
        int n = 0;
        for ( int y = 0; y < a.getHeight(); ++y ) {
            for ( int x = 0; x < a.getWidth(); ++x ) {
                if ( ( a.getRGB( x, y ) & 0xFFFFFF ) != ( b.getRGB( x, y ) & 0xFFFFFF ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [GeologicAxisRenderTest] " + msg );
        return false;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [GeologicAxisRenderTest] " + msg );
        ok[ 0 ] = false;
    }

    private GeologicAxisRenderTest() {
    }
}

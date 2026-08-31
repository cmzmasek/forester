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

import java.io.File;
import java.io.IOException;
import java.util.Locale;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;

/**
 * Renders a phylogeny to a figure file without a user ever seeing a window -- the engine behind the
 * {@code aptx_render} command-line tool.
 * <p>
 * <b>This is the public facade over a package-private world.</b> The exporters, {@link ExportSizeSpec} and the
 * display toggles all live inside this package, so a tool in {@code org.forester.application} cannot reach them;
 * everything GUI-coupled stays here, behind one entry point.
 * <p>
 * <b>It needs a display, and that is not a detail.</b> Archaeopteryx's display state IS its Swing control panel --
 * {@link TreePanel} reads it in over a hundred places -- so rendering means building the real panel. A window is
 * created but never shown. On a desktop that just works; on a headless server it needs {@code xvfb-run}. Making
 * this genuinely headless means decoupling display state from Swing, which is a refactor, not a flag.
 * <p>
 * <b>Reproducibility is the whole point, so the render must not inherit the operator's GUI preferences.</b> A tool
 * that silently picked up whatever the last person left set in Archaeopteryx would produce a different figure on
 * every machine -- which would defeat the feature. {@link #isolateSettings()} points the settings directory at a
 * throwaway location BEFORE any Archaeopteryx class loads, so the render starts from documented defaults and the
 * user's own saved settings are neither read nor written.
 */
public final class FigureRenderer {

    public enum Style {
        RECTANGULAR, CIRCULAR, UNROOTED
    }

    /** What the figure should look like. Null means "leave at the default". */
    public static final class Spec {

        public Style   style           = Style.RECTANGULAR;
        /** null = decide from the tree (a phylogram when it has branch lengths). */
        public Boolean phylogram       = null;
        public boolean support         = false;
        public boolean branch_lengths  = false;
        public String  color_by_ref    = null;
        /**
         * e.g. "170x120mm", "8x6in", "1200x900px"; null = size from {@link #width}/{@link #height}.
         * <p>
         * The default is a PHYSICAL size, not a pixel count, and deliberately so: the tree is laid out in
         * point-space, so a 12 pt font on a physically small page swamps it -- a pixel default like 1200x900 at
         * 300 dpi is only 4 x 3 inches, and figures came out cramped, with radial labels truncated to a few
         * characters. 180 x 130 mm is an ordinary journal figure and looks right with no flags at all.
         */
        public String  size            = "180x130mm";
        public int     dpi             = 300;
        public int     width           = 1200;
        public int     height          = 900;
    }

    private FigureRenderer() {
    }

    /**
     * Points the Archaeopteryx settings directory at a throwaway location. MUST be called before any
     * Archaeopteryx class is loaded -- i.e. the first statement of a tool's {@code main} -- or the saved
     * preferences will already have been read. See the class comment for why this matters.
     */
    public static void isolateSettings() {
        if ( System.getProperty( "archaeopteryx.cache.dir" ) == null ) {
            final File dir = new File( System.getProperty( "java.io.tmpdir" ),
                                       "aptx-render-" + System.nanoTime() );
            dir.mkdirs();
            System.setProperty( "archaeopteryx.cache.dir", dir.getAbsolutePath() );
            // NOT deleteOnExit(): that removes a directory only when it is EMPTY, and settings files get written
            // into this one as the JVM shuts down -- so a batch job would strew a directory per figure across the
            // temp folder. Delete the contents ourselves, then the directory.
            Runtime.getRuntime().addShutdownHook( new Thread( () -> {
                final File[] files = dir.listFiles();
                if ( files != null ) {
                    for ( final File f : files ) {
                        f.delete();
                    }
                }
                dir.delete();
            } ) );
        }
    }

    /** The export format for an output file, from its extension; null when the extension is not one we write. */
    public static GraphicsExportType formatFor( final File output ) {
        final String n = output.getName().toLowerCase( Locale.ROOT );
        final int dot = n.lastIndexOf( '.' );
        if ( dot < 0 ) {
            return null;
        }
        switch ( n.substring( dot + 1 ) ) {
            case "pdf":
                return GraphicsExportType.PDF;
            case "svg":
                return GraphicsExportType.SVG;
            case "eps":
                return GraphicsExportType.EPS;
            case "png":
                return GraphicsExportType.PNG;
            case "jpg":
            case "jpeg":
                return GraphicsExportType.JPG;
            case "tif":
            case "tiff":
                return GraphicsExportType.TIFF;
            default:
                return null;
        }
    }

    /**
     * Parses a size such as {@code 170x120mm}, {@code 8x6in} or {@code 1200x900px}.
     *
     * @return null if it does not parse, so the caller can report the user's own string back to them
     */
    static ExportSizeSpec parseSize( final String s, final int dpi ) {
        if ( ForesterUtilShim.isEmpty( s ) ) {
            return null;
        }
        final String t = s.trim().toLowerCase( Locale.ROOT );
        ExportSizeSpec.Unit unit;
        final String num;
        if ( t.endsWith( "mm" ) ) {
            unit = ExportSizeSpec.Unit.MILLIMETERS;
            num = t.substring( 0, t.length() - 2 );
        }
        else if ( t.endsWith( "in" ) ) {
            unit = ExportSizeSpec.Unit.INCHES;
            num = t.substring( 0, t.length() - 2 );
        }
        else if ( t.endsWith( "px" ) ) {
            unit = ExportSizeSpec.Unit.PIXELS;
            num = t.substring( 0, t.length() - 2 );
        }
        else {
            return null;
        }
        final int x = num.indexOf( 'x' );
        if ( x < 1 ) {
            return null;
        }
        try {
            final double w = Double.parseDouble( num.substring( 0, x ).trim() );
            final double h = Double.parseDouble( num.substring( x + 1 ).trim() );
            if ( !( w > 0 ) || !( h > 0 ) ) {
                return null;
            }
            return new ExportSizeSpec( unit, w, h, dpi );
        }
        catch ( final NumberFormatException e ) {
            return null;
        }
    }

    /**
     * Renders {@code phy} into {@code output}. Returns the export report (the same wording the GUI shows).
     *
     * @throws IOException on a write failure, an unknown output extension, or an unparseable size
     */
    public static String render( final Phylogeny phy, final File output, final Spec spec ) throws IOException {
        if ( ( phy == null ) || phy.isEmpty() ) {
            throw new IOException( "tree is empty" );
        }
        final GraphicsExportType type = formatFor( output );
        if ( type == null ) {
            throw new IOException( "cannot tell the format from the file name \"" + output.getName()
                    + "\" (use .pdf, .svg, .eps, .png, .jpg or .tiff)" );
        }
        // ALWAYS go through an ExportSizeSpec, even when no size was asked for (then: the pixel size at the
        // requested DPI). One path instead of two, and it buys two things the second path did not have: the
        // radial layouts get FITTED to the output canvas (layoutForExportSize), and -dpi actually does something
        // -- the unsized raster export used to size from the GUI's raster-export SCALE and ignore DPI entirely.
        final ExportSizeSpec size;
        if ( !ForesterUtilShim.isEmpty( spec.size ) ) {
            size = parseSize( spec.size, spec.dpi );
            if ( size == null ) {
                throw new IOException( "cannot read the size \"" + spec.size
                        + "\" (expected something like 170x120mm, 8x6in or 1200x900px)" );
            }
        }
        else {
            size = new ExportSizeSpec( ExportSizeSpec.Unit.PIXELS, spec.width, spec.height, spec.dpi );
        }
        // ON THE EDT. The frame is showing (it has to be -- see below), so the event thread is painting it; doing
        // the layout and the export from another thread lets a repaint re-run calcParametersForPainting in
        // between and quietly change the figure, which would defeat the determinism this class exists for.
        if ( javax.swing.SwingUtilities.isEventDispatchThread() ) {
            return renderOnEdt( phy, output, type, size, spec );
        }
        final String[] out = new String[ 1 ];
        final IOException[] err = new IOException[ 1 ];
        try {
            javax.swing.SwingUtilities.invokeAndWait( () -> {
                try {
                    out[ 0 ] = renderOnEdt( phy, output, type, size, spec );
                }
                catch ( final IOException e ) {
                    err[ 0 ] = e;
                }
            } );
        }
        catch ( final Exception e ) {
            throw new IOException( e.getMessage(), e );
        }
        if ( err[ 0 ] != null ) {
            throw err[ 0 ];
        }
        return out[ 0 ];
    }

    private static String renderOnEdt( final Phylogeny phy,
                                       final File output,
                                       final GraphicsExportType type,
                                       final ExportSizeSpec size,
                                       final Spec spec )
            throws IOException {
        final boolean previously_rendering_only = MainFrame.isRenderingOnly();
        MainFrame mf = null;
        try {
            MainFrame.setRenderingOnly( true ); // no load-time dialogs: a modal would hang a pipeline forever
            mf = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration() );
            // The layout measures label widths through the panel's own Graphics
            // (TreePanel.calculateLongestExtNodeInfo), and Swing hands one out only for a SHOWING component --
            // being merely displayable (pack()) is not enough. So realize it far off-screen: the measurement is
            // exactly the GUI's, and nothing appears in front of anyone. Under xvfb there is no screen anyway.
            mf.setLocation( -32000, -32000 );
            mf.setVisible( true );
            final MainPanel mp = mf.getMainPanel();
            final TreePanel tp = mp.getCurrentTreePanel();
            applyStyle( tp, mp.getControlPanel(), spec );
            // Re-fit after the style change. Switching TO a radial type invalidates the radial diameter, and the
            // GUI re-fits at that point (MainFrame.typeChanged) precisely because otherwise the first radial
            // layout uses a stale rectangular preferred size and the circle comes out off-centre.
            mp.getControlPanel().showWhole();
            final int w = size.layoutWidthPt();
            final int h = size.layoutHeightPt();
            // Size the panel so the radial paths centre on the output canvas -- but do NOT lay out here. The
            // exporters call layoutForExportSize() themselves, and laying out first at the hidden frame's size
            // CACHES label widths (calculateLongestExtNodeInfo) computed for a much smaller radial diameter:
            // circular figures then came out with a small ring and labels truncated to a few characters even
            // though the page had room for them.
            tp.setSize( w, h );
            // PDF has its OWN exporter -- neither AptxUtil entry point handles it (the GUI routes it separately
            // too), and passing PDF to the raster path silently writes nothing at all.
            if ( type == GraphicsExportType.PDF ) {
                final int[] token = tp.layoutForExportSize( w, h );
                try {
                    return PdfExporter.writePhylogenyToPdfExactSize( output.getAbsolutePath(), tp, w, h,
                                                                     tp.getOptions()
                                                                             .isGraphicsExportWhiteBackground() );
                }
                finally {
                    tp.restoreLayoutAfterExport( token );
                }
            }
            return AptxUtil.writePhylogenyToGraphicsFileAtSize( output.getAbsolutePath(), tp, type,
                                                                tp.getOptions(), size );
        }
        finally {
            if ( mf != null ) {
                mf.dispose();
            }
            // RESTORE, do not just clear: this is a process-wide switch and FigureRenderer ships in the jar, so
            // leaving it on would silently kill the load-time offers for every tree a host JVM opens afterwards
            // (including later tests in the suite).
            MainFrame.setRenderingOnly( previously_rendering_only );
        }
    }

    private static void applyStyle( final TreePanel tp, final ControlPanel cp, final Spec spec )
            throws IOException {
        switch ( spec.style ) {
            case CIRCULAR:
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                break;
            case UNROOTED:
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                break;
            default:
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                break;
        }
        if ( spec.phylogram != null ) {
            cp.setTreeDisplayType( spec.phylogram.booleanValue() ? Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM
                    : Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
        }
        if ( spec.support ) {
            cp.setCheckbox( DisplayOption.WRITE_CONFIDENCE_VALUES, true );
        }
        if ( spec.branch_lengths ) {
            cp.setCheckbox( DisplayOption.WRITE_BRANCH_LENGTH_VALUES, true );
        }
        if ( !ForesterUtilShim.isEmpty( spec.color_by_ref ) ) {
            // Refuse a ref the tree cannot colour by. Setting it regardless produces an UNCOLOURED figure, exit
            // code 0 and no message -- the worst possible outcome for an unattended pipeline, and a typo away.
            final java.util.List<String> refs = PropertyColorScheme.colorableRefs( tp.getPhylogeny() );
            if ( !refs.contains( spec.color_by_ref ) ) {
                throw new IOException( "cannot colour by \"" + spec.color_by_ref + "\"; this tree offers "
                        + ( refs.isEmpty() ? "no colourable properties" : refs.toString() ) );
            }
            cp.demoSelectColorByProperty( spec.color_by_ref );
        }
    }

    /** Tiny local helper so this class does not depend on the util package just for an emptiness test. */
    private static final class ForesterUtilShim {

        static boolean isEmpty( final String s ) {
            return ( s == null ) || ( s.trim().length() < 1 );
        }
    }
}

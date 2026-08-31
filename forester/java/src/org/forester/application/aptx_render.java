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

package org.forester.application;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.forester.archaeopteryx.FigureRenderer;
import org.forester.io.parsers.PhylogenyParser;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.factories.PhylogenyFactory;
import org.forester.util.CommandLineArguments;
import org.forester.util.ForesterUtil;

/**
 * Renders a phylogeny to a publication-ready figure from the command line -- the same renderer and the same
 * exporters Archaeopteryx uses interactively, so what a pipeline produces matches what the GUI shows.
 * <p>
 * The figure is built from documented DEFAULTS plus the flags given here; it deliberately does not inherit
 * whatever the operator last set in the GUI, or the same command would produce a different figure on a
 * different machine. See {@link FigureRenderer}.
 */
public final class aptx_render {

    final static private String PRG_NAME    = "aptx_render";
    final static private String PRG_VERSION = "1.0.0";
    final static private String PRG_DESC    = "render a phylogeny to a publication-ready figure";
    final static private String PRG_DATE    = "2026-08-30";
    final static private String E_MAIL      = "phyloxml@gmail.com";
    final static private String WWW         = "https://cmzmasek.github.io/archaeopteryx/";

    private aptx_render() {
    }

    public static void main( final String args[] ) {
        // FIRST, before any Archaeopteryx class loads: render from defaults, not from the operator's own
        // saved GUI settings -- otherwise the same command gives a different figure on a different machine.
        FigureRenderer.isolateSettings();
        ForesterUtil.printProgramInformation( PRG_NAME, PRG_DESC, PRG_VERSION, PRG_DATE, E_MAIL, WWW,
                                              ForesterUtil.getForesterLibraryInformation() );
        CommandLineArguments cla = null;
        try {
            cla = new CommandLineArguments( args );
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        if ( cla.isOptionSet( "help" ) || cla.isOptionSet( "h" ) || ( args.length == 0 ) ) {
            printHelp();
            System.exit( 0 );
        }
        if ( cla.getNumberOfNames() != 2 ) {
            printHelp();
            ForesterUtil.fatalError( PRG_NAME, "expected two arguments: <input tree> <output figure>" );
        }
        final List<String> allowed = new ArrayList<String>();
        allowed.add( "size" );
        allowed.add( "dpi" );
        allowed.add( "style" );
        allowed.add( "cladogram" );
        allowed.add( "phylogram" );
        allowed.add( "support" );
        allowed.add( "bl" );
        allowed.add( "color" );
        allowed.add( "help" );
        allowed.add( "h" );
        final String bad = cla.validateAllowedOptionsAsString( allowed );
        if ( bad.length() > 0 ) {
            ForesterUtil.fatalError( PRG_NAME, "unknown option(s): " + bad );
        }
        final File in = cla.getFile( 0 );
        final File out = cla.getFile( 1 );
        if ( !in.exists() ) {
            ForesterUtil.fatalError( PRG_NAME, "no such file: " + in );
        }
        if ( java.awt.GraphicsEnvironment.isHeadless() ) {
            ForesterUtil.fatalError( PRG_NAME,
                                     "no display available. Archaeopteryx's display settings live in its "
                                             + "control panel, so rendering needs a display even though no "
                                             + "window is shown. On a headless server, run this under a "
                                             + "virtual display, e.g.:  xvfb-run -a " + PRG_NAME + " ..." );
        }
        final FigureRenderer.Spec spec = new FigureRenderer.Spec();
        try {
            if ( cla.isOptionSet( "size" ) ) {
                spec.size = cla.getOptionValueAsCleanString( "size" );
            }
            if ( cla.isOptionSet( "dpi" ) ) {
                spec.dpi = cla.getOptionValueAsInt( "dpi" );
            }
            if ( cla.isOptionSet( "style" ) ) {
                final String s = cla.getOptionValueAsCleanString( "style" ).toLowerCase();
                if ( s.startsWith( "circ" ) ) {
                    spec.style = FigureRenderer.Style.CIRCULAR;
                }
                else if ( s.startsWith( "unroot" ) ) {
                    spec.style = FigureRenderer.Style.UNROOTED;
                }
                else if ( s.startsWith( "rect" ) ) {
                    spec.style = FigureRenderer.Style.RECTANGULAR;
                }
                else {
                    ForesterUtil.fatalError( PRG_NAME, "unknown style \"" + s
                            + "\" (expected rectangular, circular or unrooted)" );
                }
            }
            if ( cla.isOptionSet( "cladogram" ) && cla.isOptionSet( "phylogram" ) ) {
                ForesterUtil.fatalError( PRG_NAME, "-cladogram and -phylogram are mutually exclusive" );
            }
            if ( cla.isOptionSet( "cladogram" ) ) {
                spec.phylogram = Boolean.FALSE;
            }
            if ( cla.isOptionSet( "phylogram" ) ) {
                spec.phylogram = Boolean.TRUE;
            }
            spec.support = cla.isOptionSet( "support" );
            spec.branch_lengths = cla.isOptionSet( "bl" );
            if ( cla.isOptionSet( "color" ) ) {
                spec.color_by_ref = cla.getOptionValueAsCleanString( "color" );
            }
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        try {
            final PhylogenyFactory factory = ParserBasedPhylogenyFactory.getInstance();
            final PhylogenyParser parser = ParserUtils.createParserDependingOnFileType( in, true );
            final Phylogeny[] phys = factory.create( in, parser );
            if ( ( phys == null ) || ( phys.length < 1 ) ) {
                ForesterUtil.fatalError( PRG_NAME, "no tree found in: " + in );
            }
            if ( phys.length > 1 ) {
                ForesterUtil.printWarningMessage( PRG_NAME, in + " holds " + phys.length
                        + " trees; rendering the first" );
            }
            final String report = FigureRenderer.render( phys[ 0 ], out, spec );
            System.out.println();
            System.out.println( report );
            System.out.println();
        }
        catch ( final Exception e ) {
            ForesterUtil.fatalError( PRG_NAME, e.getMessage() );
        }
        // The AWT threads are not daemons and a still-displayable window keeps them alive, so a plain return
        // can leave the JVM running with nothing to do -- fatal in a pipeline. Leave deliberately.
        System.exit( 0 );
    }

    private static void printHelp() {
        System.out.println();
        System.out.println( "Usage:" );
        System.out.println();
        System.out.println( "  " + PRG_NAME + " [options] <input tree> <output figure>" );
        System.out.println();
        System.out.println( "  The output FORMAT comes from the file extension:" );
        System.out.println( "    .pdf .svg .eps        vector, for a manuscript or Illustrator" );
        System.out.println( "    .png .jpg .tiff       raster (PNG carries the true DPI)" );
        System.out.println();
        System.out.println( "Options:" );
        System.out.println( "  -size=<WxH><unit>  figure size, e.g. 170x120mm, 8x6in, 1200x900px" );
        System.out.println( "                     (default 180x130mm -- an ordinary journal figure)" );
        System.out.println( "  -dpi=<n>           dots per inch (default 300). With a mm or inch size this sets" );
        System.out.println( "                     the pixel count; with a pixel size it sets the physical size" );
        System.out.println( "  -style=<s>         rectangular (default), circular or unrooted" );
        System.out.println( "  -phylogram         draw branch lengths to scale" );
        System.out.println( "  -cladogram         ignore branch lengths" );
        System.out.println( "                     (default: a phylogram when the tree has branch lengths)" );
        System.out.println( "  -support           show confidence/support values" );
        System.out.println( "  -bl                show branch-length values" );
        System.out.println( "  -color=<ref>       colour tips by a property, e.g. data:host" );
        System.out.println( "  -help              this message" );
        System.out.println();
        System.out.println( "Examples:" );
        System.out.println();
        System.out.println( "  " + PRG_NAME + " tree.xml figure.pdf" );
        System.out.println( "  " + PRG_NAME + " -size=170x120mm -dpi=300 -support tree.xml figure.png" );
        System.out.println( "  " + PRG_NAME + " -style=circular -color=data:host tree.xml figure.svg" );
        System.out.println();
        System.out.println( "Note: a display is required (no window is shown). On a headless server:" );
        System.out.println();
        System.out.println( "  xvfb-run -a " + PRG_NAME + " tree.xml figure.pdf" );
        System.out.println();
    }
}

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
import java.io.File;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeSet;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.phylogeny.data.RenderableDomainArchitecture;
import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;

/**
 * Protein-domain colours come from the palette and nothing else. Every drawn domain name gets a palette colour when a
 * tree loads; a name first met later takes the NEXT unused palette colour, never one computed from the name's
 * characters (that let SH2 and SH3 collide); and opening another window leaves every open tree's colours alone. The
 * last was a bug: each new window replaced the shared palette with an always-empty configuration map, so every open
 * tree fell back to the character hash.
 */
public final class DomainPaletteTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "DomainPalette: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return unknownNameTakesNextPaletteColour()
                    && ( GraphicsEnvironment.isHeadless() || secondWindowKeepsColours() );
        }
        catch ( final Throwable e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    /** A name the palette has not met takes the next unused palette colour, keeps it, and collides with nothing. */
    private static boolean unknownNameTakesNextPaletteColour() {
        RenderableDomainArchitecture
                .setColorMap( AptxUtil.assignDistinctColors( new TreeSet<String>( Arrays.asList( "Pkinase", "SH2" ) ) ) );
        if ( !AptxUtil.paletteColor( 0 ).equals( RenderableDomainArchitecture.colorFor( "Pkinase" ) )
                || !AptxUtil.paletteColor( 1 ).equals( RenderableDomainArchitecture.colorFor( "SH2" ) ) ) {
            return fail( "palette names must keep their palette colours" );
        }
        final Color sh3 = RenderableDomainArchitecture.colorFor( "SH3" );
        final Color ph = RenderableDomainArchitecture.colorFor( "PH" );
        if ( !AptxUtil.paletteColor( 2 ).equals( sh3 ) || !AptxUtil.paletteColor( 3 ).equals( ph ) ) {
            return fail( "unmet names must take the next unused palette colours, got " + sh3 + " " + ph );
        }
        if ( !sh3.equals( RenderableDomainArchitecture.colorFor( "SH3" ) ) ) {
            return fail( "an unmet name must keep the colour it was given" );
        }
        if ( new TreeSet<Integer>( Arrays.asList( AptxUtil.paletteColor( 0 ).getRGB(), AptxUtil.paletteColor( 1 ).getRGB(),
                                                  sh3.getRGB(), ph.getRGB() ) ).size() != 4 ) {
            return fail( "four domain names must get four different colours" );
        }
        if ( !Color.GRAY.equals( RenderableDomainArchitecture.colorFor( null ) ) ) {
            return fail( "a nameless domain must be grey" );
        }
        // the map may start out null (no tree loaded yet): the first names still come from the palette
        RenderableDomainArchitecture.setColorMap( null );
        if ( !AptxUtil.paletteColor( 0 ).equals( RenderableDomainArchitecture.colorFor( "Ig" ) ) ) {
            return fail( "with no palette yet, the first name must take the first palette colour" );
        }
        return true;
    }

    /** Opening a second window must not change the domain colours of the tree in the first. */
    private static boolean secondWindowKeepsColours() throws Exception {
        final File f = new File( System.getProperty( "user.dir" ), "forester/demo/domain-architectures.xml" );
        final Phylogeny dom = ParserUtils.readPhylogenies( f )[ 0 ];
        final TreeSet<String> names = new TreeSet<String>();
        for( final PhylogenyNode n : dom.getExternalNodes() ) {
            if ( n.getNodeData().isHasSequence() ) {
                final DomainArchitecture da = n.getNodeData().getSequence().getDomainArchitecture();
                if ( da != null ) {
                    for( int i = 0; i < da.getNumberOfDomains(); ++i ) {
                        names.add( da.getDomain( i ).getName() );
                    }
                }
            }
        }
        if ( names.size() < 2 ) {
            return fail( "the demo tree should carry several domain names, got " + names );
        }
        final MainFrame[] mf = new MainFrame[ 2 ];
        final Map<String, Color> before = new LinkedHashMap<String, Color>();
        final Map<String, Color> after = new LinkedHashMap<String, Color>();
        try {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { dom }, new Configuration(), "palette1" ) );
            SwingUtilities.invokeAndWait( () -> {
                for( final String s : names ) {
                    before.put( s, RenderableDomainArchitecture.colorFor( s ) );
                }
            } );
            if ( !before.equals( AptxUtil.assignDistinctColors( names ) ) ) {
                return fail( "the loaded tree's domains should take the palette, got " + before );
            }
            final Phylogeny plain = Phylogeny.createInstanceFromNhxString( "((a,b),c);" );
            SwingUtilities.invokeAndWait( () -> mf[ 1 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { plain }, new Configuration(), "palette2" ) );
            // Ask in REVERSE order. A wiped map hands out palette colours in the order names are asked for, and asking
            // in sorted order reproduces the palette exactly -- a sabotage that restored the wipe passed that way.
            // The renderer asks in tree order, so a wipe really does reshuffle the colours; reverse order exposes it.
            SwingUtilities.invokeAndWait( () -> {
                for( final String s : names.descendingSet() ) {
                    after.put( s, RenderableDomainArchitecture.colorFor( s ) );
                }
            } );
            if ( !after.equals( before ) ) {
                return fail( "opening a second window changed the first window's domain colours: before " + before
                        + ", after " + after );
            }
        }
        finally {
            SwingUtilities.invokeAndWait( () -> {
                for( final MainFrame m : mf ) {
                    if ( m != null ) {
                        ( (JFrame) m ).dispose();
                    }
                }
            } );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [DomainPaletteTest] " + message );
        return false;
    }

    private DomainPaletteTest() {
        // not instantiable
    }
}

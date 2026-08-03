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
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Verifies the "Color by" palette default plumbing behind palette persistence: choosing a palette on a tree panel
 * records it as the shared default ({@link Options#getColorPaletteName()}), and a NEW tab then seeds its palette
 * from that default. That write-back + seed is what lets {@link GuiPreferences} persist the palette (it reads/writes
 * Options, not the per-panel field). Guarded to a no-op when headless (needs FlatLaf via {@code createInstance}).
 */
public final class ColorPaletteDefaultTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ColorPaletteDefault: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            // need at least 2 palettes to switch between; guard so the test fails loudly rather than passing vacuously
            final List<String> palettes = PropertyColorScheme.paletteNames();
            if ( palettes.size() < 2 ) {
                return fail( "expected at least 2 palettes to exercise the palette default" );
            }
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree() }, conf, "pal" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    final TreePanel tp1 = frame.getMainPanel().getCurrentTreePanel();
                    // choose a target DIFFERENT from the panel's current palette, so the change actually crosses
                    // setColorPaletteName's "no-op if unchanged" guard regardless of any persisted starting value
                    // (run standalone, createInstance reads the real ~/.archaeopteryx; the suite isolates it)
                    String alt = palettes.get( 0 );
                    if ( alt.equals( tp1.getColorPaletteName() ) ) {
                        alt = palettes.get( 1 );
                    }
                    // choosing a palette on the first panel must record it as the shared (Options) default
                    tp1.setColorPaletteName( alt );
                    if ( !alt.equals( frame.getOptions().getColorPaletteName() ) ) {
                        fail( ok, "setColorPaletteName must write the choice back to Options (got "
                                + frame.getOptions().getColorPaletteName() + ")" );
                    }
                    // a NEW tab must seed its palette from that shared default
                    AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { tree() }, "t2", "", conf, frame.getMainPanel() );
                    final TreePanel tp2 = frame.getMainPanel().getCurrentTreePanel();
                    if ( !alt.equals( tp2.getColorPaletteName() ) ) {
                        fail( ok, "a new tab must seed its palette from the shared default (got "
                                + tp2.getColorPaletteName() + ")" );
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

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [ColorPaletteDefaultTest] " + msg );
        ok[ 0 ] = false;
    }

    private static boolean fail( final String msg ) {
        System.out.println( "  [ColorPaletteDefaultTest] " + msg );
        return false;
    }

    private static Phylogeny tree() {
        final PhylogenyNode root = new PhylogenyNode();
        for ( int i = 0; i < 3; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private ColorPaletteDefaultTest() {
    }
}

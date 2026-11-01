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

import java.awt.Dimension;
import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.util.Arrays;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;

/**
 * The font-size bounds (Christian, 2026-10-31): every legend's font follows the small tree font within [11, 16] pt
 * ({@link TreePanel#legendFontSize}); the tree font is kept within [2, 32] pt on every way in
 * ({@link TreeFontSet#setBaseFont}), and the font chooser offers exactly that range and shows a typed size kept in it.
 */
public final class FontBoundsTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "FontBounds: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return legendSizes() & chooser() & inAWindow();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** The rule itself, pinned by literals: the floor, the band it follows the tree in, the ceiling. */
    private static boolean legendSizes() {
        boolean ok = true;
        final float[][] cases = { { 0, 11 }, { 9, 11 }, { 11, 11 }, { 12, 12 }, { 14, 14 }, { 16, 16 }, { 17, 16 },
                { 30, 16 }, { 70, 16 }, { 13.5f, 13.5f } };
        for( final float[] c : cases ) {
            ok &= ck( TreePanel.legendFontSize( c[ 0 ] ) == c[ 1 ], "a small tree font of " + c[ 0 ] + " pt: legend "
                    + c[ 1 ] + " pt, got " + TreePanel.legendFontSize( c[ 0 ] ) );
        }
        return ok;
    }

    /** Every way into the tree font lands within [2, 32]; a size inside is left alone; "bigger" at the top stays. */
    private static boolean treeFontRange( final TreeFontSet tfs ) {
        final boolean[] ok = { true };
        {
            final Font base = tfs.getBaseFont();
            tfs.setBaseFont( base.deriveFont( 72f ) ); // what the font chooser used to pass on
            ok[ 0 ] &= ck( tfs.getUserFontSize() == 32, "72 pt in: 32 pt, got " + tfs.getUserFontSize() );
            tfs.increaseUserFontSize();
            ok[ 0 ] &= ck( tfs.getUserFontSize() == 32, "'bigger' at 32 stays 32 (it used to fall from 72 to 32)" );
            tfs.decreaseUserFontSize();
            ok[ 0 ] &= ck( tfs.getUserFontSize() == 31, "'smaller' from 32: 31" );
            tfs.setBaseFont( base.deriveFont( 1f ) );
            ok[ 0 ] &= ck( tfs.getUserFontSize() == 2, "1 pt in: 2 pt, got " + tfs.getUserFontSize() );
            final Font twenty = base.deriveFont( Font.BOLD, 20f );
            tfs.setBaseFont( twenty );
            ok[ 0 ] &= ck( ( tfs.getBaseFont() == twenty ) && ( tfs.getSmallFont().getSize() == 18 ),
                           "20 pt bold in: kept as given, the small font 18" );
            tfs.setUserFontSize( 2 );
            ok[ 0 ] &= ck( tfs.getSmallFont().getSize() == 2, "at 2 pt the small font is 2 pt, not 0 (got "
                    + tfs.getSmallFont().getSize() + ")" );
            tfs.setUserFontSize( 3 );
            ok[ 0 ] &= ck( tfs.getSmallFont().getSize() == 2, "at 3 pt the small font is 2 pt, not 1" );
            ok[ 0 ] &= ck( ( TreeFontSet.clampSize( -9 ) == 2 ) && ( TreeFontSet.clampSize( 2 ) == 2 )
                           && ( TreeFontSet.clampSize( 17 ) == 17 ) && ( TreeFontSet.clampSize( 32 ) == 32 )
                           && ( TreeFontSet.clampSize( 10000 ) == 32 ), "clampSize: the one range rule" );
            ok[ 0 ] &= ck( ( TreeFontSet.withinRange( base.deriveFont( 2f ) ).getSize() == 2 )
                           && ( TreeFontSet.withinRange( base.deriveFont( 32f ) ).getSize() == 32 ),
                           "the ends of the range are inside it" );
            tfs.setBaseFont( base );
        }
        return ok[ 0 ];
    }

    /** The font chooser: its list is exactly 2..32, a typed size is kept in range and shown so, and opening it on a
     *  size and pressing OK keeps that size (a list with gaps snapped 13 to 14). */
    private static boolean chooser() throws Exception {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // a dialog needs a display
        }
        final boolean[] ok = { true };
        SwingUtilities.invokeAndWait( () -> {
            final String[] want = new String[ 31 ];
            for( int i = 0; i < 31; ++i ) {
                want[ i ] = Integer.toString( i + 2 );
            }
            final FontChooser a = new FontChooser( new Font( Font.SANS_SERIF, Font.PLAIN, 20 ) );
            ok[ 0 ] &= ck( Arrays.equals( a.sizeChoicesForTest(), want ), "the list offers 2..32: "
                    + Arrays.toString( a.sizeChoicesForTest() ) );
            ok[ 0 ] &= ck( ( a.typeSizeForTest( "72" ) == 32 ) && "32".equals( a.sizeFieldForTest() ),
                           "typed 72: 32, and the field says so" );
            ok[ 0 ] &= ck( ( a.typeSizeForTest( "1" ) == 2 ) && "2".equals( a.sizeFieldForTest() ), "typed 1: 2" );
            ok[ 0 ] &= ck( a.typeSizeForTest( "17" ) == 17, "typed 17: 17" );
            ok[ 0 ] &= ck( ( a.typeSizeForTest( "abc" ) == 17 ) && "17".equals( a.sizeFieldForTest() ),
                           "typed text: the size there was" );
            a.dispose();
            // as Settings > Tree font opens it (MainFrame.chooseFont): a default chooser, then the tree's font; the
            // real modal dialog, OK pressed while it shows
            final FontChooser b = new FontChooser();
            b.setFont( new Font( Font.SANS_SERIF, Font.PLAIN, 13 ) );
            final javax.swing.Timer press = new javax.swing.Timer( 300, e -> b.pressOkForTest() );
            press.setRepeats( false );
            press.start();
            b.showDialog( null, "font bounds" );
            ok[ 0 ] &= ck( b.getFont().getSize() == 13, "opened on 13 pt, OK: 13 pt, got " + b.getFont().getSize() );
            b.dispose();
            // a family the list does not hold (the list is read once, maybe before the bundled font was registered):
            // the chooser keeps it, never the last family of the list
            final Font unlisted = new Font( Font.SANS_SERIF, Font.PLAIN, 14 );
            final FontChooser c = new FontChooser( unlisted, new String[] { "Alpha", "Omega" } );
            c.selectCurrentForTest();
            ok[ 0 ] &= ck( unlisted.getFamily().equals( c.familyForTest() ),
                           "an unlisted family is kept: " + unlisted.getFamily() + ", got " + c.familyForTest() );
            c.pressOkForTest();
            ok[ 0 ] &= ck( unlisted.getFamily().equals( c.getFont().getFamily() ) && ( c.getFont().getSize() == 14 ),
                           "OK applies it at its size: " + c.getFont() );
            c.dispose();
            final FontChooser d = new FontChooser( unlisted, new String[] { "Alpha", unlisted.getFamily(), "Omega" } );
            d.selectCurrentForTest();
            ok[ 0 ] &= ck( unlisted.getFamily().equals( d.familyForTest() ), "control: a listed family is selected" );
            d.dispose();
        } );
        return ok[ 0 ];
    }

    /** A 20-row Color-by legend in a window as it opens: it grows with the tip labels up to 18 pt (legend 16), not
     *  past, and at the slider's 32 pt it still fits the panel. */
    private static boolean inAWindow() throws Exception {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 50; ++i ) {
            final PhylogenyNode n = new PhylogenyNode();
            n.setName( "tip_" + i );
            n.getNodeData().setProperties( new PropertiesList() );
            n.getNodeData().getProperties().addProperty( new Property( "data:group", String.format( "group %02d", i % 25 ),
                                                                       "", "xsd:string", Property.AppliesTo.NODE ) );
            root.addAsChild( n );
        }
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "font bounds" ) );
        final boolean[] ok = { true };
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                final TreeFontSet tfs = mf[ 0 ].getMainPanel().getTreeFontSet();
                ok[ 0 ] &= treeFontRange( tfs );
                tp.setColorByPropertyRef( "data:group" );
                tfs.setUserFontSize( 12 );
                ok[ 0 ] &= ck( tp.legendFontSizeForTest() == 11, "tip labels 12 pt: the legend at the 11 pt floor" );
                tfs.setUserFontSize( 16 );
                final float at16 = tp.legendFontSizeForTest();
                final Dimension d16 = tp.sharedLegendSizeForTest();
                tfs.setUserFontSize( 18 );
                final Dimension d18 = tp.sharedLegendSizeForTest();
                ok[ 0 ] &= ck( ( at16 == 14 ) && ( tp.legendFontSizeForTest() == 16 ),
                               "tip labels 16 and 18 pt: the legend follows, 14 and 16 pt" );
                ok[ 0 ] &= ck( ( d16 != null ) && ( d18 != null ) && ( d18.height > d16.height ),
                               "precondition: below the ceiling the legend grows with the font (" + d16 + " -> " + d18 + ")" );
                tfs.setUserFontSize( 32 );
                final Dimension d32 = tp.sharedLegendSizeForTest();
                ok[ 0 ] &= ck( ( tp.legendFontSizeForTest() == 16 ) && d18.equals( d32 ),
                               "tip labels 32 pt: the legend stays at 16 pt, the same box as at 18 (" + d32 + ")" );
                ok[ 0 ] &= ck( d32.height < tp.getVisibleRect().height,
                               "a 20-row legend at 32 pt tip labels fits the panel as the window opens (" + d32.height
                                       + " of " + tp.getVisibleRect().height + " px)" );
                tfs.setUserFontSize( 12 );
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    private static boolean ck( final boolean condition, final String what ) {
        if ( !condition ) {
            System.out.println( "  [FontBoundsTest] " + what );
        }
        return condition;
    }

    private FontBoundsTest() {
        // not instantiable
    }
}

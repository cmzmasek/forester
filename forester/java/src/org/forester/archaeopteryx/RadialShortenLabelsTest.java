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
import java.util.HashMap;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * A radial tip label is the text the layout measured for it (Christian, 2026-11-01: "A, shorten them in radial
 * too"): Shorten Labels shortens the node name in circular and unrooted as in rectangular (and the JS), and a name
 * that repeats the taxonomy shown is left out there too. Until 0.11.198 the radial labels drew the full name while
 * the layout measured the shortened one -- on the crowded clade-boxes demo 130 px drawn against 94 reserved. Each
 * drawn label's box is compared with the width the layout reserves for it ({@code labelTextWidth}), in both radial
 * layouts, with the option on and off. A GUI test, a no-op on a headless box.
 */
public final class RadialShortenLabelsTest {

    private static final String DEMO = "forester/demo/clade-boxes-crowded.xml";
    private static final String LONG = "Avihepatovirus isolate 03";

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Radial labels as measured (Shorten Labels, a name repeating the taxonomy): "
                + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            boolean ok = with( FigureRenderer.readTrees( new File( System.getProperty( "user.dir" ), DEMO ) )[ 0 ],
                               RadialShortenLabelsTest::shortened );
            ok &= with( duplicateTree(), RadialShortenLabelsTest::duplicate );
            return ok;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Shorten Labels on (the default) and off: every drawn label is as wide as the layout measured it. */
    private static boolean shortened( final MainFrame mf ) {
        final TreePanel tp = mf.getCurrentTreePanel();
        boolean ok = true;
        final PhylogenyNode avi = tp.getPhylogeny().getNode( LONG );
        ok &= ck( tp.getControlPanel().isShortenLabels() && !tp.nodeLabelTextForTest( avi ).contains( LONG ),
                  "precondition: Shorten Labels is on and shortens \"" + LONG + "\" (" + tp.nodeLabelTextForTest( avi )
                          + ")" );
        final Map<String, Double> on = new HashMap<String, Double>();
        for ( final PHYLOGENY_GRAPHICS_TYPE type : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
            on.put( type.name(), asMeasured( mf, tp, type, "on" ) );
            ok &= !Double.isNaN( on.get( type.name() ) );
        }
        tp.getControlPanel().setCheckbox( DisplayOption.SHORTEN_LABELS, false );
        ok &= ck( tp.nodeLabelTextForTest( avi ).contains( LONG ), "precondition: off, the full name" );
        for ( final PHYLOGENY_GRAPHICS_TYPE type : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
            final double off = asMeasured( mf, tp, type, "off" );
            ok &= ck( off > ( on.get( type.name() ) + 20 ), type + ": with Shorten Labels off, the long names are "
                    + "drawn in full, wider (the widest " + on.get( type.name() ) + " -> " + off + ")" );
        }
        tp.getControlPanel().setCheckbox( DisplayOption.SHORTEN_LABELS, true );
        return ok;
    }

    /** Lays {@code tp} out as {@code type}, paints it and checks each drawn tip label against its measure: no wider
     *  than measured (plus the one leading space a radial label without a taxonomy keeps) and no narrower. Returns
     *  the widest drawn label of a name long enough to be shortened (NaN when none was drawn or a check failed --
     *  Auto-hide keeps whichever labels fit, so which long name is drawn varies). */
    private static double asMeasured( final MainFrame mf, final TreePanel tp, final PHYLOGENY_GRAPHICS_TYPE type,
                                      final String option ) {
        tp.setPhylogenyGraphicsType( type );
        mf.showWhole();
        tp.setRecordLabelBoxesForTest( true );
        tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
        final int space = tp.getFontMetrics( mf.getMainPanel().getTreeFontSet().getLargeFont() ).stringWidth( " " );
        boolean ok = true;
        double avi = Double.NaN;
        int compared = 0;
        for ( final double[] b : tp.labelBoxesForTest() ) {
            final PhylogenyNode n = tp.getPhylogeny().getNode( ( long ) b[ 5 ] );
            if ( !n.isExternal() ) {
                continue;
            }
            ++compared;
            final double drawn = 2 * b[ 2 ];
            final int measured = tp.labelTextWidthForTest( n );
            ok &= ck( ( drawn <= ( measured + space + 1 ) ) && ( drawn >= ( measured - 1 ) ), type + ", Shorten Labels "
                    + option + ": " + n.getName() + " is drawn as wide as measured (" + drawn + " vs " + measured
                    + ")" );
            if ( n.getName().length() > AptxUtil.SHORTEN_NAME_MAX_LENGTH ) {
                avi = Double.isNaN( avi ) ? drawn : Math.max( avi, drawn );
            }
        }
        tp.setRecordLabelBoxesForTest( false );
        ok &= ck( compared >= 15, type + ", Shorten Labels " + option + ": precondition: tip labels drawn (" + compared
                + ")" );
        ok &= ck( !Double.isNaN( avi ), type + ", Shorten Labels " + option + ": precondition: a long name is drawn" );
        return ok ? avi : Double.NaN;
    }

    /** A tip named for the taxonomy it shows draws the taxonomy once, in the radial layouts as in rectangular. */
    private static boolean duplicate( final MainFrame mf ) {
        final TreePanel tp = mf.getCurrentTreePanel();
        boolean ok = true;
        final PhylogenyNode h = tp.getPhylogeny().getNode( "Homo sapiens" );
        ok &= ck( !tp.nodeLabelTextForTest( h ).contains( "Homo sapiens" ), "precondition: the measure leaves out a "
                + "name that repeats the taxonomy shown (" + tp.nodeLabelTextForTest( h ) + ")" );
        for ( final PHYLOGENY_GRAPHICS_TYPE type : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
            tp.setPhylogenyGraphicsType( type );
            mf.showWhole();
            tp.setRecordLabelBoxesForTest( true );
            tp.paintImmediately( 0, 0, tp.getWidth(), tp.getHeight() );
            double drawn = Double.NaN;
            for ( final double[] b : tp.labelBoxesForTest() ) {
                if ( ( long ) b[ 5 ] == h.getId() ) {
                    drawn = 2 * b[ 2 ];
                }
            }
            tp.setRecordLabelBoxesForTest( false );
            final int measured = tp.labelTextWidthForTest( h );
            ok &= ck( Math.abs( drawn - measured ) <= 1, type + ": \"Homo sapiens\" named for its taxonomy is drawn "
                    + "once, as wide as measured (" + drawn + " vs " + measured + ")" );
        }
        return ok;
    }

    private static Phylogeny duplicateTree() throws Exception {
        final PhylogenyNode root = new PhylogenyNode();
        final String[] names = { "Homo sapiens", "Pan troglodytes", "Gorilla gorilla", "Pongo abelii" };
        for ( final String name : names ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setName( name );
            tip.setDistanceToParent( 0.1 );
            final Taxonomy t = new Taxonomy();
            t.setScientificName( name );
            tip.getNodeData().setTaxonomy( t );
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private interface Exercise {

        boolean run( MainFrame mf );
    }

    private static boolean with( final Phylogeny phy, final Exercise ex ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "radial shorten" ) );
        final boolean[] ok = { false };
        try {
            SwingUtilities.invokeAndWait( () -> ok[ 0 ] = ex.run( mf[ 0 ] ) );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    private static boolean ck( final boolean b, final String msg ) {
        if ( !b ) {
            System.out.println( "  [RadialShortenLabelsTest] " + msg );
        }
        return b;
    }
}

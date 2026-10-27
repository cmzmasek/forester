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

import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * In the circular and unrooted layouts a tip label's text may reach as far as the ring's or fan's own reservation
 * leaves it, so the "Tree share" that splits the radius is also what decides where a label is cut (until 0.11.185
 * a fixed half of the half-diameter cut the text whatever the share: the slider moved the ring but never the cut,
 * and a tree of long scientific names opened with every long name cut -- Christian, 2026-10-16). Held: the cut
 * follows the share, in both radial layouts; ring + labels (+ the domain track) fit the canvas at every share; a
 * label shorter than its reach is never cut; the domain track rides outside the labels and is not taken from
 * their reach; the text is measured under the same cap it is drawn under.
 */
public final class RadialLabelReachTest {

    private static final int      W      = 900, H = 700;
    private static final double[] SHARES = { 0.25, 0.40, 0.80 };

    public static void main( final String[] args ) {
        System.out.println( test() ? "RadialLabelReachTest: OK." : "RadialLabelReachTest: FAILED." );
    }

    public static boolean test() {
        if ( java.awt.GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            cutFollowsTheShare( ok, PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            cutFollowsTheShare( ok, PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            shortLabelsNeverCut( ok );
            domainsRideOutsideTheLabels( ok );
        }
        catch ( final Exception e ) {
            e.printStackTrace();
            ok[ 0 ] = fail( "unexpected: " + e );
        }
        return ok[ 0 ];
    }

    /** Long labels: the cut moves with the share, the figure fits, and the ring keeps its share. */
    private static void cutFollowsTheShare( final boolean[] ok, final PHYLOGENY_GRAPHICS_TYPE type ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { star( 40, true, false ) }, new Configuration(), "reach " + type ) );
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( type );
                final int[] cap = new int[ SHARES.length ];
                final int[] text = new int[ SHARES.length ];
                for( int i = 0; i < SHARES.length; ++i ) {
                    tp.getOptions().setTreeWidthShare( SHARES[ i ] );
                    paint( tp );
                    cap[ i ] = tp.maxTipLabelTextWidthForTest();
                    text[ i ] = tp.lengthOfLongestTextOnlyForTest();
                    final int d = tp.radialDiameter();
                    if ( cap[ i ] <= 0 ) {
                        ok[ 0 ] = fail( type + " at " + SHARES[ i ] + ": the labels must get some reach, got " + cap[ i ] );
                    }
                    if ( !( text[ i ] <= cap[ i ] ) || ( text[ i ] < ( cap[ i ] - 24 ) ) ) {
                        ok[ 0 ] = fail( type + " at " + SHARES[ i ] + ": the longest label (far wider than any reach) must measure as cut "
                                + "at the reach " + cap[ i ] + ", measured " + text[ i ] );
                    }
                    if ( type == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                        final int avail = tp.circularAvailableRadiusForTest();
                        final int ring = tp.circularRadiusForTest();
                        if ( ( ring + cap[ i ] ) > avail ) {
                            ok[ 0 ] = fail( "circular at " + SHARES[ i ] + ": ring " + ring + " + label reach " + cap[ i ]
                                    + " must fit the available radius " + avail );
                        }
                        if ( ring < Math.floor( avail * SHARES[ i ] ) - 1 ) {
                            ok[ 0 ] = fail( "circular at " + SHARES[ i ] + ": the ring must keep the tree's share, " + ring + " of " + avail );
                        }
                    }
                    else {
                        // the fan reserves MOVE + labels per side within the labels' share of the half-diameter
                        final double per_side = ( d / 2.0 ) * ( 1.0 - SHARES[ i ] );
                        if ( ( cap[ i ] + 20 ) > ( per_side + 1 ) ) {
                            ok[ 0 ] = fail( "unrooted at " + SHARES[ i ] + ": margin + label reach " + ( cap[ i ] + 20 )
                                    + " must fit the labels' share of the half-diameter, " + per_side );
                        }
                    }
                }
                if ( !( ( cap[ 0 ] > cap[ 1 ] ) && ( cap[ 1 ] > cap[ 2 ] ) ) ) {
                    ok[ 0 ] = fail( type + ": the label reach must follow the share, 25/40/80% -> " + cap[ 0 ] + "/" + cap[ 1 ] + "/"
                            + cap[ 2 ] + " px" );
                }
            } );
        }
        finally {
            dispose( mf );
        }
    }

    /** Short labels in circular: never cut, at any share. */
    private static void shortLabelsNeverCut( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        final Phylogeny phy = star( 12, false, false );
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "reach short" ) );
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                int widest = 0;
                for( final PhylogenyNode tip : phy.getExternalNodes() ) {
                    widest = Math.max( widest, tp.getFontMetrics( tp.getMainPanel().getTreeFontSet().getLargeFont() )
                            .stringWidth( tip.getName() ) );
                }
                for( final double share : SHARES ) {
                    tp.getOptions().setTreeWidthShare( share );
                    paint( tp );
                    if ( !( widest < tp.maxTipLabelTextWidthForTest() ) ) {
                        ok[ 0 ] = fail( "fixture: the short names (" + widest + " px) must be under the reach at " + share + " ("
                                + tp.maxTipLabelTextWidthForTest() + ")" );
                    }
                    if ( tp.lengthOfLongestTextOnlyForTest() != widest ) {
                        ok[ 0 ] = fail( "a label under its reach is never cut: at " + share + " the widest measures "
                                + tp.lengthOfLongestTextOnlyForTest() + ", its text is " + widest );
                    }
                }
            } );
        }
        finally {
            dispose( mf );
        }
    }

    /** With domains drawn in circular, the track rides OUTSIDE the labels and is not taken from their reach: the
     *  labels come first (taking the track out of their reach cut 268 scientific names to 99 px where 182 had been,
     *  2026-10-16), the ring keeps its share, and the text is measured under the cap it is drawn under. */
    private static void domainsRideOutsideTheLabels( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { star( 20, true, true ) }, new Configuration(), "reach domains" ) );
        try {
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                tp.getOptions().setTreeWidthShare( 0.40 );
                tp.getMainPanel().getControlPanel().setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, false );
                paint( tp );
                final int without = tp.maxTipLabelTextWidthForTest();
                final int ring_without = tp.circularRadiusForTest();
                tp.getMainPanel().getControlPanel().setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, true );
                paint( tp );
                paint( tp ); // a second pass: the first measures the track
                final int with = tp.maxTipLabelTextWidthForTest();
                final int track = tp.circularDomainReachForTest();
                if ( track < 20 ) {
                    ok[ 0 ] = fail( "fixture: the domain track must be drawn and wide, got " + track + " px" );
                }
                if ( with != without ) {
                    ok[ 0 ] = fail( "the domain track is not taken from the labels' reach: " + without + " px without it, " + with + " with" );
                }
                if ( tp.lengthOfLongestTextOnlyForTest() > with ) {
                    ok[ 0 ] = fail( "the text must be measured under the cap it is drawn under: measured "
                            + tp.lengthOfLongestTextOnlyForTest() + ", cap " + with );
                }
                final int avail = tp.circularAvailableRadiusForTest();
                if ( ( tp.circularRadiusForTest() < Math.floor( avail * 0.40 ) - 1 ) || ( tp.circularRadiusForTest() > ring_without ) ) {
                    ok[ 0 ] = fail( "the ring keeps the tree's share under the track: " + tp.circularRadiusForTest() + " of " + avail
                            + " (" + ring_without + " without the track)" );
                }
            } );
        }
        finally {
            dispose( mf );
        }
    }

    // ---- fixture ----

    /** A star of {@code tips}; with {@code long_names} every tip carries a scientific name far wider than any reach,
     *  else short distinct names; with {@code domains} every tip carries a 500-residue architecture. */
    private static Phylogeny star( final int tips, final boolean long_names, final boolean domains ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < tips; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setDistanceToParent( 1.0 );
            tip.setName( "t" + i );
            if ( long_names ) {
                final Taxonomy t = new Taxonomy();
                t.setScientificName( "Pseudomonadendrogrammatophyllum megalospeciosum subspecies extralongensis variety "
                        + "sesquipedalianum forma nomenclaturalis absurdissima collectio typographica " + i );
                tip.getNodeData().addTaxonomy( t );
            }
            if ( domains ) {
                final List<PhylogenyData> ds = new ArrayList<PhylogenyData>();
                ds.add( new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 ) );
                ds.add( new ProteinDomain( "Pkinase", 185, 445, "PF00069", 1e-30 ) );
                final Sequence seq = new Sequence();
                seq.setName( tip.getName() );
                seq.setDomainArchitecture( new DomainArchitecture( ds, 500 ) );
                tip.getNodeData().addSequence( seq );
            }
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        phy.setName( "reach" );
        return phy;
    }

    private static void paint( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.fitRadialTo( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        tp.paintPhylogeny( g, false, false, W, H, 0, 0 );
        g.dispose();
    }

    private static void dispose( final MainFrame[] mf ) throws Exception {
        if ( mf[ 0 ] != null ) {
            SwingUtilities.invokeAndWait( () -> mf[ 0 ].dispose() );
        }
    }

    private static boolean fail( final String what ) {
        System.out.println( "  [RadialLabelReachTest] " + what );
        return false;
    }

    private RadialLabelReachTest() {
    }
}

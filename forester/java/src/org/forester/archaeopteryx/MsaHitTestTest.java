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
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;

import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Sequence;

/**
 * Pins the alignment HIT-TEST to what is actually PAINTED.
 * <p>
 * {@code TreePanel.msaCellAt} re-derives the alignment grid (origin, column width, window offset, per-tip row band)
 * that {@code paintMsaTrack} uses to draw it. Two copies of one geometry drift, and a hit-test that drifts is a
 * readout that confidently reports the wrong residue. So this test does not re-implement the geometry a third time:
 * it renders the panel, walks pixels, and asserts that the residue {@code msaCellAt} names for a point is the
 * residue whose COLOUR was painted at that very point.
 */
public final class MsaHitTestTest {

    private static final String[] ALIGN = { "MKQLEDPFGH-WYVAST", "MKQIEDPFGY-WYVAST", "LRQMEDANGH-WFVCST",
                                            "MK-LEDPFGH-WYVAGT" };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "MsaHitTest: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { alignedTree() }, conf, "msahit" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
                tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                tp.getOptions().setShowOverview( false );
                tp.setOvOn( false );
                tp.getOptions().setShowMsa( true );
                // BELOW the letter threshold on purpose: no glyph is drawn over the cell, so every pixel of a
                // residue cell is its class colour and can be compared with what the hit-test claims is there.
                tp.getOptions().setMsaColumnWidth( 5 );
                mf[ 0 ].showWhole();

                final int w = 1000, h = 500;
                tp.setSize( w, h );
                tp.validate();
                tp.doLayout();
                final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
                final Graphics2D g = img.createGraphics();
                g.setColor( Color.WHITE );
                g.fillRect( 0, 0, w, h );
                tp.printAll( g );
                g.dispose();

                // Per CELL, tally the pixels the hit-test attributes to it: how many carry that residue's own
                // class colour versus SOME OTHER residue's colour. Other things (the scale bar, the tree name) can
                // legitimately be drawn over the band, so a single stray pixel proves nothing -- but a hit-test that
                // is off by a row or a column would systematically attribute a cell's pixels to its neighbour's
                // colour, which this cannot miss.
                final java.util.Map<String, int[]> tally = new java.util.HashMap<>(); // key -> {right, wrong}
                for( int x = 0; x < w; ++x ) {
                    for( int y = 0; y < h; ++y ) {
                        final TreePanel.MsaCell cell = tp.msaCellAt( x, y );
                        if ( cell == null ) {
                            continue;
                        }
                        final String mol = cell._tip.getNodeData().getSequence().getMolecularSequence();
                        final Color expected = MsaColors.colorFor( mol.charAt( cell._column ), false );
                        if ( expected == null ) {
                            continue; // a gap cell has no fill to compare
                        }
                        final Color painted = new Color( img.getRGB( x, y ) );
                        final int[] t = tally.computeIfAbsent( cell._tip.getName() + "#" + cell._column,
                                                               k -> new int[ 2 ] );
                        if ( painted.equals( expected ) ) {
                            ++t[ 0 ];
                        }
                        else if ( colourAppearsIn( painted, mol ) ) {
                            ++t[ 1 ]; // another residue's colour -> the hit-test is pointing at the wrong cell
                        }
                    }
                }
                int verified = 0;
                for( final java.util.Map.Entry<String, int[]> e : tally.entrySet() ) {
                    final int right = e.getValue()[ 0 ], wrong = e.getValue()[ 1 ];
                    if ( ( right + wrong ) < 4 ) {
                        continue; // too few pixels to judge (a sliver at the window edge)
                    }
                    // ZERO tolerance, not a majority: the paint rounds each cell boundary, so a hit-test that
                    // inverts the geometry by plain division is off by one PIXEL COLUMN per cell -- which a
                    // majority vote would happily hide while the tooltip named the neighbouring residue.
                    if ( wrong > 0 ) {
                        ok[ 0 ] = false;
                        System.out.println( "  [MsaHitTestTest] cell " + e.getKey() + ": the hit-test attributes "
                                + wrong + " pixel(s) of ANOTHER residue's colour to it (" + right
                                + " of its own) -- the hit-test and the paint disagree" );
                        return;
                    }
                    ++verified;
                }
                // the tooltip WIRING: over a residue it reports that residue, and off the alignment there is none
                for( final java.util.Map.Entry<String, int[]> e : tally.entrySet() ) {
                    if ( e.getValue()[ 0 ] < 4 ) {
                        continue;
                    }
                    final String[] parts = e.getKey().split( "#" );
                    for( int x = 0; x < w; ++x ) {
                        for( int y = 0; y < h; ++y ) {
                            final TreePanel.MsaCell c = tp.msaCellAt( x, y );
                            if ( ( c == null ) || !c._tip.getName().equals( parts[ 0 ] )
                                    || ( c._column != Integer.parseInt( parts[ 1 ] ) ) ) {
                                continue;
                            }
                            final String tip_text = tp.getToolTipText( new java.awt.event.MouseEvent( tp, 0, 0, 0, x,
                                    y, 0, false ) );
                            if ( ( tip_text == null )
                                    || !tip_text.contains( "Alignment column " + ( c._column + 1 ) )
                                    || !tip_text.contains( parts[ 0 ] ) ) {
                                ok[ 0 ] = false;
                                System.out.println( "  [MsaHitTestTest] the tooltip over " + e.getKey()
                                        + " is wrong: " + tip_text );
                            }
                            x = w; // one sample is enough
                            break;
                        }
                    }
                    break;
                }
                if ( tp.getToolTipText( new java.awt.event.MouseEvent( tp, 0, 0, 0, 0, 0, 0, false ) ) != null ) {
                    ok[ 0 ] = false;
                    System.out.println( "  [MsaHitTestTest] there must be no tooltip away from the alignment" );
                }
                // The overview box FLOATS over the viewport and, in either RIGHT placement, can sit exactly where
                // the alignment window is drawn. A readout for the residue underneath would describe a cell
                // the user cannot see, so the hit-test must decline there.
                // Widen the cells so the alignment band actually reaches the viewport's right edge, where the
                // overview sits -- at the narrow width used for the colour comparison above the two do not overlap
                // and this check would prove nothing.
                tp.getOptions().setMsaColumnWidth( 24 );
                // Long tip labels push the alignment band to the right, where the overview floats -- the short
                // names used above leave the band ending well clear of it, and the check would be vacuous.
                for( final PhylogenyNode t : tp.getPhylogeny().getExternalNodes() ) {
                    t.setName( "a_rather_long_sequence_identifier_" + t.getName() );
                }
                mf[ 0 ].showWhole();
                tp.setSize( w, h );
                tp.validate();
                tp.doLayout();
                final java.awt.image.BufferedImage img2 =
                        new java.awt.image.BufferedImage( w, h, java.awt.image.BufferedImage.TYPE_INT_RGB );
                final Graphics2D g2 = img2.createGraphics();
                tp.printAll( g2 );
                g2.dispose();
                final java.util.List<int[]> before = new java.util.ArrayList<>();
                for( int x = 0; x < w; x += 3 ) {
                    for( int y = 0; y < h; y += 3 ) {
                        if ( tp.msaCellAt( x, y ) != null ) {
                            before.add( new int[] { x, y } );
                        }
                    }
                }
                tp.getOptions().setOvPlacement( Options.OVERVIEW_PLACEMENT_TYPE.LOWER_RIGHT );
                tp.getOptions().setShowOverview( true );
                tp.setOvOn( true );
                tp.updateOvSizes();
                tp.updateOvSettings(); // this is what computes the box's position from the placement
                final int[] ov = tp.floatingOverlayRectForTest();
                if ( ov == null ) {
                    ok[ 0 ] = false;
                    System.out.println( "  [MsaHitTestTest] precondition: the overview should be drawn" );
                    return;
                }
                boolean overlaps = false;
                for( final int[] pt : before ) {
                    if ( ( pt[ 0 ] >= ov[ 0 ] ) && ( pt[ 0 ] < ( ov[ 0 ] + ov[ 2 ] ) ) && ( pt[ 1 ] >= ov[ 1 ] )
                            && ( pt[ 1 ] < ( ov[ 1 ] + ov[ 3 ] ) ) ) {
                        overlaps = true;
                        break;
                    }
                }
                if ( !overlaps ) {
                    ok[ 0 ] = TestFail.here();
                    int minx = Integer.MAX_VALUE, maxx = 0, miny = Integer.MAX_VALUE, maxy = 0;
                    for( final int[] pt : before ) {
                        minx = Math.min( minx, pt[ 0 ] );
                        maxx = Math.max( maxx, pt[ 0 ] );
                        miny = Math.min( miny, pt[ 1 ] );
                        maxy = Math.max( maxy, pt[ 1 ] );
                    }
                    System.out.println( "  [MsaHitTestTest] precondition: the overview box (" + ov[ 0 ] + "," + ov[ 1 ]
                            + " " + ov[ 2 ] + "x" + ov[ 3 ] + ") must cover part of the alignment (which spans x "
                            + minx + ".." + maxx + ", y " + miny + ".." + maxy + ") for this check to mean anything" );
                    return;
                }
                int masked = 0;
                for( final int[] pt : before ) {
                    if ( tp.msaCellAt( pt[ 0 ], pt[ 1 ] ) == null ) {
                        ++masked;
                    }
                }
                if ( masked == 0 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  [MsaHitTestTest] with the overview drawn over the alignment, the cells "
                            + "beneath it must stop reporting -- none did" );
                }
                tp.getOptions().setShowOverview( false );
                tp.setOvOn( false );
                // ...and with it off again, every one of those points reports its cell once more
                for( final int[] pt : before ) {
                    if ( tp.msaCellAt( pt[ 0 ], pt[ 1 ] ) == null ) {
                        ok[ 0 ] = false;
                        System.out.println( "  [MsaHitTestTest] hiding the overview must restore the readout at ("
                                + pt[ 0 ] + "," + pt[ 1 ] + ")" );
                        break;
                    }
                }
                if ( verified < 20 ) {
                    ok[ 0 ] = false;
                    System.out.println( "  [MsaHitTestTest] only " + verified + " cells could be verified -- the "
                            + "alignment does not appear to have been drawn, so this test proved nothing" );
                }
            } );
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Whether {@code c} is the class colour of ANY residue of this row (tolerates a cell-boundary pixel). */
    private static boolean colourAppearsIn( final Color c, final String mol ) {
        for( int i = 0; i < mol.length(); ++i ) {
            if ( c.equals( MsaColors.colorFor( mol.charAt( i ), false ) ) ) {
                return true;
            }
        }
        return false;
    }

    private static Phylogeny alignedTree() {
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < ALIGN.length; ++i ) {
            root.addAsChild( alignedTip( "seq_" + i, ALIGN[ i ] ) );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode alignedTip( final String name, final String seq ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final Sequence s = new Sequence();
        s.setMolecularSequence( seq );
        s.setMolecularSequenceAligned( true );
        n.getNodeData().addSequence( s );
        return n;
    }

    private MsaHitTestTest() {
    }
}

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
import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.LayoutWidthBudget.Part;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * Tests the horizontal allocator that keeps the tree from being squeezed to a line. Pure arithmetic, so it runs
 * everywhere -- no display needed.
 */
public final class LayoutWidthBudgetTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LayoutWidthBudget: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        final boolean[] ok = { true };

        // ---- (1) the real regression: the measured wnt_deep.xml numbers must stop collapsing the tree ------------
        // Requests as the panel actually computed them at 1400 px: labels+domains 525, alignment 877, legend 249.
        // Before the allocator these summed to 1651 and the tree got -271 px (a line).
        final LayoutWidthBudget deep = new LayoutWidthBudget.Builder()
                // Text + domain track as ONE non-elastic part, the way computeWidthBudget requests them (the
                // domain track rides the label reservation). 447 px is what labelWidthCap now holds them to at
                // this width once the legend's own column has been taken off the top -- the measured 525 was the
                // pre-fix figure, when the labels were free to grow into the legend's space.
                .request( Part.LABELS, 447, 447 )
                .request( Part.MSA, 877, 120 )
                .request( Part.LEGEND, 249, 249 )
                .allocate( 1400, 40, 0.40 );
        final int usable = 1400 - 40; // DEPTH_AXIS_FIXED_MARGIN: the depth scale never gets 2 * MOVE
        final int floor40 = (int) Math.round( usable * 0.40 );
        if ( deep.treeWidth() < floor40 ) {
            ok[ 0 ] = fail( "the tree must get its 40% share (" + floor40 + " px), got " + deep.treeWidth() );
        }
        if ( deep.treeWidth() <= 0 ) {
            ok[ 0 ] = fail( "THE regression: the tree must never be squeezed to nothing, got " + deep.treeWidth() );
        }
        if ( ( deep.treeWidth() + deep.sideWidth() ) > usable ) {
            ok[ 0 ] = fail( "the allocator overspent: tree " + deep.treeWidth() + " + sides " + deep.sideWidth()
                    + " > " + usable );
        }
        // every squeezed part keeps at least its minimum, and none is granted more than it asked for
        if ( ( deep.granted( Part.MSA ) < 120 ) || ( deep.granted( Part.MSA ) > 877 ) ) {
            ok[ 0 ] = fail( "the alignment must land between its minimum and its request, got "
                    + deep.granted( Part.MSA ) );
        }
        if ( deep.granted( Part.MSA ) >= 877 ) {
            ok[ 0 ] = fail( "the alignment must actually BE squeezed here -- otherwise nothing was tested" );
        }

        // ---- (2) when everything fits, nobody is squeezed and the tree keeps the WHOLE remainder ----------------
        // A tree with a little side data must not have width confiscated just because its share is 40%.
        final LayoutWidthBudget roomy = new LayoutWidthBudget.Builder().request( Part.LABELS, 100, 40 )
                .request( Part.MSA, 200, 120 ).allocate( 1400, 40, 0.40 );
        if ( ( roomy.granted( Part.LABELS ) != 100 ) || ( roomy.granted( Part.MSA ) != 200 ) ) {
            ok[ 0 ] = fail( "requests that fit must be granted in full, got " + roomy.granted( Part.LABELS ) + "/"
                    + roomy.granted( Part.MSA ) );
        }
        if ( roomy.treeWidth() != ( usable - 300 ) ) {
            ok[ 0 ] = fail( "the tree must keep the whole remainder (" + ( usable - 300 ) + "), got "
                    + roomy.treeWidth() );
        }

        // ---- (3) the squeeze is PROPORTIONAL above the minimums -------------------------------------------------
        // Two parts, same minimum, one asking for twice as much: the bigger asker keeps the bigger grant, and the
        // gap between them narrows rather than one being starved.
        final LayoutWidthBudget prop = new LayoutWidthBudget.Builder().request( Part.MSA, 800, 100 )
                .request( Part.ANNOTATIONS, 400, 100 ).allocate( 1000, 0, 0.40 );
        if ( prop.granted( Part.MSA ) <= prop.granted( Part.ANNOTATIONS ) ) {
            ok[ 0 ] = fail( "the part that asked for more must still get more, got " + prop.granted( Part.MSA )
                    + " vs " + prop.granted( Part.ANNOTATIONS ) );
        }
        if ( ( prop.granted( Part.MSA ) < 100 ) || ( prop.granted( Part.ANNOTATIONS ) < 100 ) ) {
            ok[ 0 ] = fail( "a squeezed part must never drop below its minimum" );
        }
        // 600 px of side budget, 200 of it the minimums -> 400 headroom split 700:300 = 280:120
        if ( ( prop.granted( Part.MSA ) != 380 ) || ( prop.granted( Part.ANNOTATIONS ) != 220 ) ) {
            ok[ 0 ] = fail( "expected a 700:300 split of the headroom (380/220), got " + prop.granted( Part.MSA )
                    + "/" + prop.granted( Part.ANNOTATIONS ) );
        }

        // ---- (4) monotone: asking for more never gets you less --------------------------------------------------
        int previous = -1;
        for( int want = 200; want <= 2000; want += 50 ) {
            final int g = new LayoutWidthBudget.Builder().request( Part.MSA, want, 120 )
                    .request( Part.LABELS, 400, 60 ).allocate( 1400, 40, 0.40 ).granted( Part.MSA );
            if ( g < previous ) {
                ok[ 0 ] = fail( "asking for " + want + " granted " + g + ", less than the previous step's "
                        + previous + " -- the allocation must be monotone" );
                break;
            }
            previous = g;
        }

        // ---- (5) the share is honoured, and clamped to the shipped range ----------------------------------------
        for( final double share : new double[] { 0.25, 0.40, 0.60, 0.80 } ) {
            final LayoutWidthBudget b = new LayoutWidthBudget.Builder().request( Part.MSA, 5000, 120 )
                    .allocate( 1000, 0, share );
            final int want_floor = (int) Math.round( 1000 * share );
            if ( b.treeWidth() < want_floor ) {
                ok[ 0 ] = fail( "at share " + share + " the tree must get >= " + want_floor + ", got "
                        + b.treeWidth() );
            }
        }
        // out-of-range shares are clamped, not obeyed: a 0.0 share must NOT hand the tree nothing
        final LayoutWidthBudget clamped_low = new LayoutWidthBudget.Builder().request( Part.MSA, 5000, 120 )
                .allocate( 1000, 0, 0.0 );
        if ( clamped_low.treeWidth() < (int) Math.round( 1000 * AptxConstants.TREE_WIDTH_SHARE_MIN ) ) {
            ok[ 0 ] = fail( "a share below the minimum must clamp UP to it, got " + clamped_low.treeWidth() );
        }
        final LayoutWidthBudget clamped_high = new LayoutWidthBudget.Builder().request( Part.MSA, 5000, 120 )
                .allocate( 1000, 0, 5.0 );
        if ( clamped_high.granted( Part.MSA ) < 120 ) {
            ok[ 0 ] = fail( "a share above the maximum must clamp DOWN, leaving the alignment its minimum" );
        }

        // ---- (6) degenerate inputs must not throw, and must not hand out negative widths -------------------------
        // A window narrower than the minimums combined: the minimums win and the tree takes what is left (the one
        // case the share cannot be honoured). Nothing may go negative.
        final LayoutWidthBudget tiny = new LayoutWidthBudget.Builder().request( Part.LABELS, 300, 60 )
                .request( Part.MSA, 800, 120 ).allocate( 150, 40, 0.40 );
        if ( ( tiny.treeWidth() < 0 ) || ( tiny.granted( Part.MSA ) < 0 ) || ( tiny.granted( Part.LABELS ) < 0 ) ) {
            ok[ 0 ] = fail( "a window narrower than the minimums must not produce negative widths" );
        }
        if ( tiny.granted( Part.MSA ) != 120 ) {
            ok[ 0 ] = fail( "in the degenerate case each part keeps exactly its minimum, got "
                    + tiny.granted( Part.MSA ) );
        }
        final LayoutWidthBudget empty = new LayoutWidthBudget.Builder().allocate( 800, 40, 0.40 );
        if ( ( empty.treeWidth() != 760 ) || ( empty.sideWidth() != 0 ) ) {
            ok[ 0 ] = fail( "with no side data the tree gets everything but the fixed margin, got "
                    + empty.treeWidth() );
        }
        final LayoutWidthBudget zero = new LayoutWidthBudget.Builder().request( Part.MSA, 400, 120 )
                .allocate( 0, 40, 0.40 );
        if ( ( zero.treeWidth() != 0 ) || ( zero.granted( Part.MSA ) < 0 ) ) {
            ok[ 0 ] = fail( "a zero-width panel must yield zero, not a negative width" );
        }
        // a part that wants nothing is granted nothing, however large its stated minimum (an OFF track)
        final LayoutWidthBudget off = new LayoutWidthBudget.Builder().request( Part.MSA, 0, 120 )
                .request( Part.LABELS, 300, 60 ).allocate( 1000, 40, 0.40 );
        if ( off.granted( Part.MSA ) != 0 ) {
            ok[ 0 ] = fail( "a track that is switched off must be granted 0, got " + off.granted( Part.MSA ) );
        }
        // an unrequested part reads back as 0 rather than throwing
        if ( off.granted( Part.CLADE_BANDS ) != 0 ) {
            ok[ 0 ] = fail( "an unrequested part must read back as 0" );
        }

        // ---- (7) the allocator never overspends, over any mix ---------------------------------------------------
        // The invariant everything else rests on: tree + sides <= usable, for a wide sweep of requests.
        for( int total = 200; total <= 3000; total += 137 ) {
            for( int msa = 0; msa <= 3000; msa += 311 ) {
                for( int labels = 0; labels <= 1200; labels += 251 ) {
                    final LayoutWidthBudget b = new LayoutWidthBudget.Builder().request( Part.MSA, msa, 120 )
                            .request( Part.LABELS, labels, 60 ).request( Part.LEGEND, 249, 249 )
                            .allocate( total, 40, 0.40 );
                    final int use = Math.max( 0, total - 40 );
                    // Only where the minimums FIT. Where they do not, the allocator deliberately honours them
                    // anyway and lets the tree take what is left -- a track drawn under its minimum is unreadable,
                    // so there is nothing to gain by shaving it -- and the sum may then exceed the usable width.
                    // The Builder clamps each min into [0, want], so a part wanting nothing has a minimum of 0.
                    final boolean minimums_fit = ( Math.min( 120, msa ) + Math.min( 60, labels ) + 249 ) <= use;
                    if ( minimums_fit && ( ( b.treeWidth() + b.sideWidth() ) > use ) ) {
                        ok[ 0 ] = fail( "overspent at total=" + total + " msa=" + msa + " labels=" + labels + ": "
                                + b.treeWidth() + "+" + b.sideWidth() + " > " + use );
                        return ok[ 0 ];
                    }
                    if ( ( b.treeWidth() < 0 ) || ( b.granted( Part.MSA ) > Math.max( 0, msa ) ) ) {
                        ok[ 0 ] = fail( "bad grant at total=" + total + " msa=" + msa + " labels=" + labels );
                        return ok[ 0 ];
                    }
                }
            }
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return ok[ 0 ]; // the rest needs a real panel
        }
        try {
            onARealPanel( ok );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    /**
     * The regression itself, on a real panel: a deep tree carrying long labels, domain architectures AND an
     * alignment -- the combination that used to leave the tree a negative width and draw it as a single line.
     */
    private static void onARealPanel( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { crowdedTree() }, new Configuration(), "budget" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
        final int w = 1400, h = 900;
        final int usable = w - 40; // TreePanel.DEPTH_AXIS_FIXED_MARGIN (2 * MOVE), what the depth scale never gets

        SwingUtilities.invokeAndWait( () -> {
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
            cp.setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, true );
            cp.setCheckbox( DisplayOption.SHOW_MSA, true );
            // Keep the names at full length: "Shorten Labels" would cut them to 18 characters, and then the label
            // cap -- the thing that stops the labels eating the tree's share -- would never bind and a test of it
            // would prove nothing.
            cp.setCheckbox( DisplayOption.SHORTEN_LABELS, false );
            tp.setSize( w, h );
            tp.calcParametersForPainting( w, h );
        } );
        // (A) every track is really on -- otherwise the crowding this test is about never happened
        if ( !tp.isShowMsa() || !tp.getControlPanel().isShowDomainArchitectures() ) {
            fail( "precondition: both the alignment and the domains must be on for this to test anything" );
            ok[ 0 ] = false;
        }
        if ( tp.grantedWidthForTest( Part.MSA ) <= 0 ) {
            fail( "precondition: the alignment must be asking for width" );
            ok[ 0 ] = false;
        }
        // (B) THE regression: the tree keeps its share instead of being the leftover
        // The share is a TARGET, not an absolute floor (see LayoutWidthBudget): the parts that cannot shrink --
        // the labels once the font auto-fit has bottomed out, the clade bands, and the legend column, which keeps
        // the legend off the tracks -- are honoured first, so a crowded figure can land under it. What must never
        // happen again is the tree being left the NEGATIVE remainder, so the contract asserted here is the
        // minimum share, which is what the slider will not go below.
        final int floor = (int) Math.round( usable * AptxConstants.TREE_WIDTH_SHARE_MIN );
        if ( tp.grantedTreeWidthForTest() < floor ) {
            fail( "even crowded, the tree must keep the minimum share ("
                    + Math.round( AptxConstants.TREE_WIDTH_SHARE_MIN * 100 ) + "%, " + floor + " px) -- got "
                    + tp.grantedTreeWidthForTest() );
            ok[ 0 ] = false;
        }
        // ...and it is actually DRAWN: a zero depth scale is precisely what "there is only a line!" looked like
        if ( tp.getXcorrectionFactor() <= 0 ) {
            fail( "the depth scale collapsed to zero -- the tree would draw as a single line" );
            ok[ 0 ] = false;
        }
        // (B2) the layout must OBEY the allocation: what the reserves actually take off the depth axis cannot
        // exceed what they were granted. Without this the reserve methods could ignore their grants entirely --
        // every number above would still be right, and the tree would still be squeezed on screen.
        if ( tp.actualSideWidthForTest() > ( usable - tp.grantedTreeWidthForTest() ) ) {
            fail( "the reserves take " + tp.actualSideWidthForTest() + " px but were granted only "
                    + ( usable - tp.grantedTreeWidthForTest() ) + " -- the layout is ignoring the budget" );
            ok[ 0 ] = false;
        }

        // (C) the slider moves the split, and the alignment is what gives way
        final int[] tree_at = new int[ 2 ];
        final int[] msa_at = new int[ 2 ];
        final float[] corr_at = new float[ 2 ];
        final double[] shares = { AptxConstants.TREE_WIDTH_SHARE_MIN, 0.70 };
        for( int i = 0; i < shares.length; ++i ) {
            final int idx = i;
            SwingUtilities.invokeAndWait( () -> {
                tp.getOptions().setTreeWidthShare( shares[ idx ] );
                tp.calcParametersForPainting( w, h );
                tree_at[ idx ] = tp.grantedTreeWidthForTest();
                msa_at[ idx ] = tp.grantedWidthForTest( Part.MSA );
                corr_at[ idx ] = tp.getXcorrectionFactor();
            } );
        }
        if ( tree_at[ 1 ] <= tree_at[ 0 ] ) {
            fail( "raising the tree share must widen the tree: " + tree_at[ 0 ] + " -> " + tree_at[ 1 ] );
            ok[ 0 ] = false;
        }
        if ( msa_at[ 1 ] >= msa_at[ 0 ] ) {
            fail( "...and the alignment must be what gives way: " + msa_at[ 0 ] + " -> " + msa_at[ 1 ] );
            ok[ 0 ] = false;
        }
        // The grant has to reach the DRAWING, not just the arithmetic: the depth scale is what the tree is
        // actually drawn at, so it must widen too. Without this the reserve methods could ignore their grants
        // entirely and the numbers above would still look right.
        if ( corr_at[ 1 ] <= corr_at[ 0 ] ) {
            fail( "the drawn depth scale must follow the share (the reserves must USE their grants): "
                    + corr_at[ 0 ] + " -> " + corr_at[ 1 ] );
            ok[ 0 ] = false;
        }
        // (D) the control-panel slider is wired to the same setting, both ways
        SwingUtilities.invokeAndWait( () -> cp.setTreeSharePercentForTest( 55 ) );
        if ( Math.abs( tp.getOptions().getTreeWidthShare() - 0.55 ) > 1e-9 ) {
            fail( "dragging the Tree share slider must set the share, got " + tp.getOptions().getTreeWidthShare() );
            ok[ 0 ] = false;
        }
        // (E) Reset to Defaults puts the share back, and the slider shows it
        SwingUtilities.invokeAndWait( () -> mf[ 0 ].resetToDefaults() );
        if ( Math.abs( tp.getOptions().getTreeWidthShare() - AptxConstants.TREE_WIDTH_SHARE_DEFAULT ) > 1e-9 ) {
            fail( "Reset to Defaults must restore the default share, got " + tp.getOptions().getTreeWidthShare() );
            ok[ 0 ] = false;
        }
        if ( cp.treeSharePercentForTest() != (int) Math.round( AptxConstants.TREE_WIDTH_SHARE_DEFAULT * 100 ) ) {
            fail( "...and the slider must show it, not a stale value: " + cp.treeSharePercentForTest() );
            ok[ 0 ] = false;
        }
        // (G) the labels must not be shrunk past readability. On a window too narrow to hold everything the
        // labels OVERFLOW their cap -- the budget then squeezes the elastic tracks and, if it must, the tree --
        // rather than being drawn at an unreadable 3 px, which is where a budget-derived cap otherwise leads.
        SwingUtilities.invokeAndWait( () -> {
            tp.getOptions().setTreeWidthShare( AptxConstants.TREE_WIDTH_SHARE_DEFAULT );
            tp.setSize( 700, h );
            tp.calcParametersForPainting( 700, h );
        } );
        // Against the SHIPPED constant, deliberately -- not against autofitMinFontSize(), which is the code's own
        // answer: comparing the font to the floor the code just computed passes for every floor, including none.
        if ( tp.labelFontSizeForTest() < AptxConstants.LABEL_AUTOFIT_MIN_FONT_SIZE ) {
            fail( "the auto-fit shrank the labels to " + tp.labelFontSizeForTest() + " px, past the floor of "
                    + AptxConstants.LABEL_AUTOFIT_MIN_FONT_SIZE );
            ok[ 0 ] = false;
        }
        if ( AptxConstants.LABEL_AUTOFIT_MIN_FONT_SIZE < 6 ) {
            fail( "the shipped label floor has been lowered to " + AptxConstants.LABEL_AUTOFIT_MIN_FONT_SIZE
                    + " px, which is below readable" );
            ok[ 0 ] = false;
        }
        // ...and the tree is still a tree, not a line, even on that narrow window
        if ( tp.getXcorrectionFactor() <= 0 ) {
            fail( "even on a narrow window the tree must keep a non-zero depth scale" );
            ok[ 0 ] = false;
        }

        // (H) "Tree share" must move the tree in EVERY display type. Each layout trades a different quantity --
        // the rectangular family divides the depth WIDTH, the circular one the tip-ring RADIUS, the unrooted one
        // its fan spread -- so each is measured on its own, not on a single number that happens to exist
        // everywhere. It was inert in four of these until the budget was extended to them.
        for( final Options.TREE_ORIENTATION o : new Options.TREE_ORIENTATION[] {
                Options.TREE_ORIENTATION.ROOT_LEFT, Options.TREE_ORIENTATION.ROOT_TOP,
                Options.TREE_ORIENTATION.ROOT_BOTTOM } ) {
            final double[] got = new double[ 2 ];
            for( int i = 0; i < 2; ++i ) {
                final double share = ( i == 0 ) ? 0.25 : 0.80;
                final int idx = i;
                SwingUtilities.invokeAndWait( () -> {
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    tp.setTreeOrientation( o );
                    tp.getOptions().setTreeWidthShare( share );
                    // a square, roomy canvas: in a vertical orientation the DEPTH axis is the height, and at
                    // 900 px this fixture's labels alone exceed every share, so nothing could move
                    tp.setSize( 1600, 1600 );
                    tp.calcParametersForPainting( 1600, 1600 );
                    got[ idx ] = tp.grantedTreeWidthForTest();
                } );
            }
            if ( got[ 1 ] <= got[ 0 ] ) {
                fail( "raising the share must widen the tree in " + o + ": " + got[ 0 ] + " -> " + got[ 1 ] );
                ok[ 0 ] = false;
            }
        }
        for( final Options.PHYLOGENY_GRAPHICS_TYPE t : new Options.PHYLOGENY_GRAPHICS_TYPE[] {
                Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
            final boolean circular = ( t == Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            final double[] got = new double[ 2 ];
            for( int i = 0; i < 2; ++i ) {
                final double share = ( i == 0 ) ? 0.25 : 0.80;
                final int idx = i;
                SwingUtilities.invokeAndWait( () -> {
                    tp.setPhylogenyGraphicsType( t );
                    tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
                    tp.getOptions().setTreeWidthShare( share );
                    tp.setSize( 1600, 1600 );
                    tp.calcParametersForPainting( 1600, 1600 );
                    final java.awt.image.BufferedImage img = new java.awt.image.BufferedImage( 8, 8,
                            java.awt.image.BufferedImage.TYPE_INT_RGB );
                    final java.awt.Graphics2D g = img.createGraphics();
                    tp.paintPhylogeny( g, false, false, 1600, 1600, 0, 0 );
                    g.dispose();
                    got[ idx ] = circular ? tp.circularRadiusForTest() : tp.urtFactorForTest();
                } );
            }
            if ( got[ 1 ] <= got[ 0 ] ) {
                fail( "raising the share must grow the tree in " + t + " ("
                        + ( circular ? "tip-ring radius" : "fan spread" ) + "): " + got[ 0 ] + " -> " + got[ 1 ] );
                ok[ 0 ] = false;
            }
        }
        SwingUtilities.invokeAndWait( () -> {
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
            tp.setTreeOrientation( Options.TREE_ORIENTATION.ROOT_LEFT );
            tp.getOptions().setTreeWidthShare( AptxConstants.TREE_WIDTH_SHARE_DEFAULT );
        } );

        // (F) a radial layout has no depth column to divide -- it must not be given a budget at all
        SwingUtilities.invokeAndWait( () -> {
            tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            tp.calcParametersForPainting( w, h );
        } );
        if ( tp.grantedTreeWidthForTest() != -1 ) {
            fail( "a circular layout has no width column to divide, so it must carry no budget" );
            ok[ 0 ] = false;
        }
        SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
    }

    /** The shape that broke: many tips, long names, a domain architecture and an aligned sequence on every one. */
    private static Phylogeny crowdedTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 24; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            // deliberately far longer than any sane cap: under the old flat "labels may take 70% of the panel"
            // rule these alone would starve the tree, which is what the budget has to prevent
            tip.setName( "Homo_sapiens_WNT_family_member_" + i
                    + "_isoform_X1_preproprotein_predicted_LOC102724788_transcript_variant_7_partial_cds" );
            tip.setDistanceToParent( 0.05 + ( i * 0.01 ) );
            final Sequence seq = new Sequence();
            final List<PhylogenyData> domains = new ArrayList<PhylogenyData>();
            domains.add( new ProteinDomain( "PF0000" + ( i % 9 ), 5, 180 ) );
            seq.setDomainArchitecture( new DomainArchitecture( domains, 400 ) );
            final StringBuilder sb = new StringBuilder();
            for( int c = 0; c < 300; ++c ) {
                sb.append( "ACDEFGHIKLMNPQRSTVWY-".charAt( ( c + i ) % 21 ) );
            }
            seq.setMolecularSequence( sb.toString() );
            seq.setMolecularSequenceAligned( true );
            tip.getNodeData().addSequence( seq );
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean fail( final String message ) {
        System.out.println( "  [LayoutWidthBudgetTest] " + message );
        return false;
    }

    private LayoutWidthBudgetTest() {
        // not instantiable
    }
}

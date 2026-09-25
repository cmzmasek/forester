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
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingUtilities;

import org.forester.io.parsers.util.ParserUtils;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;
import org.forester.phylogeny.data.Taxonomy;

/**
 * "Tree share" as a control the user can actually feel.
 * <p>
 * The setting is a MINIMUM share of the width the tree is guaranteed, and the tip labels give that width up by
 * shrinking their font. Once the font reaches its readability floor the labels can give up nothing more, and the
 * share was then quietly ignored: measured on a 360-tip tree carrying domains and an alignment, the tree stalled
 * at 52% of the width whether the user asked for 60%, 70% or 80%. Past that floor the labels are now SHORTENED
 * with an ellipsis instead.
 * <p>
 * Two things have to hold, and the second is the one that bites:
 * <ol>
 * <li>raising the share past the font floor really does widen the tree, and lowering it gives the labels back;</li>
 * <li>a shortened label is DRAWN no wider than it is MEASURED. The label is painted as two segments -- taxonomy,
 * then node data -- so capping only the second one left a tip whose taxonomy alone was already over the cap being
 * drawn at full width on top of the track beside it. That is why this test reads the drawn width off the rendered
 * PIXELS rather than asking the code how wide it thinks the label is: the two disagreeing is the whole defect.</li>
 * </ol>
 * And one deliberate non-behaviour: a tree whose labels fit is never shortened, at any share.
 */
public final class TreeShareSqueezeTest {

    /** Distinct initials, so the display-time common-prefix shortening never engages on the fixture. */
    private final static String[] NAMES = { "alpha", "bravo", "cobra", "delta", "eagle", "fjord", "gamma", "hydra",
            "ibex", "jumbo", "kappa", "lemur", "mango", "nexus", "opal", "prism", "quilt", "raven", "sigma",
            "tulip" };
    private final static int W          = 1400;
    private final static int H          = 900;
    /** Antialiasing puts a faint pixel or two past the last glyph; the claims here are about tens of pixels. */
    private final static int INK_SLACK  = 4;

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TreeShareSqueeze: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            longLabels( ok );
            shortLabels( ok );
            smallGainRefused( ok );
            fittingLabelsBesideDomains( ok );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    // ---- labels too long to fit even at the smallest font: the share must still move the tree ------------------
    private static void longLabels( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { star( 40, true ) }, new Configuration(), "share-long" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
        // 40%, 80%, then 40% AGAIN: a setting that cannot be taken back is not a setting. It is also the only
        // way to catch a cap that leaks from one layout pass into the next, which would leave the labels short
        // after the user had already dragged the slider back.
        final double[] shares = { 0.40, 0.80, 0.40 };
        final int[] grant = new int[ shares.length ];
        final int[] cap = new int[ shares.length ];
        final int[] drawn = new int[ shares.length ];
        final boolean[] live = new boolean[ shares.length ];
        SwingUtilities.invokeAndWait( () -> {
            cp.setCheckbox( DisplayOption.SHOW_TAXONOMY_SCIENTIFIC_NAMES, true );
            cp.setCheckbox( DisplayOption.SHOW_SEQ_NAMES, true );
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            for( int i = 0; i < shares.length; ++i ) {
                tp.getOptions().setTreeWidthShare( shares[ i ] );
                final BufferedImage img = paint( tp );
                grant[ i ] = tp.grantedTreeWidthForTest();
                cap[ i ] = tp.maxTipLabelTextWidthForTest();
                drawn[ i ] = drawnLabelWidth( tp, img );
                live[ i ] = cp.isTreeShareSliderEnabledForTest();
            }
        } );
        // precondition: the fixture really is one the font auto-fit cannot rescue. Without this a flat result below
        // would read as "the share does nothing" when the truth is "these labels never needed shortening".
        if ( cap[ 1 ] == Integer.MAX_VALUE ) {
            fail( ok, "precondition: at an 80% share these labels must not fit whole -- no cap was applied" );
            dispose( mf );
            return;
        }
        if ( drawn[ 0 ] < 100 ) {
            fail( ok, "precondition: the 40% label must actually be long, measured " + drawn[ 0 ] + " px of ink" );
            dispose( mf );
            return;
        }
        if ( grant[ 1 ] <= grant[ 0 ] ) {
            fail( ok, "raising the share past the font floor must widen the tree: 40% gave " + grant[ 0 ]
                    + " px, 80% gave " + grant[ 1 ] );
        }
        if ( drawn[ 1 ] >= drawn[ 0 ] ) {
            fail( ok, "the labels must give the width up by getting SHORTER on screen: " + drawn[ 0 ] + " px at 40%, "
                    + drawn[ 1 ] + " px at 80%" );
        }
        // THE measure-vs-draw check, read off the pixels: every segment of the label obeys the cap, not just the
        // last one. Before the taxonomy segment was capped this overshot by the full width of a scientific name.
        if ( drawn[ 1 ] > ( cap[ 1 ] + INK_SLACK ) ) {
            fail( ok, "a capped label is drawn wider than it is measured: cap " + cap[ 1 ] + " px, drew " + drawn[ 1 ]
                    + " px -- the overflow lands on whatever is drawn beside the labels" );
        }
        // ...and the reservation the tip-aligned tracks are placed from must cover the drawn label too
        if ( drawn[ 1 ] > ( tp.getLongestExtNodeInfo() + INK_SLACK ) ) {
            fail( ok, "the label reservation (" + tp.getLongestExtNodeInfo() + " px) is smaller than the drawn label ("
                    + drawn[ 1 ] + " px)" );
        }
        if ( !live[ 0 ] || !live[ 1 ] ) {
            fail( ok, "the slider must stay live while it still changes the layout (40%: " + live[ 0 ] + ", 80%: "
                    + live[ 1 ] + ")" );
        }
        // back to 40%: the labels must come back whole, and the width with them
        if ( cap[ 2 ] != cap[ 0 ] ) {
            fail( ok, "dragging the share back must release the label cap: it was " + capText( cap[ 0 ] )
                    + " at 40%, and " + capText( cap[ 2 ] ) + " after a round trip through 80%" );
        }
        if ( ( drawn[ 2 ] != drawn[ 0 ] ) || ( grant[ 2 ] != grant[ 0 ] ) ) {
            fail( ok, "a round trip must restore the layout exactly: label " + drawn[ 0 ] + " -> " + drawn[ 2 ]
                    + " px, tree " + grant[ 0 ] + " -> " + grant[ 2 ] + " px" );
        }
        dispose( mf );
    }

    // ---- labels that fit: never shortened, and the slider says so ---------------------------------------------
    private static void shortLabels( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { star( 12, false ) }, new Configuration(), "share-short" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
        final int[] drawn = new int[ 2 ];
        final boolean[] capped = new boolean[ 2 ];
        final boolean[] live = new boolean[ 2 ];
        SwingUtilities.invokeAndWait( () -> {
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            final double[] shares = { 0.40, 0.80 };
            for( int i = 0; i < shares.length; ++i ) {
                tp.getOptions().setTreeWidthShare( shares[ i ] );
                final BufferedImage img = paint( tp );
                capped[ i ] = tp.maxTipLabelTextWidthForTest() != Integer.MAX_VALUE;
                drawn[ i ] = drawnLabelWidth( tp, img );
                live[ i ] = cp.isTreeShareSliderEnabledForTest();
            }
        } );
        if ( drawn[ 0 ] < 30 ) {
            fail( ok, "precondition: the short-label tree must still draw a label, got " + drawn[ 0 ] + " px" );
            dispose( mf );
            return;
        }
        if ( capped[ 0 ] || capped[ 1 ] ) {
            fail( ok, "a label that FITS must never be shortened -- shortening is what the font floor buys, not "
                    + "something the share does on its own (capped at 40%: " + capped[ 0 ] + ", at 80%: "
                    + capped[ 1 ] + ")" );
        }
        if ( drawn[ 1 ] != drawn[ 0 ] ) {
            fail( ok, "these labels fit at every share, so the drawn label must not change: " + drawn[ 0 ] + " px vs "
                    + drawn[ 1 ] + " px" );
        }
        // The tree already has more than the largest share would guarantee, so every position of the slider draws
        // the identical figure. Saying "no effect" is the honest alternative to a live-looking control that is not.
        if ( live[ 0 ] || live[ 1 ] ) {
            fail( ok, "with nothing competing for the width the slider must be greyed out (40%: " + live[ 0 ]
                    + ", 80%: " + live[ 1 ] + ")" );
        }
        dispose( mf );
    }

    // ---- labels that FIT, beside a wide domain track: still never shortened -----------------------------------
    /**
     * The narrow case that separates "shorten only what does not fit" from "shorten whenever it would help".
     * <p>
     * The label reservation is the text PLUS the domain track, so this is the shape where the two could come
     * apart: a track wide enough to dominate the reservation while the text itself has room to spare. The outcome
     * pinned here is that such labels keep every character.
     * <p>
     * Honest note on what this does NOT prove. Deleting the "do they actually overflow?" condition does not change
     * this result, because the cap the loop computes (label_cap minus the track) still lands outside text this
     * short, so nothing is truncated and the worth-it check discards the cap anyway. That condition is a cost
     * guard rather than a rule -- see the comment on it -- and no test can catch its removal. What this case does
     * guard is the ARITHMETIC: change how the cap is derived so that it starts biting into text that fits, and
     * this goes red.
     */
    private static void fittingLabelsBesideDomains( final boolean[] ok ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { domainStar() }, new Configuration(), "share-domains" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final int[] cap = { 0 };
        final int[] info = { 0 };
        final int[] text_only = { 0 };
        SwingUtilities.invokeAndWait( () -> {
            tp.getMainPanel().getControlPanel().setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, true );
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            tp.getOptions().setTreeWidthShare( 0.80 );
            paint( tp );
            cap[ 0 ] = tp.maxTipLabelTextWidthForTest();
            info[ 0 ] = tp.getLongestExtNodeInfo();
            text_only[ 0 ] = tp.lengthOfLongestTextOnlyForTest();
        } );
        // precondition: the track really does dominate the reservation, or the case above is not being exercised
        if ( ( info[ 0 ] - text_only[ 0 ] ) < text_only[ 0 ] ) {
            fail( ok, "precondition: the domain track must be the bulk of the reservation -- text " + text_only[ 0 ]
                    + " px of " + info[ 0 ] + " px total" );
        }
        if ( cap[ 0 ] != Integer.MAX_VALUE ) {
            fail( ok, "these labels FIT (" + text_only[ 0 ] + " px of text); only the domain track beside them is "
                    + "wide. Shortening them to " + cap[ 0 ] + " px buys the tree width by deleting text that had "
                    + "room -- the rule is 'shorten what does not fit', not 'shorten whatever helps'" );
        }
        dispose( mf );
    }

    // ---- a gain too small to pay for: the labels keep their text ----------------------------------------------
    /**
     * Shortening the labels is a trade, and a trade can be bad. On the committed {@code crowded-tracks} demo at an
     * 80% share, cutting every label by 12 px buys the tree 12 px -- 0.9% of the width, under the one-percent step
     * the slider itself moves in. The figure would look no different and every label would be missing its tail, so
     * the cap is handed back.
     * <p>
     * Measured with the threshold removed, this same tree reports a 40 px cap, labels at 232 px and a tree of 828;
     * with it, no cap, labels at 244 and a tree of 816. That is the whole behaviour, on a real file.
     */
    private static void smallGainRefused( final boolean[] ok ) throws Exception {
        final File demo = new File( System.getProperty( "user.dir" ) + File.separator + "forester" + File.separator
                + "demo" + File.separator + "crowded-tracks.xml" );
        if ( !demo.exists() ) {
            return; // the demo gallery is not where this run can see it; DemoTreesTest is what guards its presence
        }
        final Phylogeny[] p = ParserUtils.readPhylogenies( demo );
        if ( ( p == null ) || ( p.length < 1 ) ) {
            fail( ok, "precondition: could not read " + demo );
            return;
        }
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { p[ 0 ] }, new Configuration(), "share-smallgain" ) );
        final TreePanel tp = mf[ 0 ].getMainPanel().getCurrentTreePanel();
        final int[] cap = { 0 };
        final int[] tree = { 0 };
        SwingUtilities.invokeAndWait( () -> {
            tp.getOptions().setShowOverview( false );
            tp.setOvOn( false );
            tp.getOptions().setTreeWidthShare( 0.80 );
            paint( tp );
            cap[ 0 ] = tp.maxTipLabelTextWidthForTest();
            tree[ 0 ] = tp.grantedTreeWidthForTest();
        } );
        if ( cap[ 0 ] != Integer.MAX_VALUE ) {
            fail( ok, "the labels were shortened to " + cap[ 0 ] + " px to widen the tree to " + tree[ 0 ]
                    + " px -- a gain this small is below the 1% step the share slider moves in, so it is not worth "
                    + "taking the text off every label" );
        }
        dispose( mf );
    }

    // ---- fixture -----------------------------------------------------------------------------------------------

    /**
     * A star: one root, {@code tips} children all at the SAME distance, so every tip label starts at the same x and
     * the drawn width of a label is simply (rightmost ink in its row) - (tip x). The first tip carries a label far
     * too long to fit; the rest are short, so the long one is unambiguously the widest.
     */
    private static Phylogeny star( final int tips, final boolean one_very_long ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < tips; ++i ) {
            final PhylogenyNode tip = new PhylogenyNode();
            tip.setDistanceToParent( 1.0 );
            if ( one_very_long && ( i == 0 ) ) {
                tip.setName( "TIP0" );
                final Taxonomy t = new Taxonomy();
                // Long enough that the cap falls INSIDE the taxonomy segment. With a shorter one the node-data
                // segment absorbs every truncation and the taxonomy cap is never reached -- which is exactly how a
                // sabotage of it survived the first version of this test.
                t.setScientificName( "Pseudomonadendrogrammatophyllum megalospeciosum subspecies extralongensis "
                        + "variety sesquipedalianum forma nomenclaturalis absurdissima collectio typographica" );
                tip.getNodeData().addTaxonomy( t );
                final Sequence s = new Sequence();
                s.setName( "an extremely long sequence name that no sensible width could ever accommodate whole" );
                tip.getNodeData().addSequence( s );
            }
            else {
                // Deliberately NO shared prefix: names that all begin the same way are shortened on display
                // (AptxUtil.commonNamePrefix, 95% of tips), which silently turned an earlier version of this
                // fixture into twelve two-character labels and measured nothing.
                tip.setName( NAMES[ i % NAMES.length ] + ( i / NAMES.length ) );
            }
            root.addAsChild( tip );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        phy.setName( "share" );
        return phy;
    }

    /** Short names, but each tip carries a 500-residue architecture, so the domain TRACK dominates the label
     *  reservation while the text itself is nowhere near the cap. */
    private static Phylogeny domainStar() {
        final Phylogeny phy = star( 10, false );
        for( final PhylogenyNode tip : phy.getExternalNodes() ) {
            final List<PhylogenyData> ds = new ArrayList<PhylogenyData>();
            ds.add( new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 ) );
            ds.add( new ProteinDomain( "Pkinase", 185, 445, "PF00069", 1e-30 ) );
            final Sequence seq = new Sequence();
            seq.setName( tip.getName() );
            seq.setDomainArchitecture( new DomainArchitecture( ds, 500 ) );
            tip.getNodeData().addSequence( seq );
        }
        return phy;
    }

    // ---- instruments -------------------------------------------------------------------------------------------

    private static BufferedImage paint( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( tp.getTreeColorSet().getBackgroundColor() );
        g.fillRect( 0, 0, W, H );
        tp.paintPhylogeny( g, false, false, W, H, 0, 0 );
        g.dispose();
        return img;
    }

    /**
     * The drawn width of the WIDEST tip label, in pixels of actual ink: the rightmost non-background pixel in that
     * tip's own row band, minus where the label starts. Independent of anything the layout believes -- which is the
     * point, since what is being checked is that the belief and the paint agree.
     */
    private static int drawnLabelWidth( final TreePanel tp, final BufferedImage img ) {
        final PhylogenyNode tip = tp.getPhylogeny().getFirstExternalNode();
        final int bg = tp.getTreeColorSet().getBackgroundColor().getRGB();
        final int y0 = Math.max( 0, (int) tip.getYcoord() - 5 );
        final int y1 = Math.min( H - 1, (int) tip.getYcoord() + 5 );
        final int x_start = Math.max( 0, (int) tip.getXcoord() );
        int right = x_start;
        for( int y = y0; y <= y1; ++y ) {
            for( int x = W - 1; x > right; --x ) {
                if ( differs( img.getRGB( x, y ), bg ) ) {
                    right = x;
                    break;
                }
            }
        }
        return right - x_start;
    }

    /** Antialiased text fades into the background, so "ink" is a visible difference, not any difference at all. */
    private static boolean differs( final int rgb, final int bg ) {
        final Color a = new Color( rgb );
        final Color b = new Color( bg );
        return ( Math.abs( a.getRed() - b.getRed() ) + Math.abs( a.getGreen() - b.getGreen() )
                + Math.abs( a.getBlue() - b.getBlue() ) ) > 90;
    }

    private static String capText( final int cap ) {
        return ( cap == Integer.MAX_VALUE ) ? "uncapped" : ( cap + " px" );
    }

    private static void dispose( final MainFrame[] mf ) throws Exception {
        SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [TreeShareSqueezeTest] " + message );
        ok[ 0 ] = false;
    }

    private TreeShareSqueezeTest() {
        // not instantiable
    }
}

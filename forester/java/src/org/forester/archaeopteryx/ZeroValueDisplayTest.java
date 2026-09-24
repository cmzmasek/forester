// Zero is a VALUE, not an absence.
//
// A branch length of 0 and a support of 0 are facts about a tree and must be shown. Absence is carried by a
// sentinel -- PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT (-1024) for a length, and for support the presence of a
// Confidence at all -- so no display path has any business asking whether the number is non-zero.
//
// Archaeopteryx.js hit this four times in one file (2026-09-23): its Branch Lengths label skipped every
// zero-length branch, its Confidence Values label never drew a support of 0, and its node data box omitted
// "Distance to parent: 0", all from testing a number for truthiness. Nothing here pinned the equivalent
// behaviour, so this does.

package org.forester.archaeopteryx;

import java.awt.GraphicsEnvironment;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.PhylogenyDataUtil;

public class ZeroValueDisplayTest {

    /** ((A:0.0,B:0.2)[support 0]:0.3,(C:0.4,D:0.5)[support 95]:0.6) */
    private static Phylogeny fixture() throws Exception {
        final Phylogeny p = Phylogeny
                .createInstanceFromNhxString( "((A:0.0,B:0.2):0.3,(C:0.4,D:0.5):0.6)" );
        final PhylogenyNode zero_support = p.getExternalNodes().get( 0 ).getParent();
        zero_support.getBranchData().addConfidence( new Confidence( 0, "bootstrap" ) );
        final PhylogenyNode real_support = p.getExternalNodes().get( 2 ).getParent();
        real_support.getBranchData().addConfidence( new Confidence( 95, "bootstrap" ) );
        return p;
    }

    private static final int W = 800;
    private static final int H = 600;

    private static BufferedImage render( final TreePanel tp ) {
        tp.setSize( W, H );
        tp.calcParametersForPainting( W, H );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        tp.printAll( g );
        g.dispose();
        return img;
    }

    /**
     * How many pixels differ between two renders.
     *
     * Deliberately a DIFFERENCE rather than a count of non-white "ink": what colour counts as background
     * depends on the theme, and a standalone run picks up the developer's real saved preferences, so a dark
     * background would make every pixel read as ink and the comparison would say "nothing was drawn" for a
     * perfectly good build.
     */
    private static int pixelsDiffering( final BufferedImage a, final BufferedImage b ) {
        int n = 0;
        for( int x = 0; x < W; ++x ) {
            for( int y = 0; y < H; ++y ) {
                if ( a.getRGB( x, y ) != b.getRGB( x, y ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    public static boolean test() {
        try {
            // (1) the hover card names a zero distance
            final Phylogeny p1 = fixture();
            final PhylogenyNode zero_len = p1.getExternalNodes().get( 0 );
            if ( zero_len.getDistanceToParent() != 0.0 ) {
                System.out.println( "the fixture's zero-length branch is not zero: "
                        + zero_len.getDistanceToParent() );
                return false;
            }
            boolean saw = false;
            for( final NodeHoverText.Row r : NodeHoverText.rows( zero_len ) ) {
                if ( r.toString().startsWith( "Distance to parent" ) ) {
                    saw = true;
                    if ( !r.toString().equals( "Distance to parent: 0" ) ) {
                        System.out.println( "a zero distance is printed oddly: " + r );
                        return false;
                    }
                }
            }
            if ( !saw ) {
                System.out.println( "the hover card dropped a zero distance to parent" );
                return false;
            }
            // deliberate non-behaviour: a node with NO length says nothing, which is what the sentinel is for
            final PhylogenyNode no_len = p1.getExternalNodes().get( 1 );
            no_len.setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
            for( final NodeHoverText.Row r : NodeHoverText.rows( no_len ) ) {
                if ( r.toString().startsWith( "Distance to parent" ) ) {
                    System.out.println( "a node with no branch length still reported one: " + r );
                    return false;
                }
            }

            // (2) the "has any branch lengths" question is the one case where > 0 is right: a tree whose
            // lengths are all zero has none worth drawing to scale. Pinned so that someone tidying this
            // into a sentinel test is caught.
            final Phylogeny all_zero = Phylogeny.createInstanceFromNhxString( "(A:0.0,B:0.0)" );
            if ( AptxUtil.isHasAtLeastOneBranchLengthLargerThanZero( all_zero ) ) {
                System.out.println( "an all-zero tree should not count as having branch lengths" );
                return false;
            }
            if ( !AptxUtil.isHasAtLeastOneBranchLengthLargerThanZero(
                    Phylogeny.createInstanceFromNhxString( "(A:0.0,B:0.1)" ) ) ) {
                System.out.println( "a tree with one positive length should count as having them" );
                return false;
            }

            if ( GraphicsEnvironment.isHeadless() ) {
                return true; // the rest needs a display toolkit
            }
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            final Phylogeny p2 = fixture();
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { p2 }, conf, "zero" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                // dispose in a FINALLY: an early return from a failing check used to skip it, and an
                // undisposed frame keeps the AWT thread -- and so the whole JVM -- alive. The test then
                // hung instead of failing, which is worse than failing.
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().setCheckbox( DisplayOption.WRITE_BRANCH_LENGTH_VALUES, true );
                    mp.getControlPanel().setCheckbox( DisplayOption.WRITE_CONFIDENCE_VALUES, true );
                    final PhylogenyNode zl = tp.getPhylogeny().getExternalNodes().get( 0 );
                    final PhylogenyNode zs = zl.getParent();

                    // (3) the two gates say yes to a zero
                    if ( !tp.shouldWriteBranchLengthForTest( zl ) ) {
                        System.out.println( "the branch-length label is suppressed for a zero-length branch" );
                        ok[ 0 ] = false;
                        return;
                    }
                    if ( !tp.isShowConfidenceValuesForNodeForTest( zs ) ) {
                        System.out.println( "the confidence label is suppressed for a support of zero" );
                        ok[ 0 ] = false;
                        return;
                    }
                    // and the control: the ROOT has no length, so it must still be refused. Without this a
                    // gate stuck at true would pass everything above.
                    if ( tp.shouldWriteBranchLengthForTest( tp.getPhylogeny().getRoot() ) ) {
                        System.out.println( "a branch-length label was offered for the root" );
                        ok[ 0 ] = false;
                        return;
                    }

                    // (4) and each zero reaches the canvas -- checked SEPARATELY, because removing both at
                    // once passes even when only one of them was ever drawn
                    // The FIRST calcParametersForPainting settles the layout (measured: 141 pixels move
                    // between the first and second render, nothing after that), so warm up before the
                    // baseline -- otherwise that settling is counted as the label appearing or vanishing and
                    // the comparison passes whether or not anything was drawn.
                    render( tp );
                    final BufferedImage both = render( tp );
                    zl.setDistanceToParent( PhylogenyDataUtil.BRANCH_LENGTH_DEFAULT );
                    final int by_length = pixelsDiffering( both, render( tp ) );
                    zl.setDistanceToParent( 0.0 );
                    if ( by_length < 1 ) {
                        System.out.println( "removing the zero-length branch changed nothing on the canvas, "
                                + "so its label was never drawn" );
                        ok[ 0 ] = false;
                        return;
                    }
                    final BufferedImage restored = render( tp );
                    if ( pixelsDiffering( both, restored ) != 0 ) {
                        System.out.println( "restoring the zero length did not restore the drawing, so this "
                                + "comparison is measuring something else" );
                        ok[ 0 ] = false;
                        return;
                    }
                    // No render check for the confidence label: measured that printAll does not exercise
                    // paintConfidenceValues at all -- turning WRITE_CONFIDENCE_VALUES on and off with a
                    // support of 95 present changes zero pixels through this path, so a comparison here
                    // would pass whatever the code did. The gate above is the real check; the VALUE-level
                    // behaviour is checked directly below instead of pretended at.
                    if ( !"0".equals( TreePanel.confidenceLabel( List.of( new Confidence( 0, "bootstrap" ) ),
                                                                 false, 0, false, 2 ) ) ) {
                        System.out.println( "a support of zero does not produce the label \"0\": ["
                                + TreePanel.confidenceLabel( List.of( new Confidence( 0, "bootstrap" ) ),
                                                             false, 0, false, 2 )
                                + "]" );
                        ok[ 0 ] = false;
                        return;
                    }
                    // and the neighbouring case: a threshold ABOVE zero is meant to hide it, so the
                    // inclusion of zero is a real decision rather than the absence of a filter
                    if ( TreePanel.confidenceLabel( List.of( new Confidence( 0, "bootstrap" ) ),
                                                    false, 1, false, 2 ).length() != 0 ) {
                        System.out.println( "the min-confidence threshold no longer hides a support of zero" );
                        ok[ 0 ] = false;
                    }
                }
                catch ( final Throwable t ) {
                    t.printStackTrace( System.out );
                    ok[ 0 ] = false;
                }
                finally {
                    ( (JFrame) mf[ 0 ] ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Exception e ) {
            e.printStackTrace( System.out );
            return false;
        }
    }

    public static void main( final String[] args ) {
        if ( test() ) {
            System.out.println( "ZeroValueDisplayTest: OK." );
        }
        else {
            System.out.println( "ZeroValueDisplayTest: FAILED." );
        }
    }
}

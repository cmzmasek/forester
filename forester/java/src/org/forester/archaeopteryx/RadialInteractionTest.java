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
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JToggleButton;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.Options.PHYLOGENY_GRAPHICS_TYPE;
import org.forester.io.parsers.phyloxml.PhyloXmlParser;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.factories.ParserBasedPhylogenyFactory;
import org.forester.phylogeny.iterators.PhylogenyNodeIterator;

/**
 * Radial (circular/unrooted) interaction parity, increment 2: asserts that branch-click clade selection works radially
 * ({@link TreePanel#findBranch} now hit-tests the radial branch geometry), that selecting off a found radial branch
 * populates the selection, and that the "Pulse Found Nodes" halo now renders in circular AND unrooted. The findBranch
 * checks compute a point ON a node's drawn branch from that node's coords and hit-test against the same coords, so they
 * are independent of the display-density-dependent node->device transform. Also (displayTypeControlOk) asserts the
 * phylogram/cladogram P/A/C radios are ENABLED and drive the layout in circular AND unrooted (not just rectangular).
 * Headful; a green no-op when headless.
 */
public final class RadialInteractionTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RadialInteraction: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return findBranchOk() && haloOk() && radialZoomOk() && displayTypeControlOk()
                && persistedLayoutAppliedOnLoadOk() && autoHideControlOk();
    }

    /**
     * "Auto-hide Labels" is LIVE in every layout, and its state reaches both the paint and the export report.
     * Until 2026-09-27 {@code MainFrame.typeChanged} greyed the checkbox out for unrooted -- a leftover from the
     * years when that layout thinned nothing -- while the radial rule already read it (Christian: "is there a reason
     * why label hiding is not controllable via Auto-hide Labels in unrooted?"; there was none). And
     * {@code labelsDynamicallyHidden()} returned false for both radial layouts outright, so a circular or unrooted
     * export that dropped labels never carried the report's warning line. Fixture: the crowded radial-labels demo
     * (96 names at two lengths; a quarter of them hidden at the demo's 1100 x 850 window). The checkbox is driven by the
     * REAL click ({@code doClick}), not {@code setCheckbox}, because the click path is the one the user has.
     */
    private static boolean autoHideControlOk() {
        final boolean[] ok = { true };
        withFrame( "radial-labels-crowded.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            final javax.swing.JCheckBox box = cp.getDynamicallyHideData();
            if ( box == null ) {
                fail( ok, "precondition: the Auto-hide Labels checkbox exists" );
                return;
            }
            // (1) enabled in every layout, on the Type-menu path AND the tab-switch path
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.UNROOTED,
                    PHYLOGENY_GRAPHICS_TYPE.CIRCULAR, PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR } ) {
                frame.typeChanged( gt == PHYLOGENY_GRAPHICS_TYPE.UNROOTED ? frame._unrooted_type_cbmi
                        : gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ? frame._circular_type_cbmi
                                : frame._rectangular_type_cbmi );
                if ( !box.isEnabled() ) {
                    fail( ok, "Auto-hide Labels is greyed out after switching to " + gt + " via the Type menu" );
                }
                tp.setPhylogenyGraphicsType( gt );
                cp.tabChanged();
                if ( !box.isEnabled() ) {
                    fail( ok, "Auto-hide Labels is greyed out after a tab switch to " + gt );
                }
            }
            // (2) UNROOTED: the click drives the rule. On -> off -> on, reading the counters of each export paint
            //     and the report the export returns (its warning line is what the user reads).
            frame.typeChanged( frame._unrooted_type_cbmi );
            if ( !box.isSelected() ) {
                box.doClick();
            }
            //     (700 x 550, not the demo's 1100 x 850: an EXPORT shrinks the font to fit, and at the demo size that
            //     left a single collision -- a fixture that crowds by one label pins nothing worth the name)
            final String on = exportReport( tp, o, 700, 550 );
            final int hidden_on = tp.labelsHiddenForTest();
            box.doClick();
            final String off = exportReport( tp, o, 700, 550 );
            final int hidden_off = tp.labelsHiddenForTest();
            box.doClick();
            exportReport( tp, o, 700, 550 );
            final int hidden_again = tp.labelsHiddenForTest();
            if ( hidden_on < 5 ) {
                fail( ok, "UNROOTED precondition: the 700 x 550 export must crowd (>= 5 hidden), got " + hidden_on );
            }
            if ( ( hidden_off != 0 ) || ( hidden_again != hidden_on ) ) {
                fail( ok, "UNROOTED: clicking Auto-hide Labels must switch the rule off and on again (hidden with it on "
                        + hidden_on + ", off " + hidden_off + ", on again " + hidden_again + ")" );
            }
            if ( !on.contains( "Warning:" ) || !on.contains( "Auto-hide Labels" ) ) {
                fail( ok, "UNROOTED: an export that hid " + hidden_on + " labels must warn in its report; got: " + on );
            }
            if ( off.contains( "Warning:" ) ) {
                fail( ok, "UNROOTED: with Auto-hide Labels off nothing is hidden, so the report must not warn; got: " + off );
            }
            // (3) ...and a ROOMY export with the option on hides nothing, so it must not warn either: the report has
            //     to read what the paint DID, not what the checkbox says
            final String roomy = exportReport( tp, o, 3000, 3000 );
            final int hidden_roomy = tp.labelsHiddenForTest();
            if ( ( hidden_roomy != 0 ) || roomy.contains( "Warning:" ) ) {
                fail( ok, "UNROOTED: a 3000 x 3000 export hides nothing (" + hidden_roomy + " hidden) and must not warn; got: "
                        + roomy );
            }
            // (4) CIRCULAR: same rule, same warning. Cramped means SMALLER here: on the ring the 96 names sit ~26 px
            //     apart at 1100 x 850 and none collide (measured: 0 hidden); at 450 x 450 they do.
            frame.typeChanged( frame._circular_type_cbmi );
            final String circ = exportReport( tp, o, 450, 450 );
            final int hidden_circ = tp.labelsHiddenForTest();
            if ( ( hidden_circ <= 0 ) || !circ.contains( "Warning:" ) ) {
                fail( ok, "CIRCULAR: an export that hid " + hidden_circ + " labels must warn in its report; got: " + circ );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** A PNG export at {@code w x h} px, discarded; returns the export REPORT. The panel's crowding counters read
     *  right after it are those of the export paint (the restore re-lays out but paints only on the next EDT turn). */
    private static String exportReport( final TreePanel tp, final Options o, final int w, final int h ) throws Exception {
        final File png = File.createTempFile( "aptx_autohide_", ".png" );
        try {
            return AptxUtil.writePhylogenyToGraphicsFileAtSize( png.getAbsolutePath(), tp, AptxUtil.GraphicsExportType.PNG,
                    o, new ExportSizeSpec( ExportSizeSpec.Unit.PIXELS, w, h, 72 ) );
        }
        finally {
            png.delete();
        }
    }

    /** The persisted graphics type AND P/A/C display type must be APPLIED when a tree is opened, so a restart
     *  restores the layout (the user-reported gap). Exercises the real load path: with an Options as restored from
     *  GuiPreferences (circular + aligned phylogram), AptxUtil.addPhylogeniesToTabs must open the next branch-length
     *  tree CIRCULAR (a new TreePanel reads the graphics type at construction) AND ALIGNED (the display-type
     *  auto-detect lookAtSomeTreePropertiesForAptxControlSettings now reads the persisted preference); a persisted
     *  CLADOGRAM preference must open it as a cladogram (not the old hardcoded phylogram). The circular tree must
     *  also render without collapsing/throwing when its type comes from Options at construction. */
    private static boolean persistedLayoutAppliedOnLoadOk() {
        final boolean[] ok = { true };
        withFrame( "scale-axis.xml", ( frame, tp0, o ) -> {
            final MainPanel mp = frame.getMainPanel();
            final ControlPanel cp = mp.getControlPanel();
            final Configuration conf = frame.getConfiguration();
            // (1) Options as restored from preferences: circular + aligned phylogram
            o.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
            o.setPhylogenyDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM );
            AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { freshDemo( "scale-axis.xml" ) }, "a", "a", conf, mp );
            final TreePanel tpa = mp.getCurrentTreePanel();
            if ( tpa.getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
                fail( ok, "a persisted CIRCULAR graphics type must be applied to a newly opened tree (got "
                        + tpa.getPhylogenyGraphicsType() + ")" );
            }
            if ( cp.getTreeDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM ) {
                // the persisted A is applied to the tab's store AND shown: in circular every phylogram carries its
                // labels to the ring, so A is exactly what circular draws and P is the greyed flavour there
                // (Christian, 2026-09-27; for a few hours this had been the other way round)
                fail( ok, "a persisted ALIGNED_PHYLOGRAM (A) choice on a newly opened branch-length tree must show as "
                        + "A in circular (P is the greyed one there; got " + cp.getTreeDisplayType() + ")" );
            }
            if ( cp.getDisplayAsUnalignedPhylogramRb().isEnabled() || !cp.getDisplayAsAlignedPhylogramRb().isEnabled() ) {
                fail( ok, "in circular 'P' is greyed and 'A' live (P.enabled=" + cp.getDisplayAsUnalignedPhylogramRb().isEnabled()
                        + " A.enabled=" + cp.getDisplayAsAlignedPhylogramRb().isEnabled() + ")" );
            }
            // the circular tree must actually FAN OUT to a real radius (not radius-collapse to a line at the centre,
            // the 0.11.18 bug) when its type comes from Options at construction
            final int w = 500, h = 500;
            tpa.setSize( w, h );
            tpa.calcParametersForPainting( w, h );
            tpa.paintPhylogeny( new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB ).createGraphics(), false, false,
                    w, h, 0, 0 );
            if ( maxTipRadius( tpa ) < 30 ) {
                fail( ok, "a circular tree opened from a persisted CIRCULAR type must fan out to a real radius, not "
                        + "collapse to the centre (max tip radius " + maxTipRadius( tpa ) + ")" );
            }
            // (2) the P/A/C write-back is confined to a real user CLICK: a doClick updates the persisted Options
            // default, but an INTERNAL setTreeDisplayType (as the load/tab/reset paths use) must NOT -- else an
            // auto-detect on the next tree would silently overwrite the saved preference. The click goes to a LIVE
            // radio: in circular that is C (or A); a doClick on the greyed P does nothing, by Swing's own rule
            cp.getDisplayAsCladogramRb().doClick();
            if ( o.getPhylogenyDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM ) {
                fail( ok, "a P/A/C radio CLICK must update the persisted Options default (got "
                        + o.getPhylogenyDisplayType() + ")" );
            }
            cp.setTreeDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.UNALIGNED_PHYLOGRAM ); // internal call (no user gesture)
            if ( o.getPhylogenyDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM ) {
                fail( ok, "an INTERNAL setTreeDisplayType must NOT change the persisted default (only a click does); "
                        + "got " + o.getPhylogenyDisplayType() );
            }
            // (3) the node-edit re-detect path (lookAtRealBranchLengthsForAptxControlSettings) also honors the
            // persisted preference for a fully-branch-lengthed tree
            o.setPhylogenyDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM );
            AptxUtil.lookAtRealBranchLengthsForAptxControlSettings( tpa.getPhylogeny(), cp );
            if ( cp.getTreeDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.ALIGNED_PHYLOGRAM ) {
                fail( ok, "the node-edit branch-length re-detect must honor the persisted P/A/C preference, shown as "
                        + "A in circular (got " + cp.getTreeDisplayType() + ")" );
            }
            // (4) a persisted UNROOTED graphics type + CLADOGRAM preference opens a branch-length tree unrooted AND
            // as a cladogram (was: forced UNALIGNED phylogram, and UNROOTED never persisted at all)
            o.setPhylogenyGraphicsType( PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
            o.setPhylogenyDisplayType( Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM );
            AptxUtil.addPhylogeniesToTabs( new Phylogeny[] { freshDemo( "scale-axis.xml" ) }, "c", "c", conf, mp );
            if ( mp.getCurrentTreePanel().getPhylogenyGraphicsType() != PHYLOGENY_GRAPHICS_TYPE.UNROOTED ) {
                fail( ok, "a persisted UNROOTED graphics type must be applied to a newly opened tree (got "
                        + mp.getCurrentTreePanel().getPhylogenyGraphicsType() + ")" );
            }
            if ( cp.getTreeDisplayType() != Options.PHYLOGENY_DISPLAY_TYPE.CLADOGRAM ) {
                fail( ok, "a persisted CLADOGRAM (C) choice must open a branch-length tree as a cladogram (got "
                        + cp.getTreeDisplayType() + ")" );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** Max external-tip distance from the root after a radial paint -- 0 (collapsed to the centre) vs a real radius. */
    private static double maxTipRadius( final TreePanel tp ) {
        final PhylogenyNode root = tp.getPhylogeny().getRoot();
        double max = 0;
        for ( final PhylogenyNodeIterator it = tp.getPhylogeny().iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            max = Math.max( max, Math.hypot( n.getXcoord() - root.getXcoord(), n.getYcoord() - root.getYcoord() ) );
        }
        return max;
    }

    private static Phylogeny freshDemo( final String demo ) throws Exception {
        final File file = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
        return ParserBasedPhylogenyFactory.getInstance().create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
    }

    /** The phylogram/cladogram (P/A/C) buttons must WORK in CIRCULAR and UNROOTED, not just rectangular: for a
     *  branch-length tree P and C must be ENABLED in every radial layout (the pre-0.11.7 "&amp;&amp; != CIRCULAR"
     *  force-disable greyed them out in circular even though the paint responds).
     *  <p>
     *  "A" (aligned phylogram) is the one exception, and it is deliberate: aligning tip LABELS needs somewhere to
     *  pin them -- the common column of the rectangular layouts. UNROOTED has no such place, because its tips
     *  radiate in every direction; CIRCULAR has the outer ring but EVERY circular phylogram already carries its
     *  labels there (TreePanel.circularLabelsOnRing, Christian 2026-09-27), so A would change nothing. So A is
     *  DEAD in both radial layouts. Both enable-logic sites are exercised: the Type-menu path (MainFrame.typeChanged) AND the
     *  tab-switch path (ControlPanel.tabChanged) -- reverting either gate is caught. Driving the buttons (a real
     *  doClick gesture) must then CHANGE the radial layout: cladogram lays tips ~uniformly, phylogram lays them by
     *  root-distance -> the min/max tip-radius spread (a scale-independent ratio) DIFFERS between the two
     *  (abs-difference, not "clado &gt; phylo", so it is robust across circular/unrooted and to the demo's exact
     *  branch lengths). */
    private static boolean displayTypeControlOk() {
        final boolean[] ok = { true };
        withFrame( "scale-axis.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                // (1) the Type-menu path (the primary user route): drive the real MainFrame.typeChanged via the layout's
                // checkbox item -> it must ENABLE the radios in the target radial layout. Each radial layout has ONE
                // dead phylogram flavour: unrooted cannot align (A dead), circular cannot NOT align -- every circular
                // phylogram carries its labels to the ring -- so P is dead there and A is what it draws (Christian,
                // 2026-09-27: greying A in circular had told the user the opposite of the ring in front of them)
                final boolean circular = gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR;
                frame.typeChanged( circular ? frame._circular_type_cbmi : frame._unrooted_type_cbmi );
                if ( !pacEnabledAsExpected( cp, circular ) ) {
                    fail( ok, "P/A/C buttons wrong in " + gt + " via the Type menu (typeChanged): "
                            + pacRadioState( cp ) + " (expected A=" + circular + " P=" + !circular + ")" );
                    continue;
                }
                // (2) the tab-switch path: ControlPanel.tabChanged must reach the SAME conclusion for this layout.
                tp.setPhylogenyGraphicsType( gt );
                cp.tabChanged();
                if ( !pacEnabledAsExpected( cp, circular ) ) {
                    fail( ok, "P/A/C buttons wrong in " + gt + " via a tab switch (tabChanged): "
                            + pacRadioState( cp ) + " (expected A=" + circular + " P=" + !circular + ")" );
                    continue;
                }
                // (3) driving them must change the radial layout (else the radios are inert in this layout) -- the
                // LIVE phylogram button of this layout against the cladogram
                final double clado = tipRadiusSpread( tp, cp.getDisplayAsCladogramRb() );
                final double phylo = tipRadiusSpread( tp, circular ? cp.getDisplayAsAlignedPhylogramRb()
                        : cp.getDisplayAsUnalignedPhylogramRb() );
                if ( Math.abs( clado - phylo ) <= 0.05 ) {
                    fail( ok, "driving the P/C radios must change the " + gt + " layout (cladogram tip-radius spread "
                            + clado + " vs phylogram spread " + phylo + ")" );
                }
            }
        }, ok );
        return ok[ 0 ];
    }

    /** C is always live for a branch-length tree; of P and A, circular greys P (it always aligns) and unrooted
     *  greys A (it never can). */
    private static boolean pacEnabledAsExpected( final ControlPanel cp, final boolean circular ) {
        return cp.getDisplayAsCladogramRb().isEnabled()
                && ( cp.getDisplayAsAlignedPhylogramRb().isEnabled() == circular )
                && ( cp.getDisplayAsUnalignedPhylogramRb().isEnabled() == !circular );
    }

    private static String pacRadioState( final ControlPanel cp ) {
        return "P=" + cp.getDisplayAsUnalignedPhylogramRb().isEnabled() + " A="
                + cp.getDisplayAsAlignedPhylogramRb().isEnabled() + " C=" + cp.getDisplayAsCladogramRb().isEnabled();
    }

    /** doClick the given display-type radio (the real user gesture: setTreeDisplayType + showWhole), paint the radial
     *  canvas at its (square) preferred size so node coords are live device pixels, and return the external-tip
     *  min/max radius ratio about the root (1 = all tips on one ring, &lt;1 = spread by distance). */
    private static double tipRadiusSpread( final TreePanel tp, final JToggleButton display_type_rb ) {
        display_type_rb.doClick();
        final java.awt.Dimension pref = tp.getPreferredSize();
        final int w = Math.max( 50, pref.width );
        final int h = Math.max( 50, pref.height );
        tp.setSize( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        tp.paintPhylogeny( img.createGraphics(), false, false, w, h, 0, 0 );
        final PhylogenyNode root = tp.getPhylogeny().getRoot();
        double min = Double.MAX_VALUE;
        double max = 0;
        for ( final PhylogenyNodeIterator it = tp.getPhylogeny().iteratorExternalForward(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            final double r = Math.hypot( n.getXcoord() - root.getXcoord(), n.getYcoord() - root.getYcoord() );
            min = Math.min( min, r );
            max = Math.max( max, r );
        }
        return ( max > 0 ) ? ( min / max ) : 1;
    }

    /** Radial zoom is decoupled from the rectangular x/y-distance: ONE square-canvas diameter is the single knob.
     *  showWhole fits it to the viewport; BOTH the X and the Y zoom scale that same diameter (symmetric, not the old
     *  asymmetric min(w,h) behaviour); zoom-out shrinks it back. Checked via TreePanel.radialDiameter(). */
    private static boolean radialZoomOk() {
        final boolean[] ok = { true };
        withFrame( "scale-axis.xml", ( frame, tp, o ) -> {
            final ControlPanel cp = tp.getControlPanel();
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                cp.showWhole();
                final int fit = tp.radialDiameter();
                if ( fit < 50 ) {
                    fail( ok, "precondition: a sane fit diameter in " + gt + " (got " + fit + ")" );
                    return;
                }
                // fit == min(viewport); the zoom DIAMETER is the square canvas, but the PREFERRED size is PADDED to at
                // least the viewport so a tree zoomed out below fit stays centred in the window (not top-left) -- so at
                // fit the preferred size == the viewport (max(diameter, viewport) per dim).
                final int vw = frame.getMainPanel().getSizeOfViewport().width;
                final int vh = frame.getMainPanel().getSizeOfViewport().height;
                if ( Math.abs( fit - Math.min( vw, vh ) ) > 2 ) {
                    fail( ok, "showWhole must fit the radial canvas to min(viewport) in " + gt + " (fit " + fit
                            + " vs " + Math.min( vw, vh ) + ")" );
                }
                final java.awt.Dimension pref = tp.getPreferredSize();
                if ( ( pref.width != Math.max( fit, vw ) ) || ( pref.height != Math.max( fit, vh ) ) ) {
                    fail( ok, "the radial preferred size must be max(diameter, viewport) per dim in " + gt + " ("
                            + pref.width + "x" + pref.height + " vs " + Math.max( fit, vw ) + "x" + Math.max( fit, vh )
                            + ")" );
                }
                cp.zoomInX( 1.25f, 1.2f );
                final int in_x = tp.radialDiameter();
                cp.zoomInY( 1.25f );
                final int in_y = tp.radialDiameter();
                cp.zoomOutX( 0.8f, 0.83f );
                cp.zoomOutY( 0.8f );
                final int out = tp.radialDiameter();
                if ( Math.abs( in_x - ( fit * 1.25 ) ) > 3 ) {
                    fail( ok, "zoomInX must scale the radial diameter in " + gt + " (" + fit + " -> " + in_x + ")" );
                }
                if ( Math.abs( in_y - ( in_x * 1.25 ) ) > 3 ) {
                    fail( ok, "zoomInY must ALSO scale the radial diameter (symmetric) in " + gt + " (" + in_x + " -> "
                            + in_y + ")" );
                }
                if ( Math.abs( out - ( in_y * 0.64 ) ) > 3 ) {
                    fail( ok, "zoom out must shrink the radial diameter in " + gt + " (" + in_y + " -> " + out + ")" );
                }
            }
        }, ok );
        return ok[ 0 ];
    }

    /** findBranch returns the node whose incoming branch is under a point placed ON that branch -- for a LEAF and an
     *  INTERNAL clade root, in CIRCULAR and UNROOTED; and selecting off a found internal branch selects its tips. */
    private static boolean findBranchOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 820, h = 820;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            tp.getControlPanel().setActionWhenNodeClicked( ControlPanel.NodeClickAction.SELECT_NODES );
            final PhylogenyNode leaf = firstMatch( tp.getPhylogeny(), true );
            final PhylogenyNode clade = firstMatch( tp.getPhylogeny(), false );
            if ( ( leaf == null ) || ( clade == null ) ) {
                fail( ok, "precondition: need a leaf and an internal clade with a non-root parent" );
                return;
            }
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                tp.calcParametersForPainting( w, h );
                AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ); // assigns radial node coords
                final PhylogenyNode root = tp.getPhylogeny().getRoot();
                if ( ( Math.abs( leaf.getXcoord() - root.getXcoord() )
                        + Math.abs( leaf.getYcoord() - root.getYcoord() ) ) < 5 ) {
                    fail( ok, "precondition: radial layout did not assign meaningful node coords in " + gt );
                    return;
                }
                // a click ON the leaf's own incoming branch resolves to the leaf
                final double[] lp = branchPoint( leaf, root, gt );
                final PhylogenyNode leaf_hit = tp.findBranch( (int) Math.round( lp[ 0 ] ), (int) Math.round( lp[ 1 ] ) );
                if ( leaf_hit != leaf ) {
                    fail( ok, "findBranch on a leaf's branch must return that leaf in " + gt + " (got " + leaf_hit
                            + ")" );
                }
                // a click ON the clade root's branch resolves to the clade root, and selecting it selects its tips
                final double[] cp = branchPoint( clade, root, gt );
                final PhylogenyNode clade_hit = tp.findBranch( (int) Math.round( cp[ 0 ] ), (int) Math.round( cp[ 1 ] ) );
                if ( clade_hit != clade ) {
                    fail( ok, "findBranch on an internal clade's branch must return that clade in " + gt + " (got "
                            + clade_hit + ")" );
                    continue;
                }
                // hovering that branch (Select-Node(s) mode) sets the hover-preview target -- the source of the
                // translucent on-screen preview, which is now dispatched in the radial branches (was dead before).
                // Assert it resolves to THAT clade as a subtree (not just any node) -- proving the radial findBranch
                // branch-hover path, not the node-hover fallthrough (mouseMoved tries findNode first).
                tp.clearHoverPreview();
                tp.mouseMoved( new java.awt.event.MouseEvent( tp, java.awt.event.MouseEvent.MOUSE_MOVED,
                        System.currentTimeMillis(), 0, (int) Math.round( cp[ 0 ] ), (int) Math.round( cp[ 1 ] ), 0,
                        false ) );
                if ( ( tp.hoverNodeForTest() != clade ) || !tp.hoverIsSubtreeForTest() ) {
                    fail( ok, "hovering a radial branch must set that clade as the subtree hover-preview target in "
                            + gt + " (got " + tp.hoverNodeForTest() + ", subtree=" + tp.hoverIsSubtreeForTest() + ")" );
                }
                tp.clearHoverPreview();
                tp.setFoundNodes0( null );
                tp.selectSubtreeTips( clade_hit );
                final Set<Long> found = tp.getFoundNodes0();
                if ( ( found == null ) || ( found.size() != clade.getNumberOfExternalNodes() ) ) {
                    fail( ok, "branch-click selecting a clade must select all its tips in " + gt + " (selected "
                            + ( found == null ? "null" : found.size() ) + " of " + clade.getNumberOfExternalNodes()
                            + ")" );
                }
                tp.setFoundNodes0( null );
            }
        }, ok );
        return ok[ 0 ];
    }

    /** The "Pulse Found Nodes" halo now renders in circular AND unrooted: with a found node set, turning the pulse ON
     *  adds a translucent found-colour disc (a static glow in an export) -> more tinted ink than with it OFF. */
    private static boolean haloOk() {
        final boolean[] ok = { true };
        withFrame( "colorize-by-rank.xml", ( frame, tp, o ) -> {
            final int w = 780, h = 780;
            o.setGraphicsExportWhiteBackground( true );
            frame.showWhole();
            tp.setSize( w, h );
            final PhylogenyNode leaf = firstMatch( tp.getPhylogeny(), true );
            if ( leaf == null ) {
                fail( ok, "precondition: need a leaf to mark as found" );
                return;
            }
            final Set<Long> found = new HashSet<Long>();
            found.add( leaf.getId() );
            tp.setFoundNodes0( found );
            for ( final PHYLOGENY_GRAPHICS_TYPE gt : new PHYLOGENY_GRAPHICS_TYPE[] { PHYLOGENY_GRAPHICS_TYPE.CIRCULAR,
                    PHYLOGENY_GRAPHICS_TYPE.UNROOTED } ) {
                tp.setPhylogenyGraphicsType( gt );
                tp.calcParametersForPainting( w, h );
                o.setPulseFoundNodes( false );
                final int no_halo = countTinted( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                o.setPulseFoundNodes( true );
                final int with_halo = countTinted( AptxUtil.renderPhylogenyToImage( w, h, tp, o, false, 1, false ) );
                if ( with_halo <= ( no_halo + 60 ) ) {
                    fail( ok, "the found-node pulse halo must render in " + gt + " (tinted ink " + with_halo + " vs "
                            + no_halo + ")" );
                }
            }
            o.setPulseFoundNodes( false );
            tp.setFoundNodes0( null );
        }, ok );
        return ok[ 0 ];
    }

    /** A point on {@code node}'s drawn incoming branch. UNROOTED = the midpoint of the straight parent->node line,
     *  which is EXACTLY the segment paintUnrooted draws (drawLine(parent, node)) -- an independent paint-anchored
     *  oracle. CIRCULAR = the midpoint of the radial leg (node -> the point on the parent's radius at the node's
     *  angle); this mirrors both findBranch's reconstruction AND paintBranchCircular's drawn leg (verified equal:
     *  the node sits at root + r*(cos,sin), so atan2(node-root) recovers its stored angle and the inward point matches
     *  drawLine's), so it exercises dispatch + coordinate-space + tolerance + collapse-skip rather than being a fully
     *  independent geometry oracle. */
    private static double[] branchPoint( final PhylogenyNode node, final PhylogenyNode root,
            final PHYLOGENY_GRAPHICS_TYPE gt ) {
        final PhylogenyNode p = node.getParent();
        if ( gt == PHYLOGENY_GRAPHICS_TYPE.CIRCULAR ) {
            final double angle = Math.atan2( node.getYcoord() - root.getYcoord(), node.getXcoord() - root.getXcoord() );
            final double pdx = p.getXcoord() - root.getXcoord(), pdy = p.getYcoord() - root.getYcoord();
            final double parent_radius = Math.sqrt( ( pdx * pdx ) + ( pdy * pdy ) );
            final double inward_x = root.getXcoord() + ( Math.cos( angle ) * parent_radius );
            final double inward_y = root.getYcoord() + ( Math.sin( angle ) * parent_radius );
            return new double[] { ( node.getXcoord() + inward_x ) / 2.0, ( node.getYcoord() + inward_y ) / 2.0 };
        }
        return new double[] { ( node.getXcoord() + p.getXcoord() ) / 2.0, ( node.getYcoord() + p.getYcoord() ) / 2.0 };
    }

    /** First external node ({@code leaf}=true) or first internal node ({@code leaf}=false) with a non-root, non-collapsed
     *  parent, itself not collapsed; internal must have >= 2 tips (so its branch is a real clade selection). */
    private static PhylogenyNode firstMatch( final Phylogeny phy, final boolean leaf ) {
        for ( final java.util.Iterator<PhylogenyNode> it = phy.iteratorPreorder(); it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            if ( n.isRoot() || n.isCollapse() || n.getParent().isRoot() || n.getParent().isCollapse() ) {
                continue;
            }
            if ( leaf && n.isExternal() ) {
                return n;
            }
            if ( !leaf && !n.isExternal() && ( n.getNumberOfExternalNodes() >= 2 ) ) {
                return n;
            }
        }
        return null;
    }

    private interface FrameBody {
        void run( MainFrame frame, TreePanel tp, Options o ) throws Exception;
    }

    private static void withFrame( final String demo, final FrameBody body, final boolean[] ok ) {
        try {
            final File file = new File( System.getProperty( "user.dir" ), "forester/demo/" + demo );
            if ( !file.exists() ) {
                fail( ok, "demo tree missing: " + file.getAbsolutePath() );
                return;
            }
            final Phylogeny phy = ParserBasedPhylogenyFactory.getInstance()
                    .create( file, PhyloXmlParser.createPhyloXmlParser() )[ 0 ];
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, conf, "radial" ) );
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    body.run( frame, frame.getMainPanel().getCurrentTreePanel(), frame.getOptions() );
                }
                catch ( final Throwable t ) {
                    fail( ok, "unexpected: " + t );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = TestFail.here();
        }
    }

    /** Count of tinted (chromatic) pixels -- a low saturation threshold so the translucent halo disc registers. */
    private static int countTinted( final BufferedImage img ) {
        int n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                final float[] hsb = java.awt.Color.RGBtoHSB( ( rgb >> 16 ) & 0xFF, ( rgb >> 8 ) & 0xFF, rgb & 0xFF,
                        null );
                if ( ( hsb[ 1 ] >= 0.12f ) && ( hsb[ 2 ] >= 0.30f ) ) {
                    ++n;
                }
            }
        }
        return n;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [RadialInteractionTest] " + msg );
        ok[ 0 ] = false;
    }

    private RadialInteractionTest() {
    }
}

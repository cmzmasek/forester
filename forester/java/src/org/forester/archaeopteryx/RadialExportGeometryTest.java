// An export must leave the panel's layout exactly as it found it.
//
// The radial layouts centre the figure on the CANVAS (TreePanel.radialCanvasCenterX/Y), which for a full
// export is the EXPORT size rather than the panel's. Those coordinates are not local to the paint: the
// centre and radius are stored on the panel, and every node's x/y is rewritten as the tree is drawn. So a
// fixed-size PDF/PNG/SVG export of a circular or unrooted tree used to leave the panel describing the
// export canvas, and every ring hit-test -- annotation cells, the domain rollover, the legend anchor --
// answered from that geometry until something happened to repaint.
//
// The export's layout is deliberately NOT restored -- rendering at a size and reading back the coordinates it
// produced is how this codebase checks its own drawing. What the panel does instead is know that its layout is
// no longer the screen's (TreePanel.hitTestableLayout) and decline every screen hit-test until the repaint it
// schedules has laid the tree out for the panel again. See TreePanel.paintPhylogeny.

package org.forester.archaeopteryx;

import java.awt.GraphicsEnvironment;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.Sequence;

public class RadialExportGeometryTest {

    private static final int W = 900, H = 700;      // the "screen"
    private static final int EW = 2400, EH = 1800;  // a fixed-size export, deliberately a different shape

    private static Phylogeny fixture() throws Exception {
        final StringBuilder sb = new StringBuilder( "(" );
        for( int i = 0; i < 24; ++i ) {
            sb.append( i > 0 ? "," : "" ).append( "tip_" ).append( i ).append( ":0." ).append( 10 + i );
        }
        final Phylogeny p = Phylogeny.createInstanceFromNhxString( sb.append( ")" ).toString() );
        int i = 0;
        for( final PhylogenyNode tip : p.getExternalNodes() ) {
            final List<PhylogenyData> ds = new ArrayList<PhylogenyData>();
            ds.add( new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 ) );
            ds.add( new ProteinDomain( "Pkinase", 185, 445, "PF00069", 1e-30 ) );
            final Sequence seq = new Sequence();
            seq.setName( tip.getName() );
            seq.setDomainArchitecture( new DomainArchitecture( ds, 500 ) );
            tip.getNodeData().addSequence( seq );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:host", ( i % 3 == 0 ) ? "bat" : "bird", "", "xsd:string",
                                          AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:load", String.valueOf( 1 + ( i % 7 ) ), "", "xsd:decimal",
                                          AppliesTo.NODE ) );
            tip.getNodeData().setProperties( pl );
            ++i;
        }
        return p;
    }

    private static void paint( final TreePanel tp, final boolean export ) {
        final int w = export ? EW : W, h = export ? EH : H;
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        if ( export ) {
            tp.paintPhylogeny( g, false, true, EW, EH, 0, 0 ); // to_graphics_file, full canvas
        }
        else {
            tp.setSize( W, H );
            tp.calcParametersForPainting( W, H );
            tp.printAll( g );
        }
        g.dispose();
    }

    /** every node's laid-out position, in preorder */
    private static float[] coords( final Phylogeny p ) {
        final List<Float> v = new ArrayList<Float>();
        for( final org.forester.phylogeny.iterators.PhylogenyNodeIterator it = p.iteratorPreorder();
                it.hasNext(); ) {
            final PhylogenyNode n = it.next();
            v.add( Float.valueOf( n.getXcoord() ) );
            v.add( Float.valueOf( n.getYcoord() ) );
        }
        final float[] out = new float[ v.size() ];
        for( int i = 0; i < out.length; ++i ) {
            out[ i ] = v.get( i ).floatValue();
        }
        return out;
    }

    /**
     * The contract: an export may leave the tree laid out for its own canvas -- reading those coordinates back
     * is how the render tests check their drawing -- but the panel must not then ANSWER a screen point from
     * them. It declines until the repaint the export schedules has laid the tree out for the panel again.
     *
     * So this drives all three states in turn, and each one is asserted: a hit before, silence after the
     * export, and the SAME hit once the panel has repainted. The middle state is the fix; the third is what
     * stops the fix from being "return null forever", which would pass a weaker test just as well.
     */
    private static boolean declinesUntilRepainted( final TreePanel tp, final Phylogeny phy, final String what ) {
        paint( tp, false );
        final java.awt.Point centre_before = tp.circularCenterForTest();
        final int radius_before = tp.circularRadiusForTest();
        final float[] coords_before = coords( phy );
        TreePanel.DomainHit hit_before = null;
        int hx = 0, hy = 0;
        for( int y = 0; ( y < H ) && ( hit_before == null ); y += 3 ) {
            for( int x = 0; x < W; x += 3 ) {
                final TreePanel.DomainHit h = tp.domainAt( x, y );
                if ( h != null ) { hit_before = h; hx = x; hy = y; break; }
            }
        }
        if ( hit_before == null ) {
            System.out.println( what + ": nothing was hit before the export, so this pins nothing" );
            return false;
        }
        boolean ok = true;
        // Watch for the repaint the export must schedule. Without it the panel would decline for ever rather
        // than for one frame -- no rollover, no hover card, no node click -- which is WORSE than the defect
        // this fixes, so hand-driving the repaint below would leave the load-bearing half unasserted.
        javax.swing.RepaintManager.currentManager( tp ).markCompletelyClean( tp );
        paint( tp, true );   // <-- the export, at a different canvas size
        if ( javax.swing.RepaintManager.currentManager( tp ).getDirtyRegion( tp ).isEmpty() ) {
            System.out.println( what + ": the export did not ask the panel to repaint, so the layout would "
                    + "stay declined until something else happened to paint" );
            ok = false;
        }
        // The fixture only means something if the export really did move the layout out from under the panel.
        // If it did not, the silence below would prove nothing at all.
        final float[] coords_after = coords( phy );
        int moved = 0;
        for( int i = 0; i < coords_before.length; ++i ) {
            if ( Math.abs( coords_before[ i ] - coords_after[ i ] ) > 0.001f ) { ++moved; }
        }
        final boolean geom_moved = ( radius_before != tp.circularRadiusForTest() )
                || ( ( centre_before == null ) != ( tp.circularCenterForTest() == null ) )
                || ( ( centre_before != null ) && !centre_before.equals( tp.circularCenterForTest() ) );
        if ( ( moved == 0 ) && !geom_moved ) {
            System.out.println( what + ": the export left the layout untouched, so this fixture cannot show "
                    + "the defect it exists for (canvas " + EW + "x" + EH + " vs panel " + W + "x" + H + ")" );
            return false;
        }
        if ( tp.hitTestableLayout() ) {
            System.out.println( what + ": the panel still calls its layout hit-testable after an export that "
                    + "moved " + moved + " coordinates" );
            ok = false;
        }
        // Sweep the whole EXPORT canvas, not the one screen point. Asking only at the point that hit before
        // proves nothing: the export moved the tree away from there, so a hit-test with no gate at all returns
        // null there too and a removed gate survives the check. The tree is now laid out somewhere inside
        // EW x EH, so that is where an ungated hit-test would answer -- and the panel must answer nowhere.
        int answered = 0;
        String first = null;
        // Stop early once the point is made: with the gate removed EVERY one of these answers, and counting all
        // of them turns the failing run -- the run that matters -- into half a million hit-tests.
        for( int y = 0; ( y < EH ) && ( answered < 25 ); y += 5 ) {
            for( int x = 0; ( x < EW ) && ( answered < 25 ); x += 5 ) {
                final TreePanel.DomainHit d = tp.domainAt( x, y );
                final PhylogenyNode n = tp.findNode( x, y );
                final TreePanel.AnnotationCell c = tp.annotationCellAt( x, y );
                final int grab = tp.annotationColumnGrabbedAt( x, y );
                final int slot = tp.annotationColumnInsertionSlotAt( x, y );
                if ( ( d == null ) && ( n == null ) && ( c == null ) && ( grab < 0 ) && ( slot < 0 ) ) {
                    continue;
                }
                ++answered;
                if ( first == null ) {
                    first = "(" + x + "," + y + ") -> "
                            + ( d != null ? "domain " + d.domain().getName()
                                          : n != null ? "node " + n.getName()
                                          : c != null ? "a cell"
                                          : grab >= 0 ? "a GRABBED column (" + grab + ")"
                                                      : "a drop slot (" + slot + ")" );
                }
            }
        }
        if ( answered > 0 ) {
            System.out.println( what + ": the panel answered " + answered + "+ points from the export's"
                    + " geometry, e.g. " + first );
            ok = false;
        }
        paint( tp, false );  // <-- the repaint the export scheduled
        if ( !tp.hitTestableLayout() ) {
            System.out.println( what + ": the panel never took its layout back after repainting" );
            return false;
        }
        final TreePanel.DomainHit hit_again = tp.domainAt( hx, hy );
        if ( ( hit_again == null ) || ( hit_again.node() != hit_before.node() )
                || ( hit_again.domain() != hit_before.domain() ) ) {
            System.out.println( what + ": after repainting, the rollover at (" + hx + "," + hy + ") reports "
                    + ( hit_again == null ? "nothing" : hit_again.domain().getName() + " on "
                            + hit_again.node().getName() )
                    + " instead of " + hit_before.domain().getName() + " on " + hit_before.node().getName() );
            ok = false;
        }
        return ok;
    }

    /**
     * Two more ways the layout stops being the panel's, neither of which the radial sweep above can see.
     * <p>
     * The LATCH: once an export has invalidated the layout, only a SCREEN paint may declare it good again. A
     * plain assignment reads as equivalent and is not -- a second export that moves nothing would clear an
     * invalidation the first one asked a repaint to fix, and every hit-test would answer from the first
     * export's coordinates.
     * <p>
     * The FIXED-SIZE path: layoutForExportSize re-lays out EVERY display type through
     * calcParametersForPainting, so there the rectangular family is on the export's coordinates too -- the one
     * case the "rectangular is never re-centred" exemption does not cover.
     */
    private static boolean latchAndFixedSize( final MainFrame frame, final TreePanel tp, final Phylogeny phy ) {
        boolean ok = true;
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
        paint( tp, false );
        paint( tp, true );                       // invalidates: the canvas differs from the panel
        if ( tp.hitTestableLayout() ) {
            System.out.println( "latch: a radial export did not invalidate at all" );
            return false;
        }
        // a SECOND export that would not invalidate on its own must not clear the first one's invalidation
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_ARGB );
        final java.awt.Graphics2D g = img.createGraphics();
        tp.paintPhylogeny( g, false, true, W, H, 0, 0 );
        g.dispose();
        if ( tp.hitTestableLayout() ) {
            System.out.println( "latch: a second export declared the layout the screen's again, although no "
                    + "screen paint has happened since the first one moved it" );
            ok = false;
        }
        paint( tp, false );
        if ( !tp.hitTestableLayout() ) {
            System.out.println( "latch: a screen paint must be the way back" );
            return false;
        }
        // the fixed-size export path, in the RECTANGULAR family: layoutForExportSize re-lays it out
        final int[] prior = tp.layoutForExportSize( EW, EH );
        if ( tp.hitTestableLayout() ) {
            System.out.println( "fixed-size: laying the tree out for an export frame must stop the hit-tests "
                    + "answering -- calcParametersForPainting has moved every row" );
            ok = false;
        }
        tp.restoreLayoutAfterExport( prior );
        paint( tp, false );
        if ( !tp.hitTestableLayout() ) {
            System.out.println( "fixed-size: the panel never took its layout back" );
            ok = false;
        }
        return ok;
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final Phylogeny phy = fixture();
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, conf, "radialexport" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                try {
                    final MainPanel mp = mf[ 0 ].getMainPanel();
                    final TreePanel tp = mp.getCurrentTreePanel();
                    mp.getControlPanel().setCheckbox( DisplayOption.SHOW_DOMAIN_ARCHITECTURES, true );
                    mp.getOptions().setNodeLabelDirection( Options.NODE_LABEL_DIRECTION.RADIAL );
                    // a FULL export is the one that re-centres on the export canvas; visible-only keeps the
                    // panel's own width, so it could never show this and would make the test vacuous
                    mp.getOptions().setGraphicsExportVisibleOnly( false );
                    tp.setAnnotationColumns( new java.util.ArrayList<AnnotationColumns.ColumnSpec>(
                            Arrays.asList(
                                    new AnnotationColumns.ColumnSpec( "data:host",
                                                                      AnnotationColumns.Type.COLOR_STRIP ),
                                    new AnnotationColumns.ColumnSpec( "data:load",
                                                                      AnnotationColumns.Type.MATRIX ) ) ) );
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.CIRCULAR );
                    if ( !declinesUntilRepainted( tp, phy, "circular" ) ) {
                        ok[ 0 ] = false;
                    }
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.UNROOTED );
                    if ( !declinesUntilRepainted( tp, phy, "unrooted" ) ) {
                        ok[ 0 ] = false;
                    }
                    // Deliberate NON-behaviour. The rectangular family is not re-centred on the export frame
                    // -- every other use of the canvas size in the paint is chrome -- so a plain export of one
                    // must not cost the user a frame of rollover. Recorded as a decision: widening the
                    // predicate to "any export" (a plausible response to the fixed-size path, which invalidates
                    // in layoutForExportSize instead) would silently take that frame from every tree.
                    tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
                    paint( tp, false );
                    PhylogenyNode rect_hit = null;
                    int rx = 0, ry = 0;
                    for( int y = 0; ( y < H ) && ( rect_hit == null ); y += 3 ) {
                        for( int x = 0; x < W; x += 3 ) {
                            final PhylogenyNode n = tp.findNode( x, y );
                            if ( n != null ) { rect_hit = n; rx = x; ry = y; break; }
                        }
                    }
                    if ( rect_hit == null ) {
                        System.out.println( "rectangular: nothing was hit before the export, so the "
                                + "non-invalidation check pins nothing" );
                        ok[ 0 ] = false;
                    }
                    else {
                        paint( tp, true );
                        if ( !tp.hitTestableLayout() ) {
                            System.out.println( "rectangular: a plain export must not decline hit-testing -- "
                                    + "it does not re-centre the layout" );
                            ok[ 0 ] = false;
                        }
                        else if ( tp.findNode( rx, ry ) != rect_hit ) {
                            System.out.println( "rectangular: the export changed what findNode answers at ("
                                    + rx + "," + ry + ")" );
                            ok[ 0 ] = false;
                        }
                    }
                    if ( !latchAndFixedSize( mf[ 0 ], tp, phy ) ) {
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
            System.out.println( "RadialExportGeometryTest: OK." );
        }
        else {
            System.out.println( "RadialExportGeometryTest: FAILED." );
        }
    }
}

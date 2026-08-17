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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GraphicsEnvironment;
import java.awt.Point;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.util.List;

import javax.swing.JScrollPane;
import javax.swing.JViewport;

import org.forester.archaeopteryx.TanglegramColoring.Field;
import org.forester.archaeopteryx.TanglegramLinker.LinkField;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * Headful render test for {@link TanglegramPanel}: renders a reversed 4-tip pair offscreen and asserts that both
 * trees + the tip-to-tip connectors are actually drawn (ink in the central connector band, which is clear of the
 * two trees and their labels), and that the link/crossing counts are correct. A green no-op when headless.
 */
public final class TanglegramPanelRenderTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "TanglegramPanelRender: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        // constructing / rotating a panel does not paint, so those checks are headless-safe; only renderOk / hitTestOk paint
        if ( !countsOk() || !labelFallbackOk() || !rotationOk() || !autoUntangleOk() || !colorModesOk()
                || !themeColorsOk() || !phylogramLayoutOk() ) {
            return false;
        }
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        return renderOk() && hitTestOk() && colorRenderOk() && panOk() && themeBackgroundRendersOk();
    }

    /** applyThemeColors sets the panel's background + ink (matching the tree canvas) and ignores null args. */
    private static boolean themeColorsOk() {
        final TanglegramPanel panel = new TanglegramPanel( balancedABCD(), balancedABCD(), LinkField.NODE_NAME );
        panel.applyThemeColors( Color.WHITE, Color.BLACK );
        if ( !Color.WHITE.equals( panel.getBackground() ) || !Color.BLACK.equals( panel.getForeground() ) ) {
            return fail( "applyThemeColors should set the panel background/foreground" );
        }
        // a null keeps the current colour (must not clear it)
        panel.applyThemeColors( null, null );
        if ( !Color.WHITE.equals( panel.getBackground() ) || !Color.BLACK.equals( panel.getForeground() ) ) {
            return fail( "applyThemeColors(null,null) must not change the colours" );
        }
        return true;
    }

    /** The applied theme background must actually FILL the opaque panel (so a light-mode tanglegram is white, not
     *  the default panel grey) -- render with a distinctive background and sample a margin pixel. */
    private static boolean themeBackgroundRendersOk() {
        final TanglegramPanel panel = new TanglegramPanel( balancedABCD(), balancedABCD(), LinkField.NODE_NAME );
        final Color bg = new Color( 201, 211, 221 ); // distinctive, not white/black/red so it can't be tree ink
        panel.applyThemeColors( bg, Color.BLACK );
        panel.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
        final int w = 400;
        final int h = 200;
        panel.setSize( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            panel.printAll( g );
        }
        finally {
            g.dispose();
        }
        final int corner = img.getRGB( 1, 1 ) & 0xFFFFFF;
        if ( corner != ( bg.getRGB() & 0xFFFFFF ) ) {
            return fail( "the panel background should fill the canvas; corner=#" + Integer.toHexString( corner )
                    + " expected #" + Integer.toHexString( bg.getRGB() & 0xFFFFFF ) );
        }
        return true;
    }

    private static boolean panOk() {
        final TanglegramPanel panel = new TanglegramPanel( balancedABCD(), balancedABCD(), LinkField.NODE_NAME );
        final JScrollPane scroll = new JScrollPane( panel );
        scroll.getViewport().setScrollMode( JViewport.SIMPLE_SCROLL_MODE );
        scroll.setSize( 400, 140 );
        scroll.doLayout();
        scroll.getViewport().doLayout();
        final JViewport vp = scroll.getViewport();
        if ( panel.getHeight() <= vp.getHeight() ) {
            return fail( "test setup: the panel is not taller than the viewport, so it cannot pan" );
        }
        // a left-drag upward past the threshold should scroll the content (viewport y increases)
        press( panel, 100, 100 );
        drag( panel, 100, 60 );
        if ( vp.getViewPosition().y <= 0 ) {
            return fail( "a left-drag did not pan the tanglegram (viewport y = " + vp.getViewPosition().y + ")" );
        }
        if ( !panel.pannedForTest() ) {
            return fail( "a real pan should mark the gesture panned (so the trailing click does not rotate)" );
        }
        // a sub-threshold jitter must NOT pan (so a click can still rotate)
        vp.setViewPosition( new Point( 0, 0 ) );
        press( panel, 100, 100 );
        drag( panel, 101, 102 );
        if ( vp.getViewPosition().y != 0 ) {
            return fail( "a sub-threshold drag should not pan the tanglegram" );
        }
        // a big drag on a tanglegram that FITS the viewport (nothing to scroll) must NOT suppress the click-to-rotate
        final TanglegramPanel fits = new TanglegramPanel( balancedABCD(), balancedABCD(), LinkField.NODE_NAME );
        final JScrollPane big = new JScrollPane( fits );
        big.setSize( 400, 500 );
        big.doLayout();
        big.getViewport().doLayout();
        press( fits, 100, 100 );
        drag( fits, 100, 40 );
        if ( fits.pannedForTest() ) {
            return fail( "a drag that cannot scroll (tanglegram fits the window) must not suppress the click" );
        }
        return true;
    }

    private static void press( final TanglegramPanel panel, final int x, final int y ) {
        final MouseEvent e = new MouseEvent( panel, MouseEvent.MOUSE_PRESSED, 0L, InputEvent.BUTTON1_DOWN_MASK, x, y, 1,
                                             false, MouseEvent.BUTTON1 );
        for( final MouseListener l : panel.getMouseListeners() ) {
            l.mousePressed( e );
        }
    }

    private static void drag( final TanglegramPanel panel, final int x, final int y ) {
        final MouseEvent e = new MouseEvent( panel, MouseEvent.MOUSE_DRAGGED, 0L, InputEvent.BUTTON1_DOWN_MASK, x, y, 0,
                                             false );
        for( final MouseMotionListener l : panel.getMouseMotionListeners() ) {
            l.mouseDragged( e );
        }
    }

    private static boolean colorModesOk() {
        final Color warn = new Color( 220, 50, 50 );
        // a fully tangled pair: all 4 connectors cross
        final TanglegramPanel tangled = new TanglegramPanel( balancedABCD(), balancedCDAB(), LinkField.NODE_NAME );
        final Color uniform0 = tangled.connectorColorForTest( 0 );
        for( int i = 1; i < 4; i++ ) {
            if ( !uniform0.equals( tangled.connectorColorForTest( i ) ) ) {
                return fail( "uniform mode should colour every connector the same" );
            }
        }
        tangled.setCrossingColoring();
        for( int i = 0; i < 4; i++ ) {
            if ( !warn.equals( tangled.connectorColorForTest( i ) ) ) {
                return fail( "in a fully tangled pair, every connector should be the crossing colour" );
            }
        }
        // a concordant pair: no connector crosses -> none is the crossing colour
        final TanglegramPanel clean = new TanglegramPanel( balancedABCD(), balancedABCD(), LinkField.NODE_NAME );
        clean.setCrossingColoring();
        for( int i = 0; i < 4; i++ ) {
            if ( warn.equals( clean.connectorColorForTest( i ) ) ) {
                return fail( "with no crossings, no connector should be coloured as crossing" );
            }
        }
        // FIELD (taxonomy): connectors of different species get different colours
        final TanglegramPanel taxo = new TanglegramPanel( taxonTree(), taxonTree(), LinkField.SCIENTIFIC_NAME );
        final Field sci = fieldNamed( taxo, "Taxonomy: Scientific Name" );
        if ( sci == null ) {
            return fail( "expected a scientific-name colour field" );
        }
        taxo.setFieldColoring( sci );
        if ( taxo.connectorColorForTest( 0 ).equals( taxo.connectorColorForTest( 1 ) ) ) {
            return fail( "different species should give different connector colours" );
        }
        return true;
    }

    private static boolean colorRenderOk() {
        // CROSSINGS mode must actually paint the crossing connectors red
        final TanglegramPanel tangled = new TanglegramPanel( balancedABCD(), balancedCDAB(), LinkField.NODE_NAME );
        tangled.setCrossingColoring();
        if ( countReddish( renderPanel( tangled, 800, 400 ) ) <= 0 ) {
            return fail( "crossing connectors did not render in red" );
        }
        // after untangling the crossings away, the crossing highlights must recompute (no stale red)
        tangled.autoUntangle();
        if ( countReddish( renderPanel( tangled, 800, 400 ) ) != 0 ) {
            return fail( "crossing highlights should clear once the crossings are untangled" );
        }
        // FIELD mode draws a colour legend (coloured swatches) in the top-left corner
        final TanglegramPanel taxo = new TanglegramPanel( taxonTree(), taxonTree(), LinkField.SCIENTIFIC_NAME );
        taxo.setFieldColoring( fieldNamed( taxo, "Taxonomy: Scientific Name" ) );
        if ( countSaturated( renderPanel( taxo, 800, 400 ), 0, 0, 160, 160 ) <= 0 ) {
            return fail( "the field-colour legend did not render coloured swatches" );
        }
        return true;
    }

    private static Field fieldNamed( final TanglegramPanel panel, final String label ) {
        for( final Field f : panel.availableColorFields() ) {
            if ( label.equals( f.label() ) ) {
                return f;
            }
        }
        return null;
    }

    private static BufferedImage renderPanel( final TanglegramPanel panel, final int w, final int h ) {
        panel.setBackground( Color.WHITE );
        panel.setForeground( Color.BLACK );
        panel.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
        panel.setSize( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            panel.printAll( g );
        }
        finally {
            g.dispose();
        }
        return img;
    }

    private static int countReddish( final BufferedImage img ) {
        int count = 0;
        for( int x = 0; x < img.getWidth(); x++ ) {
            for( int y = 0; y < img.getHeight(); y++ ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF;
                final int gg = ( rgb >> 8 ) & 0xFF;
                final int b = rgb & 0xFF;
                if ( ( r > 150 ) && ( gg < 110 ) && ( b < 110 ) ) {
                    ++count;
                }
            }
        }
        return count;
    }

    /** Coloured (saturated, not grey/black/white) pixels in the region [x0,x1) x [y0,y1). */
    private static int countSaturated( final BufferedImage img, final int x0, final int y0, final int x1,
                                       final int y1 ) {
        int count = 0;
        for( int x = Math.max( 0, x0 ); x < Math.min( img.getWidth(), x1 ); x++ ) {
            for( int y = Math.max( 0, y0 ); y < Math.min( img.getHeight(), y1 ); y++ ) {
                final int rgb = img.getRGB( x, y );
                final int r = ( rgb >> 16 ) & 0xFF;
                final int g = ( rgb >> 8 ) & 0xFF;
                final int b = rgb & 0xFF;
                if ( ( Math.max( r, Math.max( g, b ) ) - Math.min( r, Math.min( g, b ) ) ) > 60 ) {
                    ++count;
                }
            }
        }
        return count;
    }

    private static Phylogeny taxonTree() {
        return tree( clade( taxonLeaf( "Homo sapiens" ), taxonLeaf( "Mus musculus" ) ) );
    }

    private static PhylogenyNode taxonLeaf( final String scientific_name ) {
        final PhylogenyNode n = leaf( scientific_name );
        final Taxonomy t = new Taxonomy();
        t.setScientificName( scientific_name );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static boolean rotationOk() {
        // two IDENTICAL balanced trees ((A,B),(C,D)) -> 0 crossings, tips A,B,C,D
        final TanglegramPanel panel = new TanglegramPanel( balancedABCD(), balancedABCD(), LinkField.NODE_NAME );
        if ( ( panel.getCrossingCount() != 0 ) || !"A".equals( panel.leftTipLabelForTest( 0 ) ) ) {
            return fail( "identical trees should start at 0 crossings with first tip A, got " + panel.getCrossingCount()
                    + "/" + panel.leftTipLabelForTest( 0 ) );
        }
        if ( panel.canUndo() ) {
            return fail( "no undo should be available before any rotation" );
        }
        // flip the left root -> the two halves swap -> first left tip becomes C, crossings jump to 4
        panel.rotateLeftRootForTest();
        if ( ( panel.getCrossingCount() != 4 ) || !"C".equals( panel.leftTipLabelForTest( 0 ) ) ) {
            return fail( "after flipping the left root, expected 4 crossings and first tip C, got "
                    + panel.getCrossingCount() + "/" + panel.leftTipLabelForTest( 0 ) );
        }
        if ( !panel.canUndo() || panel.canRedo() ) {
            return fail( "after a rotation, undo should be available and redo not" );
        }
        panel.undo();
        if ( ( panel.getCrossingCount() != 0 ) || !"A".equals( panel.leftTipLabelForTest( 0 ) ) ) {
            return fail( "undo should restore 0 crossings and first tip A, got " + panel.getCrossingCount() + "/"
                    + panel.leftTipLabelForTest( 0 ) );
        }
        if ( panel.canUndo() || !panel.canRedo() ) {
            return fail( "after undo, redo should be available and undo not" );
        }
        panel.redo();
        if ( panel.getCrossingCount() != 4 ) {
            return fail( "redo should re-apply the flip -> 4 crossings, got " + panel.getCrossingCount() );
        }
        return true;
    }

    private static boolean autoUntangleOk() {
        // ((A,B),(C,D)) vs ((C,D),(A,B)) -> 4 crossings; auto-untangle reaches 0 and is ONE undoable action
        final TanglegramPanel panel = new TanglegramPanel( balancedABCD(), balancedCDAB(), LinkField.NODE_NAME );
        if ( panel.getCrossingCount() != 4 ) {
            return fail( "the tangled pair should start at 4 crossings, got " + panel.getCrossingCount() );
        }
        panel.autoUntangle();
        if ( panel.getCrossingCount() != 0 ) {
            return fail( "auto-untangle should reach 0 crossings, got " + panel.getCrossingCount() );
        }
        if ( !panel.canUndo() || panel.canRedo() ) {
            return fail( "after auto-untangle: undo available, redo not" );
        }
        panel.undo();
        if ( panel.getCrossingCount() != 4 ) {
            return fail( "undo should revert the WHOLE untangle in one step (back to 4 crossings), got "
                    + panel.getCrossingCount() );
        }
        panel.redo();
        if ( panel.getCrossingCount() != 0 ) {
            return fail( "redo should re-apply the untangle -> 0 crossings, got " + panel.getCrossingCount() );
        }
        return true;
    }

    private static boolean hitTestOk() {
        final TanglegramPanel panel = new TanglegramPanel( balancedABCD(), balancedABCD(), LinkField.NODE_NAME );
        panel.setBackground( Color.WHITE );
        panel.setForeground( Color.BLACK );
        panel.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
        panel.setSize( 800, 400 );
        final BufferedImage img = new BufferedImage( 800, 400, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            panel.printAll( g );
        }
        finally {
            g.dispose();
        }
        final int[] pt = panel.leftRootBarPointForTest();
        final PhylogenyNode hit = panel.rotatableNodeAtForTest( pt[ 0 ], pt[ 1 ] );
        if ( ( hit == null ) || hit.isExternal() ) {
            return fail( "clicking the left root's vertical bar should hit an internal (rotatable) node" );
        }
        if ( panel.rotatableNodeAtForTest( -50, -50 ) != null ) {
            return fail( "clicking well outside any bar should hit nothing" );
        }
        return true;
    }

    private static boolean labelFallbackOk() {
        // tips carry ONLY a taxonomy code (no node name, no scientific name) and are linked by that code -- the
        // drawn label must still resolve to the code, not a blank row (regression for the label-fallback fix)
        try {
            final TanglegramPanel panel = new TanglegramPanel( tree( clade( codeLeaf( "FELCA" ), codeLeaf( "CANLU" ) ) ),
                                                               tree( clade( codeLeaf( "CANLU" ), codeLeaf( "FELCA" ) ) ),
                                                               LinkField.TAXONOMY_CODE );
            if ( panel.getResult().getLinks().size() != 2 ) {
                return fail( "code-linked tips should form 2 links, got " + panel.getResult().getLinks().size() );
            }
            if ( !"FELCA".equals( panel.leftTipLabelForTest( 0 ) )
                    || !"CANLU".equals( panel.leftTipLabelForTest( 1 ) ) ) {
                return fail( "code-only tips must label with their taxonomy code, got '"
                        + panel.leftTipLabelForTest( 0 ) + "'/'" + panel.leftTipLabelForTest( 1 ) + "'" );
            }
            return true;
        }
        catch ( final Exception e ) {
            return fail( "code-leaf setup failed: " + e.getMessage() );
        }
    }

    private static boolean countsOk() {
        final TanglegramPanel panel = new TanglegramPanel( treeABCD(), treeDCBA(), LinkField.NODE_NAME );
        if ( panel.getResult().getLinks().size() != 4 ) {
            return fail( "expected 4 connectors, got " + panel.getResult().getLinks().size() );
        }
        if ( panel.getUnmatchedCount() != 0 ) {
            return fail( "expected 0 unmatched, got " + panel.getUnmatchedCount() );
        }
        if ( panel.getCrossingCount() != 6 ) {
            return fail( "a fully reversed 4-tip tanglegram should have 6 crossings, got " + panel.getCrossingCount() );
        }
        return true;
    }

    private static boolean renderOk() {
        final TanglegramPanel panel = new TanglegramPanel( treeABCD(), treeDCBA(), LinkField.NODE_NAME );
        final int w = 800;
        final int h = 400;
        panel.setBackground( Color.WHITE );
        panel.setForeground( Color.BLACK );
        panel.setFont( new Font( "SansSerif", Font.PLAIN, 12 ) );
        panel.setSize( w, h );
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        final Graphics2D g = img.createGraphics();
        try {
            panel.printAll( g );
        }
        finally {
            g.dispose();
        }
        final int total_ink = countNonWhite( img, 0, w );
        if ( total_ink < 200 ) {
            return fail( "the tanglegram rendered almost no ink (" + total_ink + " px) -- trees not drawn?" );
        }
        // a narrow central band lies BETWEEN the two trees + their labels, so only connectors cross it
        final int band_ink = countNonWhite( img, ( w / 2 ) - 5, ( w / 2 ) + 5 );
        if ( band_ink <= 0 ) {
            return fail( "no ink in the central connector band -- connectors not drawn" );
        }
        return true;
    }

    private static int countNonWhite( final BufferedImage img, final int x_from, final int x_to ) {
        int count = 0;
        for( int x = Math.max( 0, x_from ); x < Math.min( img.getWidth(), x_to ); x++ ) {
            for( int y = 0; y < img.getHeight(); y++ ) {
                if ( ( img.getRGB( x, y ) & 0xFFFFFF ) != 0xFFFFFF ) {
                    ++count;
                }
            }
        }
        return count;
    }

    /** The aligned-phylogram toggle draws branches to scale: a shallow tip ends left of a deeper tip (varying x),
     *  whereas the default cladogram lines every tip up at the common column (equal x). Headless-safe (renders into a
     *  scratch raster to populate coords). */
    private static boolean phylogramLayoutOk() {
        final TanglegramPanel panel = new TanglegramPanel( unevenDepthTree(), unevenDepthTree(), LinkField.NODE_NAME );
        if ( !panel.hasBranchLengths() ) {
            return fail( "the branch-length fixture should report hasBranchLengths()" );
        }
        final int w = 600;
        final int h = 300;
        // cladogram (default): tips 0 (shallow A) and 2 (deep C) share the common column -> equal x
        panel.renderForTest( w, h );
        if ( Math.abs( panel.leftTipXForTest( 0 ) - panel.leftTipXForTest( 2 ) ) > 0.5f ) {
            return fail( "cladogram tips should line up at the common column, got " + panel.leftTipXForTest( 0 )
                    + " vs " + panel.leftTipXForTest( 2 ) );
        }
        // phylogram: the shallow tip (distance 2 of 3) ends left of the deep tip (distance 3, flush at the column)
        panel.setPhylogram( true );
        panel.renderForTest( w, h );
        final float shallow = panel.leftTipXForTest( 0 );
        final float deep = panel.leftTipXForTest( 2 );
        if ( !( shallow < ( deep - 1f ) ) ) {
            return fail( "phylogram: the shallow tip x (" + shallow + ") should be left of the deep tip x (" + deep
                    + ")" );
        }
        return true;
    }

    private static boolean fail( final String message ) {
        System.out.println( "TanglegramPanelRender test failed: " + message );
        return false;
    }

    /** ((A:1,B:1):1, C:3): tips A and B sit at distance 2 from the root, C at distance 3 -- so a phylogram layout
     *  gives them different depths while a cladogram lines them all up. */
    private static Phylogeny unevenDepthTree() {
        final PhylogenyNode ab = cladeBL( 1.0, leafBL( "A", 1.0 ), leafBL( "B", 1.0 ) );
        return tree( cladeBL( 0.0, ab, leafBL( "C", 3.0 ) ) );
    }

    private static PhylogenyNode leafBL( final String name, final double bl ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        n.setDistanceToParent( bl );
        return n;
    }

    private static PhylogenyNode cladeBL( final double bl, final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setDistanceToParent( bl );
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static Phylogeny treeABCD() {
        return tree( clade( leaf( "A" ), leaf( "B" ), leaf( "C" ), leaf( "D" ) ) );
    }

    /** A balanced binary tree ((A,B),(C,D)) -- flipping its root cleanly swaps the two halves. */
    private static Phylogeny balancedABCD() {
        return tree( clade( clade( leaf( "A" ), leaf( "B" ) ), clade( leaf( "C" ), leaf( "D" ) ) ) );
    }

    /** ((C,D),(A,B)) -- the same unordered topology as balancedABCD but with the two halves swapped (4 crossings). */
    private static Phylogeny balancedCDAB() {
        return tree( clade( clade( leaf( "C" ), leaf( "D" ) ), clade( leaf( "A" ), leaf( "B" ) ) ) );
    }

    private static Phylogeny treeDCBA() {
        return tree( clade( leaf( "D" ), leaf( "C" ), leaf( "B" ), leaf( "A" ) ) );
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    /** A tip with ONLY a taxonomy code (no node name, no scientific name). */
    private static PhylogenyNode codeLeaf( final String code ) throws Exception {
        final PhylogenyNode n = new PhylogenyNode();
        final Taxonomy t = new Taxonomy();
        t.setTaxonomyCode( code );
        n.getNodeData().setTaxonomy( t );
        return n;
    }

    private static PhylogenyNode clade( final PhylogenyNode... children ) {
        final PhylogenyNode n = new PhylogenyNode();
        for( final PhylogenyNode child : children ) {
            n.addAsChild( child );
        }
        return n;
    }

    private static Phylogeny tree( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }
}

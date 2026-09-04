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
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.phylogeny.data.PropertiesList;

/**
 * Integration test for the in-legend controls on a real {@link MainFrame}/{@link TreePanel}: a categorical
 * color-by-property legend with more distinct values than the default cap must draw a clickable sort toggle
 * and a clickable "+N more" footer; clicking them flips the sort order and expands/collapses the legend.
 * Guarded to a no-op on a headless box (needs FlatLaf via {@code createInstance}).
 */
public final class LegendControlsToolTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LegendControlsTool: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; needs a display toolkit
        }
        try {
            final Configuration conf = new Configuration();
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait(
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { tree( 35 ) }, conf,
                                                                        "legend" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
                // AUTO-COLOR (JS parity): the tree opens already colored by its best well-covered candidate
                if ( !"data:kind".equals( tp.getColorByPropertyRef() ) ) {
                    fail( ok, "a new tree should open auto-colored by its best field, got "
                            + tp.getColorByPropertyRef() );
                }
                tp.setColorByPropertyRef( "data:region" ); // 35 distinct categorical values
                frame.showWhole();
                render( tp );

                // the two controls are drawn (35 values > the default 20 cap -> "+15 more")
                if ( tp.legendSortToggleBoundsForTest() == null ) {
                    fail( ok, "sort toggle was not drawn" );
                }
                if ( tp.legendMoreBoundsForTest() == null ) {
                    fail( ok, "'+N more' control was not drawn for 35 values over the cap" );
                }
                if ( !tp.isLegendSortByCount() ) {
                    fail( ok, "default sort should be by count" );
                }

                // a static export omits the interactive chips ([by count]/[show all]) but keeps the "+N more"
                // line, so it draws fewer pixels than the on-screen legend (yet still draws the legend)
                final int screen_px = legendPixels( tp, true );
                final int export_px = legendPixels( tp, false );
                if ( export_px <= 0 ) {
                    fail( ok, "export legend drew nothing" );
                }
                if ( export_px >= screen_px ) {
                    fail( ok, "export legend should drop the interactive chips (fewer pixels than on-screen): screen="
                            + screen_px + " export=" + export_px );
                }

                // clicking the sort toggle flips count <-> A-Z (re-render between clicks: the chip label, hence
                // its bounds, changes)
                clickCenter( tp, tp.legendSortToggleBoundsForTest() );
                if ( tp.isLegendSortByCount() ) {
                    fail( ok, "clicking the sort toggle should switch to A-Z" );
                }
                render( tp );
                clickCenter( tp, tp.legendSortToggleBoundsForTest() );
                if ( !tp.isLegendSortByCount() ) {
                    fail( ok, "clicking the sort toggle again should switch back to by-count" );
                }

                // clicking "+N more" expands to show all; then a "show fewer" collapses back to the default cap
                render( tp );
                clickCenter( tp, tp.legendMoreBoundsForTest() );
                if ( tp.legendMaxEntries() != Integer.MAX_VALUE ) {
                    fail( ok, "clicking '+N more' should show all entries" );
                }
                render( tp );
                if ( tp.legendMoreBoundsForTest() == null ) {
                    fail( ok, "an expanded legend must offer a 'show fewer' control" );
                }
                clickCenter( tp, tp.legendMoreBoundsForTest() );
                if ( tp.legendMaxEntries() != 20 ) {
                    fail( ok, "clicking 'show fewer' should collapse to the default cap (20), got "
                            + tp.legendMaxEntries() );
                }

                // switching to a DIFFERENT legend subject must reset the expand state (not leak "show all")
                clickCenter( tp, tp.legendMoreBoundsForTest() ); // expand the region legend again
                render( tp );
                if ( tp.legendMaxEntries() != Integer.MAX_VALUE ) {
                    fail( ok, "test setup: legend should be expanded before the subject switch" );
                }
                tp.setColorByPropertyRef( "data:kind" ); // a different property = a different legend subject
                render( tp );
                if ( tp.legendMaxEntries() != 20 ) {
                    fail( ok, "switching legend subject must reset the expand cap to 20, got " + tp.legendMaxEntries() );
                }

                // a non-null-but-EMPTY color scheme draws no legend, so its stale bounds must not stay clickable
                final Rectangle stale = tp.getPropertyLegendBounds();
                tp.setColorByPropertyRef( "no:such:field" ); // no tip has this -> empty scheme, no legend drawn
                render( tp );
                if ( stale != null ) {
                    final MouseEvent in_old = new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(),
                                                              0, stale.x + ( stale.width / 2 ),
                                                              stale.y + ( stale.height / 2 ), 1, false );
                    if ( tp.isOnPropertyLegend( in_old ) ) {
                        fail( ok, "an empty color scheme must not leave a clickable phantom legend" );
                    }
                }

                // ---- the "no value" row (JS parity): partial coverage adds one pinned-last dashed row ----
                // data:kind (full coverage) vs data:patchy (same 4 distinct values/labels, 15 tips without a
                // value): the ONLY structural difference is the extra "no value" row, so the legend must be
                // exactly one row taller -- and carry visibly more ink (the dashed circle + "no value (15)"),
                // in the EXPORT rendering too (the row is informative, like "+N more", so exports keep it).
                tp.setColorByPropertyRef( "data:kind" );
                render( tp );
                final int h_full = tp.getPropertyLegendBounds().height;
                final int export_full = legendPixels( tp, false );
                tp.setColorByPropertyRef( "data:patchy" );
                render( tp );
                final int h_missing = tp.getPropertyLegendBounds().height;
                final int delta = h_missing - h_full;
                if ( ( delta < 8 ) || ( delta > 40 ) ) {
                    fail( ok, "partial coverage must add exactly one 'no value' row: height " + h_full + " -> "
                            + h_missing );
                }
                // the ROW ITSELF must carry ink (dashed circle + text): scan the bottom row band of the
                // on-screen legend -- a sabotage that keeps the box tall but draws a BLANK row fails here
                final Rectangle plb = tp.getPropertyLegendBounds();
                final int band_ink = inkIn( tp, true, plb.x + 4, plb.y + plb.height - 22,
                                            plb.width - 8, 18 );
                if ( band_ink < 30 ) {
                    fail( ok, "the 'no value' row band is blank (ink " + band_ink + ")" );
                }
                // exports keep the row (informative, like "+N more"): the row's ~250px of ink dwarfs the
                // title/box-size noise between the two refs, so demand a wide margin
                if ( legendPixels( tp, false ) < ( export_full + 150 ) ) {
                    fail( ok, "the export legend must keep the 'no value' row's ink" );
                }
                // the row is INERT: sweeping every y inside the legend reaches the 4 value labels and
                // nothing else -- no y resolves to a clickable "no value" entry
                final Rectangle lb = tp.getPropertyLegendBounds();
                final java.util.Set<String> hit = new java.util.HashSet<>();
                for( int yy = lb.y; yy < ( lb.y + lb.height ); ++yy ) {
                    final String v = tp.legendValueAt( new MouseEvent( tp, MouseEvent.MOUSE_CLICKED,
                                                                       System.currentTimeMillis(), 0,
                                                                       lb.x + ( lb.width / 2 ), yy, 1, false ) );
                    if ( v != null ) {
                        hit.add( v );
                    }
                }
                if ( hit.size() != 4 ) {
                    fail( ok, "expected the 4 value rows clickable, got " + hit );
                }
                if ( hit.contains( "no value" ) ) {
                    fail( ok, "the 'no value' row must be inert, not a clickable value" );
                }
                // gradient twin: data:allnum (full) vs data:num (5 tips absent) -> one extra row there too
                tp.setColorByPropertyRef( "data:allnum" );
                render( tp );
                final int gh_full = tp.getPropertyLegendBounds().height;
                tp.setColorByPropertyRef( "data:num" );
                render( tp );
                final int gdelta = tp.getPropertyLegendBounds().height - gh_full;
                if ( ( gdelta < 8 ) || ( gdelta > 40 ) ) {
                    fail( ok, "a partially-covered GRADIENT legend must add the 'no value' row: delta " + gdelta );
                }

                ( (JFrame) frame ).dispose();
            } );
            return ok[ 0 ] && identityMemoryOk() && elementSlotUiOk() && modeChipOk() && autoColorOk();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /**
     * Value-color IDENTITY through the real TreePanel wiring (JS parity): colors survive a view change
     * (collapse) instead of re-spreading; a palette switch and Reset to Defaults forget them.
     */
    private static boolean identityMemoryOk() throws Exception {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        // clade A holds 6 of the 7 "k0" tips; 5 loose "k1" tips. Collapsing A flips the visible frequency
        // order (k0: 7 -> 1, k1 stays 5), which is exactly what made the legacy rebuild swap the colors.
        final PhylogenyNode clade_a = new PhylogenyNode();
        for( int i = 0; i < 6; ++i ) {
            clade_a.addAsChild( gLeaf( "a" + i, "k0" ) );
        }
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( clade_a );
        root.addAsChild( gLeaf( "loose", "k0" ) );
        for( int i = 0; i < 5; ++i ) {
            root.addAsChild( gLeaf( "b" + i, "k1" ) );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                                    "identity" ) );
        SwingUtilities.invokeAndWait( () -> {
            final MainFrame frame = mf[ 0 ];
            final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
            tp.setColorByPropertyRef( "data:g" );
            final Color c0 = tp.getPropertyColorScheme().getValueColors().get( "K0" );
            final Color c1 = tp.getPropertyColorScheme().getValueColors().get( "K1" );
            // collapse clade A -> the rebuild sees k1 as the most frequent visible value; the colors must
            // NOT re-spread (legacy handed k1 the old k0 color here)
            clade_a.setCollapse( true );
            frame.getMainPanel().getControlPanel().displayedPhylogenyMightHaveChanged( true );
            final Map<String, Color> after = tp.getPropertyColorScheme().getValueColors();
            if ( !c0.equals( after.get( "K0" ) ) || !c1.equals( after.get( "K1" ) ) ) {
                fail( ok, "a view change must not recolor surviving values: K0 " + c0 + " -> " + after.get( "K0" )
                        + ", k1 " + c1 + " -> " + after.get( "K1" ) );
            }
            // a palette switch invalidates the remembered identities (they belong to the old palette)
            tp.setColorPaletteName( "Colorblind-friendly" );
            final Color k1_cb = tp.getPropertyColorScheme().getValueColors().get( "K1" );
            if ( c1.equals( k1_cb ) ) {
                fail( ok, "switching the palette must re-assign from the NEW palette" );
            }
            // Reset to Defaults forgets the identities: re-choosing the ref assigns fresh, by the CURRENT
            // view's frequencies -- so k1 (now most frequent) gets the default palette's FIRST color, which
            // was k0's launch color
            tp.resetColorStateToDefaults();
            tp.setColorByPropertyRef( "data:g" );
            final Color k1_fresh = tp.getPropertyColorScheme().getValueColors().get( "K1" );
            if ( !c0.equals( k1_fresh ) ) {
                fail( ok, "after Reset the memory must be forgotten (fresh frequency assignment), k1 got "
                        + k1_fresh + " expected " + c0 );
            }
            ( (JFrame) frame ).dispose();
        } );
        return ok[ 0 ];
    }

    /**
     * Element slots through the real UI (JS parity): a tree with repeating taxonomy codes offers
     * "tax:code" in the Color-by dropdown, the renderer shows the JS display label, and selecting it
     * colors the tips (scheme built, verbatim values).
     */
    private static boolean elementSlotUiOk() throws Exception {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 6; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            final org.forester.phylogeny.data.Taxonomy t = new org.forester.phylogeny.data.Taxonomy();
            try {
                t.setTaxonomyCode( ( i % 2 == 0 ) ? "MOUSE" : "CHICK" );
            }
            catch ( final Exception e ) {
                throw new RuntimeException( e );
            }
            leaf.getNodeData().setTaxonomy( t );
            root.addAsChild( leaf );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                                    "slots" ) );
        SwingUtilities.invokeAndWait( () -> {
            final MainFrame frame = mf[ 0 ];
            final ControlPanel cp = frame.getMainPanel().getControlPanel();
            final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
            cp.populateColorByPropertyBox();
            final javax.swing.JComboBox<String> box = cp.colorByPropertyBoxForTest();
            boolean offered = false;
            for( int i = 0; i < box.getItemCount(); ++i ) {
                if ( PropertyColorScheme.TAX_CODE_REF.equals( box.getItemAt( i ) ) ) {
                    offered = true;
                }
            }
            if ( !offered ) {
                fail( ok, "the Color-by dropdown should offer the taxonomy-code slot" );
            }
            // the renderer shows the JS Color-menu label, not the raw ref
            final java.awt.Component c = box.getRenderer().getListCellRendererComponent( new javax.swing.JList<>(),
                    PropertyColorScheme.TAX_CODE_REF, 0, false, false );
            if ( !( c instanceof javax.swing.JLabel )
                    || !"Taxonomy Code".equals( ( (javax.swing.JLabel) c ).getText() ) ) {
                fail( ok, "the dropdown should display 'Taxonomy Code' for tax:code" );
            }
            // selecting it colors the tips: two verbatim groups, full coverage
            tp.setColorByPropertyRef( PropertyColorScheme.TAX_CODE_REF );
            final PropertyColorScheme s = tp.getPropertyColorScheme();
            if ( ( s == null ) || ( s.getValueColors().size() != 2 )
                    || !s.getValueColors().containsKey( "MOUSE" ) || ( s.missingCount() != 0 ) ) {
                fail( ok, "coloring by tax:code should build a 2-group verbatim scheme, got "
                        + ( ( s == null ) ? "null" : s.getValueColors().keySet().toString() ) );
            }
            ( (JFrame) frame ).dispose();
        } );
        return ok[ 0 ];
    }

    /**
     * The numeric three-band system through the real TreePanel + legend: a SMALL (5-distinct) numeric field
     * defaults to individual colors with a [gradient] chip; the chip (via the real legend-click handler)
     * flips it to a gradient with a [colors] chip; a LARGE field has no chip.
     */
    private static boolean modeChipOk() throws Exception {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 25; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:code", Integer.toString( i % 5 ), "", "xsd:decimal",
                                          AppliesTo.NODE ) ); // 5 distinct numbers -> SMALL band
            pl.addProperty( new Property( "data:wide", Integer.toString( i * 7 ), "", "xsd:decimal",
                                          AppliesTo.NODE ) ); // 25 distinct numbers -> LARGE band (> 20)
            leaf.getNodeData().setProperties( pl );
            root.addAsChild( leaf );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                                    "modechip" ) );
        SwingUtilities.invokeAndWait( () -> {
            final MainFrame frame = mf[ 0 ];
            final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
            tp.setColorByPropertyRef( "data:code" );
            if ( tp.getPropertyColorScheme().isGradient() ) {
                fail( ok, "a 5-distinct numeric field should DEFAULT to individual colors (SMALL band)" );
            }
            if ( tp.getPropertyColorScheme().getValueColors().size() != 5 ) {
                fail( ok, "expected 5 individual number colors, got "
                        + tp.getPropertyColorScheme().getValueColors().size() );
            }
            if ( !tp.isColorBySwitchable() ) {
                fail( ok, "a SMALL-band field should be mode-switchable" );
            }
            render( tp );
            final Rectangle chip = tp.legendModeToggleBoundsForTest();
            if ( chip == null ) {
                fail( ok, "the [gradient] chip should be drawn on a switchable colors-mode legend" );
                ( (JFrame) frame ).dispose();
                return;
            }
            // clicking the chip through the REAL legend-click handler flips to a gradient
            clickCenter2( tp, chip );
            if ( !tp.getPropertyColorScheme().isGradient() ) {
                fail( ok, "clicking [gradient] should flip the field to a gradient" );
            }
            render( tp );
            final Rectangle back = tp.legendModeToggleBoundsForTest();
            if ( back == null ) {
                fail( ok, "the gradient legend should carry the [colors] chip" );
            }
            else {
                clickCenter2( tp, back );
                if ( tp.getPropertyColorScheme().isGradient() ) {
                    fail( ok, "clicking [colors] should flip back to individual colors" );
                }
            }
            // a LARGE (30-distinct) numeric field: gradient, not switchable, no chip
            tp.setColorByPropertyRef( "data:wide" );
            if ( tp.isColorBySwitchable() ) {
                fail( ok, "a >20-distinct numeric field must not be switchable" );
            }
            render( tp );
            if ( tp.legendModeToggleBoundsForTest() != null ) {
                fail( ok, "no mode chip on a non-switchable field" );
            }
            ( (JFrame) frame ).dispose();
        } );
        return ok[ 0 ];
    }

    /** Auto-color drives the real load hook: on, with a ref already chosen, and off. */
    private static boolean autoColorOk() throws Exception {
        final boolean[] ok = { true };
        final MainFrame[] mf = new MainFrame[ 1 ];
        final Phylogeny phy = new Phylogeny();
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < 8; ++i ) {
            final PhylogenyNode leaf = gLeaf( "t" + i, "k" + ( i % 2 ) );
            // a second, LOWER-ranked field (7-vs-1 skew -> low entropy), so "keep an already-chosen ref"
            // is provable: choosing data:h and re-running the hook must NOT hand the tree back to data:g
            leaf.getNodeData().getProperties().addProperty( new Property( "data:h",
                    ( i == 0 ) ? "y" : "x", "", "xsd:string", AppliesTo.NODE ) );
            root.addAsChild( leaf );
        }
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        SwingUtilities.invokeAndWait(
                () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { phy }, new Configuration(),
                                                                    "autocolor" ) );
        SwingUtilities.invokeAndWait( () -> {
            final MainFrame frame = mf[ 0 ];
            final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
            // the load auto-colored by the one good field
            if ( !"data:g".equals( tp.getColorByPropertyRef() ) ) {
                fail( ok, "load should auto-color by data:g, got " + tp.getColorByPropertyRef() );
            }
            // re-running the hook must NOT clobber an existing choice -- even one that is NOT the best
            // candidate (data:h ranks below data:g)
            tp.setColorByPropertyRef( "data:h" );
            AptxUtil.lookAtSomeTreePropertiesForAptxControlSettings( phy,
                    frame.getMainPanel().getControlPanel() );
            if ( !"data:h".equals( tp.getColorByPropertyRef() ) ) {
                fail( ok, "the hook must keep an already-chosen ref, got " + tp.getColorByPropertyRef() );
            }
            // with the option OFF the hook applies nothing
            tp.setColorByPropertyRef( null );
            frame.getOptions().setAutoColorNewTrees( false );
            AptxUtil.lookAtSomeTreePropertiesForAptxControlSettings( phy,
                    frame.getMainPanel().getControlPanel() );
            if ( tp.getColorByPropertyRef() != null ) {
                fail( ok, "with auto-color OFF the hook must not color, got " + tp.getColorByPropertyRef() );
            }
            frame.getOptions().setAutoColorNewTrees( true );
            ( (JFrame) frame ).dispose();
        } );
        return ok[ 0 ];
    }

    private static void clickCenter2( final TreePanel tp, final Rectangle r ) {
        tp.handleLegendClick( new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), 0,
                                              r.x + ( r.width / 2 ), r.y + ( r.height / 2 ), 1, false ) );
    }

    private static PhylogenyNode gLeaf( final String name, final String group ) {
        final PhylogenyNode leaf = new PhylogenyNode();
        leaf.setName( name );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "data:g", group, "", "xsd:string", AppliesTo.NODE ) );
        leaf.getNodeData().setProperties( pl );
        return leaf;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [LegendControlsToolTest] " + msg );
        ok[ 0 ] = false;
    }

    private static void clickCenter( final TreePanel tp, final Rectangle r ) {
        if ( r == null ) {
            return;
        }
        final int cx = r.x + ( r.width / 2 );
        final int cy = r.y + ( r.height / 2 );
        tp.handleLegendClick( new MouseEvent( tp, MouseEvent.MOUSE_CLICKED, System.currentTimeMillis(), 0, cx, cy, 1,
                                              false ) );
    }

    /** Non-white pixel count inside a sub-rectangle of the legend rendering (screen or export mode). */
    private static int inkIn( final TreePanel tp, final boolean draggable, final int rx, final int ry,
                              final int rw, final int rh ) {
        final int w = 700, h = 500;
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.drawLegendForTest( g, new Rectangle( 0, 0, w, h ), draggable );
        g.dispose();
        int n = 0;
        for( int x = Math.max( 0, rx ); x < Math.min( w, rx + rw ); ++x ) {
            for( int y = Math.max( 0, ry ); y < Math.min( h, ry + rh ); ++y ) {
                if ( ( img.getRGB( x, y ) & 0x00FFFFFF ) != 0x00FFFFFF ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Non-white pixel count of the legend alone, drawn on-screen (draggable) or in static-export mode. */
    private static int legendPixels( final TreePanel tp, final boolean draggable ) {
        final int w = 700, h = 500;
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.drawLegendForTest( g, new Rectangle( 0, 0, w, h ), draggable );
        g.dispose();
        int n = 0;
        for( int x = 0; x < w; ++x ) {
            for( int y = 0; y < h; ++y ) {
                if ( ( img.getRGB( x, y ) & 0x00FFFFFF ) != 0x00FFFFFF ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** Renders offscreen so the legend (and its control bounds) are laid out for the next hit-test. */
    private static void render( final TreePanel tp ) {
        final int w = 700, h = 500;
        tp.setSize( w, h );
        tp.validate();
        tp.doLayout();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.printAll( g );
        g.dispose();
    }

    private static Phylogeny tree( final int n ) {
        final PhylogenyNode root = new PhylogenyNode();
        for( int i = 0; i < n; ++i ) {
            final PhylogenyNode leaf = new PhylogenyNode();
            leaf.setName( "t" + i );
            final PropertiesList pl = new PropertiesList();
            pl.addProperty( new Property( "data:region", String.format( "v%02d", i ), "", "xsd:string",
                                          AppliesTo.NODE ) );
            pl.addProperty( new Property( "data:kind", "k" + ( i % 4 ), "", "xsd:string", AppliesTo.NODE ) );
            if ( i < 20 ) { // partial coverage: 15 tips carry NO value -> the legend's "no value" row
                pl.addProperty( new Property( "data:patchy", "k" + ( i % 4 ), "", "xsd:string", AppliesTo.NODE ) );
            }
            if ( i < 30 ) { // numeric, partial coverage (5 absent) -> the gradient legend's "no value" row
                pl.addProperty( new Property( "data:num", Integer.toString( i ), "", "xsd:decimal",
                                              AppliesTo.NODE ) );
            }
            pl.addProperty( new Property( "data:allnum", Integer.toString( i ), "", "xsd:decimal",
                                          AppliesTo.NODE ) );
            leaf.getNodeData().setProperties( pl );
            root.addAsChild( leaf );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }
}

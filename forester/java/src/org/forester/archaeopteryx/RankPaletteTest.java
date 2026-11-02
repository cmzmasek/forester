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
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.TreePanel.CLADE_LABEL_ANGLE;
import org.forester.archaeopteryx.TreePanel.CLADE_VIS;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * The taxonomic-rank marks' own colours (Christian, 2026-11-01: "own colors") and a saved figure's Color by off
 * (Christian, 2026-11-01: "yes"):
 * <ul>
 * <li><b>The palette.</b> {@link AptxUtil#RANK_PALETTE} pinned; each of its first five at least 16 CIEDE2000 units
 * from each of the first five of both Color-by palettes; mid-light; dealt in name order, extended past ten. The
 * distance is checked against Sharma's published CIEDE2000 test pairs first.</li>
 * <li><b>On the coronavirus demo.</b> The genus boxes deal the rank palette and share no colour with the Color-by
 * (Host) legend -- until 0.11.198 both dealt Tableau 10 from its start. A name Color by also shows (Genus, and the
 * same names in capitals) takes Color by's colour in the boxes, whichever was set first, back to the rank palette
 * with Color by off, following a value recoloured and a palette switch, at each of two hue-banded levels; a colour
 * picked on the boxes' legend wins. The bands are rebuilt only when those colours change. Colorize Subtrees by rank
 * deals the rank palette too.</li>
 * <li><b>A saved figure with Color by off</b> opens uncoloured -- the load-time auto-colour used to pick a field --
 * while a figure that does not mention Color by is auto-coloured and one that names a field opens with it.</li>
 * </ul>
 * The demo parts are GUI tests, a no-op on a headless box.
 */
public final class RankPaletteTest {

    private static final String   DEMO   = "forester/demo/clade-boxes.xml";
    // the genera in name order: the order the palette is dealt in
    private static final String[] GENERA = { "Alphacoronavirus", "Alphaletovirus", "Betacoronavirus",
            "Deltacoronavirus", "Gammacoronavirus" };
    private static final int[]    PINNED = { 0x9467BD, 0x008695, 0xE7298A, 0x999933, 0x8C785D, 0xEE8866, 0x8DA0CB,
            0x68855C, 0xA06177, 0xA6761D };

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "Rank palette (own colours, shared names, a figure's Color by off): "
                + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            boolean ok = distance();
            ok &= palette();
            ok &= lookup();
            if ( GraphicsEnvironment.isHeadless() ) {
                return ok;
            }
            ok &= withDemo( figureWith( null ), RankPaletteTest::onDemo );
            ok &= figureOff();
            return ok;
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    // ---- the palette alone -------------------------------------------------------------------------------------

    /** The CIEDE2000 implementation below, on Sharma, Wu and Dalal's published test pairs (2005): 1, 2 and 7. */
    private static boolean distance() {
        boolean ok = true;
        ok &= ck( Math.abs( de00( new double[] { 50, 2.6772, -79.7751 }, new double[] { 50, 0, -82.7485 } ) - 2.0425 )
                < 1e-4, "CIEDE2000, Sharma pair 1: 2.0425" );
        ok &= ck( Math.abs( de00( new double[] { 50, 3.1571, -77.2803 }, new double[] { 50, 0, -82.7485 } ) - 2.8615 )
                < 1e-4, "CIEDE2000, Sharma pair 2: 2.8615" );
        ok &= ck( Math.abs( de00( new double[] { 50, 0, 0 }, new double[] { 50, -1, 2 } ) - 2.3669 ) < 1e-4,
                  "CIEDE2000, Sharma pair 7: 2.3669" );
        return ok;
    }

    private static boolean palette() {
        boolean ok = true;
        final Color[] rank = AptxUtil.RANK_PALETTE;
        ok &= ck( rank.length == PINNED.length, "the rank palette has " + PINNED.length + " colours" );
        for ( int i = 0; i < Math.min( rank.length, PINNED.length ); ++i ) {
            ok &= ck( ( rank[ i ].getRGB() & 0xFFFFFF ) == PINNED[ i ], "rank colour " + i + " is "
                    + Integer.toHexString( PINNED[ i ] ) + ", got " + Integer.toHexString( rank[ i ].getRGB() & 0xFFFFFF ) );
        }
        final Color[] tableau = PropertyColorScheme.paletteForTest( "Default" );
        final Color[] okabe = PropertyColorScheme.paletteForTest( "Colorblind-friendly" );
        ok &= ck( ( ( tableau[ 0 ].getRGB() & 0xFFFFFF ) == 0x4E79A7 ) && ( ( okabe[ 0 ].getRGB() & 0xFFFFFF ) == 0xE69F00 ),
                  "precondition: the Color-by palettes are Tableau 10 and Okabe-Ito" );
        double min = Double.MAX_VALUE;
        for ( int i = 0; i < 5; ++i ) {
            for ( int j = 0; j < 5; ++j ) {
                min = Math.min( min, Math.min( de00( lab( rank[ i ] ), lab( tableau[ j ] ) ),
                                               de00( lab( rank[ i ] ), lab( okabe[ j ] ) ) ) );
            }
        }
        ok &= ck( min >= 16, "each of the rank palette's first five is >= 16 CIEDE2000 units from each of the first "
                + "five of both Color-by palettes; nearest " + min );
        for ( final Color c : rank ) {
            final double l = lab( c )[ 0 ];
            ok &= ck( ( l >= 48.5 ) && ( l <= 67.5 ), "a rank colour is mid-light (CIELAB L 48.7-67.1) for both canvases: "
                    + Integer.toHexString( c.getRGB() & 0xFFFFFF ) + " has L " + l );
        }
        // dealt in name order, extended past ten (not repeated)
        final TreeSet<String> twelve = new TreeSet<String>();
        for ( int i = 0; i < 12; ++i ) {
            twelve.add( "taxon" + ( char ) ( 'a' + i ) );
        }
        final List<Color> dealt = new ArrayList<Color>( AptxUtil.assignRankColors( twelve ).values() );
        ok &= ck( dealt.size() == 12, "every taxon gets a colour" );
        for ( int i = 0; i < 10; ++i ) {
            ok &= ck( dealt.get( i ).equals( rank[ i ] ), "the " + i + "th taxon by name takes rank colour " + i );
        }
        ok &= ck( dealt.get( 10 ).equals( AptxUtil.qualitativeColor( rank, 10 ) ) && !dealt.get( 10 ).equals( rank[ 0 ] ),
                  "past ten the palette is extended, not repeated" );
        ok &= ck( AptxUtil.assignRankColors( new TreeSet<String>() ).isEmpty()
                && AptxUtil.assignRankColors( null ).isEmpty(), "nothing to colour, no colours" );
        return ok;
    }

    /** A name's colour: picked by hand, else Color by's for that exact spelling, else Color by's case aside -- the
     *  most frequent spelling's (Color by on an element keeps "Homo sapiens" and "homo sapiens" apart). */
    private static boolean lookup() {
        boolean ok = true;
        final Color a = new Color( 0xF28E2B ), b = new Color( 0xE15759 ), mine = new Color( 0x123456 );
        final Map<String, Color> color_by = new java.util.LinkedHashMap<String, Color>(); // most frequent first
        color_by.put( "Homo sapiens", a );
        color_by.put( "homo sapiens", b );
        final Map<String, Color> none = new HashMap<String, Color>();
        ok &= ck( b.equals( TreePanel.rankNameColorForTest( none, color_by, "homo sapiens" ) )
                && a.equals( TreePanel.rankNameColorForTest( none, color_by, "Homo sapiens" ) ),
                  "each spelling Color by shows keeps its own colour" );
        ok &= ck( a.equals( TreePanel.rankNameColorForTest( none, color_by, "HOMO SAPIENS" ) ),
                  "another spelling takes the most frequent spelling's colour" );
        // ...in either order: the less frequent spelling first in the map would win a plain case-aside lookup
        final Map<String, Color> reversed = new java.util.LinkedHashMap<String, Color>();
        reversed.put( "homo sapiens", b );
        reversed.put( "Homo sapiens", a );
        ok &= ck( a.equals( TreePanel.rankNameColorForTest( none, reversed, "Homo sapiens" ) )
                && b.equals( TreePanel.rankNameColorForTest( none, reversed, "HOMO SAPIENS" ) ),
                  "the exact spelling first, then case aside by frequency" );
        final Map<String, Color> user = new HashMap<String, Color>();
        user.put( "homo sapiens", mine );
        ok &= ck( mine.equals( TreePanel.rankNameColorForTest( user, color_by, "homo sapiens" ) )
                && a.equals( TreePanel.rankNameColorForTest( user, color_by, "Homo sapiens" ) ),
                  "a colour picked by hand wins, for its own spelling" );
        ok &= ck( TreePanel.rankNameColorForTest( none, color_by, "Pan troglodytes" ) == null,
                  "a name Color by does not show has no colour here (the rank palette's stays)" );
        return ok;
    }

    // ---- on the demo ---------------------------------------------------------------------------------------------

    private static boolean onDemo( final MainFrame mf ) {
        final TreePanel tp = mf.getCurrentTreePanel();
        boolean ok = true;
        tp.setCladeBoxLabels( false );
        tp.setColorByPropertyRef( "vipr:Host" );
        ok &= ck( tp.setCladeBands( "genus", CLADE_VIS.BOXES ) == 5, "precondition: 5 genus clades" );
        // own colours: the rank palette in name order, none of them in the Color-by legend
        ok &= ck( dealtFromPalette( tp ), "the boxes deal the rank palette in name order; got " + colours( tp ) );
        final Map<String, Color> hosts = legend( tp );
        ok &= ck( hosts.size() >= 5, "precondition: the host legend has at least as many colours as there are boxes ("
                + hosts.size() + ")" );
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            ok &= ck( !hosts.containsValue( b.getColor() ), b.getTaxon() + "'s box shares no colour with the Color-by "
                    + "legend " + hosts );
        }
        // a name Color by shows keeps Color by's colour -- Color by set after the boxes...
        tp.setColorByPropertyRef( "vipr:Genus" );
        ok &= ck( legend( tp ).keySet().containsAll( Arrays.asList( GENERA ) ), "precondition: Color by Genus shows "
                + "every genus; got " + legend( tp ).keySet() );
        ok &= ck( sameAsLegend( tp ), "Color by Genus after the boxes: each box takes its genus's Color-by colour; "
                + colours( tp ) + " vs " + legend( tp ) );
        ok &= ck( !dealtFromPalette( tp ), "precondition: Color by's colours are not the rank palette's" );
        // ...and before them
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        ok &= ck( sameAsLegend( tp ), "Color by Genus before the boxes: the same" );
        // Color by off, or a gradient: the rank palette again
        tp.setColorByPropertyRef( null );
        ok &= ck( dealtFromPalette( tp ), "Color by off: the boxes go back to the rank palette; got " + colours( tp ) );
        tp.setColorByPropertyRef( "vipr:Genus" );
        tp.setColorByPropertyRef( "vipr:Genome_kb" );
        ok &= ck( tp.getPropertyColorScheme().isGradient() && dealtFromPalette( tp ), "Color by a gradient: the "
                + "boxes keep the rank palette; got " + colours( tp ) );
        // a value recoloured on the Color-by legend: the box follows
        tp.setColorByPropertyRef( "vipr:Genus" );
        final String key = tp.getPropertyColorScheme().getValueKeys().get( "Betacoronavirus" );
        final Color picked = new Color( 0x123456 );
        tp.setColorOverride( key, picked );
        ok &= ck( band( tp, "Betacoronavirus" ).getColor().equals( picked ), "a genus recoloured on the Color-by "
                + "legend recolours its box" );
        tp.clearColorOverride( key );
        ok &= ck( sameAsLegend( tp ), "...and back" );
        // a palette switch: the boxes follow Color by's new colours
        final Color before = band( tp, "Alphacoronavirus" ).getColor();
        tp.setColorPaletteName( "Colorblind-friendly" );
        ok &= ck( !band( tp, "Alphacoronavirus" ).getColor().equals( before ) && sameAsLegend( tp ),
                  "a palette switch recolours the boxes with Color by" );
        tp.setColorPaletteName( PropertyColorScheme.DEFAULT_PALETTE_NAME );
        // a colour picked on the boxes' legend wins over Color by's
        final Color mine = new Color( 0x654321 );
        tp.setRankColorOverride( "genus", "Alphacoronavirus", mine );
        ok &= ck( band( tp, "Alphacoronavirus" ).getColor().equals( mine )
                && !legend( tp ).get( "Alphacoronavirus" ).equals( mine ), "a colour picked on the boxes' legend wins "
                        + "over Color by's" );
        tp.clearRankColorOverride( "genus", "Alphacoronavirus" );
        ok &= ck( sameAsLegend( tp ), "...and cleared, Color by's again" );
        // recoloured only when Color by's colours change, and then not rebuilt (the clades are what they were)
        final int rebuilds = tp.cladeBandRebuildsForTest();
        final int recolors = tp.cladeBandRecolorsForTest();
        tp.rebuildPropertyColorScheme();
        tp.rebuildPropertyColorScheme();
        ok &= ck( ( tp.cladeBandRebuildsForTest() == rebuilds ) && ( tp.cladeBandRecolorsForTest() == recolors ),
                  "the boxes are neither rebuilt nor recoloured when Color by's colours stay ("
                          + ( tp.cladeBandRebuildsForTest() - rebuilds ) + " rebuilds, "
                          + ( tp.cladeBandRecolorsForTest() - recolors ) + " recolours)" );
        tp.setColorOverride( key, picked );
        ok &= ck( ( tp.cladeBandRecolorsForTest() > recolors ) && ( tp.cladeBandRebuildsForTest() == rebuilds ),
                  "...and recoloured, not rebuilt, when they change (" + ( tp.cladeBandRebuildsForTest() - rebuilds )
                          + " rebuilds)" );
        tp.clearColorOverride( key );
        // case aside: the same names in capitals (Color by capitalises a value's first letter, nothing else)
        for ( final PhylogenyNode tip : tp.getPhylogeny().getExternalNodes() ) {
            final String g = tip.getNodeData().getProperties().getProperties( "vipr:Genus" ).get( 0 ).getValue();
            tip.getNodeData().getProperties().addProperty( new Property( "test:upper", g.toUpperCase(), "",
                                                                         "xsd:string", AppliesTo.NODE ) );
            tip.getNodeData().getProperties().addProperty( new Property( "test:subfamily", g.equals( "Alphaletovirus" )
                    ? "Letovirinae" : "Orthocoronavirinae", "", "xsd:string", AppliesTo.NODE ) );
        }
        tp.setColorByPropertyRef( "test:upper" );
        ok &= ck( legend( tp ).containsKey( "ALPHACORONAVIRUS" ) && !legend( tp ).containsKey( "Alphacoronavirus" ),
                  "precondition: the Color-by legend spells the genera in capitals; got " + legend( tp ).keySet() );
        boolean upper = true;
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            upper &= b.getColor().equals( legend( tp ).get( b.getTaxon().toUpperCase() ) );
        }
        ok &= ck( upper, "a name Color by spells in capitals keeps Color by's colour" );
        // two levels: hue-banded, then each genus takes its Color-by colour all the same
        final List<CladeLevel.Spec> two = Arrays.asList( new CladeLevel.Spec( "genus", CLADE_LABEL_ANGLE.VERTICAL ),
                                                         new CladeLevel.Spec( "subfamily", CLADE_LABEL_ANGLE.VERTICAL ) );
        tp.setColorByPropertyRef( "vipr:Host" );
        tp.setCladeLevels( two, CLADE_VIS.BARS, false, TaxonSource.TREE );
        ok &= ck( tp.getCladeLevelSpecs().size() == 2, "precondition: two levels" );
        ok &= ck( !dealtFromPalette( tp ), "precondition: two levels are hue-banded, not dealt from the palette; "
                + colours( tp ) );
        final Map<String, Color> host_banded = colours( tp );
        tp.setColorByPropertyRef( "vipr:Genus" );
        ok &= ck( sameAsLegend( tp ) && !colours( tp ).equals( host_banded ), "two levels: each genus takes its "
                + "Color-by colour; " + colours( tp ) + " vs " + legend( tp ) );
        // collapsing a clade changes the view, not a name's colour: its box keeps the Color-by colour its collapsed
        // wedge wears, though the Color-by legend (what is on screen) loses the row
        tp.setPhylogenyGraphicsType( Options.PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR );
        tp.setColorByPropertyRef( "vipr:Genus" );
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        final Color beta = band( tp, "Betacoronavirus" ).getColor();
        final PhylogenyNode beta_root = band( tp, "Betacoronavirus" ).getRoot();
        tp.collapse( beta_root );
        ok &= ck( beta_root.isCollapse() && !legend( tp ).containsKey( "Betacoronavirus" ) && sameAsLegendBut( tp,
                "Betacoronavirus" ), "precondition: collapsed, Betacoronavirus leaves the Color-by legend; "
                        + legend( tp ).keySet() );
        ok &= ck( band( tp, "Betacoronavirus" ).getColor().equals( beta ), "a collapsed clade's box keeps its Color-by "
                + "colour; got " + band( tp, "Betacoronavirus" ).getColor() + ", was " + beta );
        tp.collapse( beta_root ); // and open again
        ok &= ck( !beta_root.isCollapse() && sameAsLegend( tp ), "...and opened again, the same" );
        // Colorize Subtrees by rank deals the rank palette too
        tp.setCladeBands( "genus", CLADE_VIS.BOXES );
        final Map<String, PhylogenyNode> roots = new HashMap<String, PhylogenyNode>();
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            roots.put( b.getTaxon(), b.getRoot() );
        }
        tp.clearCladeBands();
        ok &= ck( tp.colorByRank( "genus", TaxonSource.TREE ) == 5, "precondition: Colorize Subtrees colours 5 clades" );
        for ( int i = 0; i < GENERA.length; ++i ) {
            final Color c = AptxUtil.RANK_PALETTE[ i ];
            ok &= ck( c.equals( tp.rankLegendColorsForTest().get( GENERA[ i ] ) ), "Colorize Subtrees: " + GENERA[ i ]
                    + " takes rank colour " + i + " in the legend; got " + tp.rankLegendColorsForTest() );
            ok &= ck( c.equals( roots.get( GENERA[ i ] ).getBranchData().getBranchColor().getValue() ),
                      "Colorize Subtrees: " + GENERA[ i ] + "'s branches take rank colour " + i );
        }
        return ok;
    }

    // ---- a saved figure with Color by off ------------------------------------------------------------------------

    private static boolean figureOff() throws Exception {
        boolean ok = true;
        // the flag, as each figure sets it
        final boolean[] flags = new boolean[ 3 ];
        final boolean[] seen = { false };
        ok &= withDemo( figureWith( null ), mf -> {
            final TreePanel tp = mf.getCurrentTreePanel();
            FigureSpec.parse( "v1;colorby=" ).applyTo( tp );
            flags[ 0 ] = tp.isColorByOffByFigure();
            FigureSpec.parse( "v1;colorby=vipr:Host" ).applyTo( tp );
            flags[ 1 ] = tp.isColorByOffByFigure();
            FigureSpec.parse( "v1;colorby=" ).applyTo( tp );
            FigureSpec.parse( "v1;graphics=RECTANGULAR" ).applyTo( tp );
            flags[ 2 ] = tp.isColorByOffByFigure();
            seen[ 0 ] = FigureSpec.capture( tp ).toPropertyValue().contains( ";colorby=;" )
                    || FigureSpec.capture( tp ).toPropertyValue().endsWith( ";colorby=" );
            return true;
        } );
        ok &= ck( flags[ 0 ] && !flags[ 1 ] && !flags[ 2 ], "the flag: set by colorby=\"\", not by a field nor by a "
                + "figure that does not mention Color by; got " + Arrays.toString( flags ) );
        ok &= ck( seen[ 0 ], "a figure captured with Color by off says colorby=\"\"" );
        // opened through the real start-up, auto-colour on
        final String[] opened = new String[ 3 ];
        final boolean[] auto = { false };
        withDemo( figureWith( "" ), mf -> {
            auto[ 0 ] = mf.getOptions().isAutoColorNewTrees();
            opened[ 0 ] = mf.getCurrentTreePanel().getColorByPropertyRef();
            return true;
        } );
        withDemo( figureWith( "-" ), mf -> {
            opened[ 1 ] = mf.getCurrentTreePanel().getColorByPropertyRef();
            return true;
        } );
        withDemo( figureWith( null ), mf -> {
            opened[ 2 ] = mf.getCurrentTreePanel().getColorByPropertyRef();
            return true;
        } );
        ok &= ck( auto[ 0 ], "precondition: auto-colour is on" );
        ok &= ck( opened[ 1 ] != null, "precondition: the twin whose figure does not mention Color by is auto-coloured" );
        ok &= ck( opened[ 0 ] == null, "a figure saved with Color by off opens uncoloured; got " + opened[ 0 ] );
        ok &= ck( "vipr:Host".equals( opened[ 2 ] ), "a figure that names a field opens with it; got " + opened[ 2 ] );
        return ok;
    }

    /** The demo tree with its figure's Color by {@code null}: as saved (vipr:Host); {@code ""}: off; {@code "-"}: not
     *  mentioned. */
    private static Phylogeny figureWith( final String color_by ) throws Exception {
        final Phylogeny phy = FigureRenderer.readTrees( new File( System.getProperty( "user.dir" ), DEMO ) )[ 0 ];
        if ( color_by != null ) {
            final String saved = FigureSpec.readFrom( phy ).toPropertyValue();
            if ( !saved.contains( ";colorby=vipr:Host" ) ) {
                throw new IllegalStateException( "the demo's figure no longer colours by host: " + saved );
            }
            FigureSpec.writeToTree( phy, FigureSpec.parse( saved.replace( ";colorby=vipr:Host",
                                                                          "-".equals( color_by ) ? "" : ";colorby=" ) ) );
        }
        return phy;
    }

    // ---- helpers -------------------------------------------------------------------------------------------------

    private interface Exercise {

        boolean run( MainFrame mf );
    }

    private static boolean withDemo( final Phylogeny phy, final Exercise ex ) throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "rank palette" ) );
        final boolean[] ok = { false };
        try {
            SwingUtilities.invokeAndWait( () -> ok[ 0 ] = ex.run( mf[ 0 ] ) );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( ( JFrame ) mf[ 0 ] ).dispose() );
        }
        return ok[ 0 ];
    }

    private static boolean dealtFromPalette( final TreePanel tp ) {
        for ( int i = 0; i < GENERA.length; ++i ) {
            final CladeBand b = find( tp, GENERA[ i ] );
            if ( ( b == null ) || !b.getColor().equals( AptxUtil.RANK_PALETTE[ i ] ) ) {
                return false;
            }
        }
        return true;
    }

    private static boolean sameAsLegend( final TreePanel tp ) {
        final Map<String, Color> legend = legend( tp );
        for ( final String g : GENERA ) {
            final CladeBand b = find( tp, g );
            if ( ( b == null ) || !b.getColor().equals( legend.get( g ) ) ) {
                return false;
            }
        }
        return true;
    }

    /** {@link #sameAsLegend} for every genus but {@code but}. */
    private static boolean sameAsLegendBut( final TreePanel tp, final String but ) {
        final Map<String, Color> legend = legend( tp );
        for ( final String g : GENERA ) {
            final CladeBand b = find( tp, g );
            if ( !g.equals( but ) && ( ( b == null ) || !b.getColor().equals( legend.get( g ) ) ) ) {
                return false;
            }
        }
        return true;
    }

    private static Map<String, Color> legend( final TreePanel tp ) {
        return ( tp.getPropertyColorScheme() == null ) ? new HashMap<String, Color>()
                : tp.getPropertyColorScheme().getValueColors();
    }

    private static Map<String, Color> colours( final TreePanel tp ) {
        final Map<String, Color> m = new java.util.TreeMap<String, Color>();
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            m.put( b.getTaxon(), b.getColor() );
        }
        return m;
    }

    private static CladeBand find( final TreePanel tp, final String taxon ) {
        for ( final CladeBand b : tp.cladeBandsForTest() ) {
            if ( taxon.equals( b.getTaxon() ) ) {
                return b;
            }
        }
        return null;
    }

    private static CladeBand band( final TreePanel tp, final String taxon ) {
        final CladeBand b = find( tp, taxon );
        if ( b == null ) {
            throw new IllegalStateException( "no " + taxon + " clade" );
        }
        return b;
    }

    /** CIELAB of an sRGB colour (D65). */
    private static double[] lab( final Color c ) {
        final double r = lin( c.getRed() / 255.0 ), g = lin( c.getGreen() / 255.0 ), b = lin( c.getBlue() / 255.0 );
        final double x = ( ( 0.4124 * r ) + ( 0.3576 * g ) + ( 0.1805 * b ) ) / 0.95047;
        final double y = ( 0.2126 * r ) + ( 0.7152 * g ) + ( 0.0722 * b );
        final double z = ( ( 0.0193 * r ) + ( 0.1192 * g ) + ( 0.9505 * b ) ) / 1.08883;
        final double fx = f( x ), fy = f( y ), fz = f( z );
        return new double[] { ( 116 * fy ) - 16, 500 * ( fx - fy ), 200 * ( fy - fz ) };
    }

    private static double lin( final double c ) {
        return ( c <= 0.04045 ) ? ( c / 12.92 ) : Math.pow( ( c + 0.055 ) / 1.055, 2.4 );
    }

    private static double f( final double t ) {
        return ( t > 0.008856 ) ? Math.cbrt( t ) : ( ( 7.787 * t ) + ( 16.0 / 116.0 ) );
    }

    /** CIEDE2000 colour difference of two CIELAB colours (Sharma, Wu and Dalal 2005). */
    private static double de00( final double[] p, final double[] q ) {
        final double c1 = Math.hypot( p[ 1 ], p[ 2 ] ), c2 = Math.hypot( q[ 1 ], q[ 2 ] );
        final double cb = ( c1 + c2 ) / 2;
        final double g = 0.5 * ( 1 - Math.sqrt( Math.pow( cb, 7 ) / ( Math.pow( cb, 7 ) + Math.pow( 25, 7 ) ) ) );
        final double a1 = ( 1 + g ) * p[ 1 ], a2 = ( 1 + g ) * q[ 1 ];
        final double c1p = Math.hypot( a1, p[ 2 ] ), c2p = Math.hypot( a2, q[ 2 ] );
        final double h1p = hue( a1, p[ 2 ] ), h2p = hue( a2, q[ 2 ] );
        final double dl = q[ 0 ] - p[ 0 ], dc = c2p - c1p;
        double dh = 0;
        if ( ( c1p * c2p ) != 0 ) {
            dh = h2p - h1p;
            if ( dh > 180 ) {
                dh -= 360;
            }
            else if ( dh < -180 ) {
                dh += 360;
            }
        }
        final double dhh = 2 * Math.sqrt( c1p * c2p ) * Math.sin( Math.toRadians( dh / 2 ) );
        final double lb = ( p[ 0 ] + q[ 0 ] ) / 2, cbp = ( c1p + c2p ) / 2;
        double hb = h1p + h2p;
        if ( ( c1p * c2p ) != 0 ) {
            if ( Math.abs( h1p - h2p ) <= 180 ) {
                hb = ( h1p + h2p ) / 2;
            }
            else if ( ( h1p + h2p ) < 360 ) {
                hb = ( h1p + h2p + 360 ) / 2;
            }
            else {
                hb = ( h1p + h2p - 360 ) / 2;
            }
        }
        final double t = ( 1 - ( 0.17 * Math.cos( Math.toRadians( hb - 30 ) ) ) )
                + ( 0.24 * Math.cos( Math.toRadians( 2 * hb ) ) ) + ( 0.32 * Math.cos( Math.toRadians( ( 3 * hb ) + 6 ) ) )
                - ( 0.20 * Math.cos( Math.toRadians( ( 4 * hb ) - 63 ) ) );
        final double dth = 30 * Math.exp( -Math.pow( ( hb - 275 ) / 25, 2 ) );
        final double rc = 2 * Math.sqrt( Math.pow( cbp, 7 ) / ( Math.pow( cbp, 7 ) + Math.pow( 25, 7 ) ) );
        final double sl = 1 + ( ( 0.015 * Math.pow( lb - 50, 2 ) ) / Math.sqrt( 20 + Math.pow( lb - 50, 2 ) ) );
        final double sc = 1 + ( 0.045 * cbp ), sh = 1 + ( 0.015 * cbp * t );
        final double rt = -Math.sin( Math.toRadians( 2 * dth ) ) * rc;
        return Math.sqrt( Math.pow( dl / sl, 2 ) + Math.pow( dc / sc, 2 ) + Math.pow( dhh / sh, 2 )
                + ( rt * ( dc / sc ) * ( dhh / sh ) ) );
    }

    private static double hue( final double a, final double b ) {
        final double h = Math.toDegrees( Math.atan2( b, a ) );
        return ( h < 0 ) ? ( h + 360 ) : h;
    }

    private static boolean ck( final boolean b, final String msg ) {
        if ( !b ) {
            System.out.println( "  [RankPaletteTest] " + msg );
        }
        return b;
    }
}

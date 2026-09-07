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

import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

import org.forester.archaeopteryx.NodeHoverText.Row;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.PropertiesList;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;

/**
 * Headless tests for the node hover card: the rows {@link NodeHoverText} derives from a node (content, order,
 * headings, what is left out), and the {@link NodeHoverCard} layout/paint -- measured size, wrapping and clipping,
 * placement that keeps the card inside the viewport, and an offscreen paint whose pixels show a card with
 * transparent corners at full alpha and nothing at all at alpha 0.
 */
public final class NodeHoverCardTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "NodeHoverCard: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        try {
            return rowsOfRichNode() && rowsEdgeCases() && dateText() && layoutAndWrap() && placement() && paint();
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static boolean rowsOfRichNode() throws Exception {
        final PhylogenyNode n = NodeDataDraftTest.richNode(); // internal, one of everything
        final PhylogenyNode parent = new PhylogenyNode();
        parent.addAsChild( n );
        final List<String> got = new ArrayList<>();
        for( final Row r : NodeHoverText.rows( n ) ) {
            got.add( r.toString() );
        }
        // Every row stands either in the leading "the node itself" block or under a heading -- see NodeHoverText:
        // a row printed after a headed section reads as part of it, so "Distribution" and "Tips below" belong up
        // here with the other node facts, and the properties get a heading of their own.
        final List<String> expected = List.of( "Name: BRCA1 clade",
                                               "Distance to parent: 0.0123",
                                               "Date: 6.5 [5.0 - 8.0] mya (split)",
                                               "Distribution: San Francisco",
                                               "Depth: 1",
                                               "Tips below: 2",
                                               "Confidence [bootstrap]: 95",
                                               "Confidence [probability]: 0.98",
                                               "stdev: 0.01",
                                               "[Taxonomy]",
                                               "Id [ncbi]: 9606",
                                               "Code: HUMAN",
                                               "Scientific name: Homo sapiens",
                                               "Common name: human",
                                               "Rank: species",
                                               "[Sequence]",
                                               "Accession [UniProt]: P38398",
                                               "Symbol: BRCA1",
                                               "Name: BRCA1_HUMAN",
                                               "Gene name: BRCA1",
                                               "Location: chr17",
                                               "Type: protein",
                                               "[Sequence]",
                                               "Name: BRCA1 mRNA",
                                               "Type: rna",
                                               "[Events]",
                                               "Duplications: 1",
                                               "Losses: 2",
                                               "[Properties]",
                                               "depth: 120 METRIC:m",
                                               "habitat: coastal" );
        return eq( "rows of the rich node", expected, got );
    }

    private static boolean rowsEdgeCases() {
        final PhylogenyNode bare = new PhylogenyNode();
        boolean ok = eq( "a bare unnamed node says nothing", 0, NodeHoverText.rows( bare ).size() );
        bare.setName( "x" );
        ok = ok && eq( "a named tip: name + depth", List.of( "Name: x", "Depth: 0" ), strings( NodeHoverText.rows( bare ) ) );
        // internal aptx:* and visual-style properties are never shown; a value-less property is skipped
        final PhylogenyNode p = new PhylogenyNode();
        p.setName( "p" );
        final PropertiesList pl = new PropertiesList();
        pl.addProperty( new Property( "aptx:reimport_profile", "secret", "", "xsd:string", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:year", "2019", "", "xsd:integer", AppliesTo.NODE ) );
        pl.addProperty( new Property( "data:empty", "", "", "xsd:string", AppliesTo.NODE ) );
        p.getNodeData().setProperties( pl );
        ok = ok && eq( "properties: headed, user-visible only, namespace stripped",
                       List.of( "Name: p", "Depth: 0", "[Properties]", "year: 2019" ),
                       strings( NodeHoverText.rows( p ) ) );
        // ... and a node whose only properties are hidden ones gets no heading either
        final PhylogenyNode hp = new PhylogenyNode();
        hp.setName( "hp" );
        final PropertiesList hidden = new PropertiesList();
        hidden.addProperty( new Property( "aptx:figure", "v1;x", "", "xsd:string", AppliesTo.NODE ) );
        hp.getNodeData().setProperties( hidden );
        ok = ok && eq( "hidden properties only -> no heading", List.of( "Name: hp", "Depth: 0" ),
                       strings( NodeHoverText.rows( hp ) ) );
        // a typed event shows its type; a counts-only event does not invent one
        final PhylogenyNode e = new PhylogenyNode();
        e.setName( "e" );
        e.addAsChild( new PhylogenyNode() );
        e.getNodeData().setEvent( new Event( 0, 0, 0, "transfer" ) );
        ok = ok && eq( "typed event", List.of( "Name: e", "Depth: 0", "Tips below: 1", "[Events]", "Type: transfer" ),
                       strings( NodeHoverText.rows( e ) ) );
        e.getNodeData().setEvent( new Event( 2, 0, 0 ) );
        ok = ok && eq( "counts-only event", List.of( "Name: e", "Depth: 0", "Tips below: 1", "[Events]", "Duplications: 2" ),
                       strings( NodeHoverText.rows( e ) ) );
        // an empty taxonomy object adds no heading
        final PhylogenyNode t = new PhylogenyNode();
        t.setName( "t" );
        t.getNodeData().addTaxonomy( new org.forester.phylogeny.data.Taxonomy() );
        ok = ok && eq( "empty taxonomy -> no heading", List.of( "Name: t", "Depth: 0" ), strings( NodeHoverText.rows( t ) ) );
        return ok;
    }

    private static boolean dateText() {
        return eq( "value+range+unit", "6.5 [5 - 8] mya",
                   NodeHoverText.dateText( new Date( "", new BigDecimal( "6.5" ), new BigDecimal( "5" ), new BigDecimal( "8" ), "mya" ) ) )
                && eq( "value+unit", "90 mya", NodeHoverText.dateText( new Date( "", new BigDecimal( "90" ), null, null, "mya" ) ) )
                && eq( "range only", "[5 - 8]", NodeHoverText.dateText( new Date( "", null, new BigDecimal( "5" ), new BigDecimal( "8" ), "" ) ) )
                && eq( "desc only", "split", NodeHoverText.dateText( new Date( "split" ) ) )
                && eq( "value+desc", "90 mya (split)", NodeHoverText.dateText( new Date( "split", new BigDecimal( "90" ), null, null, "mya" ) ) )
                && eq( "nothing", "", NodeHoverText.dateText( new Date() ) );
    }

    private static Graphics2D scratch() {
        return new BufferedImage( 10, 10, BufferedImage.TYPE_INT_ARGB ).createGraphics();
    }

    private static boolean layoutAndWrap() {
        final Graphics2D g = scratch();
        final Font base = new Font( Font.SANS_SERIF, Font.PLAIN, 13 );
        final java.util.function.Function<Font, FontMetrics> fm = g::getFontMetrics;
        final List<Row> rows = List.of( Row.line( "Name", "x" ), Row.heading( "Taxonomy" ), Row.line( "Code", "HUMAN" ) );
        final NodeHoverCard small = new NodeHoverCard( rows, base, fm, false );
        boolean ok = check( "a small card is small", ( small.getWidth() < 200 ) && ( small.getHeight() < 90 ) )
                && check( "positive size", ( small.getWidth() > 0 ) && ( small.getHeight() > 0 ) );
        final String long_value = "Miki Y et al. (1994) A strong candidate for the breast and ovarian cancer susceptibility gene "
                + "BRCA1 Science 266 66 71 and then some more words to force several wrapped lines here";
        final NodeHoverCard wide = new NodeHoverCard( List.of( Row.line( "Reference", long_value ) ), base, fm, false );
        ok = ok && check( "width capped at the maximum", wide.getWidth() <= NodeHoverCard.MAX_WIDTH_PX )
                && check( "a long value wraps to more height", wide.getHeight() > small.getHeight() );
        final FontMetrics m = g.getFontMetrics( base );
        final List<String> wrapped = NodeHoverCard.wrap( "one two three four five six seven", m, m.stringWidth( "one two three" ) + 2 );
        ok = ok && check( "wraps on words: " + wrapped, ( wrapped.size() >= 3 ) && wrapped.get( 0 ).equals( "one two three" ) );
        final String clipped = NodeHoverCard.clip( "abcdefghijklmnopqrstuvwxyz", m, m.stringWidth( "abcdefgh…" ) );
        ok = ok && check( "clip ends in an ellipsis and fits: " + clipped,
                          clipped.endsWith( "…" ) && ( m.stringWidth( clipped ) <= m.stringWidth( "abcdefgh…" ) ) );
        ok = ok && eq( "a fitting value is untouched", List.of( "short" ), NodeHoverCard.wrap( "short", m, 1000 ) );
        ok = ok && eq( "empty value -> one empty line", List.of( "" ), NodeHoverCard.wrap( "", m, 100 ) );
        // the same rows at a bigger base font give a bigger card (everything scales with the font)
        final NodeHoverCard big = new NodeHoverCard( rows, base.deriveFont( 26f ), fm, false );
        ok = ok && check( "scales with the font", ( big.getWidth() > small.getWidth() ) && ( big.getHeight() > small.getHeight() ) );
        g.dispose();
        return ok;
    }

    private static boolean placement() {
        final Graphics2D g = scratch();
        final NodeHoverCard card = new NodeHoverCard( List.of( Row.line( "Name", "x" ), Row.line( "Depth", "3" ) ),
                                                      new Font( Font.SANS_SERIF, Font.PLAIN, 13 ), g::getFontMetrics, false );
        final Rectangle vis = new Rectangle( 100, 50, 600, 400 );
        final int off = NodeHoverCard.POINTER_OFFSET;
        final Rectangle a = card.placement( 200, 100, vis );
        boolean ok = eq( "below-right of the pointer by default", new Rectangle( 200 + off, 100 + off, card.getWidth(), card.getHeight() ), a );
        final Rectangle b = card.placement( 690, 100, vis );
        ok = ok && check( "flips to the LEFT near the right edge", b.x + b.width <= vis.x + vis.width )
                && eq( "left flip lands left of the pointer", 690 - off - card.getWidth(), b.x );
        final Rectangle c = card.placement( 200, 440, vis );
        ok = ok && check( "flips ABOVE near the bottom edge", c.y + c.height <= vis.y + vis.height )
                && eq( "above flip lands above the pointer", 440 - off - card.getHeight(), c.y );
        final Rectangle d = card.placement( 101, 51, vis ); // top-left corner: nothing to flip to -> clamp inside
        ok = ok && check( "always inside the visible rect", vis.contains( d ) );
        final Rectangle e = card.placement( 5000, 5000, vis ); // way outside: still clamped inside
        ok = ok && check( "clamped inside even for an absurd pointer", vis.contains( e ) );
        g.dispose();
        return ok;
    }

    private static boolean paint() {
        final Graphics2D g0 = scratch();
        final NodeHoverCard card = new NodeHoverCard( List.of( Row.line( "Name", "BRCA1_HUMAN" ), Row.heading( "Taxonomy" ),
                                                               Row.line( "Code", "HUMAN" ) ),
                                                      new Font( Font.SANS_SERIF, Font.PLAIN, 13 ), g0::getFontMetrics, false );
        g0.dispose();
        final int W = card.getWidth() + 40;
        final int H = card.getHeight() + 40;
        final BufferedImage img = new BufferedImage( W, H, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        card.paint( g, 20, 20, 1f );
        g.dispose();
        final int centre = img.getRGB( 20 + ( card.getWidth() / 2 ), 20 + ( card.getHeight() / 2 ) );
        boolean ok = check( "card centre is painted (opaque-ish white)", ( ( centre >>> 24 ) > 200 ) && ( ( centre & 0xFF ) > 200 ) )
                && check( "the exact corner pixel is (nearly) transparent: rounded corner",
                          ( ( img.getRGB( 20, 20 ) >>> 24 ) < 120 ) );
        // text was drawn: some dark pixels inside the card
        int dark = 0;
        for( int y = 20; y < 20 + card.getHeight(); ++y ) {
            for( int x = 20; x < 20 + card.getWidth(); ++x ) {
                final int p = img.getRGB( x, y );
                if ( ( ( p >>> 24 ) > 200 ) && ( ( p & 0xFF ) < 120 ) && ( ( ( p >> 8 ) & 0xFF ) < 120 ) ) {
                    ++dark;
                }
            }
        }
        ok = ok && check( "text pixels present: " + dark, dark > 30 );
        // shadow: something faint below the card, outside it
        final int below = img.getRGB( 20 + ( card.getWidth() / 2 ), 20 + card.getHeight() + 2 );
        ok = ok && check( "a faint shadow below the card", ( ( below >>> 24 ) > 0 ) && ( ( below >>> 24 ) < 120 ) );
        // alpha 0 paints nothing at all
        final BufferedImage none = new BufferedImage( W, H, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g2 = none.createGraphics();
        card.paint( g2, 20, 20, 0f );
        g2.dispose();
        boolean any = false;
        for( int y = 0; y < H && !any; ++y ) {
            for( int x = 0; x < W; ++x ) {
                if ( ( none.getRGB( x, y ) >>> 24 ) != 0 ) {
                    any = true;
                    break;
                }
            }
        }
        ok = ok && check( "alpha 0 paints nothing (the first fade frame)", !any );
        // the dark palette paints a dark card
        final NodeHoverCard darkcard = new NodeHoverCard( List.of( Row.line( "Name", "x" ) ),
                                                          new Font( Font.SANS_SERIF, Font.PLAIN, 13 ), img.createGraphics()::getFontMetrics, true );
        final BufferedImage d = new BufferedImage( W, H, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g3 = d.createGraphics();
        darkcard.paint( g3, 20, 20, 1f );
        g3.dispose();
        final int dc = d.getRGB( 20 + darkcard.getWidth() - 4, 20 + darkcard.getHeight() - 4 );
        ok = ok && check( "dark theme: dark card", ( ( dc >>> 24 ) > 200 ) && ( ( dc & 0xFF ) < 100 ) );
        return ok;
    }

    private static List<String> strings( final List<Row> rows ) {
        final List<String> out = new ArrayList<>();
        for( final Row r : rows ) {
            out.add( r.toString() );
        }
        return out;
    }

    private static boolean eq( final String what, final Object expected, final Object actual ) {
        if ( ( expected == null ) ? ( actual == null ) : expected.equals( actual ) ) {
            return true;
        }
        System.out.println( "  [NodeHoverCardTest] " + what + ": expected <" + expected + "> but got <" + actual + ">" );
        return false;
    }

    private static boolean check( final String what, final boolean condition ) {
        if ( condition ) {
            return true;
        }
        System.out.println( "  [NodeHoverCardTest] " + what );
        return false;
    }

    private NodeHoverCardTest() {
        // not instantiable
    }
}

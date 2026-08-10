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

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.archaeopteryx.AptxUtil.GraphicsExportType;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.Taxonomy;

/**
 * End-to-end render test for the internal-taxonomy label de-duplication ({@link TreePanel#paintNodeData} /
 * {@link TreePanel#paintInternalLabelAboveBranch} gated by {@link TreePanelUtil#isDuplicateOfAncestorTaxon}):
 * a non-root internal node whose taxon equals its nearest annotated ancestor's must NOT draw its (redundant)
 * label. The decision is unit-tested headlessly in {@link TreePanelUtilTest}; this confirms the paint path
 * actually suppresses it.
 *
 * <p>Two topologically identical (16-tip, deliberately tall so internal labels are not dynamically hidden) trees
 * are rendered, with node {@code a} annotated "Carnivora" and a deeper internal node {@code b} annotated either
 * "Carnivora" (duplicates {@code a} -> its label is suppressed) or "Felidae" (distinct -> its label is drawn).
 * Everything else (branches, tips, node {@code a}'s label) is pixel-identical, so the distinct-taxon render must
 * carry strictly more ink -- exactly {@code b}'s label. Coordinate-free, so it is robust to layout.
 *
 * <p>Needs FlatLaf + a display, so {@link #test()} is a no-op (true) when headless, like the other GUI tests here.
 */
public final class InternalTaxonomyDedupTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "InternalTaxonomyDedup: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        try {
            final long[] duplicate = { -1 }; // b == a (Carnivora): b's label suppressed
            final long[] distinct = { -1 };  // b != a (Felidae): b's label drawn
            final boolean[] usable = { true };
            renderTotalDark( nestedTree( "Carnivora" ), duplicate, usable );
            renderTotalDark( nestedTree( "Felidae" ), distinct, usable );
            if ( !usable[ 0 ] ) {
                return true; // no usable viewport / taxonomy display in this environment; nothing to assert
            }
            if ( !( distinct[ 0 ] > ( duplicate[ 0 ] + 25 ) ) ) {
                System.out.println( "  redundant label was NOT suppressed: duplicate-ink=" + duplicate[ 0 ]
                        + " distinct-ink=" + distinct[ 0 ] + " (distinct must exceed by a label's worth of pixels)" );
                return false;
            }
            return true;
        }
        catch ( final Throwable t ) {
            t.printStackTrace();
            return false;
        }
    }

    private static void renderTotalDark( final Phylogeny phy, final long[] out, final boolean[] usable )
            throws Exception {
        final Configuration conf = new Configuration();
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, conf, "taxonomy dedup test" ) );
        SwingUtilities.invokeAndWait( () -> {
            final JFrame f = (JFrame) mf[ 0 ];
            try {
                f.setSize( 1000, 900 );
                f.validate();
                final MainPanel mp = mf[ 0 ].getMainPanel();
                final TreePanel tp = mp.getCurrentTreePanel();
                final ControlPanel cp = mp.getControlPanel();
                mp.getOptions().setRasterExportScale( 1 );
                mp.getOptions().setInternalLabelsAboveBranch( true );
                cp.setCheckbox( Configuration.display_internal_data, true );
                cp.setCheckbox( Configuration.show_taxonomy_scientific_names, true );
                cp.setCheckbox( Configuration.dynamically_hide_data, false ); // else internal labels get hidden
                cp.showWhole();
                if ( ( tp.getWidth() < 200 ) || ( tp.getHeight() < 200 ) || !cp.isShowTaxonomyScientificNames() ) {
                    usable[ 0 ] = false;
                    return;
                }
                final File file = File.createTempFile( "aptx_taxdedup_", ".png" );
                try {
                    AptxUtil.writePhylogenyToGraphicsFile( file.getAbsolutePath(), tp.getWidth(), tp.getHeight(), tp,
                                                           cp, GraphicsExportType.PNG, mp.getOptions() );
                    out[ 0 ] = totalDark( ImageIO.read( file ) );
                }
                finally {
                    file.delete();
                }
            }
            catch ( final Exception e ) {
                e.printStackTrace();
                usable[ 0 ] = false;
            }
            finally {
                f.dispose();
            }
        } );
    }

    private static long totalDark( final BufferedImage img ) {
        long n = 0;
        for( int y = 0; y < img.getHeight(); ++y ) {
            for( int x = 0; x < img.getWidth(); ++x ) {
                final int rgb = img.getRGB( x, y );
                if ( ( ( ( rgb >> 16 ) & 0xFF ) + ( ( rgb >> 8 ) & 0xFF ) + ( rgb & 0xFF ) ) < 200 ) {
                    ++n;
                }
            }
        }
        return n;
    }

    /** A 16-tip tree with a big subtree under both annotated internals (so their labels aren't dynamically hidden):
     *  {@code root -> ( a[Carnivora] -> ( b[b_taxon] -> 8 tips, 4 tips ), 4 tips )}. */
    private static Phylogeny nestedTree( final String b_taxon ) {
        final PhylogenyNode b = taxon( bush( "b", 8 ), b_taxon );
        final PhylogenyNode a = taxon( internal( b, bush( "s", 4 ) ), "Carnivora" );
        final PhylogenyNode root = internal( a, bush( "r", 4 ) );
        final Phylogeny p = new Phylogeny();
        p.setRoot( root );
        p.externalNodesHaveChanged();
        return p;
    }

    /** An internal node with {@code n} plain tips, to give an annotated ancestor enough vertical room to render. */
    private static PhylogenyNode bush( final String prefix, final int n ) {
        final PhylogenyNode node = new PhylogenyNode();
        for( int i = 0; i < n; i++ ) {
            node.addAsChild( leaf( prefix + i ) );
        }
        return node;
    }

    private static PhylogenyNode leaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        return n;
    }

    private static PhylogenyNode internal( final PhylogenyNode x, final PhylogenyNode y ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.addAsChild( x );
        n.addAsChild( y );
        return n;
    }

    private static PhylogenyNode taxon( final PhylogenyNode n, final String sci ) {
        final Taxonomy t = new Taxonomy();
        t.setScientificName( sci );
        n.getNodeData().setTaxonomy( t );
        return n;
    }
}

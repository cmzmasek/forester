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
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * Integration test for the horizontal re-fit performed by "Extract Data from Labels…" on a real
 * {@link MainFrame}/{@link TreePanel}. Regression guard for the bug where, on a domain-bearing tree, extraction
 * shortened the names and revealed the Seq Name/Taxonomy columns but never reset the panel's preferred width --
 * so the (rightward-shifted) domain architectures were clipped at the panel edge until the user manually zoomed.
 * The fix routes {@link MainFrame#extractLabelDataAndRefit(TreePanel)} through {@link ControlPanel#fitWidth()},
 * which recalcs the longest-ext-node info AND resets the preferred size. The test asserts the extraction fired,
 * the two columns were revealed, the preferred width is left FRESH (a second {@code resetPreferredSize()} is a
 * no-op -- it would change under the old, stale path), and the painted content fits within it. Guarded to a
 * no-op on a headless box (needs FlatLaf via {@code createInstance}).
 */
public final class LabelExtractionResizeTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "LabelExtractionResize: " + ( ok ? "OK." : "FAILED." ) );
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
                    () -> mf[ 0 ] = MainFrameApplication.createInstance( new Phylogeny[] { headerTreeWithDomains() },
                                                                        conf, "extract-resize" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                final MainFrame frame = mf[ 0 ];
                try {
                    exercise( frame, ok );
                }
                finally {
                    ( (JFrame) frame ).dispose();
                }
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    private static void exercise( final MainFrame frame, final boolean[] ok ) {
        frame.setSize( 1000, 600 );
        frame.validate();
        final TreePanel tp = frame.getMainPanel().getCurrentTreePanel();
        final ControlPanel cp = tp.getControlPanel();
        cp.setCheckbox( Configuration.show_domain_architectures, true );
        cp.showWhole();

        // preconditions: the two columns the extraction reveals are off to begin with
        if ( cp.isShowSeqNames() ) {
            fail( ok, "precondition: Seq Names should start hidden" );
        }

        final int n = frame.extractLabelDataAndRefit( tp );

        // all three header tips are parsed + changed
        if ( n != 3 ) {
            fail( ok, "expected 3 labels extracted, got " + n );
        }
        // the freshly-populated columns are revealed
        if ( !cp.isShowSeqNames() || !cp.isShowTaxonomyScientificNames() ) {
            fail( ok, "extraction must reveal Seq Name + Taxonomy Scientific" );
        }
        if ( !tp.isEdited() ) {
            fail( ok, "extraction must mark the tree edited" );
        }

        // THE REGRESSION: the preferred width must be left FRESH -- consistent with the current longest-ext-node
        // info -- not stale. A second resetPreferredSize() must therefore be a no-op. Under the old path
        // (displayedPhylogenyMightHaveChanged, which recalcs the longest info but never resets the size) the
        // width stayed at its pre-extraction value, so this second reset would change it.
        final int w_after_extract = tp.getPreferredSize().width;
        tp.resetPreferredSize();
        final int w_recomputed = tp.getPreferredSize().width;
        if ( w_after_extract != w_recomputed ) {
            fail( ok, "extraction left the panel width STALE: " + w_after_extract + " -> " + w_recomputed
                    + " (domains would be clipped)" );
        }

        // and the painted content -- domains included -- is reachable within that width (not clipped off the edge)
        final int rightmost = rightmostPaintedPixel( tp );
        if ( rightmost > tp.getPreferredSize().width ) {
            fail( ok, "content extends past the preferred width (" + rightmost + " > " + tp.getPreferredSize().width
                    + "): domains clipped" );
        }
    }

    // render onto a canvas wider than the preferred width and find the rightmost non-white column
    private static int rightmostPaintedPixel( final TreePanel tp ) {
        final int w = tp.getPreferredSize().width + 500;
        final int h = Math.max( 400, tp.getPreferredSize().height + 50 );
        tp.setSize( w, h );
        tp.validate();
        tp.doLayout();
        final BufferedImage img = new BufferedImage( w, h, BufferedImage.TYPE_INT_ARGB );
        final Graphics2D g = img.createGraphics();
        g.setColor( Color.WHITE );
        g.fillRect( 0, 0, w, h );
        tp.printAll( g );
        g.dispose();
        final int white = Color.WHITE.getRGB();
        for( int x = w - 1; x >= 0; --x ) {
            for( int y = 0; y < h; ++y ) {
                if ( img.getRGB( x, y ) != white ) {
                    return x;
                }
            }
        }
        return 0;
    }

    // three UniProt-header tips, each carrying a (wide) domain architecture drawn to the right of the label
    private static Phylogeny headerTreeWithDomains() {
        final PhylogenyNode root = new PhylogenyNode();
        root.addAsChild( domainLeaf(
                "tr|A0A8H5JZG0|A0A8H5JZG0_9HYPO Radical SAM protein OS=Fusarium phyllophilum OX=47803 GN=G1 PE=4 SV=1" ) );
        root.addAsChild( domainLeaf(
                "sp|P0AEX9|MALE_ECOLI Maltose-binding periplasmic protein OS=Escherichia coli OX=83333 GN=malE PE=1 SV=1" ) );
        root.addAsChild( domainLeaf(
                "tr|Q12345|Q12345_YEAST Some other protein OS=Saccharomyces cerevisiae OX=4932 GN=G3 PE=4 SV=1" ) );
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static PhylogenyNode domainLeaf( final String name ) {
        final PhylogenyNode n = new PhylogenyNode();
        n.setName( name );
        final Sequence seq = new Sequence();
        final List<PhylogenyData> domains = new ArrayList<>();
        domains.add( new ProteinDomain( "DomA", 20, 400, 0.0001 ) );
        domains.add( new ProteinDomain( "DomB", 600, 1800, 0.0001 ) );
        seq.setDomainArchitecture( new DomainArchitecture( domains, 2000 ) );
        n.getNodeData().addSequence( seq );
        return n;
    }

    private static void fail( final boolean[] ok, final String msg ) {
        System.out.println( "  [LabelExtractionResizeTest] " + msg );
        ok[ 0 ] = false;
    }

    private LabelExtractionResizeTest() {
    }
}

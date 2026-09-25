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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;

/**
 * The Help &gt; Control Panel Cheat Sheet. Its whole point is that it is GENERATED from the live control panel --
 * each row's icon is the button's own {@code Icon} and each row's text is the description the panel already
 * attaches to that button -- so it cannot drift from what is on screen. These checks are about that property, not
 * about the wording of any one row.
 */
public final class ControlPanelCheatSheetTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "ControlPanelCheatSheet: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true;
        }
        final boolean[] ok = { true };
        try {
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { plainTree() }, new Configuration(), "cheat" ) );
            final ControlPanel cp = mf[ 0 ].getMainPanel().getControlPanel();
            final List<ControlPanelCheatSheet.Entry> es = new ArrayList<>();
            SwingUtilities.invokeAndWait( () -> es.addAll( ControlPanelCheatSheet.entries( cp ) ) );

            // (1) it describes a real panel, not a handful of rows
            if ( es.size() < 20 ) {
                fail( ok, "the sheet must cover the panel; got only " + es.size() + " rows" );
            }
            // (2) every row says something, and the icon-bearing rows really carry icons
            int with_icon = 0;
            for( final ControlPanelCheatSheet.Entry e : es ) {
                // A row must SAY something: a tooltip, or at least the control's own name (several Display Data
                // checkboxes have no tooltip because their label is the description).
                if ( ( ( e.description() == null ) || e.description().isBlank() )
                        && ( ( e.label() == null ) || e.label().isBlank() ) ) {
                    fail( ok, "a row with neither a name nor a description got in: " + e );
                    break;
                }
                if ( e.icon() != null ) {
                    ++with_icon;
                }
            }
            if ( with_icon < 8 ) {
                fail( ok, "the drawn-icon buttons must appear WITH their icons; only " + with_icon + " did" );
            }
            // (3) NO DUPLICATES. A JComboBox is itself a Container whose inner parts inherit its tooltip, so
            // walking into one listed every dropdown twice. Each (name, description) pair must appear once.
            final Set<String> seen = new HashSet<>();
            for( final ControlPanelCheatSheet.Entry e : es ) {
                final String key = e.label() + "\u0000" + e.description();
                if ( !seen.add( key ) ) {
                    fail( ok, "duplicate row (the walk descended into a control): " + e.label() + " / "
                            + e.description() );
                    break;
                }
            }
            // (4) a control with no text of its own is NAMED from the panel's own label, not left as "(dropdown)"
            boolean named = false;
            for( final ControlPanelCheatSheet.Entry e : es ) {
                if ( e.label().startsWith( "Tree share" ) || e.label().startsWith( "Color by" ) ) {
                    named = true;
                }
                if ( e.label().contains( "dropdown" ) || e.label().contains( "slider" ) ) {
                    fail( ok, "a control was left generically named: " + e.label() );
                    break;
                }
            }
            if ( !named ) {
                fail( ok, "the sliders/dropdowns must be named from the panel's labels (Tree share, Color by)" );
            }
            // (5) it renders to a real image rather than a blank one
            final BufferedImage[] img = new BufferedImage[ 1 ];
            SwingUtilities.invokeAndWait( () -> img[ 0 ] = ControlPanelCheatSheet
                    .render( ControlPanelCheatSheet.buildSheet( cp ) ) );
            if ( ( img[ 0 ] == null ) || ( img[ 0 ].getWidth() < 200 ) || ( img[ 0 ].getHeight() < 200 ) ) {
                fail( ok, "the exported sheet must be a real image, got " + img[ 0 ] );
            }
            else {
                int ink = 0;
                for( int x = 0; x < img[ 0 ].getWidth(); x += 2 ) {
                    for( int y = 0; y < img[ 0 ].getHeight(); y += 2 ) {
                        if ( ( img[ 0 ].getRGB( x, y ) & 0xffffff ) != 0xffffff ) {
                            ++ink;
                        }
                    }
                }
                if ( ink < 500 ) {
                    fail( ok, "the exported sheet is blank (" + ink + " non-white samples)" );
                }
            }
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf[ 0 ] ).dispose() );

            // (5b) EVERY control on the panel must carry a tooltip. This is the rule that keeps the cheat sheet
            // useful: the sheet is generated, so a new button needs no edit here -- but a button with no tooltip
            // is listed by bare name and explains nothing, and one with neither name nor tooltip (an icon-only
            // button) does not appear at all. Enforcing it here means the coupling cannot be forgotten, which a
            // note in the docs alone would not achieve. Driven on a deliberately RICH tree, so the data-gated
            // Display Data rows are all visible and all checked.
            final MainFrame[] mfr = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mfr[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { richTree() }, new Configuration(), "cheat-tips" ) );
            final List<String> untipped = new ArrayList<>();
            SwingUtilities.invokeAndWait(
                    () -> collectUntipped( mfr[ 0 ].getMainPanel().getControlPanel(), untipped ) );
            if ( !untipped.isEmpty() ) {
                fail( ok, untipped.size() + " control-panel control(s) carry no tooltip, so the cheat sheet cannot "
                        + "describe them: " + untipped );
            }
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mfr[ 0 ] ).dispose() );

            // (6) THE point: the sheet follows the panel. A data-gated row is on the sheet only when the tree
            // carries that data -- so the sheet describes what is on screen, not a fixed list written by hand.
            final boolean plain_has_domains = mentions( es, "Domain Architectures" );
            final MainFrame[] mf2 = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf2[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { domainTree() }, new Configuration(), "cheat2" ) );
            final ControlPanel cp2 = mf2[ 0 ].getMainPanel().getControlPanel();
            final List<ControlPanelCheatSheet.Entry> es2 = new ArrayList<>();
            SwingUtilities.invokeAndWait( () -> es2.addAll( ControlPanelCheatSheet.entries( cp2 ) ) );
            final boolean domain_has_domains = mentions( es2, "Domain Architectures" );
            if ( plain_has_domains ) {
                fail( ok, "a tree with no domains must not get a Domain Architectures row" );
            }
            if ( !domain_has_domains ) {
                fail( ok, "a tree WITH domains must get one -- the sheet is supposed to follow the panel" );
            }
            SwingUtilities.invokeAndWait( () -> ( (javax.swing.JFrame) mf2[ 0 ] ).dispose() );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            ok[ 0 ] = false;
        }
        return ok[ 0 ];
    }

    /** Every visible button / dropdown / slider on the panel, that has no tooltip. */
    private static void collectUntipped( final java.awt.Container c, final List<String> out ) {
        for( final java.awt.Component comp : c.getComponents() ) {
            if ( !comp.isVisible() ) {
                continue;
            }
            if ( ( comp instanceof javax.swing.AbstractButton ) || ( comp instanceof javax.swing.JComboBox )
                    || ( comp instanceof javax.swing.JSlider ) ) {
                final String tip = ( (javax.swing.JComponent) comp ).getToolTipText();
                if ( ( tip == null ) || tip.isBlank() ) {
                    String name = ( comp instanceof javax.swing.AbstractButton )
                            ? ( (javax.swing.AbstractButton) comp ).getText() : null;
                    if ( ( name == null ) || name.isBlank() ) {
                        name = comp.getClass().getSimpleName();
                    }
                    out.add( name );
                }
                continue; // a control is a leaf: its inner parts inherit its tooltip
            }
            if ( comp instanceof java.awt.Container ) {
                collectUntipped( (java.awt.Container) comp, out );
            }
        }
    }

    /** A tree carrying every kind of data, so each data-gated Display Data row is visible and gets checked. */
    private static Phylogeny richTree() {
        try {
            final PhylogenyNode root = new PhylogenyNode();
            final PhylogenyNode inner = new PhylogenyNode();
            inner.setDistanceToParent( 0.1 );
            inner.getBranchData().addConfidence( new org.forester.phylogeny.data.Confidence( 95, "bootstrap" ) );
            inner.getBranchData().setBranchWidth( new org.forester.phylogeny.data.BranchWidth( 3 ) );
            inner.getNodeData().setEvent( org.forester.phylogeny.data.Event.createSingleDuplicationEvent() );
            for( int i = 0; i < 3; ++i ) {
                final PhylogenyNode t = new PhylogenyNode();
                t.setName( "a_very_long_tip_name_for_shortening_" + i );
                t.setDistanceToParent( 0.2 + ( i * 0.05 ) );
                final org.forester.phylogeny.data.Taxonomy tax = new org.forester.phylogeny.data.Taxonomy();
                tax.setTaxonomyCode( "HUMAN" );
                tax.setScientificName( "Homo sapiens" );
                tax.setCommonName( "human" );
                tax.setRank( "species" );
                t.getNodeData().setTaxonomy( tax );
                final Sequence s = new Sequence();
                s.setName( "seq" + i );
                s.setGeneName( "gene" + i );
                s.setSymbol( "SYM" + i );
                s.setAccession( new org.forester.phylogeny.data.Accession( "P0000" + i, "uniprot" ) );
                s.setMolecularSequence( "MKQLEDPFGH-WYVAST" );
                s.setMolecularSequenceAligned( true );
                final List<PhylogenyData> ds = new ArrayList<>();
                ds.add( new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 ) );
                s.setDomainArchitecture( new DomainArchitecture( ds, 200 ) );
                t.getNodeData().addSequence( s );
                final org.forester.phylogeny.data.PropertiesList pl =
                        new org.forester.phylogeny.data.PropertiesList();
                pl.addProperty( new org.forester.phylogeny.data.Property( "demo:grp", "g" + i, "", "xsd:string",
                        org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
                t.getNodeData().setProperties( pl );
                final org.forester.phylogeny.data.NodeVisualData vis =
                        new org.forester.phylogeny.data.NodeVisualData();
                vis.setNodeColor( new java.awt.Color( 10 + ( i * 40 ), 90, 200 ) );
                t.getNodeData().setNodeVisualData( vis );
                inner.addAsChild( t );
            }
            root.addAsChild( inner );
            final PhylogenyNode out = new PhylogenyNode();
            out.setName( "outgroup_with_a_long_name_too" );
            out.setDistanceToParent( 0.5 );
            root.addAsChild( out );
            return wrap( root );
        }
        catch ( final Exception e ) {
            throw new RuntimeException( e );
        }
    }

    private static boolean mentions( final List<ControlPanelCheatSheet.Entry> es, final String label ) {
        for( final ControlPanelCheatSheet.Entry e : es ) {
            if ( label.equals( e.label() ) ) {
                return true;
            }
        }
        return false;
    }

    private static Phylogeny plainTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( final String n : new String[] { "a", "b", "c" } ) {
            final PhylogenyNode t = new PhylogenyNode();
            t.setName( n );
            t.setDistanceToParent( 0.2 );
            root.addAsChild( t );
        }
        return wrap( root );
    }

    private static Phylogeny domainTree() {
        final PhylogenyNode root = new PhylogenyNode();
        for( final String n : new String[] { "a", "b", "c" } ) {
            final PhylogenyNode t = new PhylogenyNode();
            t.setName( n );
            t.setDistanceToParent( 0.2 );
            final Sequence seq = new Sequence();
            final List<PhylogenyData> ds = new ArrayList<>();
            ds.add( new ProteinDomain( "SH3", 10, 60, "PF00018", 1e-6 ) );
            seq.setDomainArchitecture( new DomainArchitecture( ds, 200 ) );
            t.getNodeData().addSequence( seq );
            root.addAsChild( t );
        }
        return wrap( root );
    }

    private static Phylogeny wrap( final PhylogenyNode root ) {
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static void fail( final boolean[] ok, final String message ) {
        System.out.println( "  [ControlPanelCheatSheetTest] " + message );
        ok[ 0 ] = false;
    }

    private ControlPanelCheatSheetTest() {
        // not instantiable
    }
}

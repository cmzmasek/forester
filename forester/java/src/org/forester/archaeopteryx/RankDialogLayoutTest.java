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

import java.awt.Component;
import java.awt.Container;
import java.awt.GraphicsEnvironment;
import java.awt.Window;
import java.io.File;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.Timer;

import org.forester.phylogeny.Phylogeny;

/**
 * Tools > Annotate Clades by Rank and Tools > Colorize Subtrees via Taxonomic Rank, opened the way the menu opens them
 * (the real modal dialog): every row of the form keeps its own height ({@link MainFrame#stackedRows}) and the dialog
 * fits a laptop screen, also after the source is switched to an online one (longer entries: the dialog re-packs).
 * From 0.11.179 to 0.11.196 a GridLayout gave every row the height of the three-line "Taxa from" panel: Annotate
 * Clades by Rank was 608 x 1443 px.
 */
public final class RankDialogLayoutTest {

    public static void main( final String[] args ) {
        final boolean ok = test();
        System.out.println( "RankDialogLayout: " + ( ok ? "OK." : "FAILED." ) );
        System.exit( ok ? 0 : 1 );
    }

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // a dialog needs a display
        }
        try {
            final File f = new File( System.getProperty( "user.dir" ), "forester/demo/colorize-by-rank.xml" );
            return dialogOk( FigureRenderer.readTrees( f )[ 0 ], "Annotate Clades by Rank", 23, false )
                    & dialogOk( FigureRenderer.readTrees( f )[ 0 ], "Colorize Subtrees via Taxonomic Rank", 7, false )
                    & dialogOk( longSchemeTree(), "Annotate Clades by Rank", 23, true );
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    /** Tips whose only rank is a scheme with a long name: its entry under the tree's own source fits the dialog, its
     *  online entry ("... now; 0 to look up") is wider. */
    private static Phylogeny longSchemeTree() {
        final org.forester.phylogeny.PhylogenyNode root = new org.forester.phylogeny.PhylogenyNode();
        for( int i = 0; i < 6; ++i ) {
            final org.forester.phylogeny.PhylogenyNode n = new org.forester.phylogeny.PhylogenyNode();
            n.setName( "tip_" + i );
            n.getNodeData().setProperties( new org.forester.phylogeny.data.PropertiesList() );
            n.getNodeData().getProperties()
                    .addProperty( new org.forester.phylogeny.data.Property( "data:a_considerably_longer_custom_name_for_this_virus_genotype_lineage",
                                                                             "G" + ( i % 2 ), "", "xsd:string",
                                                                             org.forester.phylogeny.data.Property.AppliesTo.NODE ) );
            root.addAsChild( n );
        }
        final Phylogeny phy = new Phylogeny();
        phy.setRoot( root );
        phy.setRooted( true );
        phy.externalNodesHaveChanged();
        return phy;
    }

    private static boolean dialogOk( final Phylogeny phy, final String title, final int rows, final boolean switch_online )
            throws Exception {
        final MainFrame[] mf = new MainFrame[ 1 ];
        SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                .createInstance( new Phylogeny[] { phy }, new Configuration(), "rank dialogs" ) );
        final boolean[] ok = { true };
        final boolean[] seen = { false };
        try {
            SwingUtilities.invokeAndWait( () -> {
                final Timer t = new Timer( 150, null );
                final int[] ticks = { 0 };
                t.addActionListener( e -> {
                    if ( ++ticks[ 0 ] > 60 ) { // ~9 s and no such dialog: close whatever modal is up, fail, never hang
                        t.stop();
                        for( final Window w : Window.getWindows() ) {
                            if ( ( w instanceof JDialog ) && w.isShowing() ) {
                                final JOptionPane any = find( ( (JDialog) w ).getContentPane(), JOptionPane.class );
                                if ( any != null ) {
                                    any.setValue( JOptionPane.CANCEL_OPTION );
                                }
                                w.dispose();
                            }
                        }
                        return;
                    }
                    for( final Window w : Window.getWindows() ) {
                        if ( ( w instanceof JDialog ) && w.isShowing() && title.equals( ( (JDialog) w ).getTitle() ) ) {
                            t.stop();
                            seen[ 0 ] = true;
                            final JOptionPane op = find( ( (JDialog) w ).getContentPane(), JOptionPane.class );
                            final JPanel form = (JPanel) op.getMessage();
                            if ( switch_online ) {
                                final javax.swing.JComboBox<?> rank = find( form, javax.swing.JComboBox.class );
                                // it opened on the tree's own source (the one chosen last), packed to those entries
                                final int before = rank.getPreferredSize().width;
                                final int form_before = form.getWidth();
                                find( form, javax.swing.JRadioButton.class, "NCBI and UniProt for every tip" ).doClick();
                                w.validate();
                                ok[ 0 ] &= ck( ( rank.getPreferredSize().width > before )
                                        && ( rank.getPreferredSize().width > form_before ), title
                                        + ": precondition: the online entries outgrow the dialog as it opened (" + before
                                        + " -> " + rank.getPreferredSize().width + " px, the form " + form_before + ")" );
                                ok[ 0 ] &= ck( before <= form_before, title + ": precondition: the tree's own entries fit ("
                                        + before + " px in a " + form_before + " px form)" );
                                // the DIALOG must be as wide as its content asks: a form wider than the dialog is
                                // laid out at its full width anyway and simply cut off at the dialog's edge
                                ok[ 0 ] &= ck( w.getWidth() >= w.getPreferredSize().width, title
                                        + ": switched online, the dialog grows to its content (" + w.getWidth() + " of "
                                        + w.getPreferredSize().width + " px)" );
                            }
                            ok[ 0 ] &= ck( form.getComponentCount() == rows, title + ": " + rows + " rows, got "
                                    + form.getComponentCount() );
                            int label_h = 0, tallest = 0, stretched = 0, narrow = 0;
                            for( final Component c : form.getComponents() ) {
                                if ( ( c instanceof JLabel ) && ( label_h == 0 ) ) {
                                    label_h = c.getPreferredSize().height;
                                }
                                tallest = Math.max( tallest, c.getPreferredSize().height );
                                if ( c.getHeight() != c.getPreferredSize().height ) {
                                    ++stretched;
                                }
                                if ( c.getWidth() != form.getWidth() ) {
                                    ++narrow;
                                }
                            }
                            ok[ 0 ] &= ck( tallest >= ( 2 * label_h ), title
                                    + ": precondition: the rows differ in height (a label " + label_h + " px, the tallest "
                                    + tallest + ")" );
                            ok[ 0 ] &= ck( stretched == 0, title + ": every row keeps its own height (" + stretched
                                    + " stretched)" );
                            ok[ 0 ] &= ck( narrow == 0, title + ": every row spans the form's width, the combo boxes too ("
                                    + narrow + " narrower)" );
                            ok[ 0 ] &= ck( w.getHeight() <= 800, title + ": the dialog fits a laptop screen ("
                                    + w.getWidth() + " x " + w.getHeight() + " px)" );
                            op.setValue( JOptionPane.CANCEL_OPTION );
                            w.dispose();
                            return;
                        }
                    }
                } );
                t.start();
                if ( switch_online ) {
                    mf[ 0 ].setLastTaxonSourceForTest( TaxonSource.TREE );
                }
                if ( title.startsWith( "Annotate" ) ) {
                    mf[ 0 ].labelCladesByRank();
                }
                else {
                    mf[ 0 ].colorRank();
                }
            } );
        }
        finally {
            SwingUtilities.invokeAndWait( () -> ( (JFrame) mf[ 0 ] ).dispose() );
        }
        return ck( seen[ 0 ], title + ": the dialog opened" ) && ok[ 0 ];
    }

    private static javax.swing.JRadioButton find( final Component c, final Class<javax.swing.JRadioButton> k,
                                                  final String prefix ) {
        if ( ( c instanceof javax.swing.JRadioButton ) && ( (javax.swing.JRadioButton) c ).getText().startsWith( prefix ) ) {
            return (javax.swing.JRadioButton) c;
        }
        if ( c instanceof Container ) {
            for( final Component x : ( (Container) c ).getComponents() ) {
                final javax.swing.JRadioButton r = find( x, k, prefix );
                if ( r != null ) {
                    return r;
                }
            }
        }
        return null;
    }

    private static <T> T find( final Component c, final Class<T> k ) {
        if ( k.isInstance( c ) ) {
            return k.cast( c );
        }
        if ( c instanceof Container ) {
            for( final Component x : ( (Container) c ).getComponents() ) {
                final T r = find( x, k );
                if ( r != null ) {
                    return r;
                }
            }
        }
        return null;
    }

    private static boolean ck( final boolean condition, final String what ) {
        if ( !condition ) {
            System.out.println( "  [RankDialogLayoutTest] " + what );
        }
        return condition;
    }

    private RankDialogLayoutTest() {
        // not instantiable
    }
}

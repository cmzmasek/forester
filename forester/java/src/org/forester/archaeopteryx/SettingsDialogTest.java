// $Id:
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2026 Christian M. Zmasek
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// Contact: phylosoft @ gmail . com

package org.forester.archaeopteryx;

import java.awt.Component;
import java.awt.Container;
import java.awt.GraphicsEnvironment;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

/**
 * Integration test for {@link SettingsDialog}: builds the dialog over a real {@link MainFrame} and checks
 * that it has the expected tabs and that a bound checkbox drives its backing menu item (the {@code doClick}
 * binding that all the apply logic depends on). Guarded to a no-op on a headless box (it needs a toolkit).
 */
public final class SettingsDialogTest {

    public static boolean test() {
        if ( GraphicsEnvironment.isHeadless() ) {
            return true; // GUI integration test; nothing meaningful to do without a display toolkit
        }
        try {
            final Phylogeny phy = new Phylogeny();
            final PhylogenyNode root = new PhylogenyNode();
            root.addAsChild( new PhylogenyNode() );
            root.addAsChild( new PhylogenyNode() );
            phy.setRoot( root );
            phy.externalNodesHaveChanged();
            final Configuration conf = new Configuration( null, false, false, true );
            final MainFrame[] mf = new MainFrame[ 1 ];
            SwingUtilities.invokeAndWait( () -> mf[ 0 ] = MainFrameApplication
                    .createInstance( new Phylogeny[] { phy }, conf, "settings test" ) );
            final boolean[] ok = { true };
            SwingUtilities.invokeAndWait( () -> {
                // the menu bar must carry a "Settings" launcher and no longer the old "Options"/"Type" menus
                final JMenuBar bar = ( (JFrame) mf[ 0 ] ).getJMenuBar();
                boolean has_settings = false;
                boolean has_old = false;
                for ( int i = 0; i < bar.getMenuCount(); ++i ) {
                    final JMenu m = bar.getMenu( i );
                    if ( m == null ) {
                        continue;
                    }
                    if ( "Settings".equals( m.getText() ) ) {
                        has_settings = true;
                    }
                    if ( "Options".equals( m.getText() ) || "Type".equals( m.getText() ) ) {
                        has_old = true;
                    }
                }
                if ( !has_settings || has_old ) {
                    ok[ 0 ] = false;
                }
                final SettingsDialog dlg = new SettingsDialog( mf[ 0 ] );
                dlg.pack();
                final List<JTabbedPane> tabs = new ArrayList<>();
                collect( dlg.getContentPane(), JTabbedPane.class, tabs );
                if ( tabs.isEmpty() || ( tabs.get( 0 ).getTabCount() != 6 ) ) {
                    ok[ 0 ] = false;
                }
                // a dialog checkbox must drive its backing menu item (the doClick binding)
                final boolean before = mf[ 0 ]._show_scale_cbmi.isSelected();
                final JCheckBox cb = findCheckBox( dlg.getContentPane(), mf[ 0 ]._show_scale_cbmi.getText() );
                if ( cb == null ) {
                    ok[ 0 ] = false;
                }
                else {
                    cb.doClick();
                    if ( mf[ 0 ]._show_scale_cbmi.isSelected() == before ) {
                        ok[ 0 ] = false; // toggling the dialog control did not flip the menu item
                    }
                }
                dlg.dispose();
                ( (JFrame) mf[ 0 ] ).dispose();
            } );
            return ok[ 0 ];
        }
        catch ( final Throwable e ) {
            e.printStackTrace();
            return false;
        }
    }

    @SuppressWarnings( "unchecked" )
    private static <T> void collect( final Container c, final Class<T> type, final List<T> out ) {
        for ( final Component comp : c.getComponents() ) {
            if ( type.isInstance( comp ) ) {
                out.add( (T) comp );
            }
            if ( comp instanceof Container ) {
                collect( (Container) comp, type, out );
            }
        }
    }

    private static JCheckBox findCheckBox( final Container c, final String text ) {
        final List<JCheckBox> all = new ArrayList<>();
        collect( c, JCheckBox.class, all );
        for ( final JCheckBox cb : all ) {
            if ( text.equals( cb.getText() ) ) {
                return cb;
            }
        }
        return null;
    }

    public static void main( final String[] args ) {
        System.out.println( test() ? "OK" : "FAILED" );
    }

    private SettingsDialogTest() {
    }
}
